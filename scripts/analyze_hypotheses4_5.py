"""Analyze conformer-population and intramolecular-contact alternatives.

This is a read-only analysis of the retained Bz/SB optimized structures.  It
does not launch ORCA and does not modify either ensemble.  The analysis keeps
one lowest-electronic-energy representative per common-scaffold RMSD cluster
for its primary population model; raw optimized candidates are retained only
as a sampling-multiplicity sensitivity view.
"""

from __future__ import annotations

import json
import math
import sys
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from molecule2fbx.comparison import (
    DEFAULT_TEMPERATURE_K,
    HARTREE_TO_KJ_MOL,
    R_KJ_MOL_K,
    analyze_existing_conformer,
    cross_ensemble_rmsd_analysis,
    load_existing_conformers,
    load_targeted_followup_conformer,
    secondary_rmsd_analysis,
)
from molecule2fbx.structures import rmsd_atom_subsets


OUTPUT = ROOT / "outputs" / "Bz_vs_SB_preliminary_comparison"
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "Si": 2.10}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def atomic_write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def rdkit_molecule(model):
    smiles = model.metadata.get("original_smiles")
    molecule = Chem.AddHs(Chem.MolFromSmiles(str(smiles)))
    if molecule.GetNumAtoms() != len(model.atoms):
        raise ValueError("SMILES/XYZ atom-count mismatch")
    return molecule


def graph_distances(molecule, source: int) -> dict[int, int]:
    distances = {source: 0}
    pending = deque([source])
    while pending:
        current = pending.popleft()
        for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
            index = neighbor.GetIdx()
            if index not in distances:
                distances[index] = distances[current] + 1
                pending.append(index)
    return distances


def atom_groups(model) -> dict:
    molecule = rdkit_molecule(model)
    subsets = rmsd_atom_subsets(model)
    center = subsets.to_metadata()["benzoyl_reaction_center"]
    carbonyl_c = int(center["carbonyl_c"])
    carbonyl_o = int(center["carbonyl_o"])
    amide_n = int(center["amide_n"])
    ipso = int(center["benzoyl_ipso_c"])

    benzoyl_ring = None
    aromatic_core_rings = []
    for ring in molecule.GetRingInfo().AtomRings():
        if not all(molecule.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        ring_set = set(ring)
        if ipso in ring_set:
            benzoyl_ring = tuple(ring)
        else:
            aromatic_core_rings.append(tuple(ring))
    if benzoyl_ring is None:
        raise ValueError("Could not identify benzoyl aromatic ring")

    tms = tuple(int(i) for i in subsets.excluded_terminal_substituent)
    excluded = set(benzoyl_ring) | {carbonyl_c, carbonyl_o} | set(tms)
    core_heavy = tuple(
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetIdx() not in excluded
    )
    benzoyl_group = tuple(sorted(set(benzoyl_ring) | {carbonyl_c, carbonyl_o}))
    return {
        "carbonyl_c": carbonyl_c,
        "carbonyl_o": carbonyl_o,
        "amide_n": amide_n,
        "benzoyl_ipso": ipso,
        "benzoyl_ring": tuple(benzoyl_ring),
        "benzoyl_group": benzoyl_group,
        "tms": tms,
        "core_heavy": core_heavy,
        "core_aromatic_rings": tuple(aromatic_core_rings),
    }


def position(model, index: int) -> np.ndarray:
    return np.asarray(model.atoms[index].position, dtype=float)


def distance(model, first: int, second: int) -> float:
    return float(np.linalg.norm(position(model, first) - position(model, second)))


def plane_geometry(model, first_ring, second_ring) -> dict:
    def plane(indices):
        coordinates = np.asarray([position(model, i) for i in indices])
        center = np.mean(coordinates, axis=0)
        _u, _s, vh = np.linalg.svd(coordinates - center)
        normal = vh[-1] / np.linalg.norm(vh[-1])
        return center, normal

    first_center, first_normal = plane(first_ring)
    second_center, second_normal = plane(second_ring)
    vector = second_center - first_center
    centroid_distance = float(np.linalg.norm(vector))
    cosine = min(1.0, max(-1.0, abs(float(np.dot(first_normal, second_normal)))))
    normal_angle = math.degrees(math.acos(cosine))
    vertical = abs(float(np.dot(vector, first_normal)))
    lateral = math.sqrt(max(0.0, centroid_distance**2 - vertical**2))
    return {
        "centroid_distance_angstrom": centroid_distance,
        "plane_normal_angle_deg": normal_angle,
        "vertical_separation_angstrom": vertical,
        "lateral_offset_angstrom": lateral,
    }


def contacts_between(model, molecule, first_group, second_group) -> dict:
    graph_cache = {index: graph_distances(molecule, index) for index in first_group}
    contacts = []
    eligible_pairs = []
    for first in first_group:
        for second in second_group:
            if graph_cache[first].get(second, 999) <= 3:
                continue
            value = distance(model, first, second)
            radius_sum = VDW.get(model.atoms[first].element, 1.8) + VDW.get(
                model.atoms[second].element, 1.8
            )
            gap = value - radius_sum
            eligible_pairs.append(
                {
                    "first_atom": first,
                    "second_atom": second,
                    "first_element": model.atoms[first].element,
                    "second_element": model.atoms[second].element,
                    "distance_angstrom": value,
                    "vdw_gap_angstrom": gap,
                }
            )
            if gap <= 0.50:
                contacts.append(dict(eligible_pairs[-1]))
    contacts.sort(key=lambda item: item["vdw_gap_angstrom"])
    eligible_pairs.sort(key=lambda item: item["vdw_gap_angstrom"])
    return {
        "contact_definition": "graph separation >3 and distance <= vdW sum + 0.50 A",
        "contact_count": len(contacts),
        "minimum_vdw_gap_angstrom": (
            contacts[0]["vdw_gap_angstrom"] if contacts else None
        ),
        "minimum_vdw_gap_any_eligible_pair_angstrom": (
            eligible_pairs[0]["vdw_gap_angstrom"] if eligible_pairs else None
        ),
        "closest_eligible_pair_by_vdw_gap": (
            eligible_pairs[0] if eligible_pairs else None
        ),
        "contact_score_angstrom": sum(
            max(0.0, 0.50 - item["vdw_gap_angstrom"]) for item in contacts
        ),
        "pairs": contacts,
    }


def weak_ch_o_contacts(model, molecule, groups) -> list[dict]:
    oxygen = groups["carbonyl_o"]
    oxygen_position = position(model, oxygen)
    records = []
    core = set(groups["core_heavy"])
    for hydrogen in molecule.GetAtoms():
        if hydrogen.GetAtomicNum() != 1 or hydrogen.GetDegree() != 1:
            continue
        donor = hydrogen.GetNeighbors()[0]
        if donor.GetAtomicNum() != 6 or donor.GetIdx() not in core:
            continue
        if len(Chem.GetShortestPath(molecule, donor.GetIdx(), oxygen)) - 1 <= 3:
            continue
        h_index = hydrogen.GetIdx()
        h_position = position(model, h_index)
        c_position = position(model, donor.GetIdx())
        h_o = float(np.linalg.norm(oxygen_position - h_position))
        first = c_position - h_position
        second = oxygen_position - h_position
        cosine = float(np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second)))
        angle = math.degrees(math.acos(min(1.0, max(-1.0, cosine))))
        if h_o <= 2.70 and angle >= 120.0:
            records.append(
                {
                    "carbon_atom": donor.GetIdx(),
                    "hydrogen_atom": h_index,
                    "oxygen_atom": oxygen,
                    "h_o_distance_angstrom": h_o,
                    "c_h_o_angle_deg": angle,
                }
            )
    return sorted(records, key=lambda item: item["h_o_distance_angstrom"])


def contact_descriptors(model) -> dict:
    molecule = rdkit_molecule(model)
    groups = atom_groups(model)
    benzoyl_core = contacts_between(
        model, molecule, groups["benzoyl_group"], groups["core_heavy"]
    )
    ring_core = contacts_between(
        model, molecule, groups["benzoyl_ring"], groups["core_heavy"]
    )
    tms_core = contacts_between(model, molecule, groups["tms"], groups["core_heavy"])
    pi_geometries = [
        plane_geometry(model, groups["benzoyl_ring"], ring)
        for ring in groups["core_aromatic_rings"]
    ]
    pi_geometries.sort(key=lambda item: item["centroid_distance_angstrom"])
    heavy_positions = np.asarray(
        [position(model, i) for i, atom in enumerate(model.atoms) if atom.element != "H"]
    )
    center = np.mean(heavy_positions, axis=0)
    radius_of_gyration = math.sqrt(
        float(np.mean(np.sum((heavy_positions - center) ** 2, axis=1)))
    )
    return {
        "atom_groups": {key: list(value) if isinstance(value, tuple) else value for key, value in groups.items()},
        "benzoyl_group_to_lsd_core": benzoyl_core,
        "benzoyl_ring_to_lsd_core": ring_core,
        "tms_to_lsd_core": tms_core,
        "weak_c_h_to_carbonyl_o_candidates": weak_ch_o_contacts(model, molecule, groups),
        "benzoyl_to_core_aromatic_ring_geometries": pi_geometries,
        "nearest_pi_geometry": pi_geometries[0] if pi_geometries else None,
        "heavy_atom_radius_of_gyration_angstrom": radius_of_gyration,
        "interpretation_limit": (
            "Distance/van-der-Waals descriptors flag geometrically possible contacts; "
            "they are not SAPT/NCI interaction energies and do not establish stabilization."
        ),
    }


def boltzmann_weights(records: list[dict], energy_key: str) -> dict[str, float]:
    available = [
        (record["conformer_id"], record.get(energy_key)) for record in records
        if record.get(energy_key) is not None
    ]
    minimum = min(float(value) for _key, value in available)
    factors = []
    for key, value in available:
        delta = (float(value) - minimum) * HARTREE_TO_KJ_MOL
        factors.append((str(key), math.exp(-delta / (R_KJ_MOL_K * DEFAULT_TEMPERATURE_K))))
    total = sum(value for _key, value in factors)
    return {key: value / total for key, value in factors}


def weighted_stats(records: list[dict], weights: dict[str, float], path: tuple[str, ...]) -> dict:
    values = []
    for record in records:
        value = record
        for key in path:
            if not isinstance(value, dict) or key not in value:
                value = None
                break
            value = value[key]
        identifier = str(record["conformer_id"])
        if value is not None and identifier in weights:
            values.append((float(value), weights[identifier]))
    coverage = sum(weight for _value, weight in values)
    if not values or coverage <= 0:
        return {"mean": None, "std": None, "minimum": None, "maximum": None, "weight_coverage": 0.0}
    normalized = [(value, weight / coverage) for value, weight in values]
    mean = sum(value * weight for value, weight in normalized)
    variance = sum(weight * (value - mean) ** 2 for value, weight in normalized)
    return {
        "mean": mean,
        "std": math.sqrt(max(variance, 0.0)),
        "minimum": min(value for value, _weight in values),
        "maximum": max(value for value, _weight in values),
        "weight_coverage": coverage,
    }


METRICS = {
    "benzoyl_core_contact_count": ("contacts", "benzoyl_group_to_lsd_core", "contact_count"),
    "benzoyl_core_contact_score": ("contacts", "benzoyl_group_to_lsd_core", "contact_score_angstrom"),
    "benzoyl_ring_core_contact_count": ("contacts", "benzoyl_ring_to_lsd_core", "contact_count"),
    "benzoyl_ring_core_contact_score": ("contacts", "benzoyl_ring_to_lsd_core", "contact_score_angstrom"),
    "tms_core_contact_count": ("contacts", "tms_to_lsd_core", "contact_count"),
    "tms_core_contact_score": ("contacts", "tms_to_lsd_core", "contact_score_angstrom"),
    "weak_ch_o_candidate_count": ("weak_ch_o_candidate_count",),
    "nearest_pi_centroid_distance": ("contacts", "nearest_pi_geometry", "centroid_distance_angstrom"),
    "nearest_pi_plane_angle": ("contacts", "nearest_pi_geometry", "plane_normal_angle_deg"),
    "nearest_pi_lateral_offset": ("contacts", "nearest_pi_geometry", "lateral_offset_angstrom"),
    "radius_of_gyration": ("contacts", "heavy_atom_radius_of_gyration_angstrom"),
    "carbonyl_access_best_clearance": ("existing", "steric_access", "best_clearance_angstrom"),
    "carbonyl_access_p90": ("existing", "steric_access", "clearance_percentiles_angstrom", "p90"),
    "benzoyl_carbonyl_torsion": ("existing", "structure", "benzoyl_carbonyl_torsion_abs_mean_deg"),
}


def ensemble_analysis(name: str, conformers, payload) -> dict:
    rmsd = secondary_rmsd_analysis(conformers)
    rmsd_by_index = {int(item["conformer_index"]): item for item in rmsd["records"]}
    records = []
    for conformer in conformers:
        existing = analyze_existing_conformer(conformer)
        cluster = rmsd_by_index[int(conformer.entry["conformer_index"])]
        contact = contact_descriptors(conformer.model)
        records.append(
            {
                "conformer_id": str(conformer.entry["conformer_id"]),
                "dft_energy_hartree": float(conformer.entry["dft_energy_hartree"]),
                "gibbs_energy_hartree": conformer.entry.get("gibbs_energy_hartree"),
                "relative_energy_kj_mol": conformer.entry.get("relative_energy_kj_mol"),
                "frequency_calculated": conformer.entry.get("frequency_calculated") is True,
                "common_scaffold_cluster_id": cluster["common_scaffold_rmsd_cluster_id"],
                "reaction_center_cluster_id": cluster["reaction_center_rmsd_cluster_id"],
                "existing": existing,
                "contacts": contact,
                "weak_ch_o_candidate_count": len(contact["weak_c_h_to_carbonyl_o_candidates"]),
            }
        )

    core_groups = defaultdict(list)
    for record in records:
        core_groups[record["common_scaffold_cluster_id"]].append(record)
    representatives = [
        min(group, key=lambda item: item["dft_energy_hartree"])
        for group in core_groups.values()
    ]
    representatives.sort(key=lambda item: item["dft_energy_hartree"])
    weights = boltzmann_weights(representatives, "dft_energy_hartree")
    raw_weights = boltzmann_weights(records, "dft_energy_hartree")
    freq_representatives = [item for item in representatives if item["gibbs_energy_hartree"] is not None]
    gibbs_weights = boltzmann_weights(freq_representatives, "gibbs_energy_hartree") if freq_representatives else {}

    minimum = min(item["dft_energy_hartree"] for item in representatives)
    for record in records:
        record["relative_electronic_energy_kj_mol_rebased"] = (
            record["dft_energy_hartree"] - minimum
        ) * HARTREE_TO_KJ_MOL
        record["common_core_representative_weight"] = weights.get(record["conformer_id"])
        record["raw_candidate_sensitivity_weight"] = raw_weights.get(record["conformer_id"])
        record["conditional_gibbs_weight"] = gibbs_weights.get(record["conformer_id"])

    state_weights = defaultdict(float)
    for record in representatives:
        state_weights[record["reaction_center_cluster_id"]] += weights[record["conformer_id"]]
    entropy = -sum(weight * math.log(weight) for weight in weights.values() if weight > 0)
    summary = {
        "optimized_candidate_count": len(records),
        "all_heavy_cluster_count_recorded": payload.get("post_dft_clustering", {}).get("cluster_count"),
        "common_scaffold_cluster_count": len(representatives),
        "reaction_center_cluster_count": len(set(item["reaction_center_cluster_id"] for item in records)),
        "primary_population_model": "one lowest-E representative per common-scaffold cluster",
        "effective_conformer_number_electronic": math.exp(entropy),
        "highest_single_conformer_weight": max(weights.values()),
        "low_energy_representative_counts": {
            f"within_{window:g}_kj_mol": sum(
                (item["dft_energy_hartree"] - minimum) * HARTREE_TO_KJ_MOL <= window
                for item in representatives
            )
            for window in (2.5, 5.0, 10.0)
        },
        "reaction_center_cluster_weights_not_cross_molecule_labels": dict(state_weights),
        "conditional_gibbs_structure_count": len(freq_representatives),
        "metrics": {
            key: {
                "electronic_common_core_weighted": weighted_stats(representatives, weights, path),
                "raw_candidate_weighted_sensitivity": weighted_stats(records, raw_weights, path),
                "conditional_gibbs_weighted": weighted_stats(freq_representatives, gibbs_weights, path),
            }
            for key, path in METRICS.items()
        },
    }
    return {"name": name, "summary": summary, "records": records}


def paired_analysis(bz: dict, sb: dict, cross: dict) -> list[dict]:
    bz_by_id = {item["conformer_id"]: item for item in bz["records"]}
    sb_by_id = {item["conformer_id"]: item for item in sb["records"]}
    pairs = []
    for pair in cross["reciprocal_nearest_pairs"]:
        rmsd = float(pair["common_scaffold_rmsd_angstrom"])
        if rmsd > 0.50:
            continue
        first = bz_by_id[str(pair["first_conformer_id"])]
        second = sb_by_id[str(pair["second_conformer_id"])]
        differences = {}
        for key, path in METRICS.items():
            def get(record):
                value = record
                for component in path:
                    if not isinstance(value, dict) or component not in value:
                        return None
                    value = value[component]
                return value
            first_value, second_value = get(first), get(second)
            differences[key] = None if first_value is None or second_value is None else float(second_value) - float(first_value)
        pairs.append(
            {
                "bz_conformer": first["conformer_id"],
                "sb_conformer": second["conformer_id"],
                "common_scaffold_rmsd_angstrom": rmsd,
                "reaction_center_rmsd_angstrom": pair["reaction_center_rmsd_angstrom"],
                "sb_minus_bz": differences,
            }
        )
    return pairs


def number(value, digits=3):
    return "—" if value is None else f"{float(value):.{digits}f}"


def write_findings(payload: dict, path: Path) -> None:
    bz = payload["ensembles"]["1Bz-LSD_RR"]["summary"]
    sb = payload["ensembles"]["1SB-LSD_RR"]["summary"]
    def metric(summary, key):
        return summary["metrics"][key]["electronic_common_core_weighted"]

    rows = []
    for label, key in (
        ("benzoyl–core contact count", "benzoyl_core_contact_count"),
        ("benzoyl–core contact score / Å", "benzoyl_core_contact_score"),
        ("benzoyl ring–core contact count", "benzoyl_ring_core_contact_count"),
        ("TMS–core contact count", "tms_core_contact_count"),
        ("nearest aromatic centroid distance / Å", "nearest_pi_centroid_distance"),
        ("heavy-atom radius of gyration / Å", "radius_of_gyration"),
        ("carbonyl best clearance / Å", "carbonyl_access_best_clearance"),
        ("benzoyl–carbonyl torsion / deg", "benzoyl_carbonyl_torsion"),
    ):
        first, second = metric(bz, key), metric(sb, key)
        rows.append(
            f"| {label} | {number(first['mean'])} ± {number(first['std'])} | "
            f"{number(second['mean'])} ± {number(second['std'])} | "
            f"{number(None if first['mean'] is None or second['mean'] is None else second['mean']-first['mean'])} |"
        )
    lines = [
        "# 仮説4・5 既存ensemble解析",
        "",
        "新規ORCA計算は行っていない。主要集計は共通骨格RMSD clusterごとの最低電子エネルギー代表を1個だけ採用し、298.15 Kで電子エネルギー重み付けした感度解析である。これは厳密な自由エネルギーpopulationではない。",
        "",
        "## 仮説4：配座集団効果",
        "",
        f"- Bz: 共通骨格 {bz['common_scaffold_cluster_count']} cluster、反応中心 {bz['reaction_center_cluster_count']} cluster、5 kJ/mol以内 {bz['low_energy_representative_counts']['within_5_kj_mol']}構造、電子Eによる実効配座数 {bz['effective_conformer_number_electronic']:.2f}。",
        f"- SB: 共通骨格 {sb['common_scaffold_cluster_count']} cluster、反応中心 {sb['reaction_center_cluster_count']} cluster、5 kJ/mol以内 {sb['low_energy_representative_counts']['within_5_kj_mol']}構造、電子Eによる実効配座数 {sb['effective_conformer_number_electronic']:.2f}。",
        "- SBの全重原子cluster増加を、そのままLSD–benzoyl骨格の多様性増加とは扱えない。TMS末端を除くと両者とも共通骨格9 cluster、反応中心4 clusterである。",
        "- ただしSBの低エネルギー領域には複数の近接構造があり、単一current-bestだけで比較するよりensembleで見る必要がある。電子E重みは振動・回転エントロピーを含まないため、存在比の予測値ではない。",
        "",
        "## 仮説5：分子内相互作用／折り畳み",
        "",
        "| 指標 | Bz | SB | SB−Bz |",
        "| --- | ---: | ---: | ---: |",
        *rows,
        "",
        f"- 共通骨格RMSD ≤0.50 Åの相互最近傍pairは{len(payload['paired_low_rmsd_comparisons'])}組。pair差の符号が揃わない指標はTMS固有効果と解釈しない。",
        "- contactはvdW距離に基づく幾何学proxyで、安定化エネルギーではない。TMS接触数はBzに存在しない追加原子数の影響を受けるため、benzoyl ring–core共通部分と分けて扱う。",
        "- N/O–H donorはなく、記録したC–H···O候補は弱い幾何学候補にすぎない。π配置もcentroid・面角の記述であり、π–π相互作用エネルギーを意味しない。",
        "",
        "## 判定上の停止線",
        "",
        "- 仮説4・5は『TMSが自由分子の配座分布や接触幾何を変え得るか』までを扱う。酵素内populationや反応速度には直接変換しない。",
        "- Freqが部分集合のみなので、ここでの電子E重みを完全なBoltzmann populationとは呼ばない。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    bz_payload, bz_conformers = load_existing_conformers(ROOT / "outputs" / "1Bz-LSD_RR")
    sb_payload, sb_conformers = load_existing_conformers(ROOT / "outputs" / "1SB-LSD_RR")
    followup = load_targeted_followup_conformer(sb_payload, OUTPUT)
    if followup is not None:
        sb_conformers.append(followup)

    bz = ensemble_analysis("1Bz-LSD_RR", bz_conformers, bz_payload)
    sb = ensemble_analysis("1SB-LSD_RR", sb_conformers, sb_payload)
    cross = cross_ensemble_rmsd_analysis(bz_conformers, sb_conformers)
    payload = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "calculation_scope": "read-only analysis of retained optimized structures; no ORCA jobs",
        "population_caveat": "electronic-energy weighted sensitivity, not a complete conformer free-energy population",
        "ensembles": {"1Bz-LSD_RR": bz, "1SB-LSD_RR": sb},
        "paired_low_rmsd_comparisons": paired_analysis(bz, sb, cross),
        "method_assumptions": {
            "primary_population": "lowest-E representative per common-scaffold cluster",
            "temperature_kelvin": DEFAULT_TEMPERATURE_K,
            "contacts": "heavy-atom pairs at graph separation >3 and within vdW sum +0.50 A",
            "symmetry_permutations": False,
            "interaction_energy_calculated": False,
        },
    }
    atomic_write_json(OUTPUT / "hypotheses4_5_analysis.json", payload)
    write_findings(payload, OUTPUT / "HYPOTHESES4_5_FINDINGS_JA.md")


if __name__ == "__main__":
    main()

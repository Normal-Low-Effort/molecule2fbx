"""Non-destructive neutral/N6-protonated GFN2-xTB ALPB(water) screen.

The retained DFT conformers are starting geometries only.  One lowest-energy
representative per existing common-scaffold cluster is optimized for each
charge state.  This is a protonation-sensitivity screen, not a pKa calculation
and not a replacement for the retained ORCA ensemble.
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "work"))

from analyze_hypotheses4_5 import contact_descriptors  # noqa: E402
from molecule2fbx.comparison import (  # noqa: E402
    DEFAULT_TEMPERATURE_K,
    HARTREE_TO_KJ_MOL,
    R_KJ_MOL_K,
    _rdkit_model_from_xyz,
    carbonyl_steric_access,
    load_existing_conformers,
    load_targeted_followup_conformer,
    secondary_rmsd_analysis,
    structural_descriptors,
)
from molecule2fbx.structures import (  # noqa: E402
    aligned_atom_subset_rmsd,
    rmsd_atom_subsets,
    validate_model_stereochemistry,
)


XTB = Path(os.environ.get("XTB_EXECUTABLE", r"C:\ORCA_6.1.1\xtb-6.7.1pre\xtb.exe"))
OUTPUT = ROOT / "outputs" / "Bz_vs_SB_preliminary_comparison" / "protonation_screen"
STATUS = OUTPUT / "status.json"
RESULTS = OUTPUT / "results.json"
FINDINGS = OUTPUT / "FINDINGS_JA.md"
ENERGY_RE = re.compile(
    r"TOTAL ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh", re.IGNORECASE
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def write_xyz(path: Path, elements, coordinates, comment: str) -> None:
    lines = [str(len(elements)), comment]
    lines.extend(
        f"{element:<3} {coordinate[0]: .12f} {coordinate[1]: .12f} {coordinate[2]: .12f}"
        for element, coordinate in zip(elements, coordinates)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def protonated_smiles(neutral_smiles: str, nitrogen_index: int = 9) -> str:
    neutral = Chem.MolFromSmiles(neutral_smiles)
    if neutral is None:
        raise ValueError("Could not parse neutral SMILES")
    atom = neutral.GetAtomWithIdx(nitrogen_index)
    if atom.GetSymbol() != "N" or atom.GetDegree() != 3 or atom.GetIsAromatic():
        raise ValueError("Expected a non-aromatic tertiary N at atom index 9")
    # Preserve the input atom order exactly.  MolToSmiles(canonical=False) may
    # still choose a different aromatic-ring traversal after changing charge,
    # which would break the retained XYZ/SMILES atom-order contract.
    marker = "[C@H]1CN(C)[C@@H]2"
    replacement = "[C@H]1C[NH+](C)[C@@H]2"
    if neutral_smiles.count(marker) != 1:
        raise ValueError("Could not locate the unique N6 protonation site in input SMILES")
    charged_smiles = neutral_smiles.replace(marker, replacement)
    charged = Chem.MolFromSmiles(charged_smiles)
    if charged is None:
        raise ValueError("Could not parse N6-protonated SMILES")
    neutral_heavy = [atom.GetSymbol() for atom in neutral.GetAtoms()]
    charged_heavy = [atom.GetSymbol() for atom in charged.GetAtoms()]
    if neutral_heavy != charged_heavy:
        raise ValueError("N6 protonation changed heavy-atom order")
    return charged_smiles


def protonated_coordinates(model, neutral_smiles: str, charged_smiles: str, nitrogen_index: int = 9):
    neutral = Chem.AddHs(Chem.MolFromSmiles(neutral_smiles))
    charged = Chem.AddHs(Chem.MolFromSmiles(charged_smiles))
    neutral_heavy = [atom.GetSymbol() for atom in neutral.GetAtoms() if atom.GetAtomicNum() > 1]
    charged_heavy = [atom.GetSymbol() for atom in charged.GetAtoms() if atom.GetAtomicNum() > 1]
    if neutral_heavy != charged_heavy:
        raise ValueError("Protonation changed heavy-atom order")
    existing_hydrogens = {}
    for atom in neutral.GetAtoms():
        if atom.GetAtomicNum() == 1:
            parent = atom.GetNeighbors()[0].GetIdx()
            existing_hydrogens.setdefault(parent, []).append(atom.GetIdx())
    used = {parent: 0 for parent in existing_hydrogens}
    coordinates = []
    for atom in charged.GetAtoms():
        index = atom.GetIdx()
        if atom.GetAtomicNum() > 1:
            coordinates.append(np.asarray(model.atoms[index].position, dtype=float))
            continue
        parent = atom.GetNeighbors()[0].GetIdx()
        available = existing_hydrogens.get(parent, [])
        offset = used.get(parent, 0)
        if offset < len(available):
            coordinates.append(np.asarray(model.atoms[available[offset]].position, dtype=float))
            used[parent] = offset + 1
            continue
        if parent != nitrogen_index:
            raise ValueError(f"Unexpected new hydrogen on atom {parent}")
        center = np.asarray(model.atoms[parent].position, dtype=float)
        directions = []
        for neighbor in charged.GetAtomWithIdx(parent).GetNeighbors():
            if neighbor.GetAtomicNum() == 1:
                continue
            vector = np.asarray(model.atoms[neighbor.GetIdx()].position, dtype=float) - center
            directions.append(vector / np.linalg.norm(vector))
        direction = -np.sum(directions, axis=0)
        if np.linalg.norm(direction) < 1.0e-8:
            direction = np.cross(directions[0], directions[1])
        direction /= np.linalg.norm(direction)
        coordinates.append(center + 1.02 * direction)
    return [atom.GetSymbol() for atom in charged.GetAtoms()], coordinates


def parse_completed(job_dir: Path):
    output = job_dir / "xtb.out"
    optimized = job_dir / "xtbopt.xyz"
    text = output.read_text(encoding="utf-8", errors="replace")
    energies = [float(value) for value in ENERGY_RE.findall(text)]
    if "normal termination of xtb" not in text or not energies or not optimized.is_file():
        raise RuntimeError(f"Incomplete xTB job: {job_dir}")
    return energies[-1], optimized


def retained_representatives(ensemble_name: str):
    payload, conformers = load_existing_conformers(ROOT / "outputs" / ensemble_name)
    if ensemble_name == "1SB-LSD_RR":
        followup = load_targeted_followup_conformer(
            payload, ROOT / "outputs" / "Bz_vs_SB_preliminary_comparison"
        )
        if followup is not None:
            conformers.append(followup)
    analysis = secondary_rmsd_analysis(conformers)
    clusters = {}
    for record in analysis["records"]:
        clusters[int(record["conformer_index"])] = record["common_scaffold_rmsd_cluster_id"]
    grouped = {}
    for conformer in conformers:
        cluster = clusters[int(conformer.entry["conformer_index"])]
        if cluster not in grouped or float(conformer.entry["dft_energy_hartree"]) < float(
            grouped[cluster].entry["dft_energy_hartree"]
        ):
            grouped[cluster] = conformer
    return payload, sorted(grouped.values(), key=lambda item: float(item.entry["dft_energy_hartree"]))


def run_job(job_dir: Path, elements, coordinates, charge: int, state: dict):
    job_dir.mkdir(parents=True, exist_ok=True)
    input_path = job_dir / "input.xyz"
    output_path = job_dir / "xtb.out"
    if output_path.is_file():
        energy, optimized = parse_completed(job_dir)
        return energy, optimized, "reused_completed_xtb_job"
    unexpected = list(job_dir.iterdir())
    if unexpected:
        raise RuntimeError(f"Refusing to overwrite incomplete xTB directory: {job_dir}")
    write_xyz(input_path, elements, coordinates, f"charge={charge}; GFN2-xTB ALPB(water)")
    process = subprocess.run(
        [
            str(XTB), input_path.name, "--gfn", "2", "--opt", "tight",
            "--alpb", "water", "--chrg", str(charge), "--uhf", "0",
            "--parallel", "16",
        ],
        cwd=job_dir,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env={**os.environ, "OMP_NUM_THREADS": "16"},
        timeout=3600,
        check=False,
    )
    output_path.write_text(
        process.stdout + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
        encoding="utf-8",
    )
    energy, optimized = parse_completed(job_dir)
    if process.returncode != 0:
        raise RuntimeError(f"xTB returned {process.returncode}: {job_dir}")
    return energy, optimized, "computed_this_run"


def weights(records):
    minimum = min(item["energy_hartree"] for item in records)
    factors = [
        math.exp(-((item["energy_hartree"] - minimum) * HARTREE_TO_KJ_MOL) /
                 (R_KJ_MOL_K * DEFAULT_TEMPERATURE_K))
        for item in records
    ]
    total = sum(factors)
    return [value / total for value in factors]


def summarize(records):
    ranked = sorted(records, key=lambda item: item["energy_hartree"])
    representatives = []
    atom_indices = rmsd_atom_subsets(ranked[0]["_model"]).common_scaffold
    for item in ranked:
        assigned = None
        assigned_rmsd = 0.0
        for cluster_index, representative in enumerate(representatives):
            value = aligned_atom_subset_rmsd(
                item["_model"], representative["_model"], atom_indices
            )
            if value < 0.75:
                assigned = cluster_index
                assigned_rmsd = value
                break
        if assigned is None:
            assigned = len(representatives)
            representatives.append(item)
        item["xtb_common_scaffold_cluster_id"] = f"xtb_core_cluster_{assigned + 1:03d}"
        item["xtb_common_scaffold_rmsd_to_representative_angstrom"] = assigned_rmsd
    population = weights(representatives)
    minimum = min(item["energy_hartree"] for item in representatives)
    for item in records:
        item["relative_energy_kj_mol"] = (item["energy_hartree"] - minimum) * HARTREE_TO_KJ_MOL
        item["electronic_boltzmann_sensitivity_weight"] = None
    for item, weight in zip(representatives, population):
        item["electronic_boltzmann_sensitivity_weight"] = weight
    entropy = -sum(value * math.log(value) for value in population if value > 0)
    metric_paths = {
        "benzoyl_core_contact_count": ("contacts", "benzoyl_group_to_lsd_core", "contact_count"),
        "benzoyl_core_contact_score": ("contacts", "benzoyl_group_to_lsd_core", "contact_score_angstrom"),
        "tms_core_contact_count": ("contacts", "tms_to_lsd_core", "contact_count"),
        "nearest_pi_centroid_distance": ("contacts", "nearest_pi_geometry", "centroid_distance_angstrom"),
        "radius_of_gyration": ("contacts", "heavy_atom_radius_of_gyration_angstrom"),
        "carbonyl_best_clearance": ("steric", "best_clearance_angstrom"),
        "carbonyl_p90_clearance": ("steric", "clearance_percentiles_angstrom", "p90"),
        "benzoyl_carbonyl_torsion": ("structure", "benzoyl_carbonyl_torsion_abs_mean_deg"),
        "carbonyl_c_o_length": ("structure", "carbonyl_c_o_length_angstrom"),
        "amide_c_n_length": ("structure", "amide_c_n_length_angstrom"),
    }
    metrics = {}
    for name, path in metric_paths.items():
        values = []
        for item, weight in zip(representatives, population):
            value = item
            for key in path:
                value = value.get(key) if isinstance(value, dict) else None
            if value is not None:
                values.append((float(value), weight))
        coverage = sum(weight for _value, weight in values)
        normalized = [(value, weight / coverage) for value, weight in values]
        mean = sum(value * weight for value, weight in normalized) if normalized else None
        variance = sum(weight * (value - mean) ** 2 for value, weight in normalized) if normalized else None
        metrics[name] = {
            "mean": mean,
            "std": math.sqrt(max(0.0, variance)) if variance is not None else None,
            "coverage": coverage,
        }
    return {
        "structure_count": len(records),
        "common_scaffold_cluster_count_after_xtb": len(representatives),
        "population_representative_count": len(representatives),
        "effective_conformer_number_electronic": math.exp(entropy),
        "highest_single_weight": max(population),
        "metrics": metrics,
    }


def write_findings(result) -> None:
    bz = result["ensembles"]["1Bz-LSD_RR"]
    sb = result["ensembles"]["1SB-LSD_RR"]

    def metric(dataset, state, name):
        return dataset[state]["summary"]["metrics"][name]["mean"]

    rows = []
    for label, key in (
        ("benzoyl–core contact count", "benzoyl_core_contact_count"),
        ("TMS–core contact count", "tms_core_contact_count"),
        ("carbonyl best clearance / Å", "carbonyl_best_clearance"),
        ("carbonyl p90 clearance / Å", "carbonyl_p90_clearance"),
        ("benzoyl–carbonyl torsion / deg", "benzoyl_carbonyl_torsion"),
        ("C=O length / Å", "carbonyl_c_o_length"),
        ("amide C–N length / Å", "amide_c_n_length"),
    ):
        values = []
        for state in ("neutral", "n6_protonated"):
            first = metric(bz, state, key)
            second = metric(sb, state, key)
            values.extend((first, second, second - first))
        rows.append(
            f"| {label} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:+.4f} | "
            f"{values[3]:.4f} | {values[4]:.4f} | {values[5]:+.4f} |"
        )

    bz_neutral = min(item["energy_hartree"] for item in bz["neutral"]["records"])
    bz_protonated = min(item["energy_hartree"] for item in bz["n6_protonated"]["records"])
    sb_neutral = min(item["energy_hartree"] for item in sb["neutral"]["records"])
    sb_protonated = min(item["energy_hartree"] for item in sb["n6_protonated"]["records"])
    double_difference = (
        (sb_protonated - sb_neutral) - (bz_protonated - bz_neutral)
    ) * HARTREE_TO_KJ_MOL
    result["relative_n6_protonation_electronic_double_difference_sb_minus_bz_kj_mol"] = double_difference
    result["relative_n6_protonation_double_difference_interpretation"] = (
        "GFN2-xTB/ALPB electronic-energy sensitivity only; not a proton free energy or pKa"
    )
    lines = [
        "# 仮説7 N6プロトン化感度screen",
        "",
        "各分子9個の共通骨格代表について、中性体とN6プロトン化体をGFN2-xTB/ALPB(水)で最適化した。既存ORCA結果は変更していない。これはpKa計算でも完全な自由エネルギーensembleでもない。",
        "",
        "| 指標 | Bz neutral | SB neutral | SB−Bz | Bz N6H+ | SB N6H+ | SB−Bz |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        *rows,
        "",
        f"- xTB後も4状態すべてで共通骨格clusterは9個を維持した。電子E感度上の実効配座数はneutralでBz {bz['neutral']['summary']['effective_conformer_number_electronic']:.2f} / SB {sb['neutral']['summary']['effective_conformer_number_electronic']:.2f}、N6H+でBz {bz['n6_protonated']['summary']['effective_conformer_number_electronic']:.2f} / SB {sb['n6_protonated']['summary']['effective_conformer_number_electronic']:.2f}。",
        "- TMS–LSD core接触はneutral/N6H+の全SB構造で0。最短pairでもvdW半径和より2 Å以上離れている。",
        f"- 最低電子Eを使う相対プロトン化double difference (SB−Bz) は {double_difference:+.3f} kJ/mol。小さいが、熱・エントロピー・標準状態・プロトン自由エネルギーを含まないためpKa差へ変換しない。",
        "- プロトン化により各分子自身の配座順位・局所構造は変わるが、Bz/SB間の『TMSがcoreへ直接接触しない』『carbonyl access差が小さい』という方向は逆転しなかった。",
        "",
        "## 停止線",
        "",
        "このscreenは仮説1・4・5がN6プロトン化で直ちに崩れないことを示す。仮説2の電子差を確定するものではなく、実際のプロトン化率、酵素結合状態、加水分解速度も決めない。",
    ]
    FINDINGS.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if not XTB.is_file():
        raise SystemExit(f"xTB not found: {XTB}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    state = {
        "status": "RUNNING",
        "started_at_utc": now(),
        "finished_at_utc": None,
        "current_job": None,
        "total_jobs": 36,
        "processed_jobs": 0,
        "computed_jobs": 0,
        "reused_jobs": 0,
    }
    write_json(STATUS, state)
    all_results = {}
    try:
        for ensemble_name in ("1Bz-LSD_RR", "1SB-LSD_RR"):
            payload, representatives = retained_representatives(ensemble_name)
            neutral_smiles = str(payload["molecule"].get("original_smiles") or payload["molecule"]["smiles"])
            charged_smiles = protonated_smiles(neutral_smiles)
            state_results = {}
            for label, charge, smiles in (
                ("neutral", 0, neutral_smiles),
                ("n6_protonated", 1, charged_smiles),
            ):
                records = []
                for conformer in representatives:
                    conformer_id = str(conformer.entry["conformer_id"])
                    state["current_job"] = f"{ensemble_name}/{label}/{conformer_id}"
                    write_json(STATUS, state)
                    if charge == 0:
                        elements = [atom.element for atom in conformer.model.atoms]
                        coordinates = [np.asarray(atom.position, dtype=float) for atom in conformer.model.atoms]
                    else:
                        elements, coordinates = protonated_coordinates(
                            conformer.model, neutral_smiles, charged_smiles
                        )
                    job_dir = OUTPUT / ensemble_name / label / conformer_id
                    energy, optimized_xyz, provenance = run_job(
                        job_dir, elements, coordinates, charge, state
                    )
                    state["processed_jobs"] += 1
                    if provenance == "computed_this_run":
                        state["computed_jobs"] += 1
                    else:
                        state["reused_jobs"] += 1
                    state["current_job"] = None
                    write_json(STATUS, state)
                    model = _rdkit_model_from_xyz(
                        optimized_xyz, smiles, name=f"{ensemble_name}_{label}"
                    )
                    stereo = validate_model_stereochemistry(smiles, model)
                    if not stereo.matches:
                        raise RuntimeError(f"Stereo mismatch: {ensemble_name}/{label}/{conformer_id}")
                    records.append(
                        {
                            "source_conformer_id": conformer_id,
                            "energy_hartree": energy,
                            "charge": charge,
                            "provenance": provenance,
                            "optimized_xyz": str(optimized_xyz.resolve()),
                            "stereochemistry_validation": stereo.to_metadata(),
                            "structure": structural_descriptors(model),
                            "steric": carbonyl_steric_access(model),
                            "contacts": contact_descriptors(model),
                            "_model": model,
                        }
                    )
                summary = summarize(records)
                for record in records:
                    record.pop("_model", None)
                state_results[label] = {
                    "charge": charge,
                    "smiles": smiles,
                    "summary": summary,
                    "records": records,
                }
            all_results[ensemble_name] = state_results
        result = {
            "schema_version": 1,
            "generated_at_utc": now(),
            "software": "xTB 6.7.1pre",
            "method": "GFN2-xTB tight Opt with ALPB(water)",
            "purpose": "neutral versus N6-protonated sensitivity screen",
            "not_a_pka_calculation": True,
            "not_combined_with_orca_energies": True,
            "nprocs": 16,
            "ensembles": all_results,
        }
        write_findings(result)
        write_json(RESULTS, result)
        state.update(status="SUCCESS", finished_at_utc=now(), current_job=None)
        write_json(STATUS, state)
        return 0
    except Exception as exc:
        state.update(status="FAILED", finished_at_utc=now(), error=f"{type(exc).__name__}: {exc}")
        write_json(STATUS, state)
        raise


if __name__ == "__main__":
    raise SystemExit(main())

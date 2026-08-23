"""Analyze five retained para-benzoyl ensembles without starting ORCA.

The primary comparison uses only each ensemble's strict
``final_conformer_ensemble``.  Para substituent branches are excluded from the
common-scaffold RMSD, but remain present in every steric-access calculation.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

from molecule2fbx.comparison import (
    DEFAULT_TEMPERATURE_K,
    R_KJ_MOL_K,
    carbonyl_continuous_steric_access,
    cross_ensemble_rmsd_analysis,
    load_existing_conformers,
    parse_existing_orca_electronic,
    repair_ensemble_analysis_metadata,
    secondary_rmsd_analysis,
    structural_descriptors,
)
from molecule2fbx.ensemble import HARTREE_TO_KJ_MOL
from molecule2fbx.structures import rmsd_atom_subsets


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "para_substituent_series_analysis"
ENSEMBLES = (
    ("1Bz-LSD_RR", "H"),
    ("1pMeBz-LSD_RR", "Me"),
    ("1p-iPrBz-LSD_RR", "iPr"),
    ("1ptBuBz-LSD_RR", "tBu"),
    ("1SB-LSD_RR", "TMS"),
)
PROBE_KEY = "probe_radius_0.50_angstrom"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def nested(payload: Mapping[str, object], path: Sequence[str]):
    value: object = payload
    for key in path:
        if not isinstance(value, Mapping) or key not in value:
            return None
        value = value[key]
    return value


def electronic_weights(records: Sequence[Mapping[str, object]]) -> dict[str, float]:
    minimum = min(float(record["dft_energy_hartree"]) for record in records)
    factors = []
    for record in records:
        delta = (
            float(record["dft_energy_hartree"]) - minimum
        ) * HARTREE_TO_KJ_MOL
        factors.append(
            (
                str(record["conformer_id"]),
                math.exp(-delta / (R_KJ_MOL_K * DEFAULT_TEMPERATURE_K)),
            )
        )
    total = sum(value for _identifier, value in factors)
    return {identifier: value / total for identifier, value in factors}


def weighted_stats(
    records: Sequence[Mapping[str, object]],
    weights: Mapping[str, float],
    metric: str,
) -> dict[str, object]:
    values = [
        (float(record["metrics"][metric]), weights[str(record["conformer_id"])])
        for record in records
        if isinstance(record.get("metrics"), Mapping)
        and record["metrics"].get(metric) is not None
        and str(record["conformer_id"]) in weights
    ]
    coverage = sum(weight for _value, weight in values)
    if not values or coverage <= 0.0:
        return {
            "weighted_mean": None,
            "weighted_std": None,
            "minimum": None,
            "maximum": None,
            "weight_coverage": 0.0,
        }
    normalized = [(value, weight / coverage) for value, weight in values]
    mean = sum(value * weight for value, weight in normalized)
    variance = sum(weight * (value - mean) ** 2 for value, weight in normalized)
    return {
        "weighted_mean": mean,
        "weighted_std": math.sqrt(max(variance, 0.0)),
        "minimum": min(value for value, _weight in values),
        "maximum": max(value for value, _weight in values),
        "weight_coverage": coverage,
    }


def access_metrics(conformer) -> tuple[dict[str, object], dict[str, object]]:
    access = carbonyl_continuous_steric_access(
        conformer.model,
        azimuth_samples=720,
        probe_radii_angstrom=(0.0, 0.5, 1.0, 1.4),
    )
    structure = structural_descriptors(conformer.model)
    electronic = parse_existing_orca_electronic(conformer.opt_output_path)
    subsets = rmsd_atom_subsets(conformer.model)
    carbonyl_c = int(structure["atom_indices"]["benzoyl_carbonyl_c"])
    carbonyl_o = int(structure["atom_indices"]["benzoyl_carbonyl_o"])
    amide_n = int(structure["atom_indices"]["benzoyl_amide_n"])
    mayer = electronic.get("mayer_bond_orders")
    if not isinstance(mayer, Mapping):
        mayer = {}

    def electronic_charge(kind: str, atom_index: int):
        values = electronic.get(kind)
        return values.get(atom_index) if isinstance(values, Mapping) else None

    def mayer_order(first: int, second: int):
        return mayer.get(tuple(sorted((first, second))))
    substituent_distances = [
        math.dist(
            conformer.model.atoms[carbonyl_c].position,
            conformer.model.atoms[index].position,
        )
        for index in subsets.excluded_terminal_substituent
    ]
    paths = {
        "total_mean_clearance_probe_0_5_angstrom": (
            "scopes", "total", "probe_radius_sensitivity", PROBE_KEY,
            "mean_clearance_angstrom",
        ),
        "total_top_decile_clearance_probe_0_5_angstrom": (
            "scopes", "total", "probe_radius_sensitivity", PROBE_KEY,
            "top_decile_mean_clearance_angstrom",
        ),
        "total_positive_clearance_integral_probe_0_5_angstrom": (
            "scopes", "total", "probe_radius_sensitivity", PROBE_KEY,
            "positive_clearance_integral_angstrom",
        ),
        "nonlocal_mean_clearance_probe_0_5_angstrom": (
            "scopes", "nonlocal_environment", "probe_radius_sensitivity",
            PROBE_KEY, "mean_clearance_angstrom",
        ),
        "nonlocal_top_decile_clearance_probe_0_5_angstrom": (
            "scopes", "nonlocal_environment", "probe_radius_sensitivity",
            PROBE_KEY, "top_decile_mean_clearance_angstrom",
        ),
        "re_nonlocal_mean_clearance_probe_0_5_angstrom": (
            "scopes", "nonlocal_environment", "face_summary", "Re",
            "probe_radius_sensitivity", PROBE_KEY, "mean_clearance_angstrom",
        ),
        "si_nonlocal_mean_clearance_probe_0_5_angstrom": (
            "scopes", "nonlocal_environment", "face_summary", "Si",
            "probe_radius_sensitivity", PROBE_KEY, "mean_clearance_angstrom",
        ),
        "para_substituent_limiter_fraction": (
            "para_substituent_as_total_limiter_fraction",
        ),
        "direct_para_substituent_nonlocal_integral_delta_angstrom": (
            "direct_para_substituent_counterfactual", "probe_radius_sensitivity",
            PROBE_KEY, "nonlocal_environment",
            "positive_clearance_integral_with_minus_without_substituent_angstrom",
        ),
    }
    metrics = {name: nested(access, path) for name, path in paths.items()}
    for scope in ("total", "nonlocal_environment"):
        without_scope = (
            "total_without_para_substituent"
            if scope == "total"
            else "nonlocal_without_para_substituent"
        )
        with_value = nested(
            access,
            ("scopes", scope, "probe_radius_sensitivity", PROBE_KEY,
             "mean_clearance_angstrom"),
        )
        without_value = nested(
            access,
            ("scopes", without_scope, "probe_radius_sensitivity", PROBE_KEY,
             "mean_clearance_angstrom"),
        )
        metrics[
            f"direct_para_substituent_{scope}_mean_clearance_delta_angstrom"
        ] = (
            None
            if not subsets.excluded_terminal_substituent
            or with_value is None
            or without_value is None
            else float(with_value) - float(without_value)
        )
    re_value = metrics["re_nonlocal_mean_clearance_probe_0_5_angstrom"]
    si_value = metrics["si_nonlocal_mean_clearance_probe_0_5_angstrom"]
    metrics.update(
        {
            "re_minus_si_nonlocal_mean_clearance_angstrom": (
                None
                if re_value is None or si_value is None
                else float(re_value) - float(si_value)
            ),
            "minimum_para_substituent_distance_from_carbonyl_c_angstrom": (
                min(substituent_distances) if substituent_distances else None
            ),
            "carbonyl_c_o_length_angstrom": structure[
                "carbonyl_c_o_length_angstrom"
            ],
            "amide_c_n_length_angstrom": structure[
                "amide_c_n_length_angstrom"
            ],
            "mayer_bond_order_c_o": mayer_order(carbonyl_c, carbonyl_o),
            "mayer_bond_order_c_n": mayer_order(carbonyl_c, amide_n),
            "loewdin_charge_carbonyl_c": electronic_charge(
                "loewdin_charges", carbonyl_c
            ),
            "loewdin_charge_carbonyl_o": electronic_charge(
                "loewdin_charges", carbonyl_o
            ),
            "loewdin_charge_amide_n": electronic_charge(
                "loewdin_charges", amide_n
            ),
            "dipole_magnitude_debye": electronic.get(
                "dipole_magnitude_debye"
            ),
            "homo_lumo_gap_ev": electronic.get("homo_lumo_gap_ev"),
            "benzoyl_carbonyl_torsion_abs_mean_deg": structure[
                "benzoyl_carbonyl_torsion_abs_mean_deg"
            ],
        }
    )
    compact = {
        "method": access["method"],
        "parameters": access["parameters"],
        "face_assignment": access["face_assignment"],
        "para_substituent_heavy_atom_indices": access["atom_indices"][
            "para_substituent_heavy"
        ],
        "para_substituent_elements": [
            conformer.model.atoms[index].element
            for index in subsets.excluded_terminal_substituent
        ],
    }
    return metrics, compact


def analyze_ensemble(name: str, label: str, conformers, payload) -> dict[str, object]:
    rmsd = secondary_rmsd_analysis(conformers)
    by_index = {
        int(record["conformer_index"]): record for record in rmsd["records"]
    }
    records = []
    for conformer in conformers:
        index = int(conformer.entry["conformer_index"])
        metrics, access_metadata = access_metrics(conformer)
        records.append(
            {
                "conformer_id": str(conformer.entry["conformer_id"]),
                "conformer_index": index,
                "conformer_pool_index": conformer.entry.get(
                    "conformer_pool_index"
                ),
                "optimization_provenance": conformer.entry.get(
                    "optimization_provenance"
                ),
                "dft_energy_hartree": float(
                    conformer.entry["dft_energy_hartree"]
                ),
                "frequency_calculated": conformer.entry.get(
                    "frequency_calculated"
                ) is True,
                "imaginary_modes": conformer.entry.get("imaginary_modes"),
                "all_heavy_rmsd_cluster_id": by_index[index][
                    "all_heavy_rmsd_cluster_id"
                ],
                "common_scaffold_rmsd_cluster_id": by_index[index][
                    "common_scaffold_rmsd_cluster_id"
                ],
                "common_scaffold_rmsd_to_representative_angstrom": by_index[
                    index
                ]["common_scaffold_rmsd_to_representative_angstrom"],
                "reaction_center_rmsd_cluster_id": by_index[index][
                    "reaction_center_rmsd_cluster_id"
                ],
                "metrics": metrics,
                "steric_access_metadata": access_metadata,
            }
        )
    grouped = defaultdict(list)
    for record in records:
        grouped[str(record["common_scaffold_rmsd_cluster_id"])].append(record)
    representatives = [
        min(group, key=lambda item: float(item["dft_energy_hartree"]))
        for group in grouped.values()
    ]
    representatives.sort(key=lambda item: float(item["dft_energy_hartree"]))
    weights = electronic_weights(representatives)
    minimum = float(representatives[0]["dft_energy_hartree"])
    representative_ids = {str(item["conformer_id"]) for item in representatives}
    for record in records:
        identifier = str(record["conformer_id"])
        record["common_scaffold_representative"] = identifier in representative_ids
        record["electronic_weight_if_representative"] = weights.get(identifier)
        record["relative_electronic_energy_kj_mol_rebased"] = (
            float(record["dft_energy_hartree"]) - minimum
        ) * HARTREE_TO_KJ_MOL

    metric_names = tuple(records[0]["metrics"])
    dft = payload["dft"]
    search = payload["conformer_search"]
    return {
        "name": name,
        "substituent": label,
        "calculation_conditions": {
            key: dft.get(key)
            for key in (
                "software", "version", "functional", "basis", "charge",
                "multiplicity", "solvent_model", "dispersion_correction",
                "nprocs", "maxcore_mb_per_process",
            )
        },
        "conformer_search": {
            "pool_size": search.get("pool_size"),
            "random_seed": search.get("random_seed"),
            "forcefield_used": search.get("forcefield_used"),
            "energy_window_kj_mol": search.get("initial_screening", {}).get(
                "energy_window_kj_mol"
            ),
            "all_heavy_rmsd_threshold_angstrom": search.get(
                "initial_screening", {}
            ).get("rmsd_threshold_angstrom"),
        },
        "selection_protocol": {
            "dft_max_conformers": payload.get("dft_max_conformers"),
            "post_dft_all_heavy_rmsd_threshold_angstrom": dft.get(
                "post_optimization_screening", {}
            ).get("rmsd_threshold_angstrom"),
            "frequency_window_kj_mol": payload.get("frequency_window_kj"),
            "frequency_max": payload.get("frequency_max"),
            "stereochemistry_configuration": payload.get(
                "stereochemistry", {}
            ).get("configuration"),
        },
        "all_heavy_dft_candidate_count": len(records),
        "common_scaffold_cluster_count": len(representatives),
        "reaction_center_cluster_count": rmsd["reaction_center"][
            "cluster_count"
        ],
        "common_scaffold_atom_definition": rmsd["atom_subsets"],
        "representative_conformer_ids": [
            record["conformer_id"] for record in representatives
        ],
        "metric_summary": {
            metric: weighted_stats(representatives, weights, metric)
            for metric in metric_names
        },
        "records": records,
    }


def paired_to_bz(bz_conformers, other_conformers, bz, other) -> dict[str, object]:
    cross = cross_ensemble_rmsd_analysis(bz_conformers, other_conformers)
    bz_records = {record["conformer_id"]: record for record in bz["records"]}
    other_records = {
        record["conformer_id"]: record for record in other["records"]
    }
    pairs = []
    for mapping in cross["reciprocal_nearest_pairs"]:
        first = bz_records[str(mapping["first_conformer_id"])]
        second = other_records[str(mapping["second_conformer_id"])]
        deltas = {}
        for metric, first_value in first["metrics"].items():
            second_value = second["metrics"].get(metric)
            deltas[metric] = (
                None
                if first_value is None or second_value is None
                else float(second_value) - float(first_value)
            )
        pairs.append({**mapping, "other_minus_bz": deltas})
    metric_names = tuple(bz["records"][0]["metrics"])
    delta_summary = {}
    for metric in metric_names:
        values = [
            float(pair["other_minus_bz"][metric])
            for pair in pairs
            if pair["other_minus_bz"].get(metric) is not None
        ]
        delta_summary[metric] = {
            "pair_count": len(values),
            "mean_other_minus_bz": (
                sum(values) / len(values) if values else None
            ),
            "minimum_other_minus_bz": min(values) if values else None,
            "maximum_other_minus_bz": max(values) if values else None,
            "positive_pair_count": sum(value > 0.0 for value in values),
            "negative_pair_count": sum(value < 0.0 for value in values),
        }
    return {
        "atom_mapping": cross["atom_mapping"],
        "reciprocal_nearest_pair_count": len(pairs),
        "metric_delta_summary": delta_summary,
        "pairs": pairs,
    }


def number(value, digits: int = 4) -> str:
    return "—" if value is None else f"{float(value):.{digits}f}"


def write_report(payload: Mapping[str, object], path: Path) -> None:
    ensembles = payload["ensembles"]
    rows = []
    for name, label in ENSEMBLES:
        item = ensembles[name]
        summary = item["metric_summary"]
        paired_summary = payload["paired_to_1Bz_LSD_RR"].get(name, {})
        paired_delta = nested(
            paired_summary,
            (
                "metric_delta_summary",
                "nonlocal_mean_clearance_probe_0_5_angstrom",
                "mean_other_minus_bz",
            ),
        )
        total = summary["total_mean_clearance_probe_0_5_angstrom"]
        nonlocal_stats = summary["nonlocal_mean_clearance_probe_0_5_angstrom"]
        rows.append(
            "| "
            + " | ".join(
                (
                    label,
                    str(item["all_heavy_dft_candidate_count"]),
                    str(item["common_scaffold_cluster_count"]),
                    f"{number(total['weighted_mean'])} ± {number(total['weighted_std'])}",
                    f"{number(nonlocal_stats['weighted_mean'])} ± {number(nonlocal_stats['weighted_std'])}",
                    number(paired_delta),
                    number(summary["para_substituent_limiter_fraction"]["weighted_mean"]),
                    number(summary["direct_para_substituent_nonlocal_environment_mean_clearance_delta_angstrom"]["weighted_mean"]),
                    number(summary["minimum_para_substituent_distance_from_carbonyl_c_angstrom"]["weighted_mean"]),
                    number(summary["carbonyl_c_o_length_angstrom"]["weighted_mean"], 5),
                    number(summary["amide_c_n_length_angstrom"]["weighted_mean"], 5),
                    number(summary["benzoyl_carbonyl_torsion_abs_mean_deg"]["weighted_mean"], 3),
                )
            )
            + " |"
        )
    tbu = ensembles["1ptBuBz-LSD_RR"]["metric_summary"]
    tms = ensembles["1SB-LSD_RR"]["metric_summary"]
    tbu_distance = tbu[
        "minimum_para_substituent_distance_from_carbonyl_c_angstrom"
    ]["weighted_mean"]
    tms_distance = tms[
        "minimum_para_substituent_distance_from_carbonyl_c_angstrom"
    ]["weighted_mean"]
    tbu_direct = tbu[
        "direct_para_substituent_nonlocal_environment_mean_clearance_delta_angstrom"
    ]["weighted_mean"]
    tms_direct = tms[
        "direct_para_substituent_nonlocal_environment_mean_clearance_delta_angstrom"
    ]["weighted_mean"]
    tms_paired = payload["paired_to_1Bz_LSD_RR"]["1SB-LSD_RR"][
        "metric_delta_summary"
    ]
    alkyl_paired = {
        label: payload["paired_to_1Bz_LSD_RR"][name]["metric_delta_summary"]
        for name, label in ENSEMBLES[1:4]
    }
    electronic_rows = []
    for name, label in ENSEMBLES:
        summary = ensembles[name]["metric_summary"]
        paired_summary = payload["paired_to_1Bz_LSD_RR"].get(name, {})
        paired_mayer_co = nested(
            paired_summary,
            (
                "metric_delta_summary", "mayer_bond_order_c_o",
                "mean_other_minus_bz",
            ),
        )
        paired_mayer_cn = nested(
            paired_summary,
            (
                "metric_delta_summary", "mayer_bond_order_c_n",
                "mean_other_minus_bz",
            ),
        )
        electronic_rows.append(
            "| "
            + " | ".join(
                (
                    label,
                    number(summary["mayer_bond_order_c_o"]["weighted_mean"], 5),
                    number(summary["mayer_bond_order_c_n"]["weighted_mean"], 5),
                    number(paired_mayer_co, 5),
                    number(paired_mayer_cn, 5),
                    number(summary["loewdin_charge_carbonyl_c"]["weighted_mean"], 5),
                    number(summary["dipole_magnitude_debye"]["weighted_mean"], 4),
                    number(summary["homo_lumo_gap_ev"]["weighted_mean"], 4),
                )
            )
            + " |"
        )
    lines = [
        "# para置換benzoyl-LSD 5系列：共通骨格・立体アクセス解析",
        "",
        "新規ORCA計算は行っていない。各ensembleのstrict final_conformer_ensembleだけを一次集合として使用し、SB pool 109等のtargeted follow-upは混入させていない。",
        "",
        f"- 電子構造条件の一致: {payload['all_electronic_structure_conditions_identical']}",
        f"- sampling・選択protocolの一致: {payload['all_sampling_and_selection_protocols_identical']}",
        "",
        "| para置換基 | DFT構造 | 共通骨格cluster | total mean clearance / Å | nonlocal mean clearance / Å | 対応pair Δnonlocal vs H / Å | 置換基がlimiterとなる割合 | 置換基直接Δnonlocal / Å | carbonyl C–置換基最短距離 / Å | C=O / Å | N–C(O) / Å | benzoyl–carbonyl torsion / deg |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        *rows,
        "",
        "## 読み方",
        "",
        "- 共通骨格RMSDはbenzoyl para炭素より外側の重原子枝を除外する。置換基そのものは立体アクセス計算から除外していない。",
        "- 値は共通骨格clusterごとの最低電子エネルギー代表を、各分子内DFT電子エネルギーで重み付けした感度解析である。完全なGibbs/Boltzmann ensembleではない。",
        "- clearance差が配座間分布より小さい場合は、この静的モデルで検出できないと扱う。酵素反応速度や代謝経路へ直接変換しない。",
        "- 対称原子の置換によるRMSD最小化は行わず、立体化学を保持した決定論的グラフ対応を使用する。",
        "",
        "## 今回の予備判定",
        "",
        f"- TMSのcarbonyl C–置換基最短距離はtBuより {float(tms_distance)-float(tbu_distance):.4f} Å長い。Si–C結合によって置換基が外側へ位置するという幾何学的部分とは整合する。",
        f"- しかし置換基を仮想的に除いた場合との差は、TMS {float(tms_direct):.4f} Å、tBu {float(tbu_direct):.4f} Åでほぼ同程度だった。この指標ではTMSだけの特別な遮蔽緩和は検出できない。",
        f"- 対応配座でのTMS−H差はtotal mean clearanceが {float(tms_paired['total_mean_clearance_probe_0_5_angstrom']['mean_other_minus_bz']):.4f} Å、nonlocal mean clearanceが {float(tms_paired['nonlocal_mean_clearance_probe_0_5_angstrom']['mean_other_minus_bz']):.4f} Åで、いずれもアクセスをわずかに下げる方向だった。",
        "- したがって『TMSがcarbonylを大きく塞がない』とは整合するが、『長いSi–C結合がtBu等より明確に遮蔽を緩和する』という強い形の仮説1は、現データでは支持されない。",
        "",
        "## 結合幾何に見えるpara置換応答",
        "",
        f"- 対応配座でのTMS−H差はC=Oが {float(tms_paired['carbonyl_c_o_length_angstrom']['mean_other_minus_bz']):+.6f} Å、N–C(O)が {float(tms_paired['amide_c_n_length_angstrom']['mean_other_minus_bz']):+.6f} Åだった。符号はそれぞれ8/8 pair、7/8 pairで一致したが、絶対量は小さい。",
        "- C結合para置換基の対応pair平均は "
        + "; ".join(
            f"{label}: ΔC=O {float(summary['carbonyl_c_o_length_angstrom']['mean_other_minus_bz']):+.6f} Å, ΔN–C(O) {float(summary['amide_c_n_length_angstrom']['mean_other_minus_bz']):+.6f} Å"
            for label, summary in alkyl_paired.items()
        )
        + "。TMSの結合長応答はこれらより小さく、単純な置換基サイズ順でもない。",
        "- したがって結合長だけからTMS固有の強い電子効果、amide共鳴の増減、または加水分解速度の向きを決めることはできない。",
        "",
        "## 既存Opt出力の電子指標（追加single pointなし）",
        "",
        "| para置換基 | Mayer C=O | Mayer N–C(O) | 対応pair ΔMayer C=O vs H | 対応pair ΔMayer N–C(O) vs H | Loewdin carbonyl-C / e | dipole / D | HOMO–LUMO gap / eV |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        *electronic_rows,
        "",
        "- これらは全5系列で同じB3LYP/def2-SVP気相Opt出力から抽出した感度指標である。Mulliken/Loewdin電荷は分割法依存、Mayer差は小さく、canonical HOMO/LUMOは軌道同一性が保証されないため、単独では反応性の証拠にしない。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    loaded = {}
    analyses = {}
    for name, label in ENSEMBLES:
        ensemble_dir = ROOT / "outputs" / name
        repaired = repair_ensemble_analysis_metadata(ensemble_dir)
        payload, conformers = load_existing_conformers(ensemble_dir)
        if repaired["conformer_search"].get("forcefield_used") is None:
            raise ValueError(f"Could not recover forcefield_used for {name}")
        loaded[name] = (payload, conformers)
        analyses[name] = analyze_ensemble(name, label, conformers, payload)

    reference_name = ENSEMBLES[0][0]
    reference_conformers = loaded[reference_name][1]
    paired = {
        name: paired_to_bz(
            reference_conformers,
            loaded[name][1],
            analyses[reference_name],
            analyses[name],
        )
        for name, _label in ENSEMBLES[1:]
    }
    condition_signatures = {
        name: json.dumps(
            analysis["calculation_conditions"], sort_keys=True
        )
        for name, analysis in analyses.items()
    }
    protocol_signatures = {
        name: json.dumps(
            {
                "conformer_search": analysis["conformer_search"],
                "selection_protocol": analysis["selection_protocol"],
            },
            sort_keys=True,
        )
        for name, analysis in analyses.items()
    }
    payload = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "calculation_scope": (
            "read-only geometry analysis plus derived ensemble metadata refresh; "
            "no ORCA, Opt, Freq, or single-point jobs"
        ),
        "primary_analysis_set": (
            "strict final_conformer_ensemble for each molecule"
        ),
        "supplementary_results_excluded_from_primary_weights": {
            "1SB-LSD_RR": ["conf011_pool109 targeted follow-up"],
        },
        "all_electronic_structure_conditions_identical": (
            len(set(condition_signatures.values())) == 1
        ),
        "all_sampling_and_selection_protocols_identical": (
            len(set(protocol_signatures.values())) == 1
        ),
        "ensembles": analyses,
        "paired_to_1Bz_LSD_RR": paired,
        "limitations": [
            "Electronic-energy weighting is not a complete Gibbs population.",
            "Static gas-phase optimized structures do not model enzyme dynamics or reaction barriers.",
            "Equivalent symmetric atom permutations are not searched for RMSD minimization.",
            "Absolute electronic energies are not compared between different molecular formulas.",
        ],
    }
    atomic_json(OUTPUT / "analysis.json", payload)

    rows = []
    for name, _label in ENSEMBLES:
        for record in analyses[name]["records"]:
            rows.append(
                {
                    "molecule": name,
                    "conformer_id": record["conformer_id"],
                    "pool_index": record["conformer_pool_index"],
                    "common_scaffold_cluster": record[
                        "common_scaffold_rmsd_cluster_id"
                    ],
                    "common_scaffold_representative": record[
                        "common_scaffold_representative"
                    ],
                    "relative_electronic_energy_kj_mol": record[
                        "relative_electronic_energy_kj_mol_rebased"
                    ],
                    "electronic_weight_if_representative": record[
                        "electronic_weight_if_representative"
                    ],
                    **record["metrics"],
                }
            )
    with (OUTPUT / "conformer_metrics.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_report(payload, OUTPUT / "FINDINGS_JA.md")


if __name__ == "__main__":
    main()

"""Summarize paired electronic-property differences without making rate claims."""

from __future__ import annotations

import csv
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[1]
ROOT = (
    WORKSPACE
    / "outputs"
    / "Bz_vs_SB_preliminary_comparison"
    / "electronic_robustness"
)
SOURCE = ROOT / "results.json"
OUTPUT = ROOT / "paired_analysis.json"
CSV_OUTPUT = ROOT / "paired_deltas.csv"

METRICS = (
    "mbis_charge_carbonyl_c",
    "mbis_charge_carbonyl_o",
    "mbis_charge_amide_n",
    "chelpg_charge_carbonyl_c",
    "chelpg_charge_carbonyl_o",
    "chelpg_charge_amide_n",
    "loewdin_charge_carbonyl_c",
    "loewdin_charge_carbonyl_o",
    "loewdin_charge_amide_n",
    "mulliken_charge_carbonyl_c",
    "mulliken_charge_carbonyl_o",
    "mulliken_charge_amide_n",
    "mayer_bond_order_c_o",
    "mayer_bond_order_c_n",
    "homo_energy_ev",
    "lumo_energy_ev",
    "homo_lumo_gap_ev",
    "homo_loewdin_population_benzoyl_center",
    "lumo_loewdin_population_benzoyl_center",
    "dipole_magnitude_debye",
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def write_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def summarize(values: list[float]) -> dict:
    return {
        "pair_count": len(values),
        "mean_sb_minus_bz": sum(values) / len(values),
        "median_sb_minus_bz": statistics.median(values),
        "standard_deviation": statistics.stdev(values) if len(values) > 1 else 0.0,
        "minimum_sb_minus_bz": min(values),
        "maximum_sb_minus_bz": max(values),
        "positive_pair_count": sum(value > 0.0 for value in values),
        "negative_pair_count": sum(value < 0.0 for value in values),
        "same_nonzero_sign_for_all_pairs": (
            all(value > 0.0 for value in values)
            or all(value < 0.0 for value in values)
        ),
    }


def main() -> int:
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    if payload.get("status") != "SUCCESS":
        raise RuntimeError("Robustness calculations are not complete")
    by_key = {
        (item["condition"], item["molecule"], item["conformer_id"]): item
        for item in payload["results"]
    }
    rows = []
    condition_summaries = {}
    for condition in payload["conditions"]:
        condition_rows = []
        for pair_number, pair in enumerate(payload["pairs"], start=1):
            bz_key = (
                condition,
                "1Bz-LSD_RR",
                pair["first_conformer_id"],
            )
            sb_key = (
                condition,
                "1SB-LSD_RR",
                pair["second_conformer_id"],
            )
            if bz_key not in by_key or sb_key not in by_key:
                continue
            bz = by_key[bz_key]["properties"]
            sb = by_key[sb_key]["properties"]
            pair_record = {
                "condition": condition,
                "pair_id": f"pair{pair_number:03d}",
                "bz_conformer_id": pair["first_conformer_id"],
                "sb_conformer_id": pair["second_conformer_id"],
                "deltas": {},
            }
            for metric in METRICS:
                bz_value = bz.get(metric)
                sb_value = sb.get(metric)
                if bz_value is None or sb_value is None:
                    continue
                delta = float(sb_value) - float(bz_value)
                pair_record["deltas"][metric] = {
                    "bz": float(bz_value),
                    "sb": float(sb_value),
                    "sb_minus_bz": delta,
                }
                rows.append(
                    {
                        "condition": condition,
                        "pair_id": pair_record["pair_id"],
                        "bz_conformer_id": pair["first_conformer_id"],
                        "sb_conformer_id": pair["second_conformer_id"],
                        "metric": metric,
                        "bz": float(bz_value),
                        "sb": float(sb_value),
                        "sb_minus_bz": delta,
                    }
                )
            condition_rows.append(pair_record)
        condition_summaries[condition] = {
            metric: summarize(
                [
                    float(record["deltas"][metric]["sb_minus_bz"])
                    for record in condition_rows
                    if metric in record["deltas"]
                ]
            )
            for metric in METRICS
            if any(metric in record["deltas"] for record in condition_rows)
        }

    baseline = condition_summaries["baseline_b3lyp_def2_svp_gas"]
    aqueous = condition_summaries["b3lyp_d3bj_def2_tzvp_cpcm_water"]
    cross_condition = {}
    for metric in METRICS:
        if metric not in baseline or metric not in aqueous:
            continue
        baseline_by_pair = {
            str(row["pair_id"]): float(row["sb_minus_bz"])
            for row in rows
            if row["condition"] == "baseline_b3lyp_def2_svp_gas"
            and row["metric"] == metric
        }
        aqueous_by_pair = {
            str(row["pair_id"]): float(row["sb_minus_bz"])
            for row in rows
            if row["condition"] == "b3lyp_d3bj_def2_tzvp_cpcm_water"
            and row["metric"] == metric
        }
        common_pairs = sorted(set(baseline_by_pair) & set(aqueous_by_pair))
        baseline_overlap = [baseline_by_pair[pair_id] for pair_id in common_pairs]
        aqueous_overlap = [aqueous_by_pair[pair_id] for pair_id in common_pairs]
        first = sum(baseline_overlap) / len(baseline_overlap)
        second = sum(aqueous_overlap) / len(aqueous_overlap)
        cross_condition[metric] = {
            "common_pair_ids": common_pairs,
            "baseline_common_pair_mean_sb_minus_bz": first,
            "aqueous_tzvp_mean_sb_minus_bz": second,
            "same_nonzero_sign": first * second > 0.0,
            "same_sign_retained_pair_count": sum(
                baseline_value * aqueous_value > 0.0
                for baseline_value, aqueous_value in zip(
                    baseline_overlap, aqueous_overlap
                )
            ),
            "tested_pair_count": len(common_pairs),
            "aqueous_to_baseline_magnitude_ratio": (
                abs(second / first) if not math.isclose(first, 0.0) else None
            ),
        }

    primary = (
        "mbis_charge_carbonyl_c",
        "lumo_loewdin_population_benzoyl_center",
    )
    primary_direction_robust = all(
        baseline[metric]["same_nonzero_sign_for_all_pairs"]
        and aqueous[metric]["same_nonzero_sign_for_all_pairs"]
        and cross_condition[metric]["same_nonzero_sign"]
        for metric in primary
    )
    result = {
        "generated_at_utc": now(),
        "source": str(SOURCE),
        "delta_definition": "SB minus Bz on matched optimized geometries",
        "condition_summaries": condition_summaries,
        "cross_condition_robustness": cross_condition,
        "primary_direction_robust_within_tested_models": primary_direction_robust,
        "interpretation": {
            "if_true": (
                "The direction of the carbonyl-C MBIS and benzoyl-center LUMO "
                "shifts is robust across the tested conformers and two electronic "
                "conditions. This does not establish a hydrolysis rate effect."
            ),
            "scope_limit": (
                "Only B3LYP-based gas and aqueous-continuum single points on fixed "
                "geometries were tested; explicit solvent, enzyme, and activation "
                "free energies were not calculated."
            ),
        },
    }
    write_json(OUTPUT, result)
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

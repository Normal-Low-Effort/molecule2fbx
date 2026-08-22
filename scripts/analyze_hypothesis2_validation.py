"""Analyze the minimal Hypothesis-2 property calculations."""

from __future__ import annotations

import csv
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path


WORKSPACE = Path(__file__).resolve().parents[1]
ROOT = (
    WORKSPACE
    / "outputs"
    / "Bz_vs_SB_preliminary_comparison"
    / "hypothesis2_validation"
)
SOURCE = ROOT / "results.json"
OUTPUT = ROOT / "analysis.json"
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
        "count": len(values),
        "mean": sum(values) / len(values),
        "median": statistics.median(values),
        "standard_deviation": statistics.stdev(values) if len(values) > 1 else 0.0,
        "minimum": min(values),
        "maximum": max(values),
        "positive_count": sum(value > 0.0 for value in values),
        "negative_count": sum(value < 0.0 for value in values),
        "same_nonzero_sign": all(value > 0.0 for value in values)
        or all(value < 0.0 for value in values),
    }


def main() -> int:
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    if payload.get("status") != "SUCCESS":
        raise RuntimeError("Hypothesis-2 calculations are not complete")
    by_key = {
        (item["condition"], item["pair_id"], item["species"]): item
        for item in payload["results"]
    }
    rows = []
    summaries = {}
    pair_records = []
    for condition in payload["conditions"]:
        condition_deltas = {
            comparison: {metric: [] for metric in METRICS}
            for comparison in (
                "optimized_sb_minus_bz",
                "fixed_geometry_tms_minus_h",
                "geometry_residual_detms_at_sb_minus_bz",
            )
        }
        for pair_number, _pair in enumerate(payload["pairs"], start=1):
            pair_id = f"pair{pair_number:03d}"
            bz = by_key[(condition, pair_id, "bz_optimized")]["properties"]
            sb = by_key[(condition, pair_id, "sb_optimized")]["properties"]
            detms = by_key[(condition, pair_id, "sb_tms_to_h_fixed")]["properties"]
            comparisons = {
                "optimized_sb_minus_bz": (sb, bz),
                "fixed_geometry_tms_minus_h": (sb, detms),
                "geometry_residual_detms_at_sb_minus_bz": (detms, bz),
            }
            record = {"condition": condition, "pair_id": pair_id, "deltas": {}}
            for comparison, (first, second) in comparisons.items():
                record["deltas"][comparison] = {}
                for metric in METRICS:
                    first_value = first.get(metric)
                    second_value = second.get(metric)
                    if first_value is None or second_value is None:
                        continue
                    delta = float(first_value) - float(second_value)
                    condition_deltas[comparison][metric].append(delta)
                    record["deltas"][comparison][metric] = delta
                    rows.append(
                        {
                            "condition": condition,
                            "pair_id": pair_id,
                            "comparison": comparison,
                            "metric": metric,
                            "first": float(first_value),
                            "second": float(second_value),
                            "delta": delta,
                        }
                    )
            pair_records.append(record)
        summaries[condition] = {
            comparison: {
                metric: summarize(values)
                for metric, values in metric_values.items()
                if values
            }
            for comparison, metric_values in condition_deltas.items()
        }

    primary_metrics = (
        "mbis_charge_carbonyl_c",
        "lumo_loewdin_population_benzoyl_center",
    )
    cross_functional = {}
    for comparison in (
        "optimized_sb_minus_bz",
        "fixed_geometry_tms_minus_h",
    ):
        cross_functional[comparison] = {}
        for metric in METRICS:
            values = []
            signs = []
            for condition in payload["conditions"]:
                summary = summaries[condition][comparison].get(metric)
                if summary is None:
                    continue
                values.append(float(summary["mean"]))
                signs.append(
                    int(summary["positive_count"] - summary["negative_count"])
                )
            if values:
                cross_functional[comparison][metric] = {
                    "condition_means": values,
                    "same_nonzero_mean_sign": all(value > 0.0 for value in values)
                    or all(value < 0.0 for value in values),
                    "all_pairs_same_nonzero_sign_within_each_condition": all(
                        summaries[condition][comparison][metric]["same_nonzero_sign"]
                        for condition in payload["conditions"]
                        if metric in summaries[condition][comparison]
                    ),
                    "direction_votes": signs,
                }

    primary_consistent = all(
        cross_functional[comparison][metric]["same_nonzero_mean_sign"]
        and cross_functional[comparison][metric][
            "all_pairs_same_nonzero_sign_within_each_condition"
        ]
        for comparison in cross_functional
        for metric in primary_metrics
    )
    result = {
        "generated_at_utc": now(),
        "source": str(SOURCE),
        "delta_definitions": {
            "optimized_sb_minus_bz": "SB optimized geometry minus matched Bz optimized geometry",
            "fixed_geometry_tms_minus_h": (
                "SB minus the same SB common-atom geometry after SiMe3-to-H replacement"
            ),
            "geometry_residual_detms_at_sb_minus_bz": (
                "De-TMS counterfactual at SB geometry minus matched optimized Bz"
            ),
        },
        "condition_summaries": summaries,
        "pair_deltas": pair_records,
        "cross_functional_robustness": cross_functional,
        "primary_mbisc_and_lumo_directions_consistent": primary_consistent,
        "interpretation_limits": {
            "supported_if_consistent": (
                "A small ground-state electronic redistribution directly associated with TMS "
                "within the two tested functionals and two low-energy conformer pairs."
            ),
            "not_established": (
                "Hydrolysis kinetics, activation free energy, enzyme binding, metabolism, "
                "or behavior of other protonation states."
            ),
            "orbital_warning": (
                "A single canonical LUMO can change identity; MBIS carbonyl-C is the primary "
                "descriptor and the LUMO population is supporting evidence only."
            ),
        },
    }
    write_json(OUTPUT, result)
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(
        json.dumps(
            {
                "status": "SUCCESS",
                "primary_mbisc_and_lumo_directions_consistent": primary_consistent,
                "analysis": str(OUTPUT),
                "paired_deltas": str(CSV_OUTPUT),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

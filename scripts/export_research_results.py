"""Export a small, Git-friendly snapshot from the ignored calculation tree.

The files in ``outputs`` remain the canonical local calculation records.  This
script copies only explicitly approved summaries and derived tables, verifies
that every copied file is below the configured size limit, and writes SHA-256
provenance for the snapshot.  It never modifies or removes calculation data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


WORKSPACE = Path(__file__).resolve().parents[1]
DEFAULT_MAX_FILE_MIB = 5.0


# source path below outputs/, destination path below research_results/
SNAPSHOT_FILES: tuple[tuple[str, str], ...] = (
    ("1Bz-LSD_RR/ensemble.json", "ensembles/1Bz-LSD_RR/ensemble.json"),
    ("1Bz-LSD_RR/RUN_SUMMARY.md", "ensembles/1Bz-LSD_RR/RUN_SUMMARY.md"),
    ("1SB-LSD_RR/ensemble.json", "ensembles/1SB-LSD_RR/ensemble.json"),
    ("1SB-LSD_RR/RUN_SUMMARY.md", "ensembles/1SB-LSD_RR/RUN_SUMMARY.md"),
    (
        "Bz_vs_SB_preliminary_comparison/1Bz-LSD_RR_conformer_metrics.csv",
        "comparison/1Bz-LSD_RR_conformer_metrics.csv",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/1SB-LSD_RR_conformer_metrics.csv",
        "comparison/1SB-LSD_RR_conformer_metrics.csv",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/1SB-LSD_RR_regenerated_pool.json",
        "comparison/intermediate/1SB-LSD_RR_regenerated_pool.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/Bz_SB_common_scaffold_rmsd_matrix.csv",
        "comparison/Bz_SB_common_scaffold_rmsd_matrix.csv",
    ),
    ("Bz_vs_SB_preliminary_comparison/analysis.json", "comparison/analysis.json"),
    (
        "Bz_vs_SB_preliminary_comparison/PRELIMINARY_COMPARISON_JA.md",
        "comparison/PRELIMINARY_COMPARISON_JA.md",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/PRELIMINARY_COMPARISON.md",
        "comparison/PRELIMINARY_COMPARISON.md",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/paired_continuous_steric_access.json",
        "comparison/steric_access/paired_continuous_steric_access.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/paired_continuous_steric_trajectories.csv",
        "comparison/steric_access/paired_continuous_steric_trajectories.csv",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/hypotheses4_5_analysis.json",
        "comparison/hypotheses4_5/analysis.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/HYPOTHESES4_5_FINDINGS_JA.md",
        "comparison/hypotheses4_5/FINDINGS_JA.md",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/electronic_property_results.json",
        "comparison/electronic_properties/results.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/electronic_property_status.json",
        "comparison/electronic_properties/status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/electronic_robustness/results.json",
        "comparison/electronic_robustness/results.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/electronic_robustness/status.json",
        "comparison/electronic_robustness/status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/electronic_robustness/paired_analysis.json",
        "comparison/electronic_robustness/paired_analysis.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/electronic_robustness/paired_deltas.csv",
        "comparison/electronic_robustness/paired_deltas.csv",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/hypothesis2_validation/results.json",
        "comparison/hypothesis2_validation/results.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/hypothesis2_validation/status.json",
        "comparison/hypothesis2_validation/status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/hypothesis2_validation/analysis.json",
        "comparison/hypothesis2_validation/analysis.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/hypothesis2_validation/paired_deltas.csv",
        "comparison/hypothesis2_validation/paired_deltas.csv",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/hypothesis2_validation/eddb/"
        "pbe0_d3bj_def2_tzvp_cpcm_water/results.json",
        "comparison/hypothesis2_validation/eddb/results.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/hypothesis2_validation/eddb/"
        "pbe0_d3bj_def2_tzvp_cpcm_water/paired_deltas.csv",
        "comparison/hypothesis2_validation/eddb/paired_deltas.csv",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/protonation_screen/results.json",
        "comparison/protonation_screen/results.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/protonation_screen/status.json",
        "comparison/protonation_screen/status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/protonation_screen/FINDINGS_JA.md",
        "comparison/protonation_screen/FINDINGS_JA.md",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/explicit_water_screen/results.json",
        "comparison/explicit_water_screen/results.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/explicit_water_screen/status.json",
        "comparison/explicit_water_screen/status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/explicit_water_screen/FINDINGS_JA.md",
        "comparison/explicit_water_screen/FINDINGS_JA.md",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/sb_pool_gfn2_xtb_screen.json",
        "comparison/sb_pool_screen/gfn2_xtb_screen.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/sb_pool_gfn2_xtb_status.json",
        "comparison/sb_pool_screen/gfn2_xtb_status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/sb_pool_singlepoint_screen.json",
        "comparison/sb_pool_screen/singlepoint_screen.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/sb_pool_singlepoint_status.json",
        "comparison/sb_pool_screen/singlepoint_status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/sb_additional_dft_opt_result.json",
        "comparison/sb_additional_opt/result.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/sb_additional_dft_opt_status.json",
        "comparison/sb_additional_opt/status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/missing_freq_status.json",
        "provenance/missing_freq_status.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/operational_history.json",
        "provenance/operational_history.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/IMPLEMENTATION_NOTES.md",
        "provenance/IMPLEMENTATION_NOTES.md",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/ADDITIONAL_CALCULATIONS.md",
        "provenance/ADDITIONAL_CALCULATIONS.md",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/metabolic_claim_evidence_matrix.json",
        "context/metabolic_claim_evidence_matrix.json",
    ),
    (
        "Bz_vs_SB_preliminary_comparison/METABOLIC_INFERENCE_LIMITS_JA.md",
        "context/METABOLIC_INFERENCE_LIMITS_JA.md",
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def portable_snapshot_bytes(source: Path, workspace: Path) -> tuple[bytes, int]:
    """Return UTF-8 content with only the local workspace prefix redacted."""
    text = source.read_text(encoding="utf-8")
    replacements = (
        str(workspace).replace("\\", "\\\\"),
        str(workspace),
        workspace.as_posix(),
    )
    replacement_count = 0
    for local_prefix in replacements:
        count = text.count(local_prefix)
        if count:
            text = text.replace(local_prefix, "${WORKSPACE}")
            replacement_count += count
    return text.encode("utf-8"), replacement_count


def ensure_within(path: Path, parent: Path) -> None:
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError as exc:
        raise ValueError(f"Path escapes expected root: {path}") from exc


def iter_best_xyz(outputs: Path) -> Iterable[tuple[Path, Path, dict[str, str]]]:
    for name in ("1Bz-LSD_RR", "1SB-LSD_RR"):
        ensemble_path = outputs / name / "ensemble.json"
        payload = json.loads(ensemble_path.read_text(encoding="utf-8"))
        best_id = str(payload["best_conformer_id"])
        record = next(
            item
            for item in payload["conformers"]
            if str(item.get("conformer_id")) == best_id
        )
        calculation_directory = Path(record["calculation_directory"])
        ensure_within(calculation_directory, outputs)
        expected_xyz = calculation_directory / f"conformer_{int(record['conformer_index']) + 1:03d}.xyz"
        if not expected_xyz.is_file():
            candidates = sorted(
                path
                for path in calculation_directory.glob("*.xyz")
                if not path.name.endswith("_trj.xyz")
            )
            if len(candidates) != 1:
                raise FileNotFoundError(
                    f"Cannot identify optimized XYZ for {name} {best_id}: "
                    f"{calculation_directory}"
                )
            expected_xyz = candidates[0]
        destination = Path("current_best_structures") / name / expected_xyz.name
        metadata = {
            "molecule": name,
            "conformer_id": best_id,
            "calculation_source": str(record.get("source", "unknown")),
        }
        yield expected_xyz, destination, metadata


def export_snapshot(workspace: Path, max_file_mib: float, check_only: bool = False) -> list[dict[str, object]]:
    outputs = workspace / "outputs"
    destination_root = workspace / "research_results"
    max_bytes = int(max_file_mib * 1024 * 1024)
    entries: list[tuple[Path, Path, dict[str, str]]] = [
        (outputs / source, Path(destination), {})
        for source, destination in SNAPSHOT_FILES
    ]
    entries.extend(iter_best_xyz(outputs))

    manifest: list[dict[str, object]] = []
    for source, relative_destination, extra in entries:
        ensure_within(source, outputs)
        if not source.is_file():
            raise FileNotFoundError(f"Required snapshot source is missing: {source}")
        size = source.stat().st_size
        if size > max_bytes:
            raise ValueError(
                f"Snapshot source exceeds {max_file_mib:g} MiB limit: "
                f"{source} ({size} bytes)"
            )
        target = destination_root / relative_destination
        ensure_within(target, destination_root)
        snapshot_data, replacement_count = portable_snapshot_bytes(source, workspace)
        if not check_only:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(snapshot_data)
            shutil.copystat(source, target)
        entry: dict[str, object] = {
            "source": source.relative_to(workspace).as_posix(),
            "snapshot": target.relative_to(workspace).as_posix(),
            "source_bytes": size,
            "snapshot_bytes": len(snapshot_data),
            "source_sha256": sha256(source),
            "snapshot_sha256": sha256_bytes(snapshot_data),
            "workspace_path_replacements": replacement_count,
        }
        entry.update(extra)
        manifest.append(entry)

    if not check_only:
        latest_source_mtime = max(
            (workspace / str(item["source"])).stat().st_mtime for item in manifest
        )
        manifest_payload = {
            "schema_version": 1,
            "source_latest_modified_at_utc": datetime.fromtimestamp(
                latest_source_mtime, timezone.utc
            ).isoformat(),
            "source_policy": "outputs/ remains canonical and Git-ignored",
            "maximum_file_mib": max_file_mib,
            "file_count": len(manifest),
            "files": manifest,
        }
        manifest_path = destination_root / "SNAPSHOT_MANIFEST.json"
        manifest_path.write_text(
            json.dumps(manifest_payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate sources and size limits without copying files",
    )
    parser.add_argument(
        "--max-file-mib",
        type=float,
        default=DEFAULT_MAX_FILE_MIB,
        help=f"maximum copied file size (default: {DEFAULT_MAX_FILE_MIB:g} MiB)",
    )
    args = parser.parse_args()
    if args.max_file_mib <= 0:
        parser.error("--max-file-mib must be greater than zero")
    manifest = export_snapshot(WORKSPACE, args.max_file_mib, args.check)
    verb = "Validated" if args.check else "Exported"
    total = sum(int(item["source_bytes"]) for item in manifest)
    print(f"{verb} {len(manifest)} files ({total / 1024 / 1024:.2f} MiB).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Command-line interface for molecule2fbx."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .config import (
    DEFAULT_BASIS,
    DEFAULT_CHARGE,
    DEFAULT_FUNCTIONAL,
    DEFAULT_MULTIPLICITY,
    ConversionRequest,
    FrequencyOnlyRequest,
)
from .errors import Molecule2FBXError
from .frequency import run_frequency_only
from .pipeline import run_conversion
from .pubchem import validate_cid


def _write_failed_ensemble_artifacts(args: argparse.Namespace, exc: BaseException) -> None:
    """Leave a machine- and human-readable failure marker for unattended runs."""

    if not getattr(args, "ensemble", False):
        return
    output_dir = Path(args.output_dir or "output").expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    payload = {
        "schema_version": 1,
        "run_type": "conformer_ensemble_dft_screening",
        "calculation_status": "FAILED",
        "generated_at_utc": timestamp,
        "name": args.name,
        "error": {"type": type(exc).__name__, "message": str(exc)},
        "interpretation": {
            "result": "No successful ensemble result is claimed by this failure record."
        },
    }
    final_json = output_dir / "ensemble.json"
    if final_json.exists():
        final_json = output_dir / (
            "ensemble_failed_"
            + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            + ".json"
        )
    temporary = final_json.with_suffix(final_json.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(final_json)
    summary = output_dir / "RUN_SUMMARY.md"
    if summary.exists():
        summary = output_dir / (
            "RUN_SUMMARY_FAILED_"
            + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            + ".md"
        )
    summary.write_text(
        "# Ensemble Night Run\n\n"
        "## Calculation status\n\nFAILED\n\n"
        "## Error\n\n"
        f"- Type: `{type(exc).__name__}`\n"
        f"- Message: {exc}\n"
        f"- Timestamp: {timestamp}\n\n"
        "Completed ORCA directories were not deleted. Correct the reported cause and "
        "rerun the same command to resume compatible completed calculations.\n",
        encoding="utf-8",
    )


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer of 1 or greater") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be an integer of 1 or greater")
    return parsed


def _nonnegative_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer of zero or greater") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be an integer of zero or greater")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="molecule2fbx",
        description=(
            "Create provenance-aware molecular FBX models from PubChem CID or SMILES. "
            "Quantum calculations run only with --method dft/hf or explicit --ensemble."
        ),
    )
    parser.add_argument(
        "legacy_cid",
        nargs="?",
        help="PubChem CID (legacy syntax: molecule2fbx 441140)",
    )
    source = parser.add_argument_group("molecule input")
    source.add_argument("--cid", dest="cid_option", help="PubChem Compound ID")
    source.add_argument("--smiles", help="SMILES; explicit R/S and E/Z information is retained")
    source.add_argument("--name", help="Override the molecule name used for output files")
    source.add_argument(
        "--frequency-only",
        metavar="XYZ",
        help="Run ORCA Freq on an existing optimized XYZ without rerunning Opt",
    )
    source.add_argument(
        "--metadata",
        help="Existing molecule metadata JSON to update in --frequency-only mode",
    )

    calculation = parser.add_argument_group("structure and calculation")
    calculation.add_argument(
        "--ensemble",
        action="store_true",
        help=(
            "Run the staged conformer ensemble pipeline (defaults: pool 200, "
            "10 DFT candidates, 10 kJ/mol force-field window, top-3 Freq within 5 kJ/mol)"
        ),
    )
    calculation.add_argument(
        "--method",
        choices=("auto", "pubchem", "forcefield", "dft", "hf"),
        default="auto",
        help=(
            "auto: PubChem 3D then force field; pubchem: require PubChem 3D; "
            "forcefield: ETKDG+MMFF/UFF; dft/hf: ORCA geometry optimization (default: auto)"
        ),
    )
    calculation.add_argument(
        "--functional",
        help=f"ORCA DFT functional (default: {DEFAULT_FUNCTIONAL})",
    )
    calculation.add_argument(
        "--basis",
        help=f"ORCA basis-set keyword (default: {DEFAULT_BASIS})",
    )
    calculation.add_argument(
        "--charge",
        type=int,
        help=f"Total molecular charge (default: {DEFAULT_CHARGE})",
    )
    calculation.add_argument(
        "--multiplicity",
        type=int,
        help=f"Spin multiplicity 2S+1 (default: {DEFAULT_MULTIPLICITY})",
    )
    calculation.add_argument(
        "--conformers",
        type=_positive_int,
        default=None,
        help="Number of selected conformers to optimize (default: 1; --ensemble: 10)",
    )
    calculation.add_argument(
        "--conformer-pool",
        type=_positive_int,
        help=(
            "ETKDG/MMFF candidate pool before diversity selection; omitted uses "
            "--conformers for legacy behavior"
        ),
    )
    calculation.add_argument(
        "--conformer-rmsd-threshold",
        type=float,
        default=0.75,
        metavar="ANGSTROM",
        help="Minimum aligned heavy-atom RMSD for diverse pool selection (default: 0.75)",
    )
    calculation.add_argument(
        "--forcefield-energy-window-kj",
        type=float,
        metavar="KJ_MOL",
        help="Pre-DFT force-field energy window (default with --ensemble: 10)",
    )
    calculation.add_argument(
        "--dft-energy-window-kj",
        type=float,
        metavar="KJ_MOL",
        help="Optional post-DFT electronic-energy window before reclustering",
    )
    calculation.add_argument(
        "--dft-rmsd-threshold",
        type=float,
        metavar="ANGSTROM",
        help="Post-DFT heavy-atom RMSD threshold (default: initial RMSD threshold)",
    )
    calculation.add_argument(
        "--random-seed",
        type=_nonnegative_int,
        default=0xF00D,
        help="ETKDG random seed (default: 61453 / 0xF00D)",
    )
    calculation.add_argument(
        "--embedding-prune-rmsd",
        type=float,
        default=-1.0,
        metavar="ANGSTROM",
        help="RDKit embedding RMSD pruning; -1 disables it (default: -1)",
    )
    calculation.add_argument(
        "--strict-stereochemistry",
        action="store_true",
        help="Stop before 3D generation when SMILES contains unresolved stereochemistry",
    )
    calculation.add_argument(
        "--expected-stereocenters",
        type=_nonnegative_int,
        metavar="COUNT",
        help="Require exactly this many assigned tetrahedral stereocenters",
    )
    calculation.add_argument(
        "--stereochemistry-label",
        help="Human-readable configuration label stored in reports, e.g. (6aR,9R)",
    )
    calculation.add_argument(
        "--frequency",
        action="store_true",
        help="Run an expensive frequency calculation after quantum optimization",
    )
    calculation.add_argument(
        "--frequency-window-kj",
        type=float,
        metavar="KJ_MOL",
        help=(
            "With --frequency, run Freq only within this electronic-energy window "
            "above the lowest optimized conformer"
        ),
    )
    calculation.add_argument(
        "--frequency-max",
        type=_positive_int,
        metavar="COUNT",
        help="With --frequency, cap Freq jobs after energy ranking",
    )
    calculation.add_argument(
        "--frequency-include",
        type=_positive_int,
        action="append",
        default=[],
        metavar="CONFORMER",
        help="Always retain a 1-based DFT conformer for Freq (repeatable)",
    )
    calculation.add_argument(
        "--imaginary-threshold-cm1",
        type=float,
        default=-20.0,
        help="Imaginary-mode cutoff in cm^-1 (default: -20)",
    )
    calculation.add_argument(
        "--low-frequency-threshold-cm1",
        type=float,
        default=50.0,
        help="Low-frequency reporting cutoff in cm^-1 (default: 50)",
    )
    calculation.add_argument(
        "--save-all-conformers",
        action="store_true",
        help="Export every converged conformer in addition to the lowest-energy structure",
    )
    calculation.add_argument(
        "--allow-metals",
        action="store_true",
        help="Permit metal calculations only with explicit basis, charge, spin, and functional",
    )

    backend = parser.add_argument_group("external programs and resources")
    backend.add_argument(
        "--backend",
        choices=("orca",),
        default="orca",
        help="Quantum chemistry backend (default: orca)",
    )
    backend.add_argument(
        "--orca",
        help="Path to the separately installed ORCA executable",
    )
    backend.add_argument(
        "--blender",
        help="Path to Blender; otherwise BLENDER_EXECUTABLE or PATH is used",
    )
    detected_nprocs = os.cpu_count() or 1
    backend.add_argument(
        "--nprocs",
        type=_positive_int,
        default=None,
        help=(
            "ORCA logical processors; omitted uses os.cpu_count() "
            f"(currently detected: {detected_nprocs})"
        ),
    )
    backend.add_argument(
        "--maxcore",
        type=int,
        default=1000,
        help="Approximate ORCA memory per process in MB (default: 1000)",
    )
    backend.add_argument(
        "--max-opt-steps",
        type=int,
        default=100,
        help="Maximum geometry-optimization steps (default: 100)",
    )
    backend.add_argument(
        "--max-scf-iterations",
        type=int,
        default=200,
        help="Maximum SCF iterations per geometry step (default: 200)",
    )
    backend.add_argument(
        "--keep-calculation-files",
        action="store_true",
        help="Keep ORCA input, output, XYZ, GBW, and scratch artifacts under output",
    )
    backend.add_argument(
        "--reuse-calculations",
        metavar="DIR",
        help=(
            "Reuse compatible converged ORCA conformers in DIR and calculate only "
            "missing conformers; existing files are never overwritten"
        ),
    )

    output = parser.add_argument_group("output and timeouts")
    output.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory for outputs (default: ./output; for --frequency-only, "
            "the XYZ directory)"
        ),
    )
    output.add_argument("--api-timeout", type=float, default=30.0, help="PubChem timeout seconds")
    output.add_argument(
        "--blender-timeout", type=float, default=180.0, help="Blender timeout seconds"
    )
    output.add_argument(
        "--quantum-timeout",
        type=float,
        default=3600.0,
        help="Timeout per quantum conformer in seconds (default: 3600)",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def _request_from_args(args: argparse.Namespace) -> ConversionRequest:
    if args.metadata is not None:
        raise ValueError("--metadata requires --frequency-only")
    cid_values = [value for value in (args.legacy_cid, args.cid_option) if value is not None]
    if len(cid_values) > 1:
        raise ValueError("Use either positional CID or --cid, not both")
    if args.smiles is not None and cid_values:
        raise ValueError("Specify either CID or --smiles, not both")
    cid = validate_cid(cid_values[0]) if cid_values else None
    if cid is None and args.smiles is None:
        raise ValueError("A positional CID, --cid, or --smiles is required")
    method = "dft" if args.ensemble and args.method == "auto" else args.method
    conformers = args.conformers if args.conformers is not None else (10 if args.ensemble else 1)
    conformer_pool = args.conformer_pool
    if conformer_pool is None and args.ensemble:
        conformer_pool = 200
    frequency = args.frequency or args.ensemble
    frequency_window_kj = args.frequency_window_kj
    if frequency_window_kj is None and args.ensemble:
        frequency_window_kj = 5.0
    frequency_max = args.frequency_max
    if frequency_max is None and args.ensemble:
        frequency_max = 3
    forcefield_window = args.forcefield_energy_window_kj
    if forcefield_window is None and args.ensemble:
        forcefield_window = 10.0
    return ConversionRequest(
        cid=cid,
        smiles=args.smiles,
        name=args.name,
        method=method,
        output_dir=Path(args.output_dir or "output"),
        blender_executable=args.blender,
        quantum_backend=args.backend,
        quantum_executable=args.orca,
        functional=args.functional,
        basis=args.basis,
        charge=args.charge,
        multiplicity=args.multiplicity,
        conformers=conformers,
        conformer_pool=conformer_pool,
        conformer_rmsd_threshold=args.conformer_rmsd_threshold,
        frequency=frequency,
        frequency_window_kj=frequency_window_kj,
        frequency_max=frequency_max,
        save_all_conformers=args.save_all_conformers or args.ensemble,
        keep_calculation_files=args.keep_calculation_files or args.ensemble,
        reuse_calculations=(
            Path(args.reuse_calculations) if args.reuse_calculations is not None else None
        ),
        allow_metals=args.allow_metals,
        api_timeout=args.api_timeout,
        blender_timeout=args.blender_timeout,
        quantum_timeout=args.quantum_timeout,
        max_opt_steps=args.max_opt_steps,
        max_scf_iterations=args.max_scf_iterations,
        nprocs=args.nprocs if args.nprocs is not None else (os.cpu_count() or 1),
        maxcore_mb=args.maxcore,
        ensemble=args.ensemble,
        strict_stereochemistry=args.strict_stereochemistry or args.ensemble,
        random_seed=args.random_seed,
        embedding_prune_rmsd=args.embedding_prune_rmsd,
        forcefield_energy_window_kj=forcefield_window,
        dft_energy_window_kj=args.dft_energy_window_kj,
        dft_rmsd_threshold=args.dft_rmsd_threshold,
        frequency_include=tuple(index - 1 for index in args.frequency_include),
        imaginary_threshold_cm1=args.imaginary_threshold_cm1,
        low_frequency_threshold_cm1=args.low_frequency_threshold_cm1,
        expected_stereocenters=args.expected_stereocenters,
        stereochemistry_label=args.stereochemistry_label,
    )


def _frequency_request_from_args(args: argparse.Namespace) -> FrequencyOnlyRequest:
    cid_values = [value for value in (args.legacy_cid, args.cid_option) if value is not None]
    if cid_values or args.smiles is not None:
        raise ValueError("--frequency-only cannot be combined with CID or --smiles")
    if args.frequency:
        raise ValueError("Use either --frequency or --frequency-only, not both")
    if args.reuse_calculations is not None:
        raise ValueError("--reuse-calculations cannot be combined with --frequency-only")
    if args.conformer_pool is not None:
        raise ValueError("--conformer-pool cannot be combined with --frequency-only")
    if args.frequency_window_kj is not None or args.frequency_max is not None:
        raise ValueError(
            "Selective frequency options cannot be combined with --frequency-only"
        )
    if args.frequency_include:
        raise ValueError("--frequency-include cannot be combined with --frequency-only")
    if args.ensemble:
        raise ValueError("--ensemble cannot be combined with --frequency-only")
    method = None if args.method == "auto" else args.method
    return FrequencyOnlyRequest(
        xyz_path=Path(args.frequency_only),
        output_dir=Path(args.output_dir) if args.output_dir else None,
        metadata_path=Path(args.metadata) if args.metadata else None,
        method=method,
        functional=args.functional,
        basis=args.basis,
        charge=args.charge,
        multiplicity=args.multiplicity,
        quantum_executable=args.orca,
        quantum_timeout=args.quantum_timeout,
        max_scf_iterations=args.max_scf_iterations,
        nprocs=args.nprocs if args.nprocs is not None else (os.cpu_count() or 1),
        maxcore_mb=args.maxcore,
        imaginary_threshold_cm1=args.imaginary_threshold_cm1,
        low_frequency_threshold_cm1=args.low_frequency_threshold_cm1,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.frequency_only:
            request = _frequency_request_from_args(args)
            outcome = run_frequency_only(
                request, log=lambda message: print(message, file=sys.stderr)
            )
            print(f"ORCA input: {outcome.result.input_path}")
            print(f"ORCA output: {outcome.result.output_path}")
            print(f"Frequency metadata: {outcome.result_metadata_path}")
            if outcome.updated_metadata_path:
                print(f"Updated metadata: {outcome.updated_metadata_path}")
            return 0
        request = _request_from_args(args)
        outcome = run_conversion(request, log=lambda message: print(message, file=sys.stderr))
        for artifact in outcome.artifacts:
            print(f"Created: {artifact.fbx_path}")
            print(f"Metadata: {artifact.metadata_path}")
        if outcome.ensemble_report_path is not None:
            print(f"Ensemble report: {outcome.ensemble_report_path}")
        if outcome.run_summary_path is not None:
            print(f"Run summary: {outcome.run_summary_path}")
        return 0
    except Molecule2FBXError as exc:
        _write_failed_ensemble_artifacts(args, exc)
        print(str(exc), file=sys.stderr)
        return 1
    except ValueError as exc:
        _write_failed_ensemble_artifacts(args, exc)
        print(str(exc), file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover - defensive unattended-run guard
        _write_failed_ensemble_artifacts(args, exc)
        print(f"Unexpected failure: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Cancelled", file=sys.stderr)
        return 130

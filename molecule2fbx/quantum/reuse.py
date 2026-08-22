"""Safe loading of completed ORCA conformer calculations."""

from __future__ import annotations

import json
import re
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from ..errors import ConfigurationError, QuantumCalculationError
from ..frequency import ExistingOrcaSettings, find_related_metadata, read_orca_input_settings
from ..model import MoleculeModel
from ..structures import ConformerCandidate, geometry_warnings
from .base import FrequencyResult, QuantumResult, QuantumSettings
from .orca import parse_orca_output, read_xyz_coordinates


Logger = Callable[[str], None]
_OPT_TOKEN_RE = re.compile(r"(?:^|\s)Opt(?:\s|$)", re.IGNORECASE)
_FREQ_TOKEN_RE = re.compile(r"(?:^|\s)Freq(?:\s|$)", re.IGNORECASE)


def _legacy_reuse_template(model: MoleculeModel, index: int) -> MoleculeModel:
    """Drop force-field provenance that cannot be mapped to a legacy Opt job."""

    metadata = dict(model.metadata)
    for key in (
        "forcefield",
        "forcefield_energy_kcal_mol",
        "forcefield_converged",
        "forcefield_optimization_status",
        "forcefield_optimization_iterations",
        "forcefield_iteration_count_available",
        "forcefield_max_iterations",
        "conformer_pool_index",
        "conformer_pool_size",
        "initial_rmsd_cluster_id",
        "conformer_selection_rmsd_threshold_angstrom",
        "conformer_selection_threshold_relaxed",
    ):
        metadata.pop(key, None)
    metadata.update(
        {
            "conformer_index": index,
            "conformer_pool_index": None,
            "conformer_pool_index_provenance": "unavailable_for_legacy_reuse",
            "initial_structure_origin": "reused_existing_calculation",
            "forcefield_provenance": "unavailable_for_legacy_reuse",
        }
    )
    return MoleculeModel(model.cid, model.name, model.atoms, model.bonds, metadata)


def _simple_input_line(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    return next(
        (line.strip()[1:].strip() for line in text.splitlines() if line.lstrip().startswith("!")),
        "",
    )


def _same_token(actual: Optional[str], expected: Optional[str]) -> bool:
    return actual is not None and expected is not None and actual.casefold() == expected.casefold()


def _validate_settings(
    recorded: ExistingOrcaSettings,
    requested: QuantumSettings,
    *,
    input_path: Path,
) -> None:
    mismatches = []
    if recorded.method != requested.method:
        mismatches.append(f"method={recorded.method!r} (requested {requested.method!r})")
    if requested.method == "dft" and not _same_token(recorded.functional, requested.functional):
        mismatches.append(
            f"functional={recorded.functional!r} (requested {requested.functional!r})"
        )
    if not _same_token(recorded.basis, requested.basis):
        mismatches.append(f"basis={recorded.basis!r} (requested {requested.basis!r})")
    if recorded.charge != requested.charge:
        mismatches.append(f"charge={recorded.charge!r} (requested {requested.charge!r})")
    if recorded.multiplicity != requested.multiplicity:
        mismatches.append(
            f"multiplicity={recorded.multiplicity!r} (requested {requested.multiplicity!r})"
        )
    if mismatches:
        raise ConfigurationError(
            f"Cannot reuse {input_path}: ORCA settings differ: " + "; ".join(mismatches)
        )


def _read_metadata(path: Optional[Path]) -> Optional[Dict[str, object]]:
    if path is None:
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"Could not read reuse metadata JSON: {path}") from exc
    if not isinstance(data, dict):
        raise ConfigurationError(f"Reuse metadata JSON must contain an object: {path}")
    return data


def _validate_identity(
    candidate: ConformerCandidate,
    metadata: Optional[Dict[str, object]],
    *,
    metadata_path: Optional[Path],
) -> None:
    if metadata is None:
        return
    recorded_cid = metadata.get("cid")
    if candidate.model.cid is not None and recorded_cid != candidate.model.cid:
        raise ConfigurationError(
            f"Cannot reuse calculation: CID in {metadata_path} does not match the request"
        )
    atoms = metadata.get("atoms")
    if isinstance(atoms, list):
        recorded_elements = [
            atom.get("element") for atom in atoms if isinstance(atom, dict)
        ]
        expected_elements = [atom.element for atom in candidate.model.atoms]
        if recorded_elements != expected_elements:
            raise ConfigurationError(
                f"Cannot reuse calculation: atom order in {metadata_path} does not match"
            )
    nested = metadata.get("metadata")
    if isinstance(nested, dict):
        current_smiles = candidate.model.metadata.get("canonical_isomeric_smiles")
        recorded_smiles = nested.get("canonical_isomeric_smiles")
        if (
            isinstance(current_smiles, str)
            and isinstance(recorded_smiles, str)
            and current_smiles != recorded_smiles
        ):
            raise ConfigurationError(
                f"Cannot reuse calculation: canonical isomeric SMILES in {metadata_path} differs"
            )


def _frequency_result_from_files(
    input_path: Path,
    output_path: Path,
    settings: QuantumSettings,
) -> Optional[FrequencyResult]:
    if not input_path.exists() and not output_path.exists():
        return None
    if not input_path.is_file() or not output_path.is_file():
        raise ConfigurationError(
            f"Existing frequency calculation is incomplete: {input_path.parent}"
        )
    simple_line = _simple_input_line(input_path)
    if not _FREQ_TOKEN_RE.search(simple_line) or _OPT_TOKEN_RE.search(simple_line):
        raise ConfigurationError(f"Existing frequency input is not Freq-only: {input_path}")
    _validate_settings(read_orca_input_settings(input_path), settings, input_path=input_path)
    parsed = parse_orca_output(
        output_path.read_text(encoding="utf-8", errors="replace"),
        frequency_requested=True,
        require_geometry_convergence=False,
        imaginary_threshold_cm1=settings.imaginary_threshold_cm1,
    )
    if not parsed.frequencies_cm1:
        raise QuantumCalculationError(
            f"Existing frequency output contains no parseable frequencies: {output_path}"
        )
    return FrequencyResult(
        backend="ORCA",
        backend_version=parsed.version,
        method=settings.method,
        functional=settings.functional if settings.method == "dft" else None,
        basis=settings.basis,
        charge=settings.charge,
        multiplicity=settings.multiplicity,
        final_energy_hartree=parsed.energy_hartree,
        frequencies_cm1=parsed.frequencies_cm1,
        imaginary_frequencies_cm1=parsed.imaginary_frequencies_cm1,
        warnings=parsed.warnings,
        input_path=input_path,
        output_path=output_path,
        thermochemistry=parsed.thermochemistry,
    )


def _frequency_for_reused_model(
    model: MoleculeModel,
    settings: QuantumSettings,
    backend,
    calculation_root: Path,
    conformer_dir: Path,
    stem: str,
    log: Logger,
) -> FrequencyResult:
    direct = _frequency_result_from_files(
        conformer_dir / f"{stem}_freq.inp",
        conformer_dir / f"{stem}_freq.out",
        settings,
    )
    if direct is not None:
        log(f"Reusing completed frequency calculation for {stem}")
        return direct

    frequency_dir = calculation_root / "frequency_additions" / stem
    retained = _frequency_result_from_files(
        frequency_dir / f"{stem}_freq.inp",
        frequency_dir / f"{stem}_freq.out",
        settings,
    )
    if retained is not None:
        log(f"Reusing completed frequency calculation for {stem}")
        return retained
    if frequency_dir.exists() and any(frequency_dir.iterdir()):
        raise ConfigurationError(
            f"Refusing to overwrite incomplete frequency files: {frequency_dir}"
        )
    log(f"Running Freq only for optimized structure {stem}; Opt will not be rerun")
    return backend.frequency(model, settings, frequency_dir, f"{stem}_freq")


def _existing_frequency_for_reused_model(
    settings: QuantumSettings,
    calculation_root: Path,
    conformer_dir: Path,
    stem: str,
    log: Logger,
) -> Optional[FrequencyResult]:
    """Load a completed Freq result without ever launching or overwriting a job."""

    locations = (
        (
            conformer_dir / f"{stem}_freq.inp",
            conformer_dir / f"{stem}_freq.out",
        ),
        (
            calculation_root / "frequency_additions" / stem / f"{stem}_freq.inp",
            calculation_root / "frequency_additions" / stem / f"{stem}_freq.out",
        ),
    )
    for input_path, output_path in locations:
        try:
            result = _frequency_result_from_files(input_path, output_path, settings)
        except (ConfigurationError, QuantumCalculationError) as exc:
            log(
                f"Warning: retained frequency files for {stem} are not reusable: {exc}"
            )
            continue
        if result is not None:
            log(f"Reusing completed frequency calculation for {stem}")
            return result
    return None


def ensure_quantum_result_frequency(
    result: QuantumResult,
    settings: QuantumSettings,
    backend,
    calculation_root: Path,
    log: Logger,
) -> QuantumResult:
    """Attach a retained or new Freq-only result without rerunning optimization."""

    if result.frequencies_cm1:
        model = result.model.with_metadata(
            frequency_requested=True,
            frequencies_cm1=list(result.frequencies_cm1),
            imaginary_frequencies_cm1=list(result.imaginary_frequencies_cm1),
            frequency_provenance=result.model.metadata.get(
                "frequency_provenance", "reused_or_preexisting_unspecified"
            ),
        )
        return replace(result, model=model, frequency_requested=True)

    index = result.conformer_index
    stem = f"conformer_{index + 1:03d}"
    conformer_dir = result.calculation_directory or calculation_root / stem
    expected_existing_outputs = {
        (Path(conformer_dir) / f"{stem}_freq.out").resolve(),
        (
            Path(calculation_root)
            / "frequency_additions"
            / stem
            / f"{stem}_freq.out"
        ).resolve(),
    }
    existing_before = {
        path for path in expected_existing_outputs if path.is_file()
    }
    frequency = _frequency_for_reused_model(
        result.model,
        settings,
        backend,
        calculation_root,
        Path(conformer_dir),
        stem,
        log,
    )
    warnings = tuple(dict.fromkeys((*result.warnings, *frequency.warnings)))
    metadata_updates: Dict[str, object] = {
        "frequency_requested": True,
        "frequencies_cm1": list(frequency.frequencies_cm1),
        "imaginary_frequencies_cm1": list(frequency.imaginary_frequencies_cm1),
        "frequency_output": str(frequency.output_path.resolve()),
        "frequency_provenance": (
            "reused_existing"
            if frequency.output_path.resolve() in existing_before
            else "computed_this_run"
        ),
        "thermochemistry": (
            frequency.thermochemistry.to_metadata()
            if frequency.thermochemistry is not None
            else None
        ),
    }
    if warnings:
        metadata_updates["calculation_warnings"] = list(warnings)
    model = result.model.with_metadata(**metadata_updates)
    return replace(
        result,
        model=model,
        frequency_requested=True,
        frequencies_cm1=frequency.frequencies_cm1,
        imaginary_frequencies_cm1=frequency.imaginary_frequencies_cm1,
        warnings=warnings,
        thermochemistry=frequency.thermochemistry,
    )


def load_reusable_orca_result(
    candidate: ConformerCandidate,
    settings: QuantumSettings,
    backend,
    calculation_root: Path,
    log: Logger,
) -> Optional[QuantumResult]:
    """Return a verified completed conformer, or None when its directory is absent."""

    index = candidate.conformer_index
    stem = f"conformer_{index + 1:03d}"
    conformer_dir = calculation_root / stem
    if not conformer_dir.exists():
        return None
    if not conformer_dir.is_dir():
        raise ConfigurationError(f"Conformer calculation path is not a directory: {conformer_dir}")

    input_path = conformer_dir / f"{stem}.inp"
    output_path = conformer_dir / f"{stem}.out"
    xyz_path = conformer_dir / f"{stem}.xyz"
    missing = [path.name for path in (input_path, output_path, xyz_path) if not path.is_file()]
    if missing:
        raise ConfigurationError(
            f"Refusing to overwrite incomplete existing calculation {conformer_dir}; "
            f"missing: {', '.join(missing)}"
        )
    if not _OPT_TOKEN_RE.search(_simple_input_line(input_path)):
        raise ConfigurationError(f"Existing ORCA input is not an Opt job: {input_path}")
    _validate_settings(read_orca_input_settings(input_path), settings, input_path=input_path)

    metadata_path = find_related_metadata(xyz_path)
    metadata = _read_metadata(metadata_path)
    _validate_identity(candidate, metadata, metadata_path=metadata_path)
    if metadata_path is None:
        log(
            f"Warning: no associated metadata JSON found for {stem}; "
            "reuse is validated from ORCA files and atom order"
        )

    parsed = parse_orca_output(
        output_path.read_text(encoding="utf-8", errors="replace"),
        frequency_requested=False,
        require_geometry_convergence=True,
        imaginary_threshold_cm1=settings.imaginary_threshold_cm1,
    )
    coordinates = read_xyz_coordinates(
        xyz_path, [atom.element for atom in candidate.model.atoms]
    )
    optimized = candidate.model.with_coordinates(
        coordinates,
        structure_origin="computed_quantum",
        structure_source="Reused ORCA geometry optimization",
        structure_claim=(
            "Computed stationary structure; not an experimentally determined structure"
        ),
        quantum_backend="ORCA",
        quantum_backend_version=parsed.version,
        quantum_method=settings.method,
        functional=settings.functional if settings.method == "dft" else None,
        basis=settings.basis,
        charge=settings.charge,
        multiplicity=settings.multiplicity,
        final_energy_hartree=parsed.energy_hartree,
        scf_converged=True,
        geometry_converged=True,
        conformer_index=index,
        reused_quantum_result=True,
        reused_calculation_directory=str(conformer_dir.resolve()),
    )

    frequencies = parsed.frequencies_cm1
    imaginary = parsed.imaginary_frequencies_cm1
    warnings = list(parsed.warnings)
    frequency_source = output_path if frequencies else None
    thermochemistry = parsed.thermochemistry
    retained_frequency = None
    if not frequencies:
        retained_frequency = _existing_frequency_for_reused_model(
            settings,
            calculation_root,
            conformer_dir,
            stem,
            log,
        )
    if retained_frequency is not None:
        frequencies = retained_frequency.frequencies_cm1
        imaginary = retained_frequency.imaginary_frequencies_cm1
        warnings.extend(retained_frequency.warnings)
        frequency_source = retained_frequency.output_path
        thermochemistry = retained_frequency.thermochemistry
    elif settings.frequency and not frequencies:
        frequency_result = _frequency_for_reused_model(
            optimized,
            settings,
            backend,
            calculation_root,
            conformer_dir,
            stem,
            log,
        )
        frequencies = frequency_result.frequencies_cm1
        imaginary = frequency_result.imaginary_frequencies_cm1
        warnings.extend(frequency_result.warnings)
        frequency_source = frequency_result.output_path
        thermochemistry = frequency_result.thermochemistry

    warnings.extend(geometry_warnings(optimized))
    metadata_updates: Dict[str, object] = {
        "frequency_requested": bool(frequencies) or settings.frequency,
        "frequencies_cm1": list(frequencies),
        "imaginary_frequencies_cm1": list(imaginary),
        "thermochemistry": (
            thermochemistry.to_metadata() if thermochemistry is not None else None
        ),
    }
    if frequency_source is not None:
        metadata_updates["frequency_output"] = str(frequency_source.resolve())
        metadata_updates["frequency_provenance"] = "reused_existing"
    if warnings:
        metadata_updates["calculation_warnings"] = warnings
    optimized = optimized.with_metadata(**metadata_updates)

    return QuantumResult(
        model=optimized,
        backend="ORCA",
        backend_version=parsed.version,
        method=settings.method,
        functional=settings.functional if settings.method == "dft" else None,
        basis=settings.basis,
        charge=settings.charge,
        multiplicity=settings.multiplicity,
        final_energy_hartree=parsed.energy_hartree,
        geometry_converged=True,
        scf_converged=True,
        conformer_index=index,
        frequency_requested=bool(frequencies) or settings.frequency,
        frequencies_cm1=frequencies,
        imaginary_frequencies_cm1=imaginary,
        warnings=tuple(warnings),
        calculation_directory=conformer_dir.resolve(),
        thermochemistry=thermochemistry,
    )


def discover_reusable_orca_results(
    template: ConformerCandidate,
    settings: QuantumSettings,
    backend,
    calculation_root: Path,
    log: Logger,
    *,
    strict: bool,
) -> Tuple[List[QuantumResult], List[Dict[str, object]]]:
    """Discover numbered completed Opt jobs without writing to their directory."""

    calculation_root = Path(calculation_root).expanduser().resolve()
    if not calculation_root.exists():
        return [], []
    if not calculation_root.is_dir():
        raise ConfigurationError(
            f"Calculation reuse path is not a directory: {calculation_root}"
        )
    directories = []
    for path in calculation_root.iterdir():
        match = re.fullmatch(r"conformer_(\d{3,})", path.name)
        if path.is_dir() and match:
            directories.append((int(match.group(1)) - 1, path))
    results: List[QuantumResult] = []
    incomplete: List[Dict[str, object]] = []
    for index, path in sorted(directories):
        placeholder = ConformerCandidate(
            _legacy_reuse_template(template.model, index),
            index,
            None,
            "reused_existing_calculation",
            True,
        )
        try:
            loaded = load_reusable_orca_result(
                placeholder,
                settings,
                backend,
                calculation_root,
                log,
            )
        except (ConfigurationError, QuantumCalculationError) as exc:
            record = {
                "conformer_index": index,
                "directory": str(path),
                "error_type": type(exc).__name__,
                "message": str(exc),
            }
            if strict:
                raise ConfigurationError(
                    f"Existing calculation directory is not safely reusable: {path}: {exc}"
                ) from exc
            incomplete.append(record)
            log(f"Warning: incomplete resumable calculation {path.name}: {exc}")
            continue
        if loaded is not None:
            results.append(loaded)
    return results, incomplete


def archive_incomplete_calculation(directory: Path, calculation_root: Path) -> Path:
    """Move an incomplete local attempt aside before retrying; never delete it."""

    root = Path(calculation_root).resolve()
    source = Path(directory).resolve()
    if source.parent != root or not re.fullmatch(r"conformer_\d{3,}", source.name):
        raise ConfigurationError(
            f"Refusing to archive calculation outside the run directory: {source}"
        )
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive_root = root / "incomplete_attempts"
    archive_root.mkdir(parents=True, exist_ok=True)
    destination = archive_root / f"{source.name}_{timestamp}"
    suffix = 1
    while destination.exists():
        destination = archive_root / f"{source.name}_{timestamp}_{suffix:02d}"
        suffix += 1
    source.replace(destination)
    return destination

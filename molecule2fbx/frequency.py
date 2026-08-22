"""Frequency-only workflow for an existing ORCA-optimized XYZ geometry."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple, Type

from .config import DEFAULT_BASIS, DEFAULT_FUNCTIONAL, FrequencyOnlyRequest
from .electronic import validate_electronic_state
from .errors import ConfigurationError, QuantumCalculationError
from .model import Atom, MoleculeModel
from .quantum import FrequencyResult, OrcaBackend, QuantumSettings, find_orca


Logger = Callable[[str], None]
_COORDINATE_RE = re.compile(r"^\s*\*\s*xyz\s+(-?\d+)\s+(\d+)", re.IGNORECASE | re.MULTILINE)
_NPROCS_RE = re.compile(r"\bnprocs\s+(\d+)", re.IGNORECASE)
_MAXCORE_RE = re.compile(r"^\s*%maxcore\s+(\d+)", re.IGNORECASE | re.MULTILINE)


@dataclass(frozen=True)
class ExistingOrcaSettings:
    method: Optional[str] = None
    functional: Optional[str] = None
    basis: Optional[str] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    nprocs: Optional[int] = None
    maxcore_mb: Optional[int] = None


@dataclass(frozen=True)
class FrequencyOnlyOutcome:
    result: FrequencyResult
    result_metadata_path: Path
    updated_metadata_path: Optional[Path]
    source_xyz_path: Path


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def read_orca_input_settings(path: Path) -> ExistingOrcaSettings:
    """Read settings from the optimization input stored beside an ORCA XYZ."""

    if not path.is_file():
        return ExistingOrcaSettings()
    text = path.read_text(encoding="utf-8", errors="replace")
    simple_line = next(
        (line.strip()[1:].strip() for line in text.splitlines() if line.lstrip().startswith("!")),
        "",
    )
    tokens = simple_line.split()
    method = None
    functional = None
    basis = None
    if tokens:
        if tokens[0].upper() == "HF":
            method = "hf"
        else:
            method = "dft"
            functional = tokens[0]
        if len(tokens) > 1 and tokens[1].casefold() not in {
            "opt",
            "freq",
            "tightscf",
        }:
            basis = tokens[1]
    coordinate_match = _COORDINATE_RE.search(text)
    nprocs_match = _NPROCS_RE.search(text)
    maxcore_match = _MAXCORE_RE.search(text)
    return ExistingOrcaSettings(
        method=method,
        functional=functional,
        basis=basis,
        charge=int(coordinate_match.group(1)) if coordinate_match else None,
        multiplicity=int(coordinate_match.group(2)) if coordinate_match else None,
        nprocs=int(nprocs_match.group(1)) if nprocs_match else None,
        maxcore_mb=int(maxcore_match.group(1)) if maxcore_match else None,
    )


def _read_metadata(path: Optional[Path]) -> Optional[Dict[str, object]]:
    if path is None or not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"Could not read metadata JSON: {path}") from exc
    if not isinstance(value, dict):
        raise ConfigurationError(f"Metadata JSON must contain an object: {path}")
    return value


def find_related_metadata(xyz_path: Path, explicit: Optional[Path] = None) -> Optional[Path]:
    """Find the model metadata associated with a kept ORCA calculation directory."""

    if explicit is not None:
        explicit = explicit.expanduser().resolve()
        if not explicit.is_file():
            raise ConfigurationError(f"Metadata JSON not found: {explicit}")
        return explicit
    direct = xyz_path.with_suffix(".metadata.json")
    if direct.is_file():
        return direct

    calculation_root = xyz_path.parent.parent
    output_root = calculation_root.parent
    if output_root.is_dir():
        for candidate in sorted(output_root.glob("*.metadata.json")):
            data = _read_metadata(candidate)
            metadata = data.get("metadata", {}) if data else {}
            if not isinstance(metadata, dict):
                continue
            recorded = metadata.get("calculation_files_directory")
            if isinstance(recorded, str) and Path(recorded).expanduser().resolve() == calculation_root.resolve():
                return candidate
    return None


def read_xyz_model(path: Path, metadata: Optional[Dict[str, object]] = None) -> MoleculeModel:
    """Read atom order and coordinates exactly as stored in an optimized XYZ."""

    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        atom_count = int(lines[0].strip())
    except (OSError, IndexError, ValueError) as exc:
        raise QuantumCalculationError(f"Could not read optimized XYZ: {path}") from exc
    if atom_count < 1 or len(lines) < atom_count + 2:
        raise QuantumCalculationError("Optimized XYZ has an invalid atom count")
    atoms = []
    for index, line in enumerate(lines[2 : atom_count + 2]):
        fields = line.split()
        if len(fields) < 4:
            raise QuantumCalculationError(f"Invalid XYZ atom line at index {index}")
        try:
            atoms.append(
                Atom(index, fields[0], float(fields[1]), float(fields[2]), float(fields[3]))
            )
        except ValueError as exc:
            raise QuantumCalculationError(f"Invalid XYZ coordinate at atom {index}") from exc
    name = path.stem
    cid = None
    source_metadata: Dict[str, object] = {}
    if metadata:
        if isinstance(metadata.get("name"), str):
            name = str(metadata["name"])
        if isinstance(metadata.get("cid"), int):
            cid = int(metadata["cid"])
        nested = metadata.get("metadata")
        if isinstance(nested, dict):
            source_metadata.update(nested)
    source_metadata.update(
        {
            "frequency_source_xyz": str(path),
            "frequency_geometry_reoptimized": False,
        }
    )
    return MoleculeModel(cid, name, tuple(atoms), (), source_metadata)


def _metadata_settings(data: Optional[Dict[str, object]]) -> ExistingOrcaSettings:
    nested = data.get("metadata", {}) if data else {}
    if not isinstance(nested, dict):
        return ExistingOrcaSettings()

    def integer(key: str) -> Optional[int]:
        value = nested.get(key)
        return int(value) if isinstance(value, int) else None

    method = nested.get("quantum_method")
    functional = nested.get("functional")
    basis = nested.get("basis")
    return ExistingOrcaSettings(
        method=str(method).lower() if isinstance(method, str) else None,
        functional=str(functional) if isinstance(functional, str) else None,
        basis=str(basis) if isinstance(basis, str) else None,
        charge=integer("charge"),
        multiplicity=integer("multiplicity"),
    )


def _resolve_settings(
    request: FrequencyOnlyRequest,
    input_settings: ExistingOrcaSettings,
    metadata_settings: ExistingOrcaSettings,
    log: Logger,
) -> QuantumSettings:
    method = request.method or metadata_settings.method or input_settings.method or "dft"
    functional = (
        request.functional
        or metadata_settings.functional
        or input_settings.functional
        or DEFAULT_FUNCTIONAL
    )
    basis = request.basis or metadata_settings.basis or input_settings.basis or DEFAULT_BASIS
    charge = request.charge
    if charge is None:
        charge = metadata_settings.charge
    if charge is None:
        charge = input_settings.charge
    multiplicity = request.multiplicity
    if multiplicity is None:
        multiplicity = metadata_settings.multiplicity
    if multiplicity is None:
        multiplicity = input_settings.multiplicity
    if charge is None or multiplicity is None:
        raise ConfigurationError(
            "Could not recover charge/multiplicity from metadata or the matching .inp file; "
            "pass --charge and --multiplicity explicitly"
        )
    if request.charge is not None and input_settings.charge not in {None, request.charge}:
        log("Warning: --charge overrides the charge recorded in the optimization input")
    if request.multiplicity is not None and input_settings.multiplicity not in {
        None,
        request.multiplicity,
    }:
        log("Warning: --multiplicity overrides the value recorded in the optimization input")
    return QuantumSettings(
        method=method,
        functional=functional,
        basis=basis,
        charge=charge,
        multiplicity=multiplicity,
        frequency=True,
        timeout=request.quantum_timeout,
        max_scf_iterations=request.max_scf_iterations,
        nprocs=request.nprocs,
        maxcore_mb=request.maxcore_mb,
        imaginary_threshold_cm1=request.imaginary_threshold_cm1,
    )


def _atomic_write_json(path: Path, data: Dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def _analysis_record(
    result: FrequencyResult,
    source_xyz: Path,
    *,
    imaginary_threshold_cm1: float,
    low_frequency_threshold_cm1: float,
) -> Dict[str, object]:
    return {
        "job_type": "frequency_only",
        "source_xyz": str(source_xyz),
        "geometry_reoptimized": False,
        "backend": result.backend,
        "backend_version": result.backend_version,
        "method": result.method,
        "functional": result.functional,
        "basis": result.basis,
        "charge": result.charge,
        "multiplicity": result.multiplicity,
        "final_energy_hartree": result.final_energy_hartree,
        "frequencies_cm1": list(result.frequencies_cm1),
        "imaginary_frequency_threshold_cm1": imaginary_threshold_cm1,
        "imaginary_frequencies_cm1": list(result.imaginary_frequencies_cm1),
        "low_frequency_threshold_cm1": low_frequency_threshold_cm1,
        "low_frequency_modes_cm1": [
            value
            for value in result.frequencies_cm1
            if value < low_frequency_threshold_cm1
        ],
        "warnings": list(result.warnings),
        "thermochemistry": (
            result.thermochemistry.to_metadata()
            if result.thermochemistry is not None
            else None
        ),
        "completed_at_utc": _utc_now(),
        "orca_input": str(result.input_path),
        "orca_output": str(result.output_path),
    }


def run_frequency_only(
    request: FrequencyOnlyRequest,
    *,
    log: Optional[Logger] = None,
    orca_finder: Callable[[Optional[str]], str] = find_orca,
    orca_backend_class: Type[OrcaBackend] = OrcaBackend,
) -> FrequencyOnlyOutcome:
    """Run Freq only, preserving the XYZ coordinates and skipping Opt entirely."""

    request.validate()
    emit = log or (lambda _message: None)
    xyz_path = request.xyz_path.expanduser().resolve()
    if not xyz_path.is_file():
        raise ConfigurationError(f"Optimized XYZ not found: {xyz_path}")
    metadata_path = find_related_metadata(xyz_path, request.metadata_path)
    metadata_data = _read_metadata(metadata_path)
    input_settings = read_orca_input_settings(xyz_path.with_suffix(".inp"))
    settings = _resolve_settings(
        request,
        input_settings,
        _metadata_settings(metadata_data),
        emit,
    )
    model = read_xyz_model(xyz_path, metadata_data)
    validate_electronic_state(
        (atom.element for atom in model.atoms), settings.charge, settings.multiplicity
    )
    output_dir = (
        request.output_dir.expanduser().resolve()
        if request.output_dir is not None
        else xyz_path.parent
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    backend = orca_backend_class(orca_finder(request.quantum_executable))
    stem = f"{xyz_path.stem}_freq"
    emit(
        f"Running ORCA Freq only: {settings.method.upper()} {settings.basis}, "
        f"charge={settings.charge}, multiplicity={settings.multiplicity}"
    )
    result = backend.frequency(model, settings, output_dir, stem)
    analysis = _analysis_record(
        result,
        xyz_path,
        imaginary_threshold_cm1=request.imaginary_threshold_cm1,
        low_frequency_threshold_cm1=request.low_frequency_threshold_cm1,
    )
    result_metadata_path = output_dir / f"{stem}.metadata.json"
    _atomic_write_json(
        result_metadata_path,
        {
            "schema_version": 1,
            "frequency_analysis": analysis,
        },
    )
    updated_metadata_path = None
    if metadata_path is not None and metadata_data is not None:
        nested = metadata_data.setdefault("metadata", {})
        if not isinstance(nested, dict):
            raise ConfigurationError("Existing metadata 'metadata' field must be an object")
        analyses = nested.setdefault("frequency_analyses", [])
        if not isinstance(analyses, list):
            raise ConfigurationError("Existing frequency_analyses field must be an array")
        analyses.append(analysis)
        _atomic_write_json(metadata_path, metadata_data)
        updated_metadata_path = metadata_path
        emit(f"Updated existing metadata: {metadata_path}")
    return FrequencyOnlyOutcome(result, result_metadata_path, updated_metadata_path, xyz_path)

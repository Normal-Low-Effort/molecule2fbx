"""External ORCA geometry-optimization backend."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from ..errors import (
    ConfigurationError,
    GeometryNotConvergedError,
    QuantumBackendNotFoundError,
    QuantumCalculationError,
    SCFNotConvergedError,
)
from ..model import MoleculeModel, validate_model
from ..structures import geometry_warnings
from .base import FrequencyResult, QuantumResult, QuantumSettings, Thermochemistry


_SAFE_KEYWORD = re.compile(r"^[A-Za-z0-9][A-Za-z0-9+*().,_/-]*$")
_ENERGY_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)")
_FREQUENCY_RE = re.compile(r"^\s*\d+:\s+(-?\d+(?:\.\d+)?)\s+cm\*\*-1", re.MULTILINE)
_VERSION_RE = re.compile(r"Program Version\s+([0-9][\w.+-]*)", re.IGNORECASE)
_FLOAT = r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"
_SCF_FAILURE_MARKERS = (
    "SCF NOT CONVERGED",
    "SCF DID NOT CONVERGE",
    "SCF CONVERGENCE FAILURE",
    "ERROR: SCF",
)


@dataclass(frozen=True)
class ParsedOrcaOutput:
    energy_hartree: float
    version: Optional[str]
    frequencies_cm1: Tuple[float, ...]
    imaginary_frequencies_cm1: Tuple[float, ...]
    warnings: Tuple[str, ...]
    thermochemistry: Optional[Thermochemistry] = None


def _last_float(text: str, label: str) -> Optional[float]:
    matches = re.findall(
        rf"^\s*{label}\s+(?:\.{{2,}}\s*)?{_FLOAT}(?:\s+(?:Eh|K|atm))?(?:\s+.*)?$",
        text,
        re.IGNORECASE | re.MULTILINE,
    )
    return float(matches[-1]) if matches else None


def _parse_thermochemistry(text: str, electronic_energy: float) -> Optional[Thermochemistry]:
    temperature = _last_float(text, r"Temperature")
    pressure = _last_float(text, r"Pressure")
    zero_point = _last_float(text, r"Zero point energy")
    thermal_correction = _last_float(text, r"Total thermal correction")
    thermal_energy = _last_float(text, r"Total thermal energy")
    enthalpy = _last_float(text, r"Total Enthalpy")
    gibbs = _last_float(text, r"Final Gibbs free energy")
    rrho_match = re.findall(
        r"^\s*Quasi RRHO\s+(?:\.{2,}\s*)?(True|False)\s*$",
        text,
        re.IGNORECASE | re.MULTILINE,
    )
    if not any(
        value is not None
        for value in (
            temperature,
            pressure,
            zero_point,
            thermal_correction,
            thermal_energy,
            enthalpy,
            gibbs,
        )
    ):
        return None
    return Thermochemistry(
        temperature_kelvin=temperature,
        pressure_atm=pressure,
        electronic_energy_hartree=electronic_energy,
        zero_point_energy_hartree=zero_point,
        thermal_correction_hartree=thermal_correction,
        thermal_energy_hartree=thermal_energy,
        enthalpy_hartree=enthalpy,
        gibbs_free_energy_hartree=gibbs,
        quasi_rrho=(rrho_match[-1].casefold() == "true") if rrho_match else None,
    )


def find_orca(explicit: Optional[str] = None) -> str:
    """Locate ORCA without modifying PATH or installing licensed software."""

    candidates: List[str] = []
    if explicit:
        candidates.append(explicit)
    executable_env = os.environ.get("ORCA_EXECUTABLE")
    if executable_env:
        candidates.append(executable_env)
    orca_dir = os.environ.get("ORCADIR")
    if orca_dir:
        candidates.append(str(Path(orca_dir) / "orca.exe"))
        candidates.append(str(Path(orca_dir) / "orca"))
    candidates.extend(
        [
            "orca",
            "orca.exe",
            r"C:\Orca_6.1.0\orca.exe",
            r"C:\Orca_6.0.1\orca.exe",
            r"C:\orca\orca.exe",
        ]
    )
    for candidate in candidates:
        path = Path(candidate)
        if path.is_file():
            return str(path.resolve())
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise QuantumBackendNotFoundError(
        "ORCA executable not found. Install ORCA separately, then use --orca <path>, "
        "ORCA_EXECUTABLE, ORCADIR, or PATH."
    )


def _validate_keyword(value: str, option_name: str) -> str:
    if not value or not _SAFE_KEYWORD.fullmatch(value):
        raise ConfigurationError(
            f"Invalid {option_name} keyword {value!r}; use one ORCA simple-input token"
        )
    return value


def render_orca_input(model: MoleculeModel, settings: QuantumSettings) -> str:
    """Render a conservative ORCA optimization input in Angstrom units."""

    validate_model(model)
    basis = _validate_keyword(settings.basis, "basis")
    if settings.method == "hf":
        electronic_method = "HF"
    elif settings.method == "dft":
        electronic_method = _validate_keyword(settings.functional, "functional")
    else:
        raise ConfigurationError("ORCA backend supports method 'dft' or 'hf'")
    keywords = [electronic_method, basis, "Opt", "TightSCF"]
    if settings.frequency:
        keywords.append("Freq")
    lines = [
        "! " + " ".join(keywords),
        "",
        "%scf",
        f"  MaxIter {settings.max_scf_iterations}",
        "end",
        "%geom",
        f"  MaxIter {settings.max_opt_steps}",
        "end",
        "%pal",
        f"  nprocs {settings.nprocs}",
        "end",
        f"%maxcore {settings.maxcore_mb}",
        "",
        f"* xyz {settings.charge} {settings.multiplicity}",
    ]
    lines.extend(
        f"  {atom.element:<3} {atom.x: .12f} {atom.y: .12f} {atom.z: .12f}"
        for atom in model.atoms
    )
    lines.extend(["*", ""])
    return "\n".join(lines)


def render_orca_frequency_input(model: MoleculeModel, settings: QuantumSettings) -> str:
    """Render a frequency-only ORCA job without Opt or a geometry block."""

    validate_model(model)
    basis = _validate_keyword(settings.basis, "basis")
    if settings.method == "hf":
        electronic_method = "HF"
    elif settings.method == "dft":
        electronic_method = _validate_keyword(settings.functional, "functional")
    else:
        raise ConfigurationError("ORCA backend supports method 'dft' or 'hf'")
    lines = [
        f"! {electronic_method} {basis} Freq TightSCF",
        "",
        "%scf",
        f"  MaxIter {settings.max_scf_iterations}",
        "end",
        "%pal",
        f"  nprocs {settings.nprocs}",
        "end",
        f"%maxcore {settings.maxcore_mb}",
        "",
        f"* xyz {settings.charge} {settings.multiplicity}",
    ]
    lines.extend(
        f"  {atom.element:<3} {atom.x: .12f} {atom.y: .12f} {atom.z: .12f}"
        for atom in model.atoms
    )
    lines.extend(["*", ""])
    return "\n".join(lines)


def parse_orca_output(
    text: str,
    *,
    frequency_requested: bool = False,
    require_geometry_convergence: bool = True,
    imaginary_threshold_cm1: float = -20.0,
) -> ParsedOrcaOutput:
    """Validate termination, SCF, optimization, energy, and optional frequencies."""

    upper = text.upper()
    if any(marker in upper for marker in _SCF_FAILURE_MARKERS):
        raise SCFNotConvergedError("ORCA SCF did not converge")
    if "ORCA TERMINATED NORMALLY" not in upper:
        raise QuantumCalculationError("ORCA did not terminate normally")
    if require_geometry_convergence and "THE OPTIMIZATION HAS CONVERGED" not in upper:
        if "OPTIMIZATION DID NOT CONVERGE" in upper or "MAXIMUM NUMBER OF" in upper:
            raise GeometryNotConvergedError("ORCA geometry optimization did not converge")
        raise GeometryNotConvergedError(
            "ORCA output does not contain the geometry-convergence marker"
        )
    energies = [float(value) for value in _ENERGY_RE.findall(text)]
    if not energies:
        raise QuantumCalculationError("ORCA output contains no final electronic energy")
    frequencies = tuple(float(value) for value in _FREQUENCY_RE.findall(text))
    imaginary = tuple(value for value in frequencies if value < imaginary_threshold_cm1)
    warnings = []
    if frequency_requested and not frequencies:
        warnings.append("Frequency calculation was requested but no frequencies were parsed")
    if imaginary:
        warnings.append(
            f"Imaginary frequencies below {imaginary_threshold_cm1:g} cm^-1 were found; "
            "the structure may be a saddle point"
        )
    version_match = _VERSION_RE.search(text)
    return ParsedOrcaOutput(
        energy_hartree=energies[-1],
        version=version_match.group(1) if version_match else None,
        frequencies_cm1=frequencies,
        imaginary_frequencies_cm1=imaginary,
        warnings=tuple(warnings),
        thermochemistry=_parse_thermochemistry(text, energies[-1]),
    )


def read_xyz_coordinates(path: Path, expected_elements: Sequence[str]) -> Tuple[Tuple[float, float, float], ...]:
    """Read a single final ORCA XYZ geometry and verify atom ordering."""

    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        atom_count = int(lines[0].strip())
    except (OSError, IndexError, ValueError) as exc:
        raise QuantumCalculationError(f"Could not read ORCA optimized geometry: {path}") from exc
    if atom_count != len(expected_elements) or len(lines) < atom_count + 2:
        raise QuantumCalculationError("ORCA optimized geometry has an unexpected atom count")
    coordinates = []
    for index, (line, expected) in enumerate(zip(lines[2 : atom_count + 2], expected_elements)):
        fields = line.split()
        if len(fields) < 4 or fields[0].lower() != expected.lower():
            raise QuantumCalculationError(
                f"ORCA optimized geometry changed atom ordering at index {index}"
            )
        try:
            coordinates.append((float(fields[1]), float(fields[2]), float(fields[3])))
        except ValueError as exc:
            raise QuantumCalculationError("ORCA optimized geometry contains invalid coordinates") from exc
    return tuple(coordinates)


class OrcaBackend:
    name = "orca"

    def __init__(self, executable: str):
        self.executable = executable

    def optimize(
        self,
        model: MoleculeModel,
        settings: QuantumSettings,
        work_dir: Path,
        conformer_index: int,
    ) -> QuantumResult:
        work_dir = Path(work_dir).resolve()
        work_dir.mkdir(parents=True, exist_ok=True)
        stem = f"conformer_{conformer_index + 1:03d}"
        input_path = work_dir / f"{stem}.inp"
        output_path = work_dir / f"{stem}.out"
        xyz_path = work_dir / f"{stem}.xyz"
        input_path.write_text(render_orca_input(model, settings), encoding="utf-8")
        try:
            process = subprocess.run(
                [self.executable, input_path.name],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=settings.timeout,
                check=False,
            )
        except FileNotFoundError as exc:
            raise QuantumBackendNotFoundError(
                f"ORCA executable could not be started: {self.executable}"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise QuantumCalculationError(
                f"ORCA calculation timed out after {settings.timeout:g} seconds"
            ) from exc
        except OSError as exc:
            raise QuantumCalculationError(f"Could not start ORCA: {exc}") from exc

        combined_output = process.stdout
        if process.stderr:
            combined_output += "\n\n[stderr]\n" + process.stderr
        output_path.write_text(combined_output, encoding="utf-8")
        # Parse chemistry-specific failure markers before the generic exit code.
        parsed = parse_orca_output(
            combined_output,
            frequency_requested=settings.frequency,
            imaginary_threshold_cm1=settings.imaginary_threshold_cm1,
        )
        if process.returncode != 0:
            raise QuantumCalculationError(f"ORCA exited with code {process.returncode}")
        if not xyz_path.is_file():
            raise QuantumCalculationError(f"ORCA did not create optimized geometry: {xyz_path.name}")
        coordinates = read_xyz_coordinates(xyz_path, [atom.element for atom in model.atoms])
        warnings = list(parsed.warnings)
        optimized = model.with_coordinates(
            coordinates,
            structure_origin="computed_quantum",
            structure_source="ORCA geometry optimization",
            structure_claim="Computed stationary structure; not an experimentally determined structure",
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
            frequency_requested=settings.frequency,
            frequencies_cm1=list(parsed.frequencies_cm1),
            imaginary_frequencies_cm1=list(parsed.imaginary_frequencies_cm1),
            conformer_index=conformer_index,
            calculation_completed_at_utc=datetime.now(timezone.utc)
            .isoformat()
            .replace("+00:00", "Z"),
            thermochemistry=(
                parsed.thermochemistry.to_metadata()
                if parsed.thermochemistry is not None
                else None
            ),
        )
        warnings.extend(geometry_warnings(optimized))
        if warnings:
            optimized = optimized.with_metadata(calculation_warnings=warnings)
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
            conformer_index=conformer_index,
            frequency_requested=settings.frequency,
            frequencies_cm1=parsed.frequencies_cm1,
            imaginary_frequencies_cm1=parsed.imaginary_frequencies_cm1,
            warnings=tuple(warnings),
            calculation_directory=work_dir,
            thermochemistry=parsed.thermochemistry,
        )

    def frequency(
        self,
        model: MoleculeModel,
        settings: QuantumSettings,
        work_dir: Path,
        stem: str,
    ) -> FrequencyResult:
        """Run a frequency-only job using already optimized coordinates."""

        work_dir = Path(work_dir).resolve()
        work_dir.mkdir(parents=True, exist_ok=True)
        input_path = work_dir / f"{stem}.inp"
        output_path = work_dir / f"{stem}.out"
        input_path.write_text(render_orca_frequency_input(model, settings), encoding="utf-8")
        try:
            process = subprocess.run(
                [self.executable, input_path.name],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=settings.timeout,
                check=False,
            )
        except FileNotFoundError as exc:
            raise QuantumBackendNotFoundError(
                f"ORCA executable could not be started: {self.executable}"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise QuantumCalculationError(
                f"ORCA frequency calculation timed out after {settings.timeout:g} seconds"
            ) from exc
        except OSError as exc:
            raise QuantumCalculationError(f"Could not start ORCA: {exc}") from exc

        combined_output = process.stdout
        if process.stderr:
            combined_output += "\n\n[stderr]\n" + process.stderr
        output_path.write_text(combined_output, encoding="utf-8")
        parsed = parse_orca_output(
            combined_output,
            frequency_requested=True,
            require_geometry_convergence=False,
            imaginary_threshold_cm1=settings.imaginary_threshold_cm1,
        )
        if process.returncode != 0:
            raise QuantumCalculationError(f"ORCA exited with code {process.returncode}")
        if not parsed.frequencies_cm1:
            raise QuantumCalculationError("ORCA frequency job produced no parseable frequencies")
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

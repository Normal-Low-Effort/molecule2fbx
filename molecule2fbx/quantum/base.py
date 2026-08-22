"""Backend-neutral quantum chemistry interfaces."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Protocol, Tuple

from ..model import MoleculeModel


def detect_nprocs() -> int:
    """Return the host logical processor count, with a safe fallback."""

    return os.cpu_count() or 1


@dataclass(frozen=True)
class QuantumSettings:
    method: str
    functional: str
    basis: str
    charge: int
    multiplicity: int
    frequency: bool = False
    timeout: float = 3600.0
    max_opt_steps: int = 100
    max_scf_iterations: int = 200
    nprocs: int = field(default_factory=detect_nprocs)
    maxcore_mb: int = 1000
    imaginary_threshold_cm1: float = -20.0


@dataclass(frozen=True)
class Thermochemistry:
    """Thermochemical values parsed from a completed frequency calculation."""

    temperature_kelvin: Optional[float] = None
    pressure_atm: Optional[float] = None
    electronic_energy_hartree: Optional[float] = None
    zero_point_energy_hartree: Optional[float] = None
    thermal_correction_hartree: Optional[float] = None
    thermal_energy_hartree: Optional[float] = None
    enthalpy_hartree: Optional[float] = None
    gibbs_free_energy_hartree: Optional[float] = None
    quasi_rrho: Optional[bool] = None

    def to_metadata(self) -> dict:
        return {
            "temperature_kelvin": self.temperature_kelvin,
            "pressure_atm": self.pressure_atm,
            "electronic_energy_hartree": self.electronic_energy_hartree,
            "zero_point_energy_hartree": self.zero_point_energy_hartree,
            "thermal_correction_hartree": self.thermal_correction_hartree,
            "thermal_energy_hartree": self.thermal_energy_hartree,
            "enthalpy_hartree": self.enthalpy_hartree,
            "gibbs_free_energy_hartree": self.gibbs_free_energy_hartree,
            "quasi_rrho": self.quasi_rrho,
        }


@dataclass(frozen=True)
class QuantumResult:
    model: MoleculeModel
    backend: str
    backend_version: Optional[str]
    method: str
    functional: Optional[str]
    basis: str
    charge: int
    multiplicity: int
    final_energy_hartree: float
    geometry_converged: bool
    scf_converged: bool
    conformer_index: int
    frequency_requested: bool
    frequencies_cm1: Tuple[float, ...] = ()
    imaginary_frequencies_cm1: Tuple[float, ...] = ()
    warnings: Tuple[str, ...] = ()
    calculation_directory: Optional[Path] = None
    thermochemistry: Optional[Thermochemistry] = None


@dataclass(frozen=True)
class FrequencyResult:
    backend: str
    backend_version: Optional[str]
    method: str
    functional: Optional[str]
    basis: str
    charge: int
    multiplicity: int
    final_energy_hartree: float
    frequencies_cm1: Tuple[float, ...]
    imaginary_frequencies_cm1: Tuple[float, ...]
    warnings: Tuple[str, ...]
    input_path: Path
    output_path: Path
    thermochemistry: Optional[Thermochemistry] = None


class QuantumBackend(Protocol):
    name: str

    def optimize(
        self,
        model: MoleculeModel,
        settings: QuantumSettings,
        work_dir: Path,
        conformer_index: int,
    ) -> QuantumResult:
        ...

    def frequency(
        self,
        model: MoleculeModel,
        settings: QuantumSettings,
        work_dir: Path,
        stem: str,
    ) -> FrequencyResult:
        ...

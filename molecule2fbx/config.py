"""Configuration objects shared by the CLI and future user interfaces."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple


DEFAULT_FUNCTIONAL = "B3LYP"
DEFAULT_BASIS = "def2-SVP"
DEFAULT_CHARGE = 0
DEFAULT_MULTIPLICITY = 1


@dataclass(frozen=True)
class ConversionRequest:
    """A UI-neutral request for one molecule conversion."""

    cid: Optional[int] = None
    smiles: Optional[str] = None
    name: Optional[str] = None
    method: str = "auto"
    output_dir: Path = Path("output")
    blender_executable: Optional[str] = None
    quantum_backend: str = "orca"
    quantum_executable: Optional[str] = None
    functional: Optional[str] = None
    basis: Optional[str] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    conformers: int = 1
    conformer_pool: Optional[int] = None
    conformer_rmsd_threshold: float = 0.75
    frequency: bool = False
    frequency_window_kj: Optional[float] = None
    frequency_max: Optional[int] = None
    save_all_conformers: bool = False
    keep_calculation_files: bool = False
    reuse_calculations: Optional[Path] = None
    allow_metals: bool = False
    api_timeout: float = 30.0
    blender_timeout: float = 180.0
    quantum_timeout: float = 3600.0
    max_opt_steps: int = 100
    max_scf_iterations: int = 200
    nprocs: int = field(default_factory=lambda: os.cpu_count() or 1)
    maxcore_mb: int = 1000
    ensemble: bool = False
    strict_stereochemistry: bool = False
    random_seed: int = 0xF00D
    embedding_prune_rmsd: float = -1.0
    forcefield_energy_window_kj: Optional[float] = None
    dft_energy_window_kj: Optional[float] = None
    dft_rmsd_threshold: Optional[float] = None
    frequency_include: Tuple[int, ...] = ()
    imaginary_threshold_cm1: float = -20.0
    low_frequency_threshold_cm1: float = 50.0
    expected_stereocenters: Optional[int] = None
    stereochemistry_label: Optional[str] = None

    @property
    def effective_functional(self) -> str:
        return self.functional or DEFAULT_FUNCTIONAL

    @property
    def effective_basis(self) -> str:
        return self.basis or DEFAULT_BASIS

    @property
    def effective_charge(self) -> int:
        return DEFAULT_CHARGE if self.charge is None else self.charge

    @property
    def effective_multiplicity(self) -> int:
        return DEFAULT_MULTIPLICITY if self.multiplicity is None else self.multiplicity

    @property
    def effective_conformer_pool(self) -> int:
        return self.conformer_pool or self.conformers

    @property
    def selective_frequency(self) -> bool:
        return self.frequency and (
            self.frequency_window_kj is not None or self.frequency_max is not None
        )

    @property
    def effective_dft_rmsd_threshold(self) -> float:
        return self.dft_rmsd_threshold or self.conformer_rmsd_threshold

    @property
    def effective_forcefield_energy_window_kj(self) -> float:
        if self.forcefield_energy_window_kj is not None:
            return self.forcefield_energy_window_kj
        return 10.0 if self.ensemble else float("inf")

    def validate(self) -> None:
        if (self.cid is None) == (self.smiles is None):
            raise ValueError("Specify exactly one of CID or SMILES")
        if self.method not in {"auto", "pubchem", "forcefield", "dft", "hf"}:
            raise ValueError(f"Unsupported method: {self.method}")
        if self.smiles is not None and self.method == "pubchem":
            raise ValueError("--method pubchem requires a CID")
        if self.quantum_backend != "orca":
            raise ValueError(f"Unsupported quantum backend: {self.quantum_backend}")
        if self.reuse_calculations is not None and self.method not in {"dft", "hf"}:
            raise ValueError("--reuse-calculations requires --method dft or hf")
        if self.ensemble and self.method not in {"dft", "hf"}:
            raise ValueError("--ensemble requires --method dft or hf")
        if self.ensemble and self.smiles is None:
            raise ValueError("--ensemble currently requires an explicit --smiles input")
        if self.conformers < 1:
            raise ValueError("--conformers must be at least 1")
        if self.effective_conformer_pool < self.conformers:
            raise ValueError("--conformer-pool must be at least --conformers")
        if self.ensemble and not 100 <= self.effective_conformer_pool <= 500:
            raise ValueError("--ensemble requires --conformer-pool between 100 and 500")
        if self.conformer_pool is not None and self.method not in {"dft", "hf"}:
            raise ValueError("--conformer-pool requires --method dft or hf")
        if self.conformer_rmsd_threshold <= 0:
            raise ValueError("--conformer-rmsd-threshold must be greater than zero")
        if self.effective_dft_rmsd_threshold <= 0:
            raise ValueError("--dft-rmsd-threshold must be greater than zero")
        if self.random_seed < 0:
            raise ValueError("--random-seed must be zero or greater")
        if self.embedding_prune_rmsd < 0 and self.embedding_prune_rmsd != -1:
            raise ValueError("--embedding-prune-rmsd must be -1 or zero or greater")
        if (
            self.forcefield_energy_window_kj is not None
            and self.forcefield_energy_window_kj < 0
        ):
            raise ValueError("--forcefield-energy-window-kj must be zero or greater")
        if self.dft_energy_window_kj is not None and self.dft_energy_window_kj < 0:
            raise ValueError("--dft-energy-window-kj must be zero or greater")
        frequency_selection = (
            self.frequency_window_kj is not None or self.frequency_max is not None
        )
        if frequency_selection and not self.frequency:
            raise ValueError(
                "--frequency-window-kj and --frequency-max require --frequency"
            )
        if frequency_selection and self.method not in {"dft", "hf"}:
            raise ValueError("Selective frequency options require --method dft or hf")
        if self.frequency_window_kj is not None and self.frequency_window_kj < 0:
            raise ValueError("--frequency-window-kj must be zero or greater")
        if self.frequency_max is not None and self.frequency_max < 1:
            raise ValueError("--frequency-max must be at least 1")
        if self.frequency_include and not self.frequency:
            raise ValueError("--frequency-include requires --frequency")
        if any(index < 0 or index >= self.conformers for index in self.frequency_include):
            raise ValueError("--frequency-include must identify a requested conformer")
        if self.imaginary_threshold_cm1 >= 0:
            raise ValueError("--imaginary-threshold-cm1 must be negative")
        if self.low_frequency_threshold_cm1 <= 0:
            raise ValueError("--low-frequency-threshold-cm1 must be greater than zero")
        if self.expected_stereocenters is not None and self.expected_stereocenters < 0:
            raise ValueError("--expected-stereocenters must be zero or greater")
        if self.max_opt_steps < 1 or self.max_scf_iterations < 1:
            raise ValueError("Optimization and SCF iteration limits must be positive")
        if self.nprocs < 1 or self.maxcore_mb < 1:
            raise ValueError("--nprocs and --maxcore must be positive")
        if min(self.api_timeout, self.blender_timeout, self.quantum_timeout) <= 0:
            raise ValueError("Timeouts must be greater than zero")


@dataclass(frozen=True)
class FrequencyOnlyRequest:
    """Run ORCA frequencies on coordinates from an existing optimized XYZ."""

    xyz_path: Path
    output_dir: Optional[Path] = None
    metadata_path: Optional[Path] = None
    method: Optional[str] = None
    functional: Optional[str] = None
    basis: Optional[str] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    quantum_executable: Optional[str] = None
    quantum_timeout: float = 3600.0
    max_scf_iterations: int = 200
    nprocs: int = field(default_factory=lambda: os.cpu_count() or 1)
    maxcore_mb: int = 1000
    imaginary_threshold_cm1: float = -20.0
    low_frequency_threshold_cm1: float = 50.0

    def validate(self) -> None:
        if self.method not in {None, "dft", "hf"}:
            raise ValueError("--frequency-only supports --method dft or hf")
        if self.multiplicity is not None and self.multiplicity < 1:
            raise ValueError("--multiplicity must be at least 1")
        if self.nprocs < 1 or self.maxcore_mb < 1 or self.max_scf_iterations < 1:
            raise ValueError("ORCA resource and iteration values must be positive")
        if self.quantum_timeout <= 0:
            raise ValueError("--quantum-timeout must be greater than zero")
        if self.imaginary_threshold_cm1 >= 0:
            raise ValueError("--imaginary-threshold-cm1 must be negative")
        if self.low_frequency_threshold_cm1 <= 0:
            raise ValueError("--low-frequency-threshold-cm1 must be greater than zero")

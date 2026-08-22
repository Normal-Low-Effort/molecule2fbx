import subprocess

import pytest

from molecule2fbx.errors import (
    GeometryNotConvergedError,
    QuantumBackendNotFoundError,
    SCFNotConvergedError,
)
from molecule2fbx.quantum.base import QuantumSettings
from molecule2fbx.quantum.orca import (
    OrcaBackend,
    parse_orca_output,
    render_orca_frequency_input,
    render_orca_input,
)

from conftest import hydrogen_model


SUCCESS_OUTPUT = """
Program Version 6.1.0
SCF CONVERGED AFTER 8 CYCLES
FINAL SINGLE POINT ENERGY       -1.130000000000
***        THE OPTIMIZATION HAS CONVERGED     ***
   0:      -35.50 cm**-1
   1:       89.20 cm**-1
****ORCA TERMINATED NORMALLY****
"""

FREQUENCY_ONLY_OUTPUT = """
Program Version 6.1.0
SCF CONVERGED AFTER 7 CYCLES
FINAL SINGLE POINT ENERGY       -1.125000000000
   0:      -21.50 cm**-1
   1:      -10.00 cm**-1
   2:        0.00 cm**-1
   3:      750.25 cm**-1
****ORCA TERMINATED NORMALLY****
"""


def settings(frequency=False):
    return QuantumSettings(
        method="dft",
        functional="B3LYP",
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
        frequency=frequency,
    )


def test_renders_geometry_optimization_not_single_point():
    rendered = render_orca_input(hydrogen_model(), settings())
    assert "! B3LYP def2-SVP Opt TightSCF" in rendered
    assert "%geom" in rendered
    assert "* xyz 0 1" in rendered


def test_existing_frequency_mode_still_renders_opt_and_freq():
    rendered = render_orca_input(hydrogen_model(), settings(frequency=True))
    assert " Opt " in rendered
    assert "Freq" in rendered
    assert "%geom" in rendered


def test_frequency_only_input_never_contains_optimization():
    rendered = render_orca_frequency_input(hydrogen_model(), settings(frequency=True))
    assert "! B3LYP def2-SVP Freq TightSCF" in rendered
    assert " Opt " not in rendered
    assert "%geom" not in rendered
    assert "%pal\n" in rendered


def test_quantum_settings_detects_logical_processors(monkeypatch):
    monkeypatch.setattr("molecule2fbx.quantum.base.os.cpu_count", lambda: 12)
    detected = QuantumSettings(
        method="dft",
        functional="B3LYP",
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
    )
    assert detected.nprocs == 12


def test_orca_pal_uses_explicit_nprocs():
    explicit = QuantumSettings(
        method="dft",
        functional="B3LYP",
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
        nprocs=8,
    )
    rendered = render_orca_input(hydrogen_model(), explicit)
    assert "%pal\n  nprocs 8\nend" in rendered


def test_parses_converged_energy_and_imaginary_frequency():
    parsed = parse_orca_output(SUCCESS_OUTPUT, frequency_requested=True)
    assert parsed.energy_hartree == pytest.approx(-1.13)
    assert parsed.version == "6.1.0"
    assert parsed.imaginary_frequencies_cm1 == (-35.5,)
    assert any("saddle point" in warning for warning in parsed.warnings)


def test_frequency_only_parser_skips_geometry_convergence_requirement():
    parsed = parse_orca_output(
        FREQUENCY_ONLY_OUTPUT,
        frequency_requested=True,
        require_geometry_convergence=False,
    )
    assert parsed.frequencies_cm1 == (-21.5, -10.0, 0.0, 750.25)
    assert parsed.imaginary_frequencies_cm1 == (-21.5,)


def test_parse_frequency_thermochemistry_and_custom_imaginary_threshold():
    output = """
Program Version 6.1.1
FINAL SINGLE POINT ENERGY      -40.123456
  0:       -25.00 cm**-1
  1:       -10.00 cm**-1
Temperature                  ... 298.15 K
Pressure                     ... 1.00 atm
Quasi RRHO                   ... True
Zero point energy            ... 0.100000 Eh  62.75 kcal/mol
Total thermal correction         0.120000 Eh  75.30 kcal/mol
Total thermal energy             -40.003456 Eh
Total Enthalpy               ... -39.999000 Eh
Final Gibbs free energy      ... -40.050000 Eh
ORCA TERMINATED NORMALLY
"""
    parsed = parse_orca_output(
        output,
        frequency_requested=True,
        require_geometry_convergence=False,
        imaginary_threshold_cm1=-20.0,
    )
    assert parsed.imaginary_frequencies_cm1 == (-25.0,)
    assert parsed.thermochemistry is not None
    assert parsed.thermochemistry.zero_point_energy_hartree == pytest.approx(0.1)
    assert parsed.thermochemistry.gibbs_free_energy_hartree == pytest.approx(-40.05)
    assert parsed.thermochemistry.temperature_kelvin == pytest.approx(298.15)
    assert parsed.thermochemistry.quasi_rrho is True


def test_default_parser_still_requires_geometry_convergence():
    with pytest.raises(GeometryNotConvergedError):
        parse_orca_output(FREQUENCY_ONLY_OUTPUT, frequency_requested=True)


def test_reports_scf_nonconvergence():
    output = "SCF NOT CONVERGED\nORCA TERMINATED NORMALLY"
    with pytest.raises(SCFNotConvergedError):
        parse_orca_output(output)


def test_reports_geometry_nonconvergence_even_after_normal_termination():
    output = """
FINAL SINGLE POINT ENERGY -1.0
The optimization did not converge but reached the maximum number of optimization cycles.
ORCA TERMINATED NORMALLY
"""
    with pytest.raises(GeometryNotConvergedError):
        parse_orca_output(output)


def test_backend_missing_executable(tmp_path):
    backend = OrcaBackend(str(tmp_path / "missing-orca.exe"))
    with pytest.raises(QuantumBackendNotFoundError):
        backend.optimize(hydrogen_model(), settings(), tmp_path / "job", 0)


def test_mocked_orca_geometry_optimization_flow(monkeypatch, tmp_path):
    def fake_run(command, cwd, **kwargs):
        job_dir = __import__("pathlib").Path(cwd)
        (job_dir / "conformer_001.xyz").write_text(
            "2\noptimized\nH 0 0 0\nH 0.75 0 0\n", encoding="utf-8"
        )
        return subprocess.CompletedProcess(command, 0, stdout=SUCCESS_OUTPUT, stderr="")

    monkeypatch.setattr("molecule2fbx.quantum.orca.subprocess.run", fake_run)
    result = OrcaBackend("orca").optimize(
        hydrogen_model(), settings(frequency=True), tmp_path / "job", 0
    )
    assert result.geometry_converged is True
    assert result.scf_converged is True
    assert result.model.metadata["structure_origin"] == "computed_quantum"
    assert result.model.atoms[1].x == pytest.approx(0.75)


def test_mocked_orca_frequency_only_flow_does_not_expect_new_xyz(monkeypatch, tmp_path):
    def fake_run(command, cwd, **kwargs):
        return subprocess.CompletedProcess(
            command, 0, stdout=FREQUENCY_ONLY_OUTPUT, stderr=""
        )

    monkeypatch.setattr("molecule2fbx.quantum.orca.subprocess.run", fake_run)
    result = OrcaBackend("orca").frequency(
        hydrogen_model(), settings(frequency=True), tmp_path / "freq", "conformer_001_freq"
    )
    rendered = result.input_path.read_text(encoding="utf-8")
    assert " Opt " not in rendered
    assert "%geom" not in rendered
    assert result.imaginary_frequencies_cm1 == (-21.5,)

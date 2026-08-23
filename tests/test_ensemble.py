from pathlib import Path

import pytest

from molecule2fbx.cli import _request_from_args, build_parser
from molecule2fbx.config import ConversionRequest
from molecule2fbx.ensemble import (
    _forcefield_used_summary,
    build_ensemble_report,
    cluster_optimized_results,
    screen_forcefield_candidates,
)
from molecule2fbx.errors import ConfigurationError
from molecule2fbx.model import Atom, MoleculeModel
from molecule2fbx.pipeline import _quantum_initial_candidates
from molecule2fbx.pipeline import run_conversion
from molecule2fbx.quantum.base import FrequencyResult, QuantumResult, Thermochemistry
from molecule2fbx.structures import ConformerCandidate, StereochemistryValidation


def _model(points, *, index=0):
    return MoleculeModel(
        None,
        "ensemble-test",
        tuple(Atom(i, "C", *point) for i, point in enumerate(points)),
        (),
        {
            "conformer_index": index,
            "canonical_isomeric_smiles": "CCCC",
            "original_smiles": "CCCC",
            "forcefield": "MMFF94s",
            "stereochemistry": {"specified_count": 0},
            "etkdg": {"rdkit_version": "test", "embedding_failures": 0},
        },
    )


BASE = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
DISTINCT = ((0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 0.6))
TRANSLATED = tuple((x + 4.0, y - 3.0, z + 2.0) for x, y, z in BASE)


def test_forcefield_summary_falls_back_to_current_screening_records():
    assert _forcefield_used_summary(
        {},
        {
            "records": [
                {"forcefield": "MMFF94s"},
                {"forcefield": "MMFF94s"},
            ]
        },
    ) == "MMFF94s"


def _candidate(points, index, energy, converged=True):
    model = _model(points, index=index).with_metadata(
        forcefield_energy_kcal_mol=energy,
        forcefield_converged=converged,
        forcefield_optimization_status=0 if converged else 1,
    )
    return ConformerCandidate(
        model,
        index,
        energy,
        "MMFF94s",
        converged,
        None,
        None if converged else "forcefield_optimization_failure",
    )


def _result(points, index, energy, *, frequency=False, gibbs=None):
    thermo = (
        Thermochemistry(
            temperature_kelvin=298.15,
            pressure_atm=1.0,
            electronic_energy_hartree=energy,
            zero_point_energy_hartree=0.01,
            thermal_correction_hartree=0.02,
            enthalpy_hartree=energy + 0.03,
            gibbs_free_energy_hartree=gibbs,
        )
        if frequency
        else None
    )
    return QuantumResult(
        model=_model(points, index=index).with_metadata(structure_origin="computed_quantum"),
        backend="ORCA",
        backend_version="6.1.1",
        method="dft",
        functional="B3LYP",
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
        final_energy_hartree=energy,
        geometry_converged=True,
        scf_converged=True,
        conformer_index=index,
        frequency_requested=frequency,
        frequencies_cm1=(-10.0, 35.0, 100.0) if frequency else (),
        imaginary_frequencies_cm1=(),
        calculation_directory=Path(f"conf{index + 1:03d}"),
        thermochemistry=thermo,
    )


def test_ensemble_cli_applies_documented_defaults():
    args = build_parser().parse_args(["--smiles", "O", "--ensemble"])
    request = _request_from_args(args)
    request.validate()
    assert request.method == "dft"
    assert request.conformer_pool == 200
    assert request.conformers == 10
    assert request.forcefield_energy_window_kj == pytest.approx(10.0)
    assert request.frequency is True
    assert request.frequency_window_kj == pytest.approx(5.0)
    assert request.frequency_max == 3
    assert request.strict_stereochemistry is True
    assert request.keep_calculation_files is True


def test_strict_stereo_stops_before_conformer_generation():
    request = ConversionRequest(
        smiles="CC(F)Cl",
        method="dft",
        strict_stereochemistry=True,
    )
    with pytest.raises(ConfigurationError, match="unresolved stereochemistry"):
        _quantum_initial_candidates(request, lambda _message: None)


def test_forcefield_screening_records_energy_and_rmsd_exclusions():
    candidates = [
        _candidate(BASE, 0, 0.0),
        _candidate(TRANSLATED, 1, 1.0),
        _candidate(DISTINCT, 2, 3.0),
        _candidate(DISTINCT, 3, 0.5, converged=False),
    ]
    selected, report = screen_forcefield_candidates(
        candidates,
        energy_window_kj=10.0,
        rmsd_threshold=0.1,
        maximum=10,
    )
    assert [item.model.metadata["conformer_pool_index"] for item in selected] == [0]
    records = {item["pool_index"]: item for item in report["records"]}
    assert records[1]["excluded_reason"] == "rmsd_duplicate"
    assert records[2]["excluded_reason"] == "forcefield_energy_window"
    assert records[3]["excluded_reason"] == "forcefield_optimization_failure"
    assert report["threshold_relaxed"] is False


def test_post_dft_clustering_keeps_lowest_energy_representative():
    results = [
        _result(BASE, 0, -10.000),
        _result(TRANSLATED, 1, -9.999),
        _result(DISTINCT, 2, -9.998),
    ]
    selected, report = cluster_optimized_results(results, rmsd_threshold=0.1)
    assert [item.conformer_index for item in selected] == [0, 2]
    records = {item["conformer_index"]: item for item in report["records"]}
    assert records[1]["excluded_reason"] == "dft_rmsd_duplicate"
    assert records[1]["duplicate_of_conformer_index"] == 0


def test_ensemble_report_does_not_claim_unchecked_minimum():
    results = [
        _result(BASE, 0, -10.000, frequency=True, gibbs=-9.900),
        _result(DISTINCT, 2, -9.998),
    ]
    selected, post = cluster_optimized_results(results, rmsd_threshold=0.1)
    request = ConversionRequest(
        smiles="CCCC",
        method="dft",
        ensemble=True,
        conformers=2,
        conformer_pool=2,
        frequency=True,
        frequency_window_kj=5.0,
        frequency_max=1,
    )
    payload = build_ensemble_report(
        request,
        selected,
        initial_screening={"records": []},
        post_dft_screening=post,
        frequency_selection={"selected_conformer_indices": [0]},
    )
    assert "not a proven global minimum" in payload["interpretation"]["best_structure_claim"]
    entries = {item["conformer_index"]: item for item in payload["final_conformer_ensemble"]}
    assert entries[0]["local_minimum_assessment"] == "local_minimum_candidate"
    assert entries[0]["low_frequency_modes_cm1"] == [-10.0, 35.0]
    assert entries[2]["local_minimum_assessment"] == "not_evaluated"
    assert entries[2]["gibbs_energy_hartree"] is None


def test_ensemble_pipeline_writes_report_and_runs_freq_only_after_opt(monkeypatch, tmp_path):
    candidates = [
        _candidate(BASE, 0, 0.0),
        _candidate(DISTINCT, 1, 0.5),
    ]
    screening = {"records": [], "selected_for_dft": 2}
    monkeypatch.setattr(
        "molecule2fbx.pipeline._quantum_initial_candidates",
        lambda request, log: (candidates, "CCCC", screening),
    )
    monkeypatch.setattr(
        "molecule2fbx.pipeline.validate_model_stereochemistry",
        lambda smiles, model: StereochemistryValidation((), (), True),
    )

    class FakeBackend:
        optimized = []
        frequencies = []

        def __init__(self, executable):
            self.executable = executable

        def optimize(self, model, settings, work_dir, conformer_index):
            self.optimized.append(conformer_index)
            energy = -10.0 + conformer_index * 0.001
            return QuantumResult(
                model=model.with_metadata(structure_origin="computed_quantum"),
                backend="ORCA",
                backend_version="6.1.1",
                method="dft",
                functional="B3LYP",
                basis="def2-SVP",
                charge=0,
                multiplicity=1,
                final_energy_hartree=energy,
                geometry_converged=True,
                scf_converged=True,
                conformer_index=conformer_index,
                frequency_requested=False,
                calculation_directory=Path(work_dir),
            )

        def frequency(self, model, settings, work_dir, stem):
            conformer_index = int(model.metadata["conformer_index"])
            self.frequencies.append(conformer_index)
            energy = -10.0 + conformer_index * 0.001
            return FrequencyResult(
                backend="ORCA",
                backend_version="6.1.1",
                method="dft",
                functional="B3LYP",
                basis="def2-SVP",
                charge=0,
                multiplicity=1,
                final_energy_hartree=energy,
                frequencies_cm1=(25.0, 100.0),
                imaginary_frequencies_cm1=(),
                warnings=(),
                input_path=Path(work_dir) / f"{stem}.inp",
                output_path=Path(work_dir) / f"{stem}.out",
                thermochemistry=Thermochemistry(
                    temperature_kelvin=298.15,
                    pressure_atm=1.0,
                    electronic_energy_hartree=energy,
                    gibbs_free_energy_hartree=energy + 0.05,
                ),
            )

    def exporter(model, output_path, blender_executable, timeout):
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(b"FBX")
        return output_path

    request = ConversionRequest(
        smiles="CCCC",
        method="dft",
        output_dir=tmp_path,
        ensemble=True,
        strict_stereochemistry=True,
        conformers=2,
        conformer_pool=100,
        forcefield_energy_window_kj=10.0,
        frequency=True,
        frequency_window_kj=5.0,
        frequency_max=1,
        keep_calculation_files=True,
    )
    outcome = run_conversion(
        request,
        blender_finder=lambda explicit: "blender",
        orca_finder=lambda explicit: "orca",
        orca_backend_class=FakeBackend,
        exporter=exporter,
    )
    assert FakeBackend.optimized == [0, 1]
    assert FakeBackend.frequencies == [0]
    assert outcome.ensemble_report_path is not None
    assert outcome.ensemble_report_path.is_file()
    assert outcome.selected_model.metadata["ensemble_summary"]["best_conformer_id"] == "conf001"

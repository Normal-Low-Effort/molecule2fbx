import json
from pathlib import Path

import pytest

from molecule2fbx.config import ConversionRequest
from molecule2fbx.errors import ConfigurationError
from molecule2fbx.pipeline import (
    _run_quantum,
    _run_selective_frequencies,
    _select_quantum_candidates,
    run_conversion,
)
from molecule2fbx.quantum.base import FrequencyResult, QuantumResult, QuantumSettings
from molecule2fbx.quantum.orca import render_orca_frequency_input, render_orca_input
from molecule2fbx.structures import ConformerCandidate

from conftest import fake_exporter, hydrogen_model


OPT_OUTPUT = """
Program Version 6.1.0
SCF CONVERGED AFTER 8 CYCLES
FINAL SINGLE POINT ENERGY       -1.130000000000
***        THE OPTIMIZATION HAS CONVERGED     ***
****ORCA TERMINATED NORMALLY****
"""

FREQ_OUTPUT = """
Program Version 6.1.0
SCF CONVERGED AFTER 7 CYCLES
FINAL SINGLE POINT ENERGY       -1.129000000000
   0:      -25.00 cm**-1
   1:      100.00 cm**-1
****ORCA TERMINATED NORMALLY****
"""


def _candidates(count):
    return [
        ConformerCandidate(
            hydrogen_model().with_metadata(conformer_index=index),
            index,
            float(index),
            "MMFF94s",
            True,
        )
        for index in range(count)
    ]


def _request(
    root: Path,
    *,
    conformers=4,
    frequency=False,
    functional="B3LYP",
    conformer_pool=None,
    frequency_window_kj=None,
    frequency_max=None,
):
    return ConversionRequest(
        smiles="[H][H]",
        method="dft",
        conformers=conformers,
        conformer_pool=conformer_pool,
        functional=functional,
        frequency=frequency,
        frequency_window_kj=frequency_window_kj,
        frequency_max=frequency_max,
        reuse_calculations=root,
    )


def _write_completed_optimization(root: Path, candidate, *, functional="B3LYP"):
    index = candidate.conformer_index
    stem = f"conformer_{index + 1:03d}"
    directory = root / stem
    directory.mkdir(parents=True)
    settings = QuantumSettings(
        method="dft",
        functional=functional,
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
    )
    (directory / f"{stem}.inp").write_text(
        render_orca_input(candidate.model, settings), encoding="utf-8"
    )
    (directory / f"{stem}.out").write_text(OPT_OUTPUT, encoding="utf-8")
    (directory / f"{stem}.xyz").write_text(
        "2\noptimized\nH 0 0 0\nH 0.73 0 0\n", encoding="utf-8"
    )


def _optimized_result(model, settings, work_dir, index):
    energy = -1.12 - index / 1000
    optimized = model.with_metadata(
        structure_origin="computed_quantum",
        final_energy_hartree=energy,
        conformer_index=index,
    )
    return QuantumResult(
        model=optimized,
        backend="ORCA",
        backend_version="test",
        method=settings.method,
        functional=settings.functional,
        basis=settings.basis,
        charge=settings.charge,
        multiplicity=settings.multiplicity,
        final_energy_hartree=energy,
        geometry_converged=True,
        scf_converged=True,
        conformer_index=index,
        frequency_requested=settings.frequency,
        calculation_directory=work_dir,
    )


class RecordingBackend:
    def __init__(self):
        self.optimized_indices = []
        self.optimize_frequency_flags = []
        self.frequency_indices = []

    def optimize(self, model, settings, work_dir, conformer_index):
        self.optimized_indices.append(conformer_index)
        self.optimize_frequency_flags.append(settings.frequency)
        return _optimized_result(model, settings, work_dir, conformer_index)

    def frequency(self, model, settings, work_dir, stem):
        self.frequency_indices.append(int(stem.split("_")[1]) - 1)
        work_dir.mkdir(parents=True, exist_ok=True)
        return FrequencyResult(
            backend="ORCA",
            backend_version="test",
            method=settings.method,
            functional=settings.functional,
            basis=settings.basis,
            charge=settings.charge,
            multiplicity=settings.multiplicity,
            final_energy_hartree=-1.129,
            frequencies_cm1=(-25.0, 100.0),
            imaginary_frequencies_cm1=(-25.0,),
            warnings=("imaginary frequency",),
            input_path=work_dir / f"{stem}.inp",
            output_path=work_dir / f"{stem}.out",
        )


def test_reuses_conformer_001_and_optimizes_only_three_additions(tmp_path):
    candidates = _candidates(4)
    _write_completed_optimization(tmp_path, candidates[0])
    backend = RecordingBackend()

    results, failures = _run_quantum(
        _request(tmp_path), candidates, backend, tmp_path, lambda _message: None
    )

    assert [result.conformer_index for result in results] == [0, 1, 2, 3]
    assert backend.optimized_indices == [1, 2, 3]
    assert results[0].model.atoms[1].x == pytest.approx(0.73)
    assert results[0].model.metadata["reused_quantum_result"] is True
    assert failures == []


def test_absent_reuse_results_preserves_four_optimization_behavior(tmp_path):
    candidates = _candidates(4)
    backend = RecordingBackend()

    results, _ = _run_quantum(
        _request(tmp_path), candidates, backend, tmp_path, lambda _message: None
    )

    assert len(results) == 4
    assert backend.optimized_indices == [0, 1, 2, 3]


def test_incompatible_existing_result_is_not_overwritten(tmp_path):
    candidates = _candidates(1)
    _write_completed_optimization(tmp_path, candidates[0], functional="PBE0")
    input_path = tmp_path / "conformer_001" / "conformer_001.inp"
    original = input_path.read_text(encoding="utf-8")
    backend = RecordingBackend()

    with pytest.raises(ConfigurationError, match="settings differ"):
        _run_quantum(
            _request(tmp_path, conformers=1),
            candidates,
            backend,
            tmp_path,
            lambda _message: None,
        )

    assert backend.optimized_indices == []
    assert input_path.read_text(encoding="utf-8") == original


def test_frequency_runs_without_reoptimizing_reused_geometry(tmp_path):
    candidates = _candidates(1)
    _write_completed_optimization(tmp_path, candidates[0])
    backend = RecordingBackend()

    results, _ = _run_quantum(
        _request(tmp_path, conformers=1, frequency=True),
        candidates,
        backend,
        tmp_path,
        lambda _message: None,
    )

    assert backend.optimized_indices == []
    assert backend.frequency_indices == [0]
    assert results[0].frequencies_cm1 == (-25.0, 100.0)
    assert results[0].imaginary_frequencies_cm1 == (-25.0,)
    assert (tmp_path / "conformer_001" / "conformer_001.out").read_text(
        encoding="utf-8"
    ) == OPT_OUTPUT


def test_incomplete_existing_directory_is_never_overwritten(tmp_path):
    candidates = _candidates(1)
    directory = tmp_path / "conformer_001"
    directory.mkdir()
    marker = directory / "conformer_001.out"
    marker.write_text("partial output", encoding="utf-8")
    backend = RecordingBackend()

    with pytest.raises(ConfigurationError, match="incomplete existing calculation"):
        _run_quantum(
            _request(tmp_path, conformers=1),
            candidates,
            backend,
            tmp_path,
            lambda _message: None,
        )

    assert backend.optimized_indices == []
    assert marker.read_text(encoding="utf-8") == "partial output"


def test_existing_frequency_result_is_reused_too(tmp_path):
    candidates = _candidates(1)
    _write_completed_optimization(tmp_path, candidates[0])
    settings = QuantumSettings(
        method="dft",
        functional="B3LYP",
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
        frequency=True,
    )
    directory = tmp_path / "conformer_001"
    (directory / "conformer_001_freq.inp").write_text(
        render_orca_frequency_input(candidates[0].model, settings), encoding="utf-8"
    )
    (directory / "conformer_001_freq.out").write_text(FREQ_OUTPUT, encoding="utf-8")
    backend = RecordingBackend()

    results, _ = _run_quantum(
        _request(tmp_path, conformers=1, frequency=True),
        candidates,
        backend,
        tmp_path,
        lambda _message: None,
    )

    assert backend.optimized_indices == []
    assert backend.frequency_indices == []
    assert results[0].frequencies_cm1 == (-25.0, 100.0)


def test_related_metadata_atom_order_must_match(tmp_path):
    root = tmp_path / "saved_calculations"
    candidates = _candidates(1)
    _write_completed_optimization(root, candidates[0])
    (tmp_path / "wrong.metadata.json").write_text(
        json.dumps(
            {
                "cid": None,
                "atoms": [{"element": "O"}, {"element": "H"}],
                "metadata": {"calculation_files_directory": str(root.resolve())},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="atom order"):
        _run_quantum(
            _request(root, conformers=1),
            candidates,
            RecordingBackend(),
            root,
            lambda _message: None,
        )


def test_conversion_metadata_records_reused_conformer(monkeypatch, tmp_path):
    root = tmp_path / "saved_calculations"
    candidates = _candidates(4)
    _write_completed_optimization(root, candidates[0])
    monkeypatch.setattr(
        "molecule2fbx.pipeline._quantum_initial_candidates",
        lambda request, log: (candidates, "[H][H]"),
    )
    backend = RecordingBackend()

    outcome = run_conversion(
        ConversionRequest(
            smiles="[H][H]",
            method="dft",
            conformers=4,
            output_dir=tmp_path / "output",
            reuse_calculations=root,
        ),
        blender_finder=lambda explicit: "blender",
        orca_finder=lambda explicit: "orca",
        orca_backend_class=lambda executable: backend,
        exporter=fake_exporter,
    )

    search = outcome.selected_model.metadata["conformer_search"]
    assert search["converged_conformers"] == 4
    assert search["reused_conformers"] == [0]
    assert outcome.selected_model.metadata["calculation_files_directory"] == str(
        root.resolve()
    )
    assert backend.optimized_indices == [1, 2, 3]


def test_pool_selection_keeps_four_existing_jobs_before_adding_candidates(tmp_path):
    candidates = _candidates(10)
    for index in range(4):
        (tmp_path / f"conformer_{index + 1:03d}").mkdir()
    request = _request(tmp_path, conformers=6, conformer_pool=10)

    selected = _select_quantum_candidates(
        request, candidates, lambda _message: None
    )

    assert [
        candidate.model.metadata["conformer_pool_index"]
        for candidate in selected[:4]
    ] == [0, 1, 2, 3]
    assert [candidate.conformer_index for candidate in selected] == list(range(6))


def test_selective_frequency_runs_only_inside_window(tmp_path):
    candidates = _candidates(3)
    results = []
    for candidate in candidates:
        directory = tmp_path / f"conformer_{candidate.conformer_index + 1:03d}"
        directory.mkdir()
        results.append(
            _optimized_result(
                candidate.model,
                QuantumSettings("dft", "B3LYP", "def2-SVP", 0, 1),
                directory,
                candidate.conformer_index,
            )
        )
    backend = RecordingBackend()
    request = _request(
        tmp_path,
        conformers=3,
        frequency=True,
        frequency_window_kj=5.0,
        frequency_max=3,
    )

    updated, selection = _run_selective_frequencies(
        request, results, backend, tmp_path, lambda _message: None
    )

    assert backend.frequency_indices == [2, 1]
    assert selection["selected_conformer_indices"] == [1, 2]
    assert selection["skipped_conformer_indices"] == [0]
    by_index = {result.conformer_index: result for result in updated}
    assert by_index[0].frequency_requested is False
    assert by_index[1].frequency_requested is True
    assert by_index[2].frequency_requested is True


def test_conversion_selective_frequency_separates_opt_and_freq(monkeypatch, tmp_path):
    candidates = _candidates(3)
    monkeypatch.setattr(
        "molecule2fbx.pipeline._quantum_initial_candidates",
        lambda request, log: (candidates, "[H][H]"),
    )
    backend = RecordingBackend()

    outcome = run_conversion(
        ConversionRequest(
            smiles="[H][H]",
            method="dft",
            conformers=3,
            frequency=True,
            frequency_window_kj=10.0,
            frequency_max=1,
            output_dir=tmp_path / "output",
        ),
        blender_finder=lambda explicit: "blender",
        orca_finder=lambda explicit: "orca",
        orca_backend_class=lambda executable: backend,
        exporter=fake_exporter,
    )

    assert backend.optimized_indices == [0, 1, 2]
    assert backend.optimize_frequency_flags == [False, False, False]
    assert backend.frequency_indices == [2]
    selection = outcome.selected_model.metadata["conformer_search"][
        "frequency_selection"
    ]
    assert selection["selected_conformer_indices"] == [2]
    assert selection["skipped_conformer_indices"] == [0, 1]
    assert outcome.selected_model.metadata["frequency_requested"] is True

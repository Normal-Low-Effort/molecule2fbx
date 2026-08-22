from pathlib import Path

from molecule2fbx.config import ConversionRequest
from molecule2fbx.pipeline import run_conversion
from molecule2fbx.quantum.base import QuantumResult
from molecule2fbx.structures import ConformerCandidate

from conftest import fake_exporter, hydrogen_model


def test_existing_cid_pubchem_to_fbx(monkeypatch, tmp_path):
    pubchem = hydrogen_model("PubChemHydrogen").with_metadata(
        structure_origin="pubchem_3d"
    )
    pubchem = type(pubchem)(7, pubchem.name, pubchem.atoms, pubchem.bonds, pubchem.metadata)
    monkeypatch.setattr("molecule2fbx.pipeline.download_and_parse", lambda cid, timeout: pubchem)
    outcome = run_conversion(
        ConversionRequest(cid=7, output_dir=tmp_path),
        blender_finder=lambda explicit: "blender",
        exporter=fake_exporter,
    )
    assert outcome.artifacts[0].fbx_path.name == "PubChemHydrogen.fbx"
    assert outcome.artifacts[0].metadata_path.is_file()
    assert outcome.selected_model.metadata["structure_origin"] == "pubchem_3d"


def test_smiles_forcefield_to_fbx(monkeypatch, tmp_path):
    candidate = ConformerCandidate(
        hydrogen_model("Hydrogen").with_metadata(
            structure_origin="computed_forcefield", conformer_index=0
        ),
        0,
        0.1,
        "MMFF94s",
        True,
    )
    monkeypatch.setattr(
        "molecule2fbx.pipeline._forcefield_candidates",
        lambda request, count=None, index_offset=0: ([candidate], "[H][H]"),
    )
    outcome = run_conversion(
        ConversionRequest(smiles="[H][H]", method="forcefield", output_dir=tmp_path),
        blender_finder=lambda explicit: "blender",
        exporter=fake_exporter,
    )
    assert outcome.artifacts[0].fbx_path.name == "Hydrogen_forcefield.fbx"
    assert outcome.selected_model.metadata["structure_origin"] == "computed_forcefield"


def test_smiles_quantum_to_optimized_fbx(monkeypatch, tmp_path):
    candidate = ConformerCandidate(
        hydrogen_model("Hydrogen").with_metadata(
            structure_origin="computed_forcefield", conformer_index=0
        ),
        0,
        0.1,
        "MMFF94s",
        True,
    )
    monkeypatch.setattr(
        "molecule2fbx.pipeline._quantum_initial_candidates",
        lambda request, log: ([candidate], "[H][H]"),
    )

    class FakeOrcaBackend:
        def __init__(self, executable):
            self.executable = executable

        def optimize(self, model, settings, work_dir, conformer_index):
            optimized = model.with_coordinates(
                ((0.0, 0.0, 0.0), (0.73, 0.0, 0.0)),
                structure_origin="computed_quantum",
                final_energy_hartree=-1.12,
                conformer_index=conformer_index,
            )
            return QuantumResult(
                model=optimized,
                backend="ORCA",
                backend_version="test",
                method="dft",
                functional="B3LYP",
                basis="def2-SVP",
                charge=0,
                multiplicity=1,
                final_energy_hartree=-1.12,
                geometry_converged=True,
                scf_converged=True,
                conformer_index=conformer_index,
                frequency_requested=False,
            )

    outcome = run_conversion(
        ConversionRequest(smiles="[H][H]", method="dft", output_dir=tmp_path),
        blender_finder=lambda explicit: "blender",
        orca_finder=lambda explicit: "orca",
        orca_backend_class=FakeOrcaBackend,
        exporter=fake_exporter,
    )
    assert outcome.artifacts[0].fbx_path.name == "Hydrogen_dft.fbx"
    assert outcome.selected_model.metadata["structure_origin"] == "computed_quantum"
    assert outcome.selected_model.atoms[1].x == 0.73

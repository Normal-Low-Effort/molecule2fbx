import pytest

pytest.importorskip("rdkit")

from molecule2fbx.errors import InvalidSMILESError
from molecule2fbx.model import Atom, MoleculeModel
from molecule2fbx.pubchem import parse_sdf
from molecule2fbx.structures import (
    ConformerCandidate,
    aligned_heavy_atom_rmsd,
    generate_forcefield_conformers,
    inspect_smiles,
    select_diverse_conformers,
)


def test_invalid_smiles():
    with pytest.raises(InvalidSMILESError):
        inspect_smiles("this is not smiles")


def test_small_smiles_forcefield_conformer():
    candidates, inspection = generate_forcefield_conformers("CCO", 1, name="Ethanol")
    assert len(candidates) == 1
    assert candidates[0].model.metadata["structure_origin"] == "computed_forcefield"
    assert len(candidates[0].model.atoms) == 9
    assert inspection.canonical_isomeric_smiles == "CCO"


def test_unspecified_stereochemistry_is_recorded_without_claiming_assignment():
    inspection = inspect_smiles("CC(O)F")
    assert inspection.unspecified_stereocenters == (1,)
    assert "do not establish" in inspection.to_metadata()["stereochemistry"]["interpretation"]


def test_explicit_ez_stereochemistry_is_preserved_in_canonical_smiles():
    inspection = inspect_smiles("F/C=C/F")
    assert inspection.canonical_isomeric_smiles == "F/C=C/F"
    assert inspection.specified_stereo_count == 1


def test_numeric_pubchem_sdf_title_falls_back_to_compound_title(monkeypatch):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("O"))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    mol.SetProp("_Name", "123")
    sdf = Chem.MolToMolBlock(mol)
    monkeypatch.setattr("molecule2fbx.pubchem._fallback_title", lambda cid, timeout: "Water")
    model = parse_sdf(sdf, cid=123)
    assert model.name == "Water"


def _four_atom_model(coordinates):
    return MoleculeModel(
        None,
        "test",
        tuple(
            Atom(index, "C", *coordinate)
            for index, coordinate in enumerate(coordinates)
        ),
        (),
    )


def test_aligned_heavy_atom_rmsd_removes_rotation_and_translation():
    first = _four_atom_model(
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (2.0, 1.0, 1.0))
    )
    rotated = _four_atom_model(
        ((5.0, -2.0, 3.0), (5.0, -1.0, 3.0), (4.0, -1.0, 3.0), (4.0, 0.0, 4.0))
    )

    assert aligned_heavy_atom_rmsd(first, rotated) == pytest.approx(0.0, abs=1e-12)


def test_diverse_selection_keeps_required_reuse_prefix_and_reindexes():
    candidates = []
    for index, height in enumerate((1.0, 1.2, 2.0, 3.0, 4.0)):
        model = _four_atom_model(
            (
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0),
                (2.0, 1.0, height),
            )
        ).with_metadata(conformer_index=index)
        candidates.append(ConformerCandidate(model, index, float(index), "MMFF94s", True))

    selected = select_diverse_conformers(
        candidates, 3, 0.25, required_indices=(0, 1)
    )

    assert [candidate.conformer_index for candidate in selected] == [0, 1, 2]
    assert [candidate.model.metadata["conformer_pool_index"] for candidate in selected[:2]] == [
        0,
        1,
    ]
    assert selected[2].model.metadata["conformer_pool_size"] == 5

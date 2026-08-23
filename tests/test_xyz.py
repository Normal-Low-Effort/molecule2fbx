import json

import pytest

from molecule2fbx.config import XYZConversionRequest
from molecule2fbx.errors import RDKitError
from molecule2fbx.pipeline import run_xyz_conversion
from molecule2fbx.xyz import load_xyz_for_export, read_xyz_geometry

from conftest import fake_exporter


def _write_water(path):
    path.write_text(
        "3\noptimized water\nO 0.0000 0.0000 0.0000\n"
        "H 0.9572 0.0000 0.0000\nH -0.2390 0.9270 0.0000\n",
        encoding="utf-8",
    )


def test_load_xyz_preserves_coordinates_and_infers_bonds(tmp_path):
    xyz = tmp_path / "water.xyz"
    _write_water(xyz)
    model = load_xyz_for_export(xyz)
    assert model.name == "water"
    assert [atom.element for atom in model.atoms] == ["O", "H", "H"]
    assert model.atoms[1].x == pytest.approx(0.9572)
    assert {(bond.begin, bond.end, bond.order) for bond in model.bonds} == {
        (0, 1, 1),
        (0, 2, 1),
    }
    assert model.metadata["structure_origin"] == "imported_xyz"
    assert model.metadata["xyz_bond_inference"]["assumed_total_charge"] == 0


def test_xyz_rejects_multiple_frames(tmp_path):
    xyz = tmp_path / "trajectory.xyz"
    xyz.write_text("1\nfirst\nH 0 0 0\n1\nsecond\nH 1 0 0\n", encoding="utf-8")
    with pytest.raises(RDKitError, match="more than one frame"):
        read_xyz_geometry(xyz)


def test_xyz_rejects_nonfinite_coordinates(tmp_path):
    xyz = tmp_path / "invalid.xyz"
    xyz.write_text("1\ninvalid\nH nan 0 0\n", encoding="utf-8")
    with pytest.raises(RDKitError, match="Invalid XYZ coordinate"):
        read_xyz_geometry(xyz)


def test_xyz_pipeline_exports_without_recalculation(tmp_path):
    source_dir = tmp_path / "snapshot"
    source_dir.mkdir()
    xyz = source_dir / "best.xyz"
    _write_water(xyz)
    output_dir = tmp_path / "fbx"
    outcome = run_xyz_conversion(
        XYZConversionRequest(xyz_path=xyz, output_dir=output_dir, name="Best Water"),
        blender_finder=lambda explicit: "blender",
        exporter=fake_exporter,
    )
    artifact = outcome.artifacts[0]
    assert artifact.fbx_path.name == "Best Water.fbx"
    assert artifact.metadata_path.name == "Best Water.fbx.metadata.json"
    assert artifact.fbx_path.is_file()
    metadata = json.loads(artifact.metadata_path.read_text(encoding="utf-8"))
    assert metadata["metadata"]["structure_origin"] == "imported_xyz"
    assert metadata["metadata"]["xyz_atom_order_preserved"] is True
    assert len(metadata["bonds"]) == 2

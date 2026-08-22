import pytest

from molecule2fbx.blender_export import find_blender
from molecule2fbx.errors import BlenderNotFoundError


def test_blender_not_found(monkeypatch, tmp_path):
    monkeypatch.delenv("BLENDER_EXECUTABLE", raising=False)
    monkeypatch.setattr("molecule2fbx.blender_export.shutil.which", lambda candidate: None)
    with pytest.raises(BlenderNotFoundError):
        find_blender(str(tmp_path / "missing-blender.exe"))

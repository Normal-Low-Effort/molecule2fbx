from pathlib import Path

from molecule2fbx.model import Atom, Bond, MoleculeModel


def hydrogen_model(name: str = "Hydrogen") -> MoleculeModel:
    return MoleculeModel(
        cid=None,
        name=name,
        atoms=(
            Atom(0, "H", 0.0, 0.0, 0.0),
            Atom(1, "H", 0.74, 0.0, 0.0),
        ),
        bonds=(Bond(0, 1, 1),),
        metadata={"formal_charge": 0},
    )


def fake_exporter(model, output_path: Path, blender_executable: str, timeout: float):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(b"Kaydara FBX Binary mock")
    return output_path.resolve()

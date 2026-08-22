from molecule2fbx.model import (
    Atom,
    Bond,
    MoleculeModel,
    normalize_bond_order,
    safe_filename,
    validate_model,
)


def test_safe_filename_removes_windows_unsafe_characters():
    assert safe_filename('A:B/C*D?', 'fallback') == "A_B_C_D_"
    assert safe_filename("   ", "fallback") == "fallback"
    assert safe_filename("CON", "fallback") == "_CON"


def test_bond_order_mapping():
    assert normalize_bond_order(1.0) == 1
    assert normalize_bond_order(1.5, is_aromatic=True) == 1
    assert normalize_bond_order(2.0) == 2
    assert normalize_bond_order(3.0) == 3


def test_model_serializes_atom_and_bond_data():
    model = MoleculeModel(
        cid=1,
        name="Hydrogen",
        atoms=(
            Atom(0, "H", 0.0, 0.0, 0.0),
            Atom(1, "H", 0.7, 0.0, 0.0),
        ),
        bonds=(Bond(0, 1, 1),),
    )
    validate_model(model)
    serialized = model.to_dict()
    assert serialized["name"] == "Hydrogen"
    assert serialized["atoms"][1]["position"] == [0.7, 0.0, 0.0]
    assert serialized["bonds"] == [{"begin": 0, "end": 1, "order": 1}]

"""Read a single XYZ geometry and infer display bonds with RDKit."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Optional, Tuple

from .errors import RDKitError
from .model import Atom, Bond, MoleculeModel, normalize_bond_order


def _rdkit_modules():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RDKitError(
            "The RDKit package is required to infer bonds from XYZ coordinates"
        ) from exc
    return Chem, rdDetermineBonds


def read_xyz_geometry(path: Path, *, name: Optional[str] = None) -> MoleculeModel:
    """Read exactly one XYZ frame while preserving atom order and coordinates."""

    path = Path(path)
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        atom_count = int(lines[0].strip())
    except (OSError, IndexError, ValueError) as exc:
        raise RDKitError(f"Could not read XYZ: {path}") from exc
    if atom_count < 1 or len(lines) < atom_count + 2:
        raise RDKitError("XYZ has an invalid atom count")
    if any(line.strip() for line in lines[atom_count + 2 :]):
        raise RDKitError("XYZ contains more than one frame or trailing atom data")

    Chem, _ = _rdkit_modules()
    periodic_table = Chem.GetPeriodicTable()
    atoms = []
    for index, line in enumerate(lines[2 : atom_count + 2]):
        fields = line.split()
        if len(fields) < 4:
            raise RDKitError(f"Invalid XYZ atom line at index {index}")
        try:
            atomic_number = int(periodic_table.GetAtomicNumber(fields[0]))
        except (RuntimeError, ValueError) as exc:
            raise RDKitError(
                f"Unknown element symbol in XYZ at atom {index}: {fields[0]}"
            ) from exc
        if atomic_number < 1:
            raise RDKitError(
                f"Unknown element symbol in XYZ at atom {index}: {fields[0]}"
            )
        try:
            coordinates = tuple(float(value) for value in fields[1:4])
        except ValueError as exc:
            raise RDKitError(f"Invalid XYZ coordinate at atom {index}") from exc
        if not all(math.isfinite(value) for value in coordinates):
            raise RDKitError(f"Invalid XYZ coordinate at atom {index}")
        atoms.append(
            Atom(
                index,
                periodic_table.GetElementSymbol(atomic_number),
                coordinates[0],
                coordinates[1],
                coordinates[2],
            )
        )

    return MoleculeModel(
        cid=None,
        name=name or path.stem,
        atoms=tuple(atoms),
        bonds=(),
        metadata={
            "structure_origin": "imported_xyz",
            "structure_source": str(path.expanduser().resolve()),
            "structure_claim": (
                "Coordinates imported from XYZ without recalculation; their experimental "
                "or computational provenance was not independently verified."
            ),
            "xyz_atom_order_preserved": True,
        },
    )


def _xyz_block(model: MoleculeModel) -> str:
    lines = [str(len(model.atoms)), model.name]
    lines.extend(
        f"{atom.element} {atom.x:.16g} {atom.y:.16g} {atom.z:.16g}"
        for atom in model.atoms
    )
    return "\n".join(lines) + "\n"


def infer_xyz_bonds(
    model: MoleculeModel,
    *,
    charge: int = 0,
) -> Tuple[Tuple[Bond, ...], Dict[str, object]]:
    """Infer connectivity and, when possible, bond orders for FBX rendering."""

    Chem, rdDetermineBonds = _rdkit_modules()
    molecule = Chem.MolFromXYZBlock(_xyz_block(model))
    if molecule is None:
        raise RDKitError("RDKit could not parse the normalized XYZ geometry")

    status = "bond_orders_assigned"
    warning = None
    try:
        rdDetermineBonds.DetermineBonds(molecule, charge=int(charge))
    except Exception as exc:  # RDKit exception types vary by release
        molecule = Chem.MolFromXYZBlock(_xyz_block(model))
        if molecule is None:  # pragma: no cover - guarded by the first parse
            raise RDKitError("RDKit could not parse the normalized XYZ geometry") from exc
        try:
            rdDetermineBonds.DetermineConnectivity(molecule, charge=int(charge))
        except Exception as connectivity_exc:
            raise RDKitError(
                "RDKit could not infer molecular connectivity from the XYZ geometry"
            ) from connectivity_exc
        status = "connectivity_only"
        warning = (
            "RDKit could infer connectivity but not chemically consistent bond orders; "
            "all displayed bonds were rendered as single bonds."
        )

    bonds = tuple(
        Bond(
            min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            (
                1
                if status == "connectivity_only"
                else normalize_bond_order(
                    float(bond.GetBondTypeAsDouble()), bond.GetIsAromatic()
                )
            ),
        )
        for bond in molecule.GetBonds()
    )
    metadata: Dict[str, object] = {
        "xyz_bond_inference": {
            "software": "RDKit rdDetermineBonds",
            "status": status,
            "assumed_total_charge": int(charge),
            "bond_count": len(bonds),
            "warning": warning,
        }
    }
    return bonds, metadata


def load_xyz_for_export(
    path: Path,
    *,
    name: Optional[str] = None,
    charge: int = 0,
) -> MoleculeModel:
    """Build an FBX-ready model from one XYZ file without altering coordinates."""

    model = read_xyz_geometry(path, name=name)
    bonds, metadata = infer_xyz_bonds(model, charge=charge)
    merged = dict(model.metadata)
    merged.update(metadata)
    return MoleculeModel(model.cid, model.name, model.atoms, bonds, merged)

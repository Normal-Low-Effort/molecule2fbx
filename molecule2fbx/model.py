"""Data structures and small pure helpers for molecular geometry."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from dataclasses import field
from typing import Dict, Mapping, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Atom:
    """An atom and its 3D coordinate in the SDF coordinate system."""

    index: int
    element: str
    x: float
    y: float
    z: float

    @property
    def position(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass(frozen=True)
class Bond:
    """A bond connecting two atom indices."""

    begin: int
    end: int
    order: int


@dataclass(frozen=True)
class MoleculeModel:
    """Serializable molecule data passed from Python to Blender."""

    cid: Optional[int]
    name: str
    atoms: Tuple[Atom, ...]
    bonds: Tuple[Bond, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        return {
            "cid": self.cid,
            "name": self.name,
            "metadata": dict(self.metadata),
            "atoms": [
                {
                    "index": atom.index,
                    "element": atom.element,
                    "position": list(atom.position),
                }
                for atom in self.atoms
            ],
            "bonds": [
                {
                    "begin": bond.begin,
                    "end": bond.end,
                    "order": bond.order,
                }
                for bond in self.bonds
            ],
        }

    def with_metadata(self, **updates: object) -> "MoleculeModel":
        metadata = dict(self.metadata)
        metadata.update(updates)
        return MoleculeModel(
            cid=self.cid,
            name=self.name,
            atoms=self.atoms,
            bonds=self.bonds,
            metadata=metadata,
        )

    def with_coordinates(
        self,
        coordinates: Sequence[Tuple[float, float, float]],
        **metadata_updates: object,
    ) -> "MoleculeModel":
        if len(coordinates) != len(self.atoms):
            raise ValueError("Optimized coordinate count does not match the molecule")
        atoms = tuple(
            Atom(atom.index, atom.element, float(x), float(y), float(z))
            for atom, (x, y, z) in zip(self.atoms, coordinates)
        )
        metadata = dict(self.metadata)
        metadata.update(metadata_updates)
        return MoleculeModel(self.cid, self.name, atoms, self.bonds, metadata)


def normalize_bond_order(order: float, is_aromatic: bool = False) -> int:
    """Map an RDKit bond order to the three visual styles supported by the tool."""

    # Aromatic bonds are represented by one cylinder.  This avoids drawing
    # every aromatic ring as a double bond while retaining a useful model.
    if is_aromatic:
        return 1
    if not math.isfinite(order):
        return 1
    if order <= 1.25:
        return 1
    if order <= 2.5:
        return 2
    return 3


_INVALID_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
_WINDOWS_RESERVED_NAMES = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
}


def safe_filename(name: str, fallback: str) -> str:
    """Return a filesystem-safe filename stem without changing its extension."""

    cleaned = _INVALID_FILENAME_CHARS.sub("_", (name or "").strip())
    cleaned = cleaned.rstrip(" .")
    if not cleaned:
        cleaned = fallback
    if cleaned.upper() in _WINDOWS_RESERVED_NAMES:
        cleaned = f"_{cleaned}"
    return cleaned


def validate_model(model: MoleculeModel) -> None:
    """Validate data before handing it to Blender."""

    if not model.atoms:
        raise ValueError("The molecule contains no atoms")
    atom_indices = {atom.index for atom in model.atoms}
    if atom_indices != set(range(len(model.atoms))):
        raise ValueError("Atom indices must be contiguous and zero-based")
    for atom in model.atoms:
        if not all(math.isfinite(value) for value in atom.position):
            raise ValueError(f"Atom {atom.index} has invalid coordinates")
    for bond in model.bonds:
        if bond.begin not in atom_indices or bond.end not in atom_indices:
            raise ValueError("Bond references an unknown atom")
        if bond.begin == bond.end:
            raise ValueError("A bond cannot connect an atom to itself")
        if bond.order not in (1, 2, 3):
            raise ValueError("Bond order must be 1, 2, or 3")

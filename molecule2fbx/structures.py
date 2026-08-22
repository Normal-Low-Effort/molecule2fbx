"""SMILES inspection, RDKit conformer generation, and geometry checks."""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .errors import ForceFieldError, InvalidSMILESError, RDKitError
from .model import Atom, Bond, MoleculeModel, normalize_bond_order


@dataclass(frozen=True)
class SmilesInspection:
    original_smiles: str
    canonical_isomeric_smiles: str
    formal_charge: int
    radical_electrons: int
    unspecified_stereocenters: Tuple[int, ...]
    unspecified_double_bonds: Tuple[int, ...]
    specified_stereo_count: int
    specified_stereocenters: Tuple[Tuple[int, str], ...]

    @property
    def has_unspecified_stereo(self) -> bool:
        return bool(self.unspecified_stereocenters or self.unspecified_double_bonds)

    def to_metadata(self) -> Dict[str, object]:
        return {
            "original_smiles": self.original_smiles,
            "canonical_isomeric_smiles": self.canonical_isomeric_smiles,
            "formal_charge": self.formal_charge,
            "radical_electrons": self.radical_electrons,
            "stereochemistry": {
                "specified_count": self.specified_stereo_count,
                "specified_stereocenters": [
                    {"atom_index": index, "cip": label}
                    for index, label in self.specified_stereocenters
                ],
                "unspecified_stereocenters": list(self.unspecified_stereocenters),
                "unspecified_double_bonds": list(self.unspecified_double_bonds),
                "interpretation": (
                    "Generated coordinates are one possible conformer and do not establish "
                    "unspecified stereochemistry."
                    if self.has_unspecified_stereo
                    else "Explicit SMILES stereochemistry was retained during conformer generation."
                ),
            },
        }


@dataclass(frozen=True)
class ConformerCandidate:
    model: MoleculeModel
    conformer_index: int
    energy_kcal_mol: Optional[float]
    forcefield: str
    converged: bool
    optimization_iterations: Optional[int] = None
    failure_status: Optional[str] = None


@dataclass(frozen=True)
class StereochemistryValidation:
    expected_centers: Tuple[Tuple[int, str], ...]
    observed_centers: Tuple[Tuple[int, str], ...]
    atom_order_matches: bool

    @property
    def matches(self) -> bool:
        return self.atom_order_matches and self.expected_centers == self.observed_centers

    def to_metadata(self) -> Dict[str, object]:
        return {
            "matches": self.matches,
            "atom_order_matches": self.atom_order_matches,
            "expected_centers": [
                {"atom_index": index, "cip": label}
                for index, label in self.expected_centers
            ],
            "observed_centers": [
                {"atom_index": index, "cip": label}
                for index, label in self.observed_centers
            ],
        }


@dataclass(frozen=True)
class RMSDAtomSubsets:
    """Fixed atom-order subsets used for complementary conformer RMSDs.

    The subsets deliberately do not replace the historical all-heavy-atom
    metric.  They add views that ignore a terminal trimethylsilyl group and
    that focus on the N-benzoyl reaction centre.
    """

    all_heavy: Tuple[int, ...]
    common_scaffold: Tuple[int, ...]
    reaction_center: Tuple[int, ...]
    excluded_terminal_substituent: Tuple[int, ...]
    benzoyl_carbonyl_c: Optional[int]
    benzoyl_carbonyl_o: Optional[int]
    benzoyl_amide_n: Optional[int]
    benzoyl_ipso_c: Optional[int]

    def to_metadata(self) -> Dict[str, object]:
        return {
            "all_heavy_atom_indices": list(self.all_heavy),
            "common_scaffold_atom_indices": list(self.common_scaffold),
            "reaction_center_atom_indices": list(self.reaction_center),
            "excluded_terminal_substituent_atom_indices": list(
                self.excluded_terminal_substituent
            ),
            "benzoyl_reaction_center": {
                "carbonyl_c": self.benzoyl_carbonyl_c,
                "carbonyl_o": self.benzoyl_carbonyl_o,
                "amide_n": self.benzoyl_amide_n,
                "benzoyl_ipso_c": self.benzoyl_ipso_c,
            },
            "atom_correspondence": "fixed_input_atom_order",
            "symmetry_permutations": False,
            "assumption": (
                "Equivalent atoms are not permuted. The common-scaffold view excludes "
                "a terminal aryl-Si(CH3)3 branch so methyl rotation does not create a "
                "separate scaffold cluster."
            ),
        }


@dataclass(frozen=True)
class CrossMoleculeAtomMapping:
    """Deterministic Bz/TMS-Bz scaffold correspondence."""

    first_indices: Tuple[int, ...]
    second_indices: Tuple[int, ...]
    first_reaction_indices: Tuple[int, ...]
    second_reaction_indices: Tuple[int, ...]

    def to_metadata(self) -> Dict[str, object]:
        return {
            "first_atom_indices": list(self.first_indices),
            "second_atom_indices": list(self.second_indices),
            "first_reaction_center_atom_indices": list(
                self.first_reaction_indices
            ),
            "second_reaction_center_atom_indices": list(
                self.second_reaction_indices
            ),
            "method": "terminal_TMS_removal_then_graph_isomorphism",
            "symmetry_permutations_for_rmsd_minimization": False,
            "isomorphism_choice": (
                "benzoyl reaction-centre anchors followed by minimum retained "
                "input-order displacement"
            ),
        }


def _rdkit_modules():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RDKitError("The RDKit package is required; install the project dependencies first") from exc
    return Chem, AllChem


def _parse_smiles(smiles: str):
    Chem, _ = _rdkit_modules()
    if not isinstance(smiles, str) or not smiles.strip():
        raise InvalidSMILESError("SMILES must be a non-empty string")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise InvalidSMILESError("Invalid SMILES")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem, mol


def inspect_smiles(smiles: str) -> SmilesInspection:
    Chem, mol = _parse_smiles(smiles)
    chiral = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    unspecified_centers = tuple(index for index, label in chiral if label == "?")
    specified_count = sum(1 for _, label in chiral if label != "?")
    unspecified_bonds: List[int] = []
    try:
        for info in Chem.FindPotentialStereo(mol):
            info_type = str(info.type)
            specified = str(info.specified)
            if "Bond_Double" in info_type:
                if "Unspecified" in specified:
                    unspecified_bonds.append(int(info.centeredOn))
                else:
                    specified_count += 1
    except (AttributeError, RuntimeError):
        # Older RDKit versions may not expose FindPotentialStereo. Explicit
        # atom stereochemistry is still preserved by ETKDG in those versions.
        pass
    return SmilesInspection(
        original_smiles=smiles,
        canonical_isomeric_smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
        formal_charge=int(Chem.GetFormalCharge(mol)),
        radical_electrons=int(sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())),
        unspecified_stereocenters=unspecified_centers,
        unspecified_double_bonds=tuple(unspecified_bonds),
        specified_stereo_count=specified_count,
        specified_stereocenters=tuple(
            (index, label) for index, label in chiral if label != "?"
        ),
    )


def validate_model_stereochemistry(
    smiles: str, model: MoleculeModel
) -> StereochemistryValidation:
    """Reassign CIP labels from 3D coordinates and compare explicit SMILES centers."""

    Chem, _ = _rdkit_modules()
    _, base_mol = _parse_smiles(smiles)
    expected = tuple(
        (index, label)
        for index, label in Chem.FindMolChiralCenters(
            base_mol, includeUnassigned=True, useLegacyImplementation=False
        )
        if label != "?"
    )
    mol = Chem.AddHs(base_mol)
    expected_elements = [atom.GetSymbol().casefold() for atom in mol.GetAtoms()]
    actual_elements = [atom.element.casefold() for atom in model.atoms]
    if expected_elements != actual_elements:
        return StereochemistryValidation(expected, (), False)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for atom in model.atoms:
        conformer.SetAtomPosition(atom.index, atom.position)
    mol.RemoveAllConformers()
    mol.AddConformer(conformer, assignId=True)
    Chem.RemoveStereochemistry(mol)
    Chem.AssignAtomChiralTagsFromStructure(
        mol, confId=0, replaceExistingTags=True
    )
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    observed_all = dict(
        Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False
        )
    )
    observed = tuple((index, observed_all.get(index, "?")) for index, _ in expected)
    return StereochemistryValidation(expected, observed, True)


def default_smiles_name(smiles: str) -> str:
    digest = hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:10]
    return f"SMILES_{digest}"


def _model_from_rdkit(
    mol,
    conformer_id: int,
    *,
    cid: Optional[int],
    name: str,
    metadata: Dict[str, object],
) -> MoleculeModel:
    conformer = mol.GetConformer(conformer_id)
    atoms = []
    for atom in mol.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        atoms.append(
            Atom(atom.GetIdx(), atom.GetSymbol(), float(position.x), float(position.y), float(position.z))
        )
    bonds = tuple(
        Bond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            normalize_bond_order(float(bond.GetBondTypeAsDouble()), bond.GetIsAromatic()),
        )
        for bond in mol.GetBonds()
    )
    return MoleculeModel(cid, name, tuple(atoms), bonds, metadata)


def generate_forcefield_conformers(
    smiles: str,
    count: int,
    *,
    cid: Optional[int] = None,
    name: Optional[str] = None,
    index_offset: int = 0,
    max_iterations: int = 500,
    random_seed: int = 0xF00D,
    prune_rms_threshold: float = -1.0,
) -> Tuple[List[ConformerCandidate], SmilesInspection]:
    """Generate ETKDG conformers and pre-optimize each with MMFF or UFF."""

    if count < 1:
        raise ValueError("Conformer count must be positive")
    Chem, AllChem = _rdkit_modules()
    inspection = inspect_smiles(smiles)
    _, base_mol = _parse_smiles(smiles)
    mol = Chem.AddHs(base_mol)
    try:
        params = AllChem.ETKDGv3()
    except AttributeError:  # pragma: no cover - legacy RDKit
        params = AllChem.ETKDGv2()
    params.randomSeed = random_seed
    params.enforceChirality = True
    params.pruneRmsThresh = prune_rms_threshold
    params.numThreads = 1
    conformer_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=count, params=params))
    if not conformer_ids:
        raise ForceFieldError(
            f"RDKit generated 0 of {count} requested conformers"
        )

    if AllChem.MMFFHasAllMoleculeParams(mol):
        forcefield = "MMFF94s"
        optimization = list(
            AllChem.MMFFOptimizeMoleculeConfs(
                mol,
                numThreads=1,
                maxIters=max_iterations,
                mmffVariant="MMFF94s",
            )
        )
    elif AllChem.UFFHasAllMoleculeParams(mol):
        forcefield = "UFF"
        optimization = list(
            AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=max_iterations)
        )
    else:
        raise ForceFieldError("Neither MMFF nor UFF has parameters for all atoms in this molecule")

    molecule_name = name or default_smiles_name(inspection.canonical_isomeric_smiles)
    candidates = []
    stereo_metadata = inspection.to_metadata()
    try:
        from rdkit import rdBase

        rdkit_version = rdBase.rdkitVersion
    except (ImportError, AttributeError):  # pragma: no cover - legacy RDKit
        rdkit_version = None
    for local_index, (conformer_id, result) in enumerate(zip(conformer_ids, optimization)):
        status, energy = int(result[0]), float(result[1])
        conformer_index = index_offset + local_index
        metadata = dict(stereo_metadata)
        metadata.update(
            {
                "structure_origin": "computed_forcefield",
                "structure_source": f"RDKit ETKDG + {forcefield}",
                "forcefield": forcefield,
                "forcefield_energy_kcal_mol": energy,
                "forcefield_converged": status == 0,
                "forcefield_optimization_status": status,
                "forcefield_optimization_iterations": None,
                "forcefield_iteration_count_available": False,
                "forcefield_max_iterations": max_iterations,
                "conformer_index": conformer_index,
                "etkdg": {
                    "requested_conformers": count,
                    "generated_conformers": len(conformer_ids),
                    "embedding_failures": count - len(conformer_ids),
                    "random_seed": random_seed,
                    "prune_rms_threshold_angstrom": prune_rms_threshold,
                    "enforce_chirality": True,
                    "rdkit_version": rdkit_version,
                },
                "structure_claim": "Computed conformer; not an experimentally determined structure",
            }
        )
        model = _model_from_rdkit(
            mol,
            conformer_id,
            cid=cid,
            name=molecule_name,
            metadata=metadata,
        )
        candidates.append(
            ConformerCandidate(
                model,
                conformer_index,
                energy,
                forcefield,
                status == 0,
                None,
                None if status == 0 else "forcefield_optimization_failure",
            )
        )
    return candidates, inspection


def aligned_atom_subset_rmsd(
    first: MoleculeModel,
    second: MoleculeModel,
    atom_indices: Sequence[int],
) -> float:
    """Return fixed-order Kabsch RMSD for an explicit atom subset."""

    first_elements = [atom.element.casefold() for atom in first.atoms]
    second_elements = [atom.element.casefold() for atom in second.atoms]
    if first_elements != second_elements:
        raise ValueError("Conformer atom order does not match")
    indices = [int(index) for index in atom_indices]
    if len(set(indices)) != len(indices):
        raise ValueError("RMSD atom indices must be unique")
    if any(index < 0 or index >= len(first.atoms) for index in indices):
        raise ValueError("RMSD atom index is outside the molecule")
    if not indices:
        raise ValueError("Cannot calculate RMSD for an empty molecule")
    first_coordinates = np.asarray(
        [first.atoms[index].position for index in indices], dtype=float
    )
    second_coordinates = np.asarray(
        [second.atoms[index].position for index in indices], dtype=float
    )
    first_coordinates -= first_coordinates.mean(axis=0)
    second_coordinates -= second_coordinates.mean(axis=0)
    left, _singular_values, right_transpose = np.linalg.svd(
        first_coordinates.T @ second_coordinates
    )
    handedness = np.linalg.det(left @ right_transpose)
    rotation = left @ np.diag((1.0, 1.0, handedness)) @ right_transpose
    difference = first_coordinates @ rotation - second_coordinates
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def aligned_mapped_atom_rmsd(
    first: MoleculeModel,
    second: MoleculeModel,
    first_indices: Sequence[int],
    second_indices: Sequence[int],
) -> float:
    """Return Kabsch RMSD for an explicit cross-molecule atom mapping."""

    first_selected = [int(index) for index in first_indices]
    second_selected = [int(index) for index in second_indices]
    if len(first_selected) != len(second_selected) or not first_selected:
        raise ValueError("Mapped RMSD atom lists must have the same non-zero length")
    if len(set(first_selected)) != len(first_selected) or len(
        set(second_selected)
    ) != len(second_selected):
        raise ValueError("Mapped RMSD atom indices must be unique")
    if any(index < 0 or index >= len(first.atoms) for index in first_selected):
        raise ValueError("First mapped RMSD atom index is outside the molecule")
    if any(index < 0 or index >= len(second.atoms) for index in second_selected):
        raise ValueError("Second mapped RMSD atom index is outside the molecule")
    for first_index, second_index in zip(first_selected, second_selected):
        if (
            first.atoms[first_index].element.casefold()
            != second.atoms[second_index].element.casefold()
        ):
            raise ValueError("Mapped RMSD atom elements do not match")
    first_coordinates = np.asarray(
        [first.atoms[index].position for index in first_selected], dtype=float
    )
    second_coordinates = np.asarray(
        [second.atoms[index].position for index in second_selected], dtype=float
    )
    first_coordinates -= first_coordinates.mean(axis=0)
    second_coordinates -= second_coordinates.mean(axis=0)
    left, _singular_values, right_transpose = np.linalg.svd(
        first_coordinates.T @ second_coordinates
    )
    handedness = np.linalg.det(left @ right_transpose)
    rotation = left @ np.diag((1.0, 1.0, handedness)) @ right_transpose
    difference = first_coordinates @ rotation - second_coordinates
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def aligned_heavy_atom_rmsd(first: MoleculeModel, second: MoleculeModel) -> float:
    """Return fixed-order Kabsch RMSD over matching heavy atoms."""

    indices = [
        atom.index for atom in first.atoms if atom.element.casefold() != "h"
    ]
    if len(indices) < 3:
        indices = list(range(len(first.atoms)))
    return aligned_atom_subset_rmsd(first, second, indices)


def _analysis_rdkit_molecule(model: MoleculeModel):
    """Recover the explicit-H RDKit graph while preserving model atom order."""

    # RDKit canonicalization may reorder atoms.  The original SMILES is the
    # atom-order contract used by ETKDG, ORCA XYZ, and the retained metadata.
    smiles = model.metadata.get("original_smiles")
    if not isinstance(smiles, str) or not smiles:
        smiles = model.metadata.get("canonical_isomeric_smiles")
    if not isinstance(smiles, str) or not smiles:
        raise ValueError("Canonical or original SMILES is required for RMSD atom subsets")
    Chem, _ = _rdkit_modules()
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    expected = [atom.GetSymbol().casefold() for atom in mol.GetAtoms()]
    observed = [atom.element.casefold() for atom in model.atoms]
    heavy_expected = [
        atom.GetSymbol().casefold() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1
    ]
    if expected != observed and heavy_expected != observed:
        raise ValueError("SMILES atom order does not match the optimized model")
    return mol


def _benzoyl_reaction_center(mol) -> Optional[Tuple[int, int, int, int]]:
    """Find aromatic-N-C(=O)-aryl, distinguishing it from the side-chain amide."""

    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygen = None
        aromatic_nitrogen = None
        aromatic_carbon = None
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            order = float(bond.GetBondTypeAsDouble())
            if neighbor.GetAtomicNum() == 8 and order >= 1.75:
                oxygen = neighbor.GetIdx()
            elif neighbor.GetAtomicNum() == 7 and neighbor.GetIsAromatic():
                aromatic_nitrogen = neighbor.GetIdx()
            elif neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                aromatic_carbon = neighbor.GetIdx()
        if oxygen is not None and aromatic_nitrogen is not None and aromatic_carbon is not None:
            matches.append(
                (atom.GetIdx(), oxygen, aromatic_nitrogen, aromatic_carbon)
            )
    return matches[0] if len(matches) == 1 else None


def _terminal_tms_heavy_atoms(mol) -> Tuple[int, ...]:
    excluded = set()
    for silicon in mol.GetAtoms():
        if silicon.GetAtomicNum() != 14:
            continue
        aromatic_anchors = [
            neighbor
            for neighbor in silicon.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic()
        ]
        carbon_branches = [
            neighbor
            for neighbor in silicon.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic()
        ]
        if len(aromatic_anchors) != 1 or len(carbon_branches) != 3:
            continue
        excluded.add(silicon.GetIdx())
        stack = [atom.GetIdx() for atom in carbon_branches]
        anchor_index = aromatic_anchors[0].GetIdx()
        while stack:
            index = stack.pop()
            if index == anchor_index or index in excluded:
                continue
            atom = mol.GetAtomWithIdx(index)
            if atom.GetAtomicNum() == 1:
                continue
            excluded.add(index)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in excluded and neighbor.GetIdx() != anchor_index:
                    stack.append(neighbor.GetIdx())
    return tuple(sorted(excluded))


def rmsd_atom_subsets(model: MoleculeModel) -> RMSDAtomSubsets:
    """Derive all-heavy, common-scaffold, and local reaction-centre subsets."""

    Chem, _ = _rdkit_modules()
    mol = _analysis_rdkit_molecule(model)
    all_heavy = tuple(
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1
    )
    excluded = _terminal_tms_heavy_atoms(mol)
    common = tuple(index for index in all_heavy if index not in set(excluded))
    center = _benzoyl_reaction_center(mol)
    reaction: Tuple[int, ...] = ()
    if center is not None:
        carbonyl_c, _oxygen, _nitrogen, _ipso = center
        distances = Chem.GetDistanceMatrix(mol)
        reaction = tuple(
            index
            for index in all_heavy
            if distances[carbonyl_c, index] <= 2.0 + 1.0e-9
        )
    return RMSDAtomSubsets(
        all_heavy=all_heavy,
        common_scaffold=common,
        reaction_center=reaction,
        excluded_terminal_substituent=excluded,
        benzoyl_carbonyl_c=center[0] if center else None,
        benzoyl_carbonyl_o=center[1] if center else None,
        benzoyl_amide_n=center[2] if center else None,
        benzoyl_ipso_c=center[3] if center else None,
    )


def common_scaffold_atom_mapping(
    first: MoleculeModel, second: MoleculeModel
) -> CrossMoleculeAtomMapping:
    """Map a Bz scaffold to its TMS-Bz counterpart without reordering XYZ."""

    Chem, _ = _rdkit_modules()
    first_mol = Chem.RemoveHs(_analysis_rdkit_molecule(first))
    second_mol = Chem.RemoveHs(_analysis_rdkit_molecule(second))
    first_subsets = rmsd_atom_subsets(first)
    second_subsets = rmsd_atom_subsets(second)

    def stripped(mol, excluded: Sequence[int]):
        for atom in mol.GetAtoms():
            atom.SetIntProp("_m2f_original_index", atom.GetIdx())
        editable = Chem.RWMol(mol)
        for index in sorted(set(int(value) for value in excluded), reverse=True):
            editable.RemoveAtom(index)
        result = editable.GetMol()
        Chem.SanitizeMol(result)
        return result

    first_core = stripped(first_mol, first_subsets.excluded_terminal_substituent)
    second_core = stripped(second_mol, second_subsets.excluded_terminal_substituent)
    if first_core.GetNumAtoms() != second_core.GetNumAtoms():
        raise ValueError("Molecules do not share the same TMS-excluded heavy scaffold")
    matches = second_core.GetSubstructMatches(
        first_core,
        useChirality=True,
        uniquify=False,
        maxMatches=10000,
    )
    if not matches:
        raise ValueError("No stereochemistry-preserving common-scaffold mapping found")
    first_original = [
        atom.GetIntProp("_m2f_original_index") for atom in first_core.GetAtoms()
    ]
    first_rank = {value: rank for rank, value in enumerate(first_original)}
    second_original_by_core = {
        atom.GetIdx(): atom.GetIntProp("_m2f_original_index")
        for atom in second_core.GetAtoms()
    }
    second_rank = {
        value: rank
        for rank, value in enumerate(
            sorted(second_original_by_core.values())
        )
    }
    anchors = (
        (first_subsets.benzoyl_carbonyl_c, second_subsets.benzoyl_carbonyl_c),
        (first_subsets.benzoyl_carbonyl_o, second_subsets.benzoyl_carbonyl_o),
        (first_subsets.benzoyl_amide_n, second_subsets.benzoyl_amide_n),
        (first_subsets.benzoyl_ipso_c, second_subsets.benzoyl_ipso_c),
    )

    def expanded(match: Sequence[int]) -> Tuple[int, ...]:
        return tuple(second_original_by_core[int(index)] for index in match)

    expanded_matches = [expanded(match) for match in matches]
    anchored = [
        match
        for match in expanded_matches
        if all(
            first_anchor is None
            or second_anchor is None
            or match[first_original.index(int(first_anchor))] == int(second_anchor)
            for first_anchor, second_anchor in anchors
        )
    ]
    choices = anchored or expanded_matches

    def order_score(match: Sequence[int]):
        return (
            sum(
                abs(first_rank[first_index] - second_rank[second_index])
                for first_index, second_index in zip(first_original, match)
            ),
            tuple(match),
        )

    second_selected = min(choices, key=order_score)
    mapping = dict(zip(first_original, second_selected))
    first_reaction = tuple(first_subsets.reaction_center)
    second_reaction = tuple(mapping[index] for index in first_reaction)
    if set(second_reaction) != set(second_subsets.reaction_center):
        raise ValueError("Common-scaffold mapping did not preserve the reaction centre")
    return CrossMoleculeAtomMapping(
        tuple(first_original),
        tuple(second_selected),
        first_reaction,
        second_reaction,
    )


def select_diverse_conformers(
    candidates: Sequence[ConformerCandidate],
    count: int,
    rmsd_threshold: float,
    *,
    required_indices: Sequence[int] = (),
) -> List[ConformerCandidate]:
    """Select low-energy, RMSD-diverse candidates and assign contiguous job indices."""

    if count < 1:
        raise ValueError("Selected conformer count must be positive")
    if rmsd_threshold <= 0:
        raise ValueError("RMSD threshold must be greater than zero")
    converged = [candidate for candidate in candidates if candidate.converged]
    if len(converged) < count:
        raise ForceFieldError(
            f"Only {len(converged)} converged force-field conformers are available; "
            f"{count} are required"
        )
    by_index = {candidate.conformer_index: candidate for candidate in converged}
    required = []
    for index in sorted(set(required_indices)):
        if index not in by_index:
            raise ForceFieldError(
                f"Reusable conformer {index + 1} has no matching force-field candidate"
            )
        if len(required) < count:
            required.append(by_index[index])

    def energy_key(candidate: ConformerCandidate) -> Tuple[float, int]:
        energy = candidate.energy_kcal_mol
        return (
            float("inf") if energy is None else energy,
            candidate.conformer_index,
        )

    selected = list(required)
    selected_pool_indices = {candidate.conformer_index for candidate in selected}
    for candidate in sorted(converged, key=energy_key):
        if candidate.conformer_index in selected_pool_indices:
            continue
        if all(
            aligned_heavy_atom_rmsd(candidate.model, retained.model) >= rmsd_threshold
            for retained in selected
        ):
            selected.append(candidate)
            selected_pool_indices.add(candidate.conformer_index)
            if len(selected) == count:
                break

    threshold_relaxed = set()
    if len(selected) < count:
        remaining = [
            candidate
            for candidate in sorted(converged, key=energy_key)
            if candidate.conformer_index not in selected_pool_indices
        ]
        for candidate in remaining:
            selected.append(candidate)
            selected_pool_indices.add(candidate.conformer_index)
            threshold_relaxed.add(candidate.conformer_index)
            if len(selected) == count:
                break

    pool_size = len(candidates)
    reindexed = []
    for job_index, candidate in enumerate(selected):
        pool_index = candidate.conformer_index
        model = candidate.model.with_metadata(
            conformer_index=job_index,
            conformer_pool_index=pool_index,
            conformer_pool_index_provenance="recorded_from_current_etkdg_pool",
            conformer_pool_size=pool_size,
            conformer_selection_rmsd_threshold_angstrom=rmsd_threshold,
            conformer_selection_threshold_relaxed=pool_index in threshold_relaxed,
        )
        reindexed.append(
            ConformerCandidate(
                model,
                job_index,
                candidate.energy_kcal_mol,
                candidate.forcefield,
                candidate.converged,
                candidate.optimization_iterations,
                candidate.failure_status,
            )
        )
    return reindexed


def geometry_warnings(model: MoleculeModel) -> Tuple[str, ...]:
    """Detect gross coordinate problems without claiming chemical validation."""

    warnings = []
    positions = [atom.position for atom in model.atoms]
    bonded = {tuple(sorted((bond.begin, bond.end))) for bond in model.bonds}
    for first in range(len(positions)):
        for second in range(first + 1, len(positions)):
            distance = math.dist(positions[first], positions[second])
            pair = (first, second)
            if distance < 0.45:
                warnings.append(
                    f"Abnormally short atom distance: {first}-{second} = {distance:.3f} A"
                )
            elif pair in bonded and distance > 3.20:
                warnings.append(
                    f"Unusually long modeled bond: {first}-{second} = {distance:.3f} A"
                )
            elif pair not in bonded and distance < 0.60:
                warnings.append(
                    f"Close nonbonded atoms: {first}-{second} = {distance:.3f} A"
                )
    return tuple(warnings)

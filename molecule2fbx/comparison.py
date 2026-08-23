"""Read-only-first analysis of completed Bz/TMS-Bz ORCA ensembles."""

from __future__ import annotations

import copy
import csv
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .ensemble import HARTREE_TO_KJ_MOL, write_run_summary
from .frequency import read_xyz_model
from .model import Atom, Bond, MoleculeModel, normalize_bond_order, validate_model
from .quantum.orca import parse_orca_output
from .structures import (
    aligned_atom_subset_rmsd,
    aligned_mapped_atom_rmsd,
    common_scaffold_atom_mapping,
    generate_forcefield_conformers,
    rmsd_atom_subsets,
)


R_KJ_MOL_K = 0.00831446261815324
DEFAULT_TEMPERATURE_K = 298.15


@dataclass(frozen=True)
class ExistingConformer:
    entry: Dict[str, object]
    model: MoleculeModel
    xyz_path: Path
    opt_output_path: Path
    frequency_output_path: Optional[Path]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def _load_json(path: Path) -> Dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return value


def _synchronize_ensemble_conformer_aliases(
    payload: Dict[str, object],
) -> None:
    """Keep the legacy ``conformers`` alias identical to the canonical array.

    JSON serialization breaks the shared-list identity used by
    :func:`build_ensemble_report`.  Any later metadata repair must therefore
    update both keys explicitly or consumers of the legacy alias can observe
    stale Freq/provenance values.
    """

    entries = payload.get("final_conformer_ensemble")
    if isinstance(entries, list):
        payload["conformers"] = copy.deepcopy(entries)


def _rdkit_model_from_xyz(
    xyz_path: Path,
    smiles: str,
    *,
    name: str,
) -> MoleculeModel:
    from rdkit import Chem

    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError("Could not parse ensemble SMILES")
    mol = Chem.AddHs(base)
    xyz_model = read_xyz_model(xyz_path)
    expected = [atom.GetSymbol().casefold() for atom in mol.GetAtoms()]
    observed = [atom.element.casefold() for atom in xyz_model.atoms]
    if expected != observed:
        raise ValueError(f"SMILES/XYZ atom order mismatch: {xyz_path}")
    bonds = tuple(
        Bond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            normalize_bond_order(
                float(bond.GetBondTypeAsDouble()), bond.GetIsAromatic()
            ),
        )
        for bond in mol.GetBonds()
    )
    return MoleculeModel(
        None,
        name,
        xyz_model.atoms,
        bonds,
        {
            "canonical_isomeric_smiles": Chem.MolToSmiles(
                base, isomericSmiles=True
            ),
            "original_smiles": smiles,
        },
    )


def _conformer_stem(entry: Mapping[str, object]) -> str:
    index = int(entry["conformer_index"])
    return f"conformer_{index + 1:03d}"


def _locate_frequency_output(
    ensemble_dir: Path,
    calculation_dir: Path,
    stem: str,
) -> Optional[Path]:
    candidates = (
        calculation_dir / f"{stem}_freq.out",
        calculation_dir.parent / "frequency_additions" / stem / f"{stem}_freq.out",
        ensemble_dir / "conformers" / "frequency_additions" / stem / f"{stem}_freq.out",
    )
    for path in candidates:
        if path.is_file():
            return path.resolve()
    return None


def load_existing_conformers(
    ensemble_dir: Path,
) -> Tuple[Dict[str, object], List[ExistingConformer]]:
    ensemble_dir = Path(ensemble_dir).resolve()
    payload = _load_json(ensemble_dir / "ensemble.json")
    molecule = payload.get("molecule", {})
    if not isinstance(molecule, dict) or not isinstance(molecule.get("smiles"), str):
        raise ValueError("ensemble.json does not contain a usable molecule SMILES")
    smiles = str(molecule.get("original_smiles") or molecule["smiles"])
    name = str(payload.get("name") or molecule.get("name") or ensemble_dir.name)
    entries = payload.get("final_conformer_ensemble", [])
    if not isinstance(entries, list):
        raise ValueError("final_conformer_ensemble must be an array")
    conformers = []
    for raw in entries:
        if not isinstance(raw, dict):
            continue
        entry = dict(raw)
        calculation = entry.get("calculation_directory")
        if not isinstance(calculation, str):
            raise ValueError(f"Missing calculation directory for {entry.get('conformer_id')}")
        calculation_dir = Path(calculation).resolve()
        stem = _conformer_stem(entry)
        xyz_path = calculation_dir / f"{stem}.xyz"
        opt_output = calculation_dir / f"{stem}.out"
        if not xyz_path.is_file() or not opt_output.is_file():
            raise ValueError(f"Missing retained Opt files for {entry.get('conformer_id')}")
        model = _rdkit_model_from_xyz(xyz_path, smiles, name=name)
        conformers.append(
            ExistingConformer(
                entry,
                model,
                xyz_path.resolve(),
                opt_output.resolve(),
                _locate_frequency_output(ensemble_dir, calculation_dir, stem),
            )
        )
    return payload, conformers


def load_targeted_followup_conformer(
    ensemble_payload: Mapping[str, object], output_dir: Path
) -> Optional[ExistingConformer]:
    """Load the optional pool-109 follow-up without mutating the original run."""

    result_path = output_dir / "sb_additional_dft_opt_result.json"
    if not result_path.is_file():
        return None
    result = _load_json(result_path)
    calculation_dir = Path(str(result["calculation_directory"])).resolve()
    xyz_path = Path(str(result["optimized_xyz"])).resolve()
    output_path = calculation_dir / "conformer_011.out"
    if not xyz_path.is_file() or not output_path.is_file():
        raise ValueError("Targeted follow-up Opt files are incomplete")
    molecule = ensemble_payload.get("molecule", {})
    if not isinstance(molecule, Mapping):
        raise ValueError("SB ensemble molecule metadata is unavailable")
    smiles = molecule.get("original_smiles") or molecule.get("smiles")
    if not isinstance(smiles, str):
        raise ValueError("SB ensemble SMILES is unavailable")
    model = _rdkit_model_from_xyz(xyz_path, smiles, name="1SB-LSD_RR")
    entry = {
        "conformer_id": "conf011_pool109",
        "conformer_index": 10,
        "dft_energy_hartree": float(result["energy_hartree"]),
        "relative_energy_kj_mol": float(
            result["relative_to_existing_current_best_kj_mol"]
        ),
        "gibbs_energy_hartree": None,
        "relative_gibbs_energy_kj_mol": None,
        "frequency_calculated": False,
        "frequency_provenance": None,
        "number_of_modes": None,
        "imaginary_modes": None,
        "imaginary_frequencies_cm1": [],
        "low_frequency_modes_cm1": [],
        "local_minimum_assessment": "not_evaluated",
        "optimization_provenance": "computed_targeted_followup",
        "source": "generated_followup",
        "conformer_pool_index": 109,
        "conformer_pool_index_provenance": "regenerated_exact_seed_pool",
        "calculation_directory": str(calculation_dir),
        "post_dft_stereochemistry_validation": result.get(
            "stereochemistry_validation"
        ),
    }
    return ExistingConformer(entry, model, xyz_path, output_path, None)


def _greedy_subset_clusters(
    conformers: Sequence[ExistingConformer],
    atom_indices: Sequence[int],
    threshold: float,
    prefix: str,
) -> Tuple[Dict[int, str], Dict[int, float]]:
    ranked = sorted(
        conformers, key=lambda item: float(item.entry["dft_energy_hartree"])
    )
    representatives: List[ExistingConformer] = []
    assignments: Dict[int, str] = {}
    rmsds: Dict[int, float] = {}
    for conformer in ranked:
        assigned = None
        assigned_rmsd = 0.0
        for cluster_index, representative in enumerate(representatives):
            value = aligned_atom_subset_rmsd(
                conformer.model, representative.model, atom_indices
            )
            if value < threshold:
                assigned = cluster_index
                assigned_rmsd = value
                break
        if assigned is None:
            assigned = len(representatives)
            representatives.append(conformer)
        index = int(conformer.entry["conformer_index"])
        assignments[index] = f"{prefix}_{assigned + 1:03d}"
        rmsds[index] = assigned_rmsd
    return assignments, rmsds


def secondary_rmsd_analysis(
    conformers: Sequence[ExistingConformer],
    *,
    common_threshold: float = 0.75,
    reaction_threshold: float = 0.25,
) -> Dict[str, object]:
    if not conformers:
        return {}
    subsets = rmsd_atom_subsets(conformers[0].model)
    common_ids, common_rmsds = _greedy_subset_clusters(
        conformers,
        subsets.common_scaffold,
        common_threshold,
        "dft_core_cluster",
    )
    reaction_ids, reaction_rmsds = _greedy_subset_clusters(
        conformers,
        subsets.reaction_center,
        reaction_threshold,
        "dft_site_cluster",
    )
    records = []
    for conformer in sorted(
        conformers, key=lambda item: int(item.entry["conformer_index"])
    ):
        index = int(conformer.entry["conformer_index"])
        records.append(
            {
                "conformer_index": index,
                "conformer_id": conformer.entry["conformer_id"],
                "all_heavy_rmsd_cluster_id": conformer.entry.get(
                    "all_heavy_rmsd_cluster_id",
                    conformer.entry.get("rmsd_cluster_id"),
                ),
                "common_scaffold_rmsd_cluster_id": common_ids[index],
                "common_scaffold_rmsd_to_representative_angstrom": common_rmsds[
                    index
                ],
                "reaction_center_rmsd_cluster_id": reaction_ids[index],
                "reaction_center_rmsd_to_representative_angstrom": reaction_rmsds[
                    index
                ],
            }
        )
    return {
        "method": "greedy_energy_ordered_fixed_atom_order_kabsch",
        "atom_subsets": subsets.to_metadata(),
        "common_scaffold": {
            "threshold_angstrom": common_threshold,
            "cluster_count": len(set(common_ids.values())),
        },
        "reaction_center": {
            "threshold_angstrom": reaction_threshold,
            "cluster_count": len(set(reaction_ids.values())),
        },
        "records": records,
    }


def cross_ensemble_rmsd_analysis(
    first_conformers: Sequence[ExistingConformer],
    second_conformers: Sequence[ExistingConformer],
) -> Dict[str, object]:
    """Compare Bz and TMS-Bz using a graph-derived common atom mapping."""

    if not first_conformers or not second_conformers:
        return {"records": [], "reciprocal_nearest_pairs": []}
    mapping = common_scaffold_atom_mapping(
        first_conformers[0].model, second_conformers[0].model
    )
    records = []
    for first in first_conformers:
        for second in second_conformers:
            records.append(
                {
                    "first_conformer_id": first.entry["conformer_id"],
                    "second_conformer_id": second.entry["conformer_id"],
                    "common_scaffold_rmsd_angstrom": aligned_mapped_atom_rmsd(
                        first.model,
                        second.model,
                        mapping.first_indices,
                        mapping.second_indices,
                    ),
                    "reaction_center_rmsd_angstrom": aligned_mapped_atom_rmsd(
                        first.model,
                        second.model,
                        mapping.first_reaction_indices,
                        mapping.second_reaction_indices,
                    ),
                }
            )
    first_nearest = {}
    second_nearest = {}
    for record in records:
        first_id = str(record["first_conformer_id"])
        second_id = str(record["second_conformer_id"])
        if (
            first_id not in first_nearest
            or float(record["common_scaffold_rmsd_angstrom"])
            < float(first_nearest[first_id]["common_scaffold_rmsd_angstrom"])
        ):
            first_nearest[first_id] = record
        if (
            second_id not in second_nearest
            or float(record["common_scaffold_rmsd_angstrom"])
            < float(second_nearest[second_id]["common_scaffold_rmsd_angstrom"])
        ):
            second_nearest[second_id] = record
    reciprocal = [
        record
        for record in first_nearest.values()
        if second_nearest[str(record["second_conformer_id"])] is record
    ]
    return {
        "first_name": first_conformers[0].model.name,
        "second_name": second_conformers[0].model.name,
        "atom_mapping": mapping.to_metadata(),
        "records": records,
        "first_to_nearest_second": list(first_nearest.values()),
        "second_to_nearest_first": list(second_nearest.values()),
        "reciprocal_nearest_pairs": reciprocal,
        "interpretation": (
            "RMSD is minimized only by rigid Kabsch alignment for one deterministic "
            "graph mapping; symmetry-equivalent atom permutations are not searched."
        ),
    }


_ATOMIC_CHARGE_RE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]+)\s*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
    re.MULTILINE,
)
_MAYER_BOND_RE = re.compile(
    r"B\(\s*(\d+)-[A-Za-z]+\s*,\s*(\d+)-[A-Za-z]+\s*\)\s*:\s*"
    r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"
)
_ORBITAL_RE = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+"
    r"(-?\d+(?:\.\d+)?)\s*$",
    re.MULTILINE,
)


def _last_section(text: str, marker: str) -> str:
    start = text.rfind(marker)
    return text[start:] if start >= 0 else ""


def _atomic_charges(text: str, marker: str) -> Optional[Dict[int, float]]:
    section = _last_section(text, marker)
    if not section:
        return None
    values = {}
    started = False
    for line in section.splitlines()[1:]:
        match = _ATOMIC_CHARGE_RE.match(line)
        if match:
            started = True
            values[int(match.group(1))] = float(match.group(3))
        elif started and not line.strip():
            break
    return values or None


def _mbis_charges(text: str) -> Optional[Dict[int, float]]:
    section = _last_section(text, "MBIS ANALYSIS")
    if not section:
        return None
    values = {}
    row = re.compile(
        r"^\s*(\d+)\s+[A-Za-z]+\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+",
        re.MULTILINE,
    )
    started = False
    for line in section.splitlines()[1:]:
        match = row.match(line)
        if match:
            started = True
            values[int(match.group(1))] = float(match.group(2))
        elif started and not line.strip():
            break
    return values or None


def _chelpg_charges(text: str) -> Optional[Dict[int, float]]:
    section = _last_section(text, "CHELPG Charges")
    if not section:
        section = _last_section(text, "CHELPG CHARGES")
    if not section:
        return None
    values = {}
    row = re.compile(
        r"^\s*(\d+)\s+[A-Za-z]+\s*:?[ \t]+"
        r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
        re.MULTILINE,
    )
    started = False
    for line in section.splitlines()[1:]:
        match = row.match(line)
        if match:
            started = True
            values[int(match.group(1))] = float(match.group(2))
        elif started and not line.strip():
            break
    return values or None


def _frontier_mo_populations(text: str) -> Optional[Dict[str, object]]:
    section = _last_section(text, "FRONTIER MOLECULAR ORBITAL POPULATION ANALYSIS")
    if not section:
        return None
    orbital_match = re.search(r"HOMO=\s*(\d+)\s+LUMO=\s*(\d+)", section)
    row = re.compile(
        r"^\s*(\d+)-([A-Za-z]+)\s+"
        r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+"
        r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$",
        re.MULTILINE,
    )
    atoms = {
        int(index): {
            "element": element,
            "homo_mulliken": float(homo_mulliken),
            "homo_loewdin": float(homo_loewdin),
            "lumo_mulliken": float(lumo_mulliken),
            "lumo_loewdin": float(lumo_loewdin),
        }
        for index, element, homo_mulliken, homo_loewdin, lumo_mulliken, lumo_loewdin in row.findall(
            section[:30000]
        )
    }
    if not atoms:
        return None
    return {
        "homo_index": int(orbital_match.group(1)) if orbital_match else None,
        "lumo_index": int(orbital_match.group(2)) if orbital_match else None,
        "atom_populations": atoms,
    }


def parse_existing_orca_electronic(path: Path) -> Dict[str, object]:
    """Extract only properties already printed by a completed ORCA Opt job."""

    text = path.read_text(encoding="utf-8", errors="replace")
    mulliken = _atomic_charges(text, "MULLIKEN ATOMIC CHARGES")
    loewdin = _atomic_charges(text, "LOEWDIN ATOMIC CHARGES")
    mbis = _mbis_charges(text)
    chelpg = _chelpg_charges(text)
    frontier = _frontier_mo_populations(text)

    mayer_section = _last_section(text, "MAYER POPULATION ANALYSIS")
    mayer = {}
    if mayer_section:
        for first, second, value in _MAYER_BOND_RE.findall(
            mayer_section[:30000]
        ):
            mayer[tuple(sorted((int(first), int(second))))] = float(value)

    orbital_section = _last_section(text, "ORBITAL ENERGIES")
    orbitals = []
    if orbital_section:
        for index, occupation, energy_h, energy_ev in _ORBITAL_RE.findall(
            orbital_section[:20000]
        ):
            orbitals.append(
                {
                    "index": int(index),
                    "occupation": float(occupation),
                    "energy_hartree": float(energy_h),
                    "energy_ev": float(energy_ev),
                }
            )
    occupied = [item for item in orbitals if item["occupation"] > 1.0e-8]
    virtual = [item for item in orbitals if item["occupation"] <= 1.0e-8]
    homo = occupied[-1] if occupied else None
    lumo = virtual[0] if virtual else None

    dipole_section = _last_section(text, "DIPOLE MOMENT")
    dipole_match = re.search(
        r"Magnitude \(Debye\)\s*:\s*(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
        dipole_section,
    )
    return {
        "mulliken_charges": mulliken,
        "loewdin_charges": loewdin,
        "mbis_charges": mbis,
        "chelpg_charges": chelpg,
        "mayer_bond_orders": mayer,
        "dipole_magnitude_debye": (
            float(dipole_match.group(1)) if dipole_match else None
        ),
        "homo": homo,
        "lumo": lumo,
        "homo_lumo_gap_ev": (
            float(lumo["energy_ev"]) - float(homo["energy_ev"])
            if homo is not None and lumo is not None
            else None
        ),
        "frontier_orbital_local_contributions": frontier,
        "frontier_orbital_local_contributions_reason": (
            None
            if frontier is not None
            else "The output does not print per-MO atom populations; a new single-point "
            "property job or GBW post-processing is required."
        ),
    }


def _distance(model: MoleculeModel, first: int, second: int) -> float:
    return math.dist(model.atoms[first].position, model.atoms[second].position)


def _dihedral_degrees(
    model: MoleculeModel, first: int, second: int, third: int, fourth: int
) -> float:
    points = np.asarray(
        [model.atoms[index].position for index in (first, second, third, fourth)],
        dtype=float,
    )
    b0 = -(points[1] - points[0])
    b1 = points[2] - points[1]
    b2 = points[3] - points[2]
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _folded_planarity_angle(value: float) -> float:
    value = abs(value) % 180.0
    return min(value, 180.0 - value)


def structural_descriptors(model: MoleculeModel) -> Dict[str, object]:
    subsets = rmsd_atom_subsets(model)
    if subsets.benzoyl_carbonyl_c is None:
        raise ValueError("Could not identify the aromatic N-benzoyl carbonyl")
    c = int(subsets.benzoyl_carbonyl_c)
    o = int(subsets.benzoyl_carbonyl_o)  # type: ignore[arg-type]
    n = int(subsets.benzoyl_amide_n)  # type: ignore[arg-type]
    ipso = int(subsets.benzoyl_ipso_c)  # type: ignore[arg-type]
    neighbors: Dict[int, List[int]] = {index: [] for index in range(len(model.atoms))}
    for bond in model.bonds:
        neighbors[bond.begin].append(bond.end)
        neighbors[bond.end].append(bond.begin)
    aryl_neighbors = [
        index
        for index in neighbors[ipso]
        if index != c and model.atoms[index].element.casefold() == "c"
    ]
    n_neighbors = [
        index
        for index in neighbors[n]
        if index != c and model.atoms[index].element.casefold() != "h"
    ]
    aryl_torsions = [
        _folded_planarity_angle(_dihedral_degrees(model, n, c, ipso, index))
        for index in aryl_neighbors
    ]
    amide_torsions = [
        _folded_planarity_angle(_dihedral_degrees(model, index, n, c, o))
        for index in n_neighbors
    ]
    c_si_lengths = []
    for index in subsets.excluded_terminal_substituent:
        if model.atoms[index].element.casefold() != "si":
            continue
        for neighbor in neighbors[index]:
            if neighbor not in subsets.excluded_terminal_substituent:
                c_si_lengths.append(_distance(model, index, neighbor))
    return {
        "atom_indices": {
            "benzoyl_carbonyl_c": c,
            "benzoyl_carbonyl_o": o,
            "benzoyl_amide_n": n,
            "benzoyl_ipso_c": ipso,
            "benzoyl_ortho_c": aryl_neighbors,
            "amide_n_ring_neighbors": n_neighbors,
        },
        "carbonyl_c_o_length_angstrom": _distance(model, c, o),
        "amide_c_n_length_angstrom": _distance(model, c, n),
        "benzoyl_carbonyl_torsions_deg": aryl_torsions,
        "benzoyl_carbonyl_torsion_abs_mean_deg": (
            sum(aryl_torsions) / len(aryl_torsions) if aryl_torsions else None
        ),
        "amide_n_planarity_torsions_deg": amide_torsions,
        "amide_n_planarity_torsion_abs_mean_deg": (
            sum(amide_torsions) / len(amide_torsions) if amide_torsions else None
        ),
        "aryl_c_si_length_angstrom": c_si_lengths[0] if c_si_lengths else None,
    }


def replace_terminal_tms_with_hydrogen(
    model: MoleculeModel,
    *,
    aromatic_c_h_length_angstrom: float = 1.085,
) -> MoleculeModel:
    """Return a fixed-scaffold counterfactual with terminal SiMe3 replaced by H.

    The common-atom coordinates are retained exactly.  The replacement hydrogen
    is placed along the original aryl-C--Si vector, so the model can be used to
    separate the direct electronic effect of TMS from geometry relaxation.  The
    function deliberately accepts only one unambiguous terminal trimethylsilyl
    group and does not attempt to edit other organosilicon topologies.
    """

    if not math.isfinite(aromatic_c_h_length_angstrom) or aromatic_c_h_length_angstrom <= 0:
        raise ValueError("aromatic_c_h_length_angstrom must be positive")

    neighbors: Dict[int, List[int]] = {index: [] for index in range(len(model.atoms))}
    for bond in model.bonds:
        neighbors[bond.begin].append(bond.end)
        neighbors[bond.end].append(bond.begin)

    candidates: List[Tuple[int, int]] = []
    for atom in model.atoms:
        if atom.element.casefold() != "si":
            continue
        carbon_neighbors = [
            index
            for index in neighbors[atom.index]
            if model.atoms[index].element.casefold() == "c"
        ]
        if len(carbon_neighbors) != 4:
            continue
        anchors = []
        terminal_methyls = []
        for index in carbon_neighbors:
            heavy_degree = sum(
                model.atoms[neighbor].element.casefold() != "h"
                for neighbor in neighbors[index]
            )
            if heavy_degree > 1:
                anchors.append(index)
            else:
                terminal_methyls.append(index)
        if len(anchors) == 1 and len(terminal_methyls) == 3:
            candidates.append((atom.index, anchors[0]))
    if len(candidates) != 1:
        raise ValueError(
            "Expected exactly one unambiguous terminal trimethylsilyl group, "
            f"found {len(candidates)}"
        )

    silicon, anchor = candidates[0]
    removed = {silicon}
    pending = [index for index in neighbors[silicon] if index != anchor]
    while pending:
        index = pending.pop()
        if index == anchor or index in removed:
            continue
        removed.add(index)
        pending.extend(
            neighbor
            for neighbor in neighbors[index]
            if neighbor != anchor and neighbor not in removed
        )

    anchor_atom = model.atoms[anchor]
    silicon_atom = model.atoms[silicon]
    vector = np.asarray(silicon_atom.position) - np.asarray(anchor_atom.position)
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-12:
        raise ValueError("Aryl-C--Si vector has zero length")
    hydrogen_position = (
        np.asarray(anchor_atom.position)
        + vector * (aromatic_c_h_length_angstrom / norm)
    )
    # Match RDKit AddHs ordering: retained heavy atoms first, followed by each
    # heavy atom's hydrogens in heavy-atom order.  This keeps the edited SMILES,
    # XYZ, and population-analysis atom indices under one explicit contract.
    retained_heavy = [
        atom.index
        for atom in model.atoms
        if atom.index not in removed and atom.element.casefold() != "h"
    ]
    ordered_source_indices = list(retained_heavy)
    replacement_after_source_position = None
    for heavy_index in retained_heavy:
        ordered_source_indices.extend(
            sorted(
                index
                for index in neighbors[heavy_index]
                if index not in removed
                and model.atoms[index].element.casefold() == "h"
            )
        )
        if heavy_index == anchor:
            replacement_after_source_position = len(ordered_source_indices)
    retained_source_indices = {
        atom.index for atom in model.atoms if atom.index not in removed
    }
    if set(ordered_source_indices) != retained_source_indices:
        raise ValueError("Could not establish explicit-H counterfactual atom order")
    if replacement_after_source_position is None:
        raise ValueError("Could not place the replacement aromatic hydrogen")

    old_to_new: Dict[int, int] = {}
    retained_atoms: List[Atom] = []
    hydrogen_index = -1
    for source_position in range(len(ordered_source_indices) + 1):
        if source_position == replacement_after_source_position:
            hydrogen_index = len(retained_atoms)
            retained_atoms.append(
                Atom(
                    hydrogen_index,
                    "H",
                    float(hydrogen_position[0]),
                    float(hydrogen_position[1]),
                    float(hydrogen_position[2]),
                )
            )
        if source_position == len(ordered_source_indices):
            continue
        old_index = ordered_source_indices[source_position]
        atom = model.atoms[old_index]
        new_index = len(retained_atoms)
        old_to_new[old_index] = new_index
        retained_atoms.append(Atom(new_index, atom.element, atom.x, atom.y, atom.z))
    if hydrogen_index < 0:
        raise ValueError("Replacement hydrogen was not added")

    retained_bonds = [
        Bond(old_to_new[bond.begin], old_to_new[bond.end], bond.order)
        for bond in model.bonds
        if bond.begin not in removed and bond.end not in removed
    ]
    retained_bonds.append(Bond(old_to_new[anchor], hydrogen_index, 1))
    metadata = dict(model.metadata)
    source_smiles = metadata.get("original_smiles")
    if not isinstance(source_smiles, str) or not source_smiles:
        source_smiles = metadata.get("canonical_isomeric_smiles")
    if isinstance(source_smiles, str) and source_smiles:
        from rdkit import Chem

        source_mol = Chem.MolFromSmiles(source_smiles)
        if source_mol is None:
            raise ValueError("Could not parse source SMILES for TMS-to-H editing")
        heavy_elements = [
            atom.element.casefold()
            for atom in model.atoms
            if atom.element.casefold() != "h"
        ]
        if [atom.GetSymbol().casefold() for atom in source_mol.GetAtoms()] != heavy_elements:
            raise ValueError("Source SMILES heavy-atom order does not match the model")
        editable = Chem.RWMol(source_mol)
        for index in sorted(
            (index for index in removed if model.atoms[index].element.casefold() != "h"),
            reverse=True,
        ):
            editable.RemoveAtom(index)
        edited = editable.GetMol()
        Chem.SanitizeMol(edited)
        metadata["original_smiles"] = Chem.MolToSmiles(
            edited, canonical=False, isomericSmiles=True
        )
        metadata["canonical_isomeric_smiles"] = Chem.MolToSmiles(
            edited, canonical=True, isomericSmiles=True
        )
    metadata.update(
        {
            "counterfactual": "terminal_tms_replaced_by_hydrogen",
            "counterfactual_geometry_optimized": False,
            "common_atom_coordinates_retained": True,
            "replacement_c_h_length_angstrom": aromatic_c_h_length_angstrom,
            "source_silicon_atom_index": silicon,
            "source_aryl_anchor_atom_index": anchor,
            "removed_source_atom_indices": sorted(removed),
            "source_to_counterfactual_atom_indices": old_to_new,
            "replacement_hydrogen_atom_index": hydrogen_index,
        }
    )
    result = MoleculeModel(
        model.cid,
        f"{model.name}_TMS_to_H_fixed_geometry",
        tuple(retained_atoms),
        tuple(retained_bonds),
        metadata,
    )
    validate_model(result)
    return result


_VDW_RADII = {
    "h": 1.20,
    "c": 1.70,
    "n": 1.55,
    "o": 1.52,
    "f": 1.47,
    "si": 2.10,
    "p": 1.80,
    "s": 1.80,
    "cl": 1.75,
}


def carbonyl_steric_access(
    model: MoleculeModel,
    *,
    burgi_dunitz_angle_deg: float = 107.0,
    attack_distance_angstrom: float = 2.2,
    probe_radius_angstrom: float = 0.5,
    azimuth_samples: int = 72,
) -> Dict[str, object]:
    """Sample hard-sphere clearance around a Bürgi-Dunitz attack cone.

    This is a geometry descriptor, not a reaction barrier or an enzyme-pocket
    model.  The two faces are kept as signed geometric faces because assigning
    Re/Si labels without a dedicated stereochemical face implementation would
    overstate what is encoded here.
    """

    subsets = rmsd_atom_subsets(model)
    if subsets.benzoyl_carbonyl_c is None:
        raise ValueError("Could not identify the aromatic N-benzoyl carbonyl")
    c = int(subsets.benzoyl_carbonyl_c)
    o = int(subsets.benzoyl_carbonyl_o)  # type: ignore[arg-type]
    n = int(subsets.benzoyl_amide_n)  # type: ignore[arg-type]
    center = np.asarray(model.atoms[c].position, dtype=float)
    co = np.asarray(model.atoms[o].position, dtype=float) - center
    cn = np.asarray(model.atoms[n].position, dtype=float) - center
    axis = co / np.linalg.norm(co)
    in_plane = cn - np.dot(cn, axis) * axis
    in_plane /= np.linalg.norm(in_plane)
    normal = np.cross(axis, in_plane)
    normal /= np.linalg.norm(normal)
    theta = math.radians(burgi_dunitz_angle_deg)
    records = []
    for sample in range(azimuth_samples):
        phi = 2.0 * math.pi * sample / azimuth_samples
        radial = math.cos(phi) * in_plane + math.sin(phi) * normal
        direction = math.cos(theta) * axis + math.sin(theta) * radial
        probe = center + attack_distance_angstrom * direction
        clearances = []
        limiting_atom = None
        for atom in model.atoms:
            if atom.index == c:
                continue
            radius = _VDW_RADII.get(atom.element.casefold(), 1.80)
            clearance = (
                float(np.linalg.norm(probe - np.asarray(atom.position, dtype=float)))
                - radius
                - probe_radius_angstrom
            )
            clearances.append(clearance)
            if limiting_atom is None or clearance < limiting_atom["clearance_angstrom"]:
                limiting_atom = {
                    "atom_index": atom.index,
                    "element": atom.element,
                    "clearance_angstrom": clearance,
                }
        signed_face = (
            "face_positive"
            if sample < azimuth_samples // 2
            else "face_negative"
        )
        records.append(
            {
                "azimuth_deg": math.degrees(phi),
                "face": signed_face,
                "clearance_angstrom": min(clearances),
                "accessible": min(clearances) >= 0.0,
                "limiting_atom": limiting_atom,
            }
        )

    face_summary = {}
    for face in ("face_positive", "face_negative"):
        face_records = [record for record in records if record["face"] == face]
        face_summary[face] = {
            "accessible_fraction": sum(
                1 for record in face_records if record["accessible"]
            )
            / len(face_records),
            "best_clearance_angstrom": max(
                float(record["clearance_angstrom"]) for record in face_records
            ),
            "mean_clearance_angstrom": sum(
                float(record["clearance_angstrom"]) for record in face_records
            )
            / len(face_records),
        }

    clearance_values = sorted(
        float(record["clearance_angstrom"]) for record in records
    )

    def percentile(values: Sequence[float], fraction: float) -> float:
        if len(values) == 1:
            return values[0]
        position = fraction * (len(values) - 1)
        lower = int(math.floor(position))
        upper = int(math.ceil(position))
        if lower == upper:
            return values[lower]
        weight = position - lower
        return values[lower] * (1.0 - weight) + values[upper] * weight

    sensitivity = {}
    for radius in (0.0, 0.5, 0.8, 1.0, 1.2):
        shifted = [
            float(record["clearance_angstrom"])
            + probe_radius_angstrom
            - radius
            for record in records
        ]
        sensitivity[f"probe_radius_{radius:.1f}_angstrom"] = {
            "accessible_fraction": sum(value >= 0.0 for value in shifted)
            / len(shifted),
            "best_clearance_angstrom": max(shifted),
        }

    heavy_distances = []
    tms_indices = set(subsets.excluded_terminal_substituent)
    for atom in model.atoms:
        if atom.index == c or atom.element.casefold() == "h":
            continue
        heavy_distances.append(
            {
                "atom_index": atom.index,
                "element": atom.element,
                "distance_angstrom": _distance(model, c, atom.index),
                "terminal_tms_atom": atom.index in tms_indices,
            }
        )
    heavy_distances.sort(key=lambda item: float(item["distance_angstrom"]))
    terminal_tms_distances = [
        item for item in heavy_distances if item["terminal_tms_atom"]
    ]
    terminal_si_distance = next(
        (
            float(item["distance_angstrom"])
            for item in terminal_tms_distances
            if str(item["element"]).casefold() == "si"
        ),
        None,
    )
    return {
        "method": "sampled_vdw_probe_on_burgi_dunitz_cone",
        "parameters": {
            "burgi_dunitz_o_c_probe_angle_deg": burgi_dunitz_angle_deg,
            "carbonyl_c_probe_distance_angstrom": attack_distance_angstrom,
            "probe_radius_angstrom": probe_radius_angstrom,
            "azimuth_samples": azimuth_samples,
            "vdw_radii_angstrom": dict(_VDW_RADII),
        },
        "face_assignment": {
            "labels": ["face_positive", "face_negative"],
            "re_si_assigned": False,
            "reason": "Signed faces are reported without an unvalidated Re/Si assignment.",
        },
        "overall_accessible_fraction": sum(
            1 for record in records if record["accessible"]
        )
        / len(records),
        "best_clearance_angstrom": max(
            float(record["clearance_angstrom"]) for record in records
        ),
        "mean_clearance_angstrom": sum(clearance_values) / len(clearance_values),
        "clearance_percentiles_angstrom": {
            "p10": percentile(clearance_values, 0.10),
            "p50": percentile(clearance_values, 0.50),
            "p90": percentile(clearance_values, 0.90),
        },
        "face_summary": face_summary,
        "probe_radius_sensitivity": sensitivity,
        "nearest_heavy_atoms_from_carbonyl_c": heavy_distances[:12],
        "terminal_tms_distances_from_carbonyl_c": terminal_tms_distances,
        "minimum_terminal_tms_distance_from_carbonyl_c_angstrom": (
            min(float(item["distance_angstrom"]) for item in terminal_tms_distances)
            if terminal_tms_distances
            else None
        ),
        "silicon_distance_from_carbonyl_c_angstrom": terminal_si_distance,
        "terminal_tms_heavy_atoms_within_5_angstrom": sum(
            float(item["distance_angstrom"]) <= 5.0
            for item in terminal_tms_distances
        ),
        "interpretation_limit": (
            "Hard-sphere access is a directional reactant-geometry descriptor. It does "
            "not include solvent, protein motion, nucleophile identity, or activation energy."
        ),
    }


def _percentile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("Cannot calculate a percentile for an empty sequence")
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _clearance_statistics(values: Sequence[float]) -> Dict[str, object]:
    numeric = [float(value) for value in values]
    if not numeric:
        return {
            "sample_count": 0,
            "mean_clearance_angstrom": None,
            "best_clearance_angstrom": None,
            "minimum_clearance_angstrom": None,
            "clearance_percentiles_angstrom": None,
            "positive_clearance_fraction": None,
            "positive_clearance_integral_angstrom": None,
        }
    positive = [max(value, 0.0) for value in numeric]
    top_count = max(1, int(math.ceil(0.10 * len(numeric))))
    return {
        "sample_count": len(numeric),
        "mean_clearance_angstrom": sum(numeric) / len(numeric),
        "best_clearance_angstrom": max(numeric),
        "minimum_clearance_angstrom": min(numeric),
        "clearance_percentiles_angstrom": {
            "p10": _percentile(numeric, 0.10),
            "p50": _percentile(numeric, 0.50),
            "p90": _percentile(numeric, 0.90),
        },
        "top_decile_mean_clearance_angstrom": (
            sum(sorted(numeric, reverse=True)[:top_count]) / top_count
        ),
        "positive_clearance_fraction": (
            sum(value >= 0.0 for value in numeric) / len(numeric)
        ),
        "positive_clearance_integral_angstrom": sum(positive) / len(positive),
    }


def _point_to_segment_distance(
    point: np.ndarray, first: np.ndarray, second: np.ndarray
) -> Tuple[float, float]:
    segment = second - first
    denominator = float(np.dot(segment, segment))
    if denominator <= 0.0:
        return float(np.linalg.norm(point - first)), 0.0
    fraction = float(np.dot(point - first, segment) / denominator)
    fraction = min(1.0, max(0.0, fraction))
    closest = first + fraction * segment
    return float(np.linalg.norm(point - closest)), fraction


def _probe_sensitivity(
    base_clearances: Sequence[float], probe_radii: Sequence[float]
) -> Dict[str, object]:
    return {
        f"probe_radius_{float(radius):.2f}_angstrom": _clearance_statistics(
            [float(value) - float(radius) for value in base_clearances]
        )
        for radius in probe_radii
    }


def carbonyl_continuous_steric_access(
    model: MoleculeModel,
    *,
    burgi_dunitz_angle_deg: float = 107.0,
    near_distance_angstrom: float = 2.0,
    far_distance_angstrom: float = 8.0,
    probe_radii_angstrom: Sequence[float] = (0.0, 0.5, 1.0, 1.4),
    azimuth_samples: int = 720,
    include_trajectory: bool = False,
) -> Dict[str, object]:
    """Measure continuous hard-sphere clearance along carbonyl attack paths.

    Each azimuth defines a straight approach corridor from ``far_distance`` to
    ``near_distance`` at the requested Buerghi-Dunitz angle.  Clearance is the
    minimum distance between that corridor and an atomic van der Waals sphere.
    ``total`` includes all atoms except the attacked carbonyl carbon, while
    ``nonlocal_environment`` excludes the carbonyl O, amide N, and benzoyl ipso
    carbon as well.  The latter prevents the shared local trigonal geometry from
    hiding differences in remote steric shielding.

    This remains a static gas-phase geometry descriptor.  It is not an enzyme
    pocket model and does not calculate a hydrolysis activation barrier.
    """

    if azimuth_samples < 8:
        raise ValueError("azimuth_samples must be at least 8")
    if near_distance_angstrom <= 0.0:
        raise ValueError("near_distance_angstrom must be positive")
    if far_distance_angstrom <= near_distance_angstrom:
        raise ValueError("far_distance_angstrom must exceed near_distance_angstrom")
    radii = tuple(float(radius) for radius in probe_radii_angstrom)
    if not radii or any(radius < 0.0 for radius in radii):
        raise ValueError("probe radii must contain non-negative values")

    subsets = rmsd_atom_subsets(model)
    if subsets.benzoyl_carbonyl_c is None:
        raise ValueError("Could not identify the aromatic N-benzoyl carbonyl")
    c = int(subsets.benzoyl_carbonyl_c)
    o = int(subsets.benzoyl_carbonyl_o)  # type: ignore[arg-type]
    n = int(subsets.benzoyl_amide_n)  # type: ignore[arg-type]

    neighbors: Dict[int, List[int]] = {atom.index: [] for atom in model.atoms}
    for bond in model.bonds:
        neighbors[bond.begin].append(bond.end)
        neighbors[bond.end].append(bond.begin)
    ipso_candidates = [
        index
        for index in neighbors[c]
        if index not in {o, n}
        and model.atoms[index].element.casefold() == "c"
    ]
    if len(ipso_candidates) != 1:
        raise ValueError("Could not identify the benzoyl ipso carbon")
    ipso = ipso_candidates[0]

    terminal_tms_heavy = set(subsets.excluded_terminal_substituent)
    terminal_tms_group = set(terminal_tms_heavy)
    for index in tuple(terminal_tms_heavy):
        terminal_tms_group.update(
            neighbor
            for neighbor in neighbors[index]
            if model.atoms[neighbor].element.casefold() == "h"
        )

    center = np.asarray(model.atoms[c].position, dtype=float)
    priority_vectors = [
        np.asarray(model.atoms[index].position, dtype=float) - center
        for index in (o, n, ipso)
    ]
    co = priority_vectors[0]
    cn = priority_vectors[1]
    axis = co / np.linalg.norm(co)
    in_plane = cn - np.dot(cn, axis) * axis
    in_plane /= np.linalg.norm(in_plane)
    normal = np.cross(axis, in_plane)
    normal /= np.linalg.norm(normal)
    theta = math.radians(burgi_dunitz_angle_deg)
    local_environment = {c, o, n, ipso}

    records: List[Dict[str, object]] = []
    for sample in range(azimuth_samples):
        phi = 2.0 * math.pi * sample / azimuth_samples
        radial = math.cos(phi) * in_plane + math.sin(phi) * normal
        direction = math.cos(theta) * axis + math.sin(theta) * radial
        near = center + near_distance_angstrom * direction
        far = center + far_distance_angstrom * direction

        orientation = float(
            np.dot(
                direction,
                np.cross(
                    priority_vectors[1] - priority_vectors[0],
                    priority_vectors[2] - priority_vectors[0],
                ),
            )
        )
        if abs(orientation) <= 1.0e-10:
            face = "coplanar"
        else:
            # Viewed from the approaching nucleophile toward the carbonyl C,
            # clockwise O > N > aryl is Re and gives a negative determinant.
            face = "Re" if orientation < 0.0 else "Si"

        scope_values: Dict[str, Optional[float]] = {
            "total": None,
            "total_without_terminal_tms": None,
            "nonlocal_environment": None,
            "nonlocal_without_terminal_tms": None,
            "terminal_tms": None,
        }
        scope_limiters: Dict[str, Optional[Dict[str, object]]] = {
            key: None for key in scope_values
        }
        for atom in model.atoms:
            if atom.index == c:
                continue
            atom_position = np.asarray(atom.position, dtype=float)
            distance, path_fraction = _point_to_segment_distance(
                atom_position, near, far
            )
            clearance = distance - _VDW_RADII.get(
                atom.element.casefold(), 1.80
            )
            scopes = ["total"]
            if atom.index not in terminal_tms_group:
                scopes.append("total_without_terminal_tms")
            if atom.index not in local_environment:
                scopes.append("nonlocal_environment")
                if atom.index not in terminal_tms_group:
                    scopes.append("nonlocal_without_terminal_tms")
            if atom.index in terminal_tms_group:
                scopes.append("terminal_tms")
            for scope in scopes:
                current = scope_values[scope]
                if current is None or clearance < current:
                    scope_values[scope] = clearance
                    scope_limiters[scope] = {
                        "atom_index": atom.index,
                        "element": atom.element,
                        "clearance_without_probe_angstrom": clearance,
                        "path_fraction_near_to_far": path_fraction,
                    }
        records.append(
            {
                "azimuth_deg": math.degrees(phi),
                "face": face,
                "clearance_without_probe_angstrom": scope_values,
                "limiting_atoms": scope_limiters,
            }
        )

    scopes = {}
    for scope in (
        "total",
        "total_without_terminal_tms",
        "nonlocal_environment",
        "nonlocal_without_terminal_tms",
        "terminal_tms",
    ):
        values = [
            float(value)
            for record in records
            for value in [record["clearance_without_probe_angstrom"][scope]]  # type: ignore[index]
            if value is not None
        ]
        face_summary = {}
        for face in ("Re", "Si"):
            face_values = [
                float(record["clearance_without_probe_angstrom"][scope])  # type: ignore[index]
                for record in records
                if record["face"] == face
                and record["clearance_without_probe_angstrom"][scope] is not None  # type: ignore[index]
            ]
            face_summary[face] = {
                "clearance_without_probe": _clearance_statistics(face_values),
                "probe_radius_sensitivity": _probe_sensitivity(face_values, radii),
            }
        scopes[scope] = {
            "clearance_without_probe": _clearance_statistics(values),
            "probe_radius_sensitivity": _probe_sensitivity(values, radii),
            "face_summary": face_summary,
        }

    tms_limiter_count = sum(
        record["limiting_atoms"]["total"] is not None  # type: ignore[index]
        and int(record["limiting_atoms"]["total"]["atom_index"])  # type: ignore[index]
        in terminal_tms_group
        for record in records
    )
    direct_tms_effect: Dict[str, object] = {
        "total_clearance_changed_fraction": None,
        "nonlocal_clearance_changed_fraction": None,
        "probe_radius_sensitivity": {},
    }
    if terminal_tms_group:
        for scope, without_scope, output_key in (
            ("total", "total_without_terminal_tms", "total_clearance_changed_fraction"),
            (
                "nonlocal_environment",
                "nonlocal_without_terminal_tms",
                "nonlocal_clearance_changed_fraction",
            ),
        ):
            direct_tms_effect[output_key] = sum(
                float(record["clearance_without_probe_angstrom"][scope])  # type: ignore[index]
                < float(
                    record["clearance_without_probe_angstrom"][without_scope]  # type: ignore[index]
                )
                - 1.0e-12
                for record in records
            ) / len(records)
        for radius in radii:
            key = f"probe_radius_{radius:.2f}_angstrom"
            direct_tms_effect["probe_radius_sensitivity"][key] = {}
            for scope, without_scope in (
                ("total", "total_without_terminal_tms"),
                ("nonlocal_environment", "nonlocal_without_terminal_tms"),
            ):
                with_probe = scopes[scope]["probe_radius_sensitivity"][key]
                without_probe = scopes[without_scope][
                    "probe_radius_sensitivity"
                ][key]
                direct_tms_effect["probe_radius_sensitivity"][key][scope] = {
                    "positive_clearance_fraction_with_minus_without_tms": (
                        float(with_probe["positive_clearance_fraction"])
                        - float(without_probe["positive_clearance_fraction"])
                    ),
                    "positive_clearance_integral_with_minus_without_tms_angstrom": (
                        float(with_probe["positive_clearance_integral_angstrom"])
                        - float(
                            without_probe[
                                "positive_clearance_integral_angstrom"
                            ]
                        )
                    ),
                }
    result = {
        "method": "continuous_vdw_corridor_on_burgi_dunitz_cone",
        "parameters": {
            "burgi_dunitz_o_c_approach_angle_deg": burgi_dunitz_angle_deg,
            "near_distance_from_carbonyl_c_angstrom": near_distance_angstrom,
            "far_distance_from_carbonyl_c_angstrom": far_distance_angstrom,
            "probe_radii_angstrom": list(radii),
            "azimuth_samples": azimuth_samples,
            "vdw_radii_angstrom": dict(_VDW_RADII),
        },
        "atom_indices": {
            "benzoyl_carbonyl_c": c,
            "benzoyl_carbonyl_o_priority_1": o,
            "benzoyl_amide_n_priority_2": n,
            "benzoyl_ipso_c_priority_3": ipso,
            "terminal_tms_heavy": sorted(terminal_tms_heavy),
            "terminal_tms_group_including_h": sorted(terminal_tms_group),
        },
        "face_assignment": {
            "re_si_assigned": True,
            "priority_order": ["carbonyl O", "amide N", "benzoyl ipso C"],
            "convention": (
                "Viewed from the approaching point toward carbonyl C, clockwise "
                "priority 1>2>3 is Re; determinant-negative paths are Re."
            ),
            "coplanar_sample_count": sum(
                record["face"] == "coplanar" for record in records
            ),
        },
        "scopes": scopes,
        "terminal_tms_as_total_limiter_fraction": (
            tms_limiter_count / len(records) if terminal_tms_group else None
        ),
        "direct_terminal_tms_counterfactual": direct_tms_effect,
        "interpretation_limit": (
            "This is a static hard-sphere approach-corridor descriptor. It omits "
            "enzyme geometry, solvent, molecular dynamics, nucleophile-specific "
            "interactions, and activation free energy."
        ),
    }
    if include_trajectory:
        result["trajectory"] = records
    return result


def _nested_float(payload: Mapping[str, object], path: Sequence[str]) -> float:
    value: object = payload
    for key in path:
        if not isinstance(value, Mapping):
            raise KeyError(".".join(path))
        value = value[key]
    return float(value)


def _paired_delta_summary(values: Sequence[float]) -> Dict[str, object]:
    numeric = [float(value) for value in values]
    return {
        "pair_count": len(numeric),
        "mean_sb_minus_bz": sum(numeric) / len(numeric),
        "median_sb_minus_bz": _percentile(numeric, 0.50),
        "minimum_sb_minus_bz": min(numeric),
        "maximum_sb_minus_bz": max(numeric),
        "sb_greater_pair_count": sum(value > 0.0 for value in numeric),
        "bz_greater_pair_count": sum(value < 0.0 for value in numeric),
        "equal_pair_count": sum(abs(value) <= 1.0e-12 for value in numeric),
    }


def _pointwise_delta_summary(values: Sequence[float]) -> Dict[str, object]:
    numeric = [float(value) for value in values]
    absolute = [abs(value) for value in numeric]
    summary = _paired_delta_summary(numeric)
    summary.update(
        {
            "mean_absolute_sb_minus_bz": sum(absolute) / len(absolute),
            "rms_sb_minus_bz": math.sqrt(
                sum(value * value for value in numeric) / len(numeric)
            ),
            "p95_absolute_sb_minus_bz": _percentile(absolute, 0.95),
        }
    )
    return summary


def run_paired_continuous_steric_analysis(
    bz_ensemble_dir: Path,
    sb_ensemble_dir: Path,
    output_dir: Path,
    *,
    maximum_common_scaffold_rmsd_angstrom: float = 0.80,
    maximum_reaction_center_rmsd_angstrom: float = 0.05,
    azimuth_samples: int = 720,
    probe_radii_angstrom: Sequence[float] = (0.0, 0.5, 1.0, 1.4),
) -> Dict[str, object]:
    """Compare continuous carbonyl access only for matched Bz/SB conformers."""

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    bz_payload, bz_conformers = load_existing_conformers(bz_ensemble_dir)
    sb_payload, sb_conformers = load_existing_conformers(sb_ensemble_dir)
    targeted = load_targeted_followup_conformer(sb_payload, output_dir)
    if targeted is not None:
        sb_conformers.append(targeted)
    cross = cross_ensemble_rmsd_analysis(bz_conformers, sb_conformers)
    pairs = [
        record
        for record in cross.get("reciprocal_nearest_pairs", [])
        if float(record["common_scaffold_rmsd_angstrom"])
        <= maximum_common_scaffold_rmsd_angstrom
        and float(record["reaction_center_rmsd_angstrom"])
        <= maximum_reaction_center_rmsd_angstrom
    ]
    bz_by_id = {str(item.entry["conformer_id"]): item for item in bz_conformers}
    sb_by_id = {str(item.entry["conformer_id"]): item for item in sb_conformers}

    metric_paths = {
        "total_mean_clearance_probe_0.50_angstrom": (
            "scopes", "total", "probe_radius_sensitivity",
            "probe_radius_0.50_angstrom", "mean_clearance_angstrom",
        ),
        "total_top_decile_clearance_probe_0.50_angstrom": (
            "scopes", "total", "probe_radius_sensitivity",
            "probe_radius_0.50_angstrom", "top_decile_mean_clearance_angstrom",
        ),
        "total_positive_clearance_integral_probe_0.50_angstrom": (
            "scopes", "total", "probe_radius_sensitivity",
            "probe_radius_0.50_angstrom", "positive_clearance_integral_angstrom",
        ),
        "nonlocal_mean_clearance_probe_0.50_angstrom": (
            "scopes", "nonlocal_environment", "probe_radius_sensitivity",
            "probe_radius_0.50_angstrom", "mean_clearance_angstrom",
        ),
        "nonlocal_top_decile_clearance_probe_0.50_angstrom": (
            "scopes", "nonlocal_environment", "probe_radius_sensitivity",
            "probe_radius_0.50_angstrom", "top_decile_mean_clearance_angstrom",
        ),
        "nonlocal_positive_clearance_integral_probe_0.50_angstrom": (
            "scopes", "nonlocal_environment", "probe_radius_sensitivity",
            "probe_radius_0.50_angstrom", "positive_clearance_integral_angstrom",
        ),
        "re_nonlocal_mean_clearance_probe_0.50_angstrom": (
            "scopes", "nonlocal_environment", "face_summary", "Re",
            "probe_radius_sensitivity", "probe_radius_0.50_angstrom",
            "mean_clearance_angstrom",
        ),
        "si_nonlocal_mean_clearance_probe_0.50_angstrom": (
            "scopes", "nonlocal_environment", "face_summary", "Si",
            "probe_radius_sensitivity", "probe_radius_0.50_angstrom",
            "mean_clearance_angstrom",
        ),
    }
    pair_records = []
    trajectory_rows = []
    for number, pair in enumerate(pairs, start=1):
        bz = bz_by_id[str(pair["first_conformer_id"])]
        sb = sb_by_id[str(pair["second_conformer_id"])]
        bz_access = carbonyl_continuous_steric_access(
            bz.model,
            azimuth_samples=azimuth_samples,
            probe_radii_angstrom=probe_radii_angstrom,
            include_trajectory=True,
        )
        sb_access = carbonyl_continuous_steric_access(
            sb.model,
            azimuth_samples=azimuth_samples,
            probe_radii_angstrom=probe_radii_angstrom,
            include_trajectory=True,
        )
        bz_trajectory = list(bz_access["trajectory"])
        sb_trajectory = list(sb_access["trajectory"])
        metric_values = {}
        for name, path in metric_paths.items():
            bz_value = _nested_float(bz_access, path)
            sb_value = _nested_float(sb_access, path)
            metric_values[name] = {
                "bz": bz_value,
                "sb": sb_value,
                "sb_minus_bz": sb_value - bz_value,
            }
        pointwise_deltas = {}
        for scope in ("total", "nonlocal_environment"):
            values = [
                float(sb_row["clearance_without_probe_angstrom"][scope])
                - float(bz_row["clearance_without_probe_angstrom"][scope])
                for bz_row, sb_row in zip(bz_trajectory, sb_trajectory)
            ]
            pointwise_deltas[scope] = _pointwise_delta_summary(values)
        pair_id = f"pair{number:03d}"
        bz_access.pop("trajectory")
        sb_access.pop("trajectory")
        pair_records.append(
            {
                "pair_id": pair_id,
                **pair,
                "metrics": metric_values,
                "pointwise_azimuth_delta_summary": pointwise_deltas,
                "bz_access_summary": bz_access,
                "sb_access_summary": sb_access,
                "bz_terminal_tms_as_total_limiter_fraction": bz_access[
                    "terminal_tms_as_total_limiter_fraction"
                ],
                "sb_terminal_tms_as_total_limiter_fraction": sb_access[
                    "terminal_tms_as_total_limiter_fraction"
                ],
            }
        )
        for label, conformer, access in (
            ("Bz", bz, bz_trajectory),
            ("SB", sb, sb_trajectory),
        ):
            for row in access:
                trajectory_rows.append(
                    {
                        "pair_id": pair_id,
                        "molecule": label,
                        "conformer_id": conformer.entry["conformer_id"],
                        "azimuth_deg": row["azimuth_deg"],
                        "face": row["face"],
                        "total_clearance_without_probe_angstrom": row[
                            "clearance_without_probe_angstrom"
                        ]["total"],
                        "nonlocal_clearance_without_probe_angstrom": row[
                            "clearance_without_probe_angstrom"
                        ]["nonlocal_environment"],
                        "terminal_tms_clearance_without_probe_angstrom": row[
                            "clearance_without_probe_angstrom"
                        ]["terminal_tms"],
                    }
                )

    aggregate = {
        name: _paired_delta_summary(
            [
                float(record["metrics"][name]["sb_minus_bz"])
                for record in pair_records
            ]
        )
        for name in metric_paths
    }
    probe_sensitivity = {}
    for radius in probe_radii_angstrom:
        key = f"probe_radius_{float(radius):.2f}_angstrom"
        for scope in ("total", "nonlocal_environment"):
            deltas = []
            for record in pair_records:
                path = (
                    "scopes", scope, "probe_radius_sensitivity", key,
                    "positive_clearance_integral_angstrom",
                )
                deltas.append(
                    _nested_float(record["sb_access_summary"], path)
                    - _nested_float(record["bz_access_summary"], path)
                )
            probe_sensitivity[f"{scope}.{key}"] = _paired_delta_summary(deltas)

    result = {
        "schema_version": 1,
        "generated_at_utc": _utc_now(),
        "purpose": "Paired continuous carbonyl steric-access sensitivity analysis",
        "source_ensembles": {
            "bz": str(Path(bz_ensemble_dir).resolve()),
            "sb": str(Path(sb_ensemble_dir).resolve()),
            "bz_name": bz_payload.get("name"),
            "sb_name": sb_payload.get("name"),
        },
        "pair_selection": {
            "method": "reciprocal nearest common-scaffold RMSD",
            "maximum_common_scaffold_rmsd_angstrom": (
                maximum_common_scaffold_rmsd_angstrom
            ),
            "maximum_reaction_center_rmsd_angstrom": (
                maximum_reaction_center_rmsd_angstrom
            ),
            "selected_pair_count": len(pair_records),
            "selected_pairs": [
                {
                    "first_conformer_id": pair["first_conformer_id"],
                    "second_conformer_id": pair["second_conformer_id"],
                    "common_scaffold_rmsd_angstrom": pair[
                        "common_scaffold_rmsd_angstrom"
                    ],
                    "reaction_center_rmsd_angstrom": pair[
                        "reaction_center_rmsd_angstrom"
                    ],
                }
                for pair in pairs
            ],
        },
        "analysis_parameters": {
            "azimuth_samples": azimuth_samples,
            "probe_radii_angstrom": list(probe_radii_angstrom),
            "corridor_near_distance_angstrom": 2.0,
            "corridor_far_distance_angstrom": 8.0,
            "burgi_dunitz_angle_deg": 107.0,
        },
        "pairs": pair_records,
        "paired_delta_summary": aggregate,
        "probe_radius_sensitivity": probe_sensitivity,
        "interpretation_limit": (
            "A paired static-geometry clearance difference can support or oppose "
            "direct steric shielding at this model level, but cannot determine an "
            "enzyme hydrolysis rate or establish a metabolic pathway."
        ),
    }
    _atomic_write_json(output_dir / "paired_continuous_steric_access.json", result)
    if trajectory_rows:
        with (output_dir / "paired_continuous_steric_trajectories.csv").open(
            "w", newline="", encoding="utf-8-sig"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(trajectory_rows[0]))
            writer.writeheader()
            writer.writerows(trajectory_rows)
    return result


def _parse_hess_modes(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    def marker_index(marker: str) -> int:
        for index, line in enumerate(lines):
            if line.strip() == marker:
                return index
        raise ValueError(f"Missing {marker} in {path}")

    frequency_start = marker_index("$vibrational_frequencies")
    dimension = int(lines[frequency_start + 1].strip())
    frequencies = np.asarray(
        [
            float(lines[frequency_start + 2 + index].split()[1])
            for index in range(dimension)
        ],
        dtype=float,
    )
    mode_start = marker_index("$normal_modes")
    rows, columns = (int(value) for value in lines[mode_start + 1].split()[:2])
    if rows != dimension or columns != dimension:
        raise ValueError("ORCA Hessian mode dimensions are inconsistent")
    matrix = np.zeros((rows, columns), dtype=float)
    cursor = mode_start + 2
    while cursor < len(lines) and columns:
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        header = [int(value) for value in lines[cursor].split()]
        cursor += 1
        if not header:
            break
        for expected_row in range(rows):
            fields = lines[cursor].split()
            cursor += 1
            if int(fields[0]) != expected_row or len(fields) != len(header) + 1:
                raise ValueError("Malformed ORCA normal-mode matrix")
            for offset, column in enumerate(header):
                matrix[expected_row, column] = float(fields[offset + 1])
        if header[-1] == columns - 1:
            break
    return frequencies, matrix


def carbonyl_stretch_mode(
    frequency_output_path: Optional[Path],
    model: MoleculeModel,
) -> Dict[str, object]:
    if frequency_output_path is None:
        return {
            "available": False,
            "reason": "frequency_not_calculated",
        }
    hess_path = frequency_output_path.with_suffix(".hess")
    if not hess_path.is_file():
        return {
            "available": False,
            "reason": "orca_hessian_not_retained",
            "frequency_output": str(frequency_output_path),
        }
    subsets = rmsd_atom_subsets(model)
    c = int(subsets.benzoyl_carbonyl_c)  # type: ignore[arg-type]
    o = int(subsets.benzoyl_carbonyl_o)  # type: ignore[arg-type]
    frequencies, modes = _parse_hess_modes(hess_path)
    bond = np.asarray(model.atoms[o].position) - np.asarray(model.atoms[c].position)
    bond /= np.linalg.norm(bond)
    candidates = []
    for mode_index, frequency in enumerate(frequencies):
        if frequency < 500.0:
            continue
        carbon = modes[3 * c : 3 * c + 3, mode_index]
        oxygen = modes[3 * o : 3 * o + 3, mode_index]
        projection = abs(float(np.dot(oxygen - carbon, bond)))
        candidates.append((projection, mode_index, float(frequency)))
    if not candidates:
        return {"available": False, "reason": "no_vibrational_mode_above_500_cm1"}
    projection, mode_index, frequency = max(candidates)
    ranked = sorted(candidates, reverse=True)[:5]
    return {
        "available": True,
        "frequency_cm1": frequency,
        "mode_index": mode_index,
        "relative_bond_stretch_projection": projection,
        "top_projection_candidates": [
            {
                "frequency_cm1": item[2],
                "mode_index": item[1],
                "relative_projection": item[0],
            }
            for item in ranked
        ],
        "method": (
            "Largest absolute projection of the C/O relative normal-mode displacement "
            "onto the optimized benzoyl C=O bond, considering modes above 500 cm^-1."
        ),
        "hessian_path": str(hess_path),
    }


def _bond_order(
    values: Mapping[Tuple[int, int], float], first: int, second: int
) -> Optional[float]:
    return values.get(tuple(sorted((first, second))))


def analyze_existing_conformer(conformer: ExistingConformer) -> Dict[str, object]:
    structure = structural_descriptors(conformer.model)
    atoms = structure["atom_indices"]
    c = int(atoms["benzoyl_carbonyl_c"])
    o = int(atoms["benzoyl_carbonyl_o"])
    n = int(atoms["benzoyl_amide_n"])
    parsed = parse_existing_orca_electronic(conformer.opt_output_path)

    def charge(kind: str, index: int) -> Optional[float]:
        values = parsed.get(kind)
        return values.get(index) if isinstance(values, dict) else None

    mayer = parsed.get("mayer_bond_orders")
    if not isinstance(mayer, dict):
        mayer = {}
    electronic = {
        "mulliken_charge_carbonyl_c": charge("mulliken_charges", c),
        "mulliken_charge_carbonyl_o": charge("mulliken_charges", o),
        "mulliken_charge_amide_n": charge("mulliken_charges", n),
        "loewdin_charge_carbonyl_c": charge("loewdin_charges", c),
        "loewdin_charge_carbonyl_o": charge("loewdin_charges", o),
        "loewdin_charge_amide_n": charge("loewdin_charges", n),
        "mbis_charge_carbonyl_c": charge("mbis_charges", c),
        "mbis_charge_carbonyl_o": charge("mbis_charges", o),
        "mbis_charge_amide_n": charge("mbis_charges", n),
        "chelpg_charge_carbonyl_c": charge("chelpg_charges", c),
        "chelpg_charge_carbonyl_o": charge("chelpg_charges", o),
        "chelpg_charge_amide_n": charge("chelpg_charges", n),
        "mayer_bond_order_c_o": _bond_order(mayer, c, o),
        "mayer_bond_order_c_n": _bond_order(mayer, c, n),
        "dipole_magnitude_debye": parsed["dipole_magnitude_debye"],
        "homo_energy_ev": (
            parsed["homo"]["energy_ev"] if parsed["homo"] is not None else None
        ),
        "lumo_energy_ev": (
            parsed["lumo"]["energy_ev"] if parsed["lumo"] is not None else None
        ),
        "homo_lumo_gap_ev": parsed["homo_lumo_gap_ev"],
        "frontier_orbital_local_contributions": parsed[
            "frontier_orbital_local_contributions"
        ],
        "frontier_orbital_local_contributions_reason": parsed[
            "frontier_orbital_local_contributions_reason"
        ],
        "charge_model_note": (
            "MBIS/CHELPG are null when they were not requested in the retained Opt job. "
            "Mulliken and Loewdin values are extracted but are not treated as observables."
        ),
    }
    return {
        "conformer_id": conformer.entry["conformer_id"],
        "conformer_index": conformer.entry["conformer_index"],
        "dft_energy_hartree": conformer.entry["dft_energy_hartree"],
        "relative_energy_kj_mol": conformer.entry["relative_energy_kj_mol"],
        "gibbs_energy_hartree": conformer.entry.get("gibbs_energy_hartree"),
        "relative_gibbs_energy_kj_mol": conformer.entry.get(
            "relative_gibbs_energy_kj_mol"
        ),
        "frequency_calculated": conformer.entry.get("frequency_calculated", False),
        "frequency_provenance": conformer.entry.get("frequency_provenance"),
        "imaginary_modes": conformer.entry.get("imaginary_modes"),
        "imaginary_frequencies_cm1": conformer.entry.get(
            "imaginary_frequencies_cm1", []
        ),
        "low_frequency_modes_cm1": conformer.entry.get(
            "low_frequency_modes_cm1", []
        ),
        "local_minimum_assessment": conformer.entry.get(
            "local_minimum_assessment", "not_evaluated"
        ),
        "optimization_provenance": conformer.entry.get(
            "optimization_provenance", conformer.entry.get("source")
        ),
        "all_heavy_rmsd_cluster_id": conformer.entry.get(
            "all_heavy_rmsd_cluster_id", conformer.entry.get("rmsd_cluster_id")
        ),
        "common_scaffold_rmsd_cluster_id": conformer.entry.get(
            "common_scaffold_rmsd_cluster_id"
        ),
        "reaction_center_rmsd_cluster_id": conformer.entry.get(
            "reaction_center_rmsd_cluster_id"
        ),
        "calculation_directory": conformer.entry.get("calculation_directory"),
        "structure": structure,
        "steric_access": carbonyl_steric_access(conformer.model),
        "electronic": electronic,
        "carbonyl_stretch": carbonyl_stretch_mode(
            conformer.frequency_output_path, conformer.model
        ),
    }


def merge_property_single_points(
    molecule_name: str,
    records: Sequence[Dict[str, object]],
    property_summary_path: Path,
) -> Dict[str, object]:
    """Attach optional property-only results without changing Opt energies."""

    if not property_summary_path.is_file():
        return {
            "available": False,
            "matched_conformer_count": 0,
            "reason": "No completed property-only summary was found.",
        }
    payload = _load_json(property_summary_path)
    by_id = {str(record["conformer_id"]): record for record in records}
    matched = []
    for item in payload.get("results", []):
        if not isinstance(item, dict) or item.get("molecule") != molecule_name:
            continue
        identifier = str(item.get("conformer_id"))
        record = by_id.get(identifier)
        output = item.get("output")
        if record is None or not isinstance(output, str) or not Path(output).is_file():
            continue
        parsed = parse_existing_orca_electronic(Path(output))
        structure = record.get("structure", {})
        atom_indices = (
            structure.get("atom_indices", {})
            if isinstance(structure, dict)
            else {}
        )
        if not isinstance(atom_indices, dict):
            continue
        carbonyl_c = int(atom_indices["benzoyl_carbonyl_c"])
        carbonyl_o = int(atom_indices["benzoyl_carbonyl_o"])
        amide_n = int(atom_indices["benzoyl_amide_n"])
        electronic = record.get("electronic", {})
        if not isinstance(electronic, dict):
            continue

        def value(kind: str, atom_index: int) -> Optional[float]:
            values = parsed.get(kind)
            return values.get(atom_index) if isinstance(values, dict) else None

        electronic.update(
            {
                "mbis_charge_carbonyl_c": value("mbis_charges", carbonyl_c),
                "mbis_charge_carbonyl_o": value("mbis_charges", carbonyl_o),
                "mbis_charge_amide_n": value("mbis_charges", amide_n),
                "chelpg_charge_carbonyl_c": value("chelpg_charges", carbonyl_c),
                "chelpg_charge_carbonyl_o": value("chelpg_charges", carbonyl_o),
                "chelpg_charge_amide_n": value("chelpg_charges", amide_n),
                "frontier_orbital_local_contributions": parsed.get(
                    "frontier_orbital_local_contributions"
                ),
                "frontier_orbital_local_contributions_reason": parsed.get(
                    "frontier_orbital_local_contributions_reason"
                ),
                "property_single_point_output": str(Path(output).resolve()),
                "property_single_point_provenance": item.get("provenance"),
                "charge_model_note": (
                    "MBIS/CHELPG were evaluated in a property-only single point on "
                    "the retained optimized geometry. Mulliken/Loewdin are retained "
                    "for reference but are not treated as observables."
                ),
            }
        )
        frontier = parsed.get("frontier_orbital_local_contributions")
        populations = (
            frontier.get("atom_populations", {})
            if isinstance(frontier, dict)
            else {}
        )
        if isinstance(populations, dict):
            region = (carbonyl_c, carbonyl_o, amide_n, int(atom_indices["benzoyl_ipso_c"]))
            electronic["homo_loewdin_population_benzoyl_center"] = sum(
                float(populations[index]["homo_loewdin"])
                for index in region
                if index in populations
            )
            electronic["lumo_loewdin_population_benzoyl_center"] = sum(
                float(populations[index]["lumo_loewdin"])
                for index in region
                if index in populations
            )
        matched.append(identifier)
    return {
        "available": bool(matched),
        "matched_conformer_count": len(matched),
        "matched_conformer_ids": matched,
        "selection": payload.get("selection"),
        "calculation": payload.get("calculation"),
        "geometry_optimization_performed": payload.get(
            "geometry_optimization_performed"
        ),
    }


def _frequency_entry_fields(
    path: Path,
    *,
    imaginary_threshold_cm1: float,
    low_frequency_threshold_cm1: float,
) -> Dict[str, object]:
    parsed = parse_orca_output(
        path.read_text(encoding="utf-8", errors="replace"),
        frequency_requested=True,
        require_geometry_convergence=False,
        imaginary_threshold_cm1=imaginary_threshold_cm1,
    )
    frequencies = list(parsed.frequencies_cm1)
    thermo = parsed.thermochemistry
    return {
        "frequency_calculated": bool(frequencies),
        "number_of_modes": len(frequencies) if frequencies else None,
        "lowest_frequencies_cm1": sorted(frequencies)[:10] if frequencies else [],
        "imaginary_modes": (
            len(parsed.imaginary_frequencies_cm1) if frequencies else None
        ),
        "imaginary_frequencies_cm1": list(parsed.imaginary_frequencies_cm1),
        "translation_rotation_near_zero_modes_cm1": [
            value for value in frequencies if abs(value) < 1.0
        ],
        "low_frequency_modes_cm1": [
            value
            for value in frequencies
            if abs(value) >= 1.0 and value < low_frequency_threshold_cm1
        ],
        "local_minimum_assessment": (
            "not_evaluated"
            if not frequencies
            else (
                "not_a_confirmed_local_minimum"
                if parsed.imaginary_frequencies_cm1
                else "local_minimum_candidate"
            )
        ),
        "gibbs_energy_hartree": (
            thermo.gibbs_free_energy_hartree if thermo is not None else None
        ),
        "thermochemistry": thermo.to_metadata() if thermo is not None else None,
        "frequency_output": str(path),
    }


def _recorded_pool_index(
    ensemble_dir: Path,
    conformer_id: str,
    calculation_dir: Optional[Path] = None,
) -> Optional[int]:
    matches = list(
        (ensemble_dir / "structures").glob(f"*_dft_{conformer_id}.metadata.json")
    )
    if len(matches) == 1:
        data = _load_json(matches[0])
        metadata = data.get("metadata", {})
        if isinstance(metadata, dict):
            value = metadata.get("conformer_pool_index")
            if isinstance(value, int):
                return int(value)

    if calculation_dir is None:
        return None
    calculation_dir = Path(calculation_dir).resolve()
    match = re.fullmatch(r"conf(\d{3,})", conformer_id)
    if match is None:
        return None
    conformer_index = int(match.group(1)) - 1
    calculation_root = calculation_dir.parent
    metadata_candidates = list(calculation_dir.glob("*.metadata.json"))
    metadata_candidates.extend(calculation_root.parent.glob("*.metadata.json"))
    for path in sorted(set(metadata_candidates)):
        data = _load_json(path)
        metadata = data.get("metadata", {})
        if not isinstance(metadata, dict):
            continue
        if metadata.get("conformer_index") != conformer_index:
            continue
        recorded_root = metadata.get("calculation_files_directory")
        if (
            isinstance(recorded_root, str)
            and Path(recorded_root).expanduser().resolve() != calculation_root
        ):
            continue
        value = metadata.get("conformer_pool_index")
        if isinstance(value, int):
            return int(value)
    return None


def _successful_selection_repair_pool_indices(ensemble_dir: Path) -> set[int]:
    """Return pool members added by a completed strict-selection repair."""

    path = Path(ensemble_dir) / "strict_selection_repair_plan.json"
    if not path.is_file():
        return set()
    data = _load_json(path)
    if str(data.get("status", "")).upper() != "SUCCESS":
        return set()
    target = data.get("target", {})
    result = data.get("result", {})
    if not isinstance(target, dict) or not isinstance(result, dict):
        return set()
    pool_index = target.get("pool_index")
    if not isinstance(pool_index, int):
        return set()
    if result.get("geometry_optimization_converged") is not True:
        return set()
    return {int(pool_index)}


def refresh_existing_ensemble_metadata(
    ensemble_dir: Path,
    *,
    write: bool = True,
) -> Tuple[Dict[str, object], List[ExistingConformer], Dict[str, object]]:
    """Repair provenance/summary semantics and attach retained Freq results."""

    ensemble_dir = Path(ensemble_dir).resolve()
    payload, conformers = load_existing_conformers(ensemble_dir)
    frequency = payload.get("frequency", {})
    if not isinstance(frequency, dict):
        frequency = {}
        payload["frequency"] = frequency
    imaginary_threshold = float(frequency.get("imaginary_threshold_cm1", -20.0))
    low_threshold = float(frequency.get("low_frequency_threshold_cm1", 50.0))
    frequency["translation_rotation_near_zero_threshold_abs_cm1"] = 1.0
    frequency["low_frequency_modes_exclude_near_zero_translation_rotation"] = True
    selection = frequency.get("selection", {})
    if not isinstance(selection, dict):
        selection = {}
        frequency["selection"] = selection
    originally_frequency = {
        int(entry["conformer_index"])
        for entry in payload.get("final_conformer_ensemble", [])
        if isinstance(entry, dict) and entry.get("frequency_calculated") is True
    }
    selected_original = {
        int(value) for value in selection.get("selected_conformer_indices", [])
    }
    supplemental = []
    rmsd_analysis = secondary_rmsd_analysis(conformers)
    rmsd_records = {
        int(record["conformer_index"]): record
        for record in rmsd_analysis.get("records", [])
    }
    strict_selection_repairs = _successful_selection_repair_pool_indices(
        ensemble_dir
    )

    entries = payload.get("final_conformer_ensemble", [])
    by_index = {int(item.entry["conformer_index"]): item for item in conformers}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        index = int(entry["conformer_index"])
        conformer = by_index[index]
        was_frequency = entry.get("frequency_calculated") is True
        retained_frequency_provenance = entry.get("frequency_provenance")
        if conformer.frequency_output_path is not None:
            entry.update(
                _frequency_entry_fields(
                    conformer.frequency_output_path,
                    imaginary_threshold_cm1=imaginary_threshold,
                    low_frequency_threshold_cm1=low_threshold,
                )
            )
            if not was_frequency:
                supplemental.append(index)
            if isinstance(retained_frequency_provenance, str) and retained_frequency_provenance:
                frequency_provenance = retained_frequency_provenance
            elif index in originally_frequency:
                frequency_provenance = (
                    "reused_external_read_only"
                    if not conformer.frequency_output_path.is_relative_to(ensemble_dir)
                    else "computed_original_ensemble_run"
                )
            else:
                frequency_provenance = "computed_supplemental_frequency_only"
            entry["frequency_provenance"] = frequency_provenance
        else:
            entry.update(
                {
                    "frequency_calculated": False,
                    "frequency_provenance": None,
                    "number_of_modes": None,
                    "imaginary_modes": None,
                    "imaginary_frequencies_cm1": [],
                    "translation_rotation_near_zero_modes_cm1": [],
                    "lowest_frequencies_cm1": [],
                    "low_frequency_modes_cm1": [],
                    "local_minimum_assessment": "not_evaluated",
                    "gibbs_energy_hartree": None,
                    "relative_gibbs_energy_kj_mol": None,
                    "thermochemistry": None,
                }
            )
        calculation = Path(str(entry["calculation_directory"])).resolve()
        reused = not calculation.is_relative_to(ensemble_dir / "conformers")
        pool_index = _recorded_pool_index(
            ensemble_dir,
            str(entry["conformer_id"]),
            calculation,
        )
        if pool_index in strict_selection_repairs:
            entry["optimization_provenance"] = "computed_strict_selection_repair"
            entry["source"] = "generated"
        else:
            entry["optimization_provenance"] = (
                "reused_external_read_only"
                if reused
                else "computed_original_ensemble_run"
            )
            entry["source"] = "reused" if reused else "generated"
        entry["conformer_pool_index"] = pool_index
        entry["conformer_pool_index_provenance"] = (
            "recorded_from_reuse_metadata"
            if reused and pool_index is not None
            else (
                "recorded_from_current_etkdg_pool"
                if pool_index is not None
                else (
                    "unavailable_for_legacy_reuse"
                    if reused
                    else "unavailable_in_retained_metadata"
                )
            )
        )
        rmsd = rmsd_records[index]
        entry.update(
            {
                "all_heavy_rmsd_cluster_id": rmsd[
                    "all_heavy_rmsd_cluster_id"
                ],
                "common_scaffold_rmsd_cluster_id": rmsd[
                    "common_scaffold_rmsd_cluster_id"
                ],
                "common_scaffold_rmsd_to_representative_angstrom": rmsd[
                    "common_scaffold_rmsd_to_representative_angstrom"
                ],
                "reaction_center_rmsd_cluster_id": rmsd[
                    "reaction_center_rmsd_cluster_id"
                ],
                "reaction_center_rmsd_to_representative_angstrom": rmsd[
                    "reaction_center_rmsd_to_representative_angstrom"
                ],
            }
        )

    gibbs_values = [
        float(entry["gibbs_energy_hartree"])
        for entry in entries
        if isinstance(entry, dict) and entry.get("gibbs_energy_hartree") is not None
    ]
    gibbs_minimum = min(gibbs_values) if gibbs_values else None
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        gibbs = entry.get("gibbs_energy_hartree")
        entry["relative_gibbs_energy_kj_mol"] = (
            (float(gibbs) - gibbs_minimum) * HARTREE_TO_KJ_MOL
            if gibbs is not None and gibbs_minimum is not None
            else None
        )

    dft = payload.get("dft", {})
    if isinstance(dft, dict):
        legacy_failures = dft.get("failures", [])
        if not isinstance(legacy_failures, list):
            legacy_failures = []
        final_failures = [
            item
            for item in legacy_failures
            if isinstance(item, dict)
            and item.get("phase") in {"dft_optimization", "post_dft_validation"}
        ]
        recovery = [
            item
            for item in legacy_failures
            if isinstance(item, dict) and item.get("phase") == "resume_recovery"
        ]
        recovery.extend(
            item
            for item in dft.get("recovery_history", [])
            if isinstance(item, dict) and item not in recovery
        )
        dft["failures"] = final_failures
        dft["recovery_history"] = recovery
        post = dft.get("post_optimization_screening")
        if isinstance(post, dict):
            post["secondary_rmsd_analysis"] = rmsd_analysis

    available = {
        int(entry["conformer_index"])
        for entry in entries
        if isinstance(entry, dict) and entry.get("frequency_calculated") is True
    }
    selection["supplemental_completed_conformer_indices"] = sorted(
        set(selection.get("supplemental_completed_conformer_indices", []))
        | set(supplemental)
    )
    summary = payload.get("summary", {})
    if isinstance(summary, dict):
        summary.update(
            {
                "frequency_candidates": len(selected_original),
                "frequency_completed": len(selected_original & available),
                "frequency_selected_this_run": len(selected_original),
                "frequency_completed_for_this_run_selection": len(
                    selected_original & available
                ),
                "frequency_available_total": len(available),
                "frequency_preexisting_or_reused_total": sum(
                    1
                    for entry in entries
                    if isinstance(entry, dict)
                    and entry.get("frequency_calculated") is True
                    and entry.get("frequency_provenance")
                    != "computed_supplemental_frequency_only"
                ),
                "frequency_computed_this_run": sum(
                    1
                    for entry in entries
                    if isinstance(entry, dict)
                    and entry.get("frequency_calculated") is True
                    and entry.get("frequency_provenance")
                    == "computed_supplemental_frequency_only"
                ),
                "frequency_computed_supplemental_analysis_run": sum(
                    1
                    for entry in entries
                    if isinstance(entry, dict)
                    and entry.get("frequency_provenance")
                    == "computed_supplemental_frequency_only"
                ),
                "frequency_preexisting_local_original_ensemble": sum(
                    1
                    for entry in entries
                    if isinstance(entry, dict)
                    and entry.get("frequency_provenance")
                    == "computed_original_ensemble_run"
                ),
                "frequency_reused_external_read_only": sum(
                    1
                    for entry in entries
                    if isinstance(entry, dict)
                    and entry.get("frequency_provenance")
                    == "reused_external_read_only"
                ),
                "frequency_supplemental_completed": len(
                    selection["supplemental_completed_conformer_indices"]
                ),
                "common_scaffold_dft_clusters": rmsd_analysis["common_scaffold"][
                    "cluster_count"
                ],
                "reaction_center_dft_clusters": rmsd_analysis["reaction_center"][
                    "cluster_count"
                ],
            }
        )
    payload["schema_version"] = max(int(payload.get("schema_version", 1)), 2)
    revisions = payload.setdefault("analysis_revisions", [])
    if not isinstance(revisions, list):
        revisions = []
        payload["analysis_revisions"] = revisions
    revision = {
        "timestamp_utc": _utc_now(),
        "type": "preliminary_bz_sb_comparison_metadata_repair",
        "changes": [
            "null imaginary_modes for structures without Freq",
            "separated final failures from resume recovery history",
            "separated optimization and frequency provenance",
            "added common-scaffold and reaction-centre RMSD analyses",
            "attached non-destructive supplemental frequency results when present",
            "synchronized conformers with final_conformer_ensemble",
            "preserved successful strict-selection repair provenance",
        ],
    }
    if not any(
        isinstance(item, dict) and item.get("type") == revision["type"]
        for item in revisions
    ):
        revisions.append(revision)
    _synchronize_ensemble_conformer_aliases(payload)
    if write:
        path = ensemble_dir / "ensemble.json"
        backup = ensemble_dir / "ensemble.before_preliminary_repair.json"
        if not backup.exists():
            backup.write_bytes(path.read_bytes())
        _atomic_write_json(path, payload)
        write_run_summary(payload, ensemble_dir)
    return payload, conformers, revision


def _nested_value(record: Mapping[str, object], path: str) -> Optional[float]:
    value: object = record
    for key in path.split("."):
        if not isinstance(value, Mapping) or key not in value:
            return None
        value = value[key]
    return float(value) if isinstance(value, (int, float)) else None


def _boltzmann_weights(
    records: Sequence[Mapping[str, object]],
    energy_path: str,
    *,
    temperature_kelvin: float = DEFAULT_TEMPERATURE_K,
) -> Dict[str, float]:
    available = [
        (str(record["conformer_id"]), _nested_value(record, energy_path))
        for record in records
    ]
    available = [(key, value) for key, value in available if value is not None]
    if not available:
        return {}
    minimum = min(float(value) for _key, value in available)
    if energy_path.endswith("hartree"):
        deltas = [
            (key, (float(value) - minimum) * HARTREE_TO_KJ_MOL)
            for key, value in available
        ]
    else:
        deltas = [(key, float(value) - minimum) for key, value in available]
    factors = [
        (key, math.exp(-delta / (R_KJ_MOL_K * temperature_kelvin)))
        for key, delta in deltas
    ]
    total = sum(value for _key, value in factors)
    return {key: value / total for key, value in factors}


def _weighted_summary(
    records: Sequence[Mapping[str, object]],
    weights: Mapping[str, float],
    metric_path: str,
) -> Dict[str, object]:
    values = []
    for record in records:
        identifier = str(record["conformer_id"])
        value = _nested_value(record, metric_path)
        if value is not None and identifier in weights:
            values.append((identifier, value, float(weights[identifier])))
    if not values:
        return {
            "weighted_mean": None,
            "weighted_std": None,
            "minimum": None,
            "maximum": None,
            "conformer_count": 0,
            "weight_coverage": 0.0,
        }
    coverage = sum(weight for _identifier, _value, weight in values)
    normalized = [
        (identifier, value, weight / coverage)
        for identifier, value, weight in values
    ]
    mean = sum(value * weight for _identifier, value, weight in normalized)
    variance = sum(
        weight * (value - mean) ** 2
        for _identifier, value, weight in normalized
    )
    return {
        "weighted_mean": mean,
        "weighted_std": math.sqrt(max(variance, 0.0)),
        "minimum": min(value for _identifier, value, _weight in values),
        "maximum": max(value for _identifier, value, _weight in values),
        "conformer_count": len(values),
        "weight_coverage": coverage,
    }


_COMPARISON_METRICS = (
    "structure.carbonyl_c_o_length_angstrom",
    "structure.amide_c_n_length_angstrom",
    "structure.benzoyl_carbonyl_torsion_abs_mean_deg",
    "structure.amide_n_planarity_torsion_abs_mean_deg",
    "structure.aryl_c_si_length_angstrom",
    "steric_access.overall_accessible_fraction",
    "steric_access.best_clearance_angstrom",
    "steric_access.mean_clearance_angstrom",
    "steric_access.clearance_percentiles_angstrom.p90",
    "steric_access.minimum_terminal_tms_distance_from_carbonyl_c_angstrom",
    "steric_access.silicon_distance_from_carbonyl_c_angstrom",
    "steric_access.terminal_tms_heavy_atoms_within_5_angstrom",
    "steric_access.face_summary.face_positive.accessible_fraction",
    "steric_access.face_summary.face_negative.accessible_fraction",
    "electronic.mulliken_charge_carbonyl_c",
    "electronic.loewdin_charge_carbonyl_c",
    "electronic.loewdin_charge_carbonyl_o",
    "electronic.mbis_charge_carbonyl_c",
    "electronic.mbis_charge_carbonyl_o",
    "electronic.chelpg_charge_carbonyl_c",
    "electronic.chelpg_charge_carbonyl_o",
    "electronic.mayer_bond_order_c_o",
    "electronic.mayer_bond_order_c_n",
    "electronic.dipole_magnitude_debye",
    "electronic.homo_energy_ev",
    "electronic.lumo_energy_ev",
    "electronic.homo_lumo_gap_ev",
    "electronic.homo_loewdin_population_benzoyl_center",
    "electronic.lumo_loewdin_population_benzoyl_center",
    "carbonyl_stretch.frequency_cm1",
)


def ensemble_descriptor_summary(
    payload: Mapping[str, object],
    conformer_metrics: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    electronic_weights = _boltzmann_weights(
        conformer_metrics, "dft_energy_hartree"
    )
    gibbs_weights = _boltzmann_weights(
        conformer_metrics, "gibbs_energy_hartree"
    )
    common_clusters: Dict[str, List[Dict[str, object]]] = {}
    for record in conformer_metrics:
        cluster = str(
            record.get("common_scaffold_rmsd_cluster_id")
            or record["conformer_id"]
        )
        common_clusters.setdefault(cluster, []).append(record)
    common_unique = [
        min(items, key=lambda item: float(item["dft_energy_hartree"]))
        for items in common_clusters.values()
    ]
    common_unique_electronic_weights = _boltzmann_weights(
        common_unique, "dft_energy_hartree"
    )
    common_unique_gibbs_candidates = []
    for items in common_clusters.values():
        available = [
            item for item in items if item.get("gibbs_energy_hartree") is not None
        ]
        if available:
            common_unique_gibbs_candidates.append(
                min(available, key=lambda item: float(item["gibbs_energy_hartree"]))
            )
    common_unique_gibbs_weights = _boltzmann_weights(
        common_unique_gibbs_candidates, "gibbs_energy_hartree"
    )
    for record in conformer_metrics:
        identifier = str(record["conformer_id"])
        record["electronic_boltzmann_weight"] = electronic_weights.get(identifier)
        record["conditional_gibbs_boltzmann_weight"] = gibbs_weights.get(identifier)
        record["common_scaffold_unique_electronic_weight"] = (
            common_unique_electronic_weights.get(identifier)
        )
        record["common_scaffold_unique_conditional_gibbs_weight"] = (
            common_unique_gibbs_weights.get(identifier)
        )
    entries = payload.get("final_conformer_ensemble", [])
    final_count = len(entries) if isinstance(entries, list) else len(conformer_metrics)
    frequency_count = sum(
        1 for record in conformer_metrics if record.get("frequency_calculated") is True
    )
    payload_summary = payload.get("summary", {})
    if not isinstance(payload_summary, Mapping):
        payload_summary = {}
    metrics = {}
    common_unique_metrics = {}
    for metric in _COMPARISON_METRICS:
        metrics[metric] = {
            "electronic_energy_weighted": _weighted_summary(
                conformer_metrics, electronic_weights, metric
            ),
            "conditional_gibbs_weighted": _weighted_summary(
                conformer_metrics, gibbs_weights, metric
            ),
        }
        common_unique_metrics[metric] = {
            "electronic_energy_weighted": _weighted_summary(
                common_unique, common_unique_electronic_weights, metric
            ),
            "conditional_gibbs_weighted": _weighted_summary(
                common_unique_gibbs_candidates,
                common_unique_gibbs_weights,
                metric,
            ),
        }
    return {
        "name": payload.get("name"),
        "temperature_kelvin": DEFAULT_TEMPERATURE_K,
        "dft_representative_count": final_count,
        "frequency_completed_count": frequency_count,
        "frequency_bookkeeping": {
            "original_selection_candidates": payload_summary.get(
                "frequency_candidates"
            ),
            "original_selection_completed": payload_summary.get(
                "frequency_completed"
            ),
            "available_total": payload_summary.get(
                "frequency_available_total", frequency_count
            ),
            "supplemental_completed": payload_summary.get(
                "frequency_supplemental_completed", 0
            ),
            "preexisting_or_reused_total": payload_summary.get(
                "frequency_preexisting_or_reused_total"
            ),
        },
        "electronic_weight_scope": (
            "All final DFT representatives in this run, weighted by electronic energy; "
            "not a complete molecular conformer ensemble."
        ),
        "gibbs_weight_scope": (
            "All final DFT representatives"
            if frequency_count == final_count
            else "Conditional on conformers with completed frequency calculations"
        ),
        "gibbs_weight_is_complete_for_final_dft_set": frequency_count == final_count,
        "electronic_weights": electronic_weights,
        "conditional_gibbs_weights": gibbs_weights,
        "metrics": metrics,
        "common_scaffold_unique": {
            "representative_count": len(common_unique),
            "representative_conformer_ids": [
                record["conformer_id"] for record in common_unique
            ],
            "electronic_weights": common_unique_electronic_weights,
            "conditional_gibbs_weights": common_unique_gibbs_weights,
            "metrics": common_unique_metrics,
            "scope": (
                "Lowest-electronic-energy representative per common-scaffold RMSD "
                "cluster. This avoids counting TMS-only rotations as independent "
                "scaffold conformers, but does not supply rigorous rotor degeneracy "
                "or configurational entropy."
            ),
        },
    }


def _flatten_metric_row(record: Mapping[str, object]) -> Dict[str, object]:
    row: Dict[str, object] = {
        "conformer_id": record["conformer_id"],
        "conformer_index": record["conformer_index"],
        "dft_energy_hartree": record["dft_energy_hartree"],
        "relative_energy_kj_mol": record["relative_energy_kj_mol"],
        "gibbs_energy_hartree": record.get("gibbs_energy_hartree"),
        "relative_gibbs_energy_kj_mol": record.get(
            "relative_gibbs_energy_kj_mol"
        ),
        "frequency_calculated": record["frequency_calculated"],
        "electronic_boltzmann_weight": record.get("electronic_boltzmann_weight"),
        "conditional_gibbs_boltzmann_weight": record.get(
            "conditional_gibbs_boltzmann_weight"
        ),
        "common_scaffold_unique_electronic_weight": record.get(
            "common_scaffold_unique_electronic_weight"
        ),
        "common_scaffold_unique_conditional_gibbs_weight": record.get(
            "common_scaffold_unique_conditional_gibbs_weight"
        ),
    }
    for metric in _COMPARISON_METRICS:
        row[metric] = _nested_value(record, metric)
    return row


def _write_metric_csv(path: Path, records: Sequence[Mapping[str, object]]) -> None:
    rows = [_flatten_metric_row(record) for record in records]
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _mean(
    summary: Mapping[str, object],
    metric: str,
    weighting: str,
    *,
    common_scaffold_unique: bool = False,
) -> Optional[float]:
    source: Mapping[str, object] = summary
    if common_scaffold_unique:
        candidate = summary.get("common_scaffold_unique", {})
        if isinstance(candidate, Mapping):
            source = candidate
    metrics = source.get("metrics", {})
    if not isinstance(metrics, Mapping):
        return None
    item = metrics.get(metric, {})
    if not isinstance(item, Mapping):
        return None
    value = item.get(weighting, {})
    if not isinstance(value, Mapping):
        return None
    mean = value.get("weighted_mean")
    return float(mean) if isinstance(mean, (int, float)) else None


def _format_number(value: Optional[float], digits: int = 5) -> str:
    return "not available" if value is None else f"{value:.{digits}f}"


def write_preliminary_comparison_report(
    output_path: Path,
    bz_summary: Mapping[str, object],
    sb_summary: Mapping[str, object],
    bz_records: Sequence[Mapping[str, object]],
    sb_records: Sequence[Mapping[str, object]],
    *,
    sb_pool_review: Optional[Mapping[str, object]] = None,
) -> None:
    weighting = "electronic_energy_weighted"
    chosen = (
        ("C=O length (A)", "structure.carbonyl_c_o_length_angstrom", 5),
        ("Amide C-N length (A)", "structure.amide_c_n_length_angstrom", 5),
        (
            "Benzoyl-carbonyl torsion (deg)",
            "structure.benzoyl_carbonyl_torsion_abs_mean_deg",
            3,
        ),
        (
            "BD-cone accessible fraction",
            "steric_access.overall_accessible_fraction",
            4,
        ),
        (
            "Best BD-cone clearance (A)",
            "steric_access.best_clearance_angstrom",
            4,
        ),
        (
            "MBIS carbonyl-C charge",
            "electronic.mbis_charge_carbonyl_c",
            5,
        ),
        ("MBIS carbonyl-O charge", "electronic.mbis_charge_carbonyl_o", 5),
        (
            "CHELPG carbonyl-C charge",
            "electronic.chelpg_charge_carbonyl_c",
            5,
        ),
        ("CHELPG carbonyl-O charge", "electronic.chelpg_charge_carbonyl_o", 5),
        ("Mayer C=O bond order", "electronic.mayer_bond_order_c_o", 5),
        ("Mayer C-N bond order", "electronic.mayer_bond_order_c_n", 5),
        ("Dipole magnitude (D)", "electronic.dipole_magnitude_debye", 4),
        (
            "LUMO Loewdin population at benzoyl center",
            "electronic.lumo_loewdin_population_benzoyl_center",
            5,
        ),
        ("Carbonyl stretch (cm-1)", "carbonyl_stretch.frequency_cm1", 2),
    )
    table = ["| Descriptor | 1Bz-LSD_RR | 1SB-LSD_RR | SB - Bz |", "| --- | ---: | ---: | ---: |"]
    for label, metric, digits in chosen:
        bz = _mean(
            bz_summary, metric, weighting, common_scaffold_unique=True
        )
        sb = _mean(
            sb_summary, metric, weighting, common_scaffold_unique=True
        )
        difference = sb - bz if bz is not None and sb is not None else None
        table.append(
            f"| {label} | {_format_number(bz, digits)} | {_format_number(sb, digits)} | "
            f"{_format_number(difference, digits)} |"
        )

    bz_freq = int(bz_summary["frequency_completed_count"])
    sb_freq = int(sb_summary["frequency_completed_count"])
    bz_count = int(bz_summary["dft_representative_count"])
    sb_count = int(sb_summary["dft_representative_count"])
    pool_text = "not performed"
    if sb_pool_review is not None:
        pool_text = (
            f"{sb_pool_review.get('common_scaffold_cluster_count')} common-scaffold clusters "
            f"within the force-field window; {len(sb_pool_review.get('unrepresented_candidates', []))} "
            "clusters lacked a selected DFT representative"
        )
    bz_frequency = bz_summary.get("frequency_bookkeeping", {})
    sb_frequency = sb_summary.get("frequency_bookkeeping", {})
    bz_frequency = bz_frequency if isinstance(bz_frequency, Mapping) else {}
    sb_frequency = sb_frequency if isinstance(sb_frequency, Mapping) else {}
    lines = [
        "# 1Bz-LSD_RR vs 1SB-LSD_RR preliminary comparison",
        "",
        "Generated from retained B3LYP/def2-SVP gas-phase structures. These are computed candidate structures, not experimental structures.",
        "",
        "## Scope",
        "",
        f"- Bz: {bz_count} final DFT representatives; original Freq selection "
        f"{bz_frequency.get('original_selection_completed')}/"
        f"{bz_frequency.get('original_selection_candidates')} completed, "
        f"{bz_frequency.get('supplemental_completed')} supplemental, {bz_freq} total available.",
        f"- SB sensitivity set: {sb_count} DFT representatives including one targeted "
        f"follow-up; original Freq selection {sb_frequency.get('original_selection_completed')}/"
        f"{sb_frequency.get('original_selection_candidates')} completed, "
        f"{sb_frequency.get('supplemental_completed')} supplemental, {sb_freq} total available.",
        f"- SB force-field pool review: {pool_text}.",
        "- Values in the main table use one lowest-energy representative per common-scaffold cluster, then electronic-energy weighting. This is not claimed to be a complete Boltzmann ensemble or a rotor-degeneracy treatment.",
        "- Gibbs-weighted values in analysis.json are conditional on structures with Freq unless every final representative has Freq.",
        "",
        "## Electronically weighted descriptor comparison",
        "",
        *table,
        "",
        "## Evidence supporting hypothesis 1 (steric relief from the longer aryl-Si bond)",
        "",
        "- SB retains essentially the same sampled directional carbonyl access as Bz. The terminal para-TMS atoms remain remote from the carbonyl center. This is consistent with absence of direct static occlusion at the carbonyl, but it does not prove that the longer Si-C bond is the cause.",
        "",
        "## Evidence against hypothesis 1",
        "",
        "- SB does not show an increased best clearance; the mean difference is much smaller than the conformer spread. The fixed-cone accessible fraction is also quantized and identical for every retained structure, so it is not sufficiently discriminating on its own.",
        "- The probe omits enzyme-pocket reorganization, water, and an actual nucleophile trajectory. Comparable static access therefore cannot be converted into a hydrolysis-rate claim.",
        "",
        "## Evidence supporting hypothesis 2 (electronic effect)",
        "",
        "- Property-only single points provide MBIS/CHELPG charges and atom-resolved frontier populations for the low-energy subset. SB shows a coherent lower carbonyl-C MBIS charge and lower benzoyl-center LUMO population; smaller shifts in C=O length, Mayer C=O order, and carbonyl stretch point in a compatible direction.",
        "",
        "## Evidence against hypothesis 2",
        "",
        "- CHELPG, bond-length, Mayer-order, dipole, and stretch differences are small relative to conformer spread or expected model sensitivity. The property subset covers about 93% of the electronic weight, not every structure. The present calculation therefore detects an electronic-structure signal but does not establish its magnitude in solution or its kinetic consequence.",
        "",
        "## Not decidable from the current data",
        "",
        "- Two molecules change steric and electronic factors simultaneously; causality cannot be separated without a control or reaction-energy calculation.",
        "- A static reactant structure does not determine hydrolysis rate or an enzymatic activation barrier.",
        "- Differences comparable to conformer spread or method uncertainty must be reported as not detected at this level, not as absent.",
        "- Re/Si face labels are intentionally not assigned by the geometric probe implementation; signed faces are reported instead.",
        "",
        "## Method limitations",
        "",
        "- B3LYP/def2-SVP, gas phase, no explicit dispersion correction and neutral singlet only.",
        "- Common-scaffold RMSD uses fixed input atom order. It does not permute symmetric atoms; the entire terminal aryl-Si(CH3)3 branch is excluded to prevent methyl rotation from defining scaffold clusters.",
        "- Low-frequency modes make sub-kJ/mol Gibbs rankings sensitive to the thermochemistry treatment.",
        "",
        "## Highest-value next calculation",
        "",
        "No further long calculation is required for this preliminary comparison. Pool 109 was already optimized as the single targeted follow-up and remained 4.03 kJ/mol above SB conf002; pool 46 is the only uncovered common-scaffold cluster and was 6.48 kJ/mol above the best retained initial geometry in the DFT single-point screen. The highest-value next checks are a finer continuous steric trajectory descriptor and matched higher-basis/solvent property single points. Pool 46 Opt is optional if closing the residual search asymmetry becomes more important than model validation. Reaction TS and substituted controls remain later-stage work.",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_japanese_preliminary_report(
    output_path: Path,
    bz_summary: Mapping[str, object],
    sb_summary: Mapping[str, object],
    bz_records: Sequence[Mapping[str, object]],
    sb_records: Sequence[Mapping[str, object]],
    *,
    sb_pool_review: Mapping[str, object],
    cross_rmsd: Mapping[str, object],
    bz_property_scope: Mapping[str, object],
    sb_property_scope: Mapping[str, object],
) -> None:
    """Write the human-facing preliminary assessment in Japanese."""

    def stat(summary: Mapping[str, object], metric: str) -> Mapping[str, object]:
        unique = summary.get("common_scaffold_unique", {})
        source = unique if isinstance(unique, Mapping) else summary
        metrics = source.get("metrics", {})
        if not isinstance(metrics, Mapping):
            return {}
        item = metrics.get(metric, {})
        if not isinstance(item, Mapping):
            return {}
        value = item.get("electronic_energy_weighted", {})
        return value if isinstance(value, Mapping) else {}

    def number(value: object, digits: int = 4) -> str:
        return "—" if not isinstance(value, (int, float)) else f"{value:.{digits}f}"

    def mean_std(summary: Mapping[str, object], metric: str, digits: int) -> str:
        item = stat(summary, metric)
        mean = item.get("weighted_mean")
        std = item.get("weighted_std")
        if not isinstance(mean, (int, float)):
            return "—"
        return f"{mean:.{digits}f} ± {float(std or 0.0):.{digits}f}"

    descriptor_rows = []
    descriptors = (
        ("C=O結合長 / Å", "structure.carbonyl_c_o_length_angstrom", 5),
        ("amide C–N結合長 / Å", "structure.amide_c_n_length_angstrom", 5),
        (
            "benzoyl–carbonyl二面角 / deg",
            "structure.benzoyl_carbonyl_torsion_abs_mean_deg",
            3,
        ),
        (
            "攻撃円錐accessible fraction",
            "steric_access.overall_accessible_fraction",
            4,
        ),
        (
            "攻撃円錐best clearance / Å",
            "steric_access.best_clearance_angstrom",
            4,
        ),
        (
            "攻撃円錐clearance p90 / Å",
            "steric_access.clearance_percentiles_angstrom.p90",
            4,
        ),
        ("MBIS carbonyl-C", "electronic.mbis_charge_carbonyl_c", 5),
        ("MBIS carbonyl-O", "electronic.mbis_charge_carbonyl_o", 5),
        ("CHELPG carbonyl-C", "electronic.chelpg_charge_carbonyl_c", 5),
        ("CHELPG carbonyl-O", "electronic.chelpg_charge_carbonyl_o", 5),
        ("Mayer C=O", "electronic.mayer_bond_order_c_o", 5),
        ("Mayer C–N", "electronic.mayer_bond_order_c_n", 5),
        ("dipole / D", "electronic.dipole_magnitude_debye", 4),
        (
            "HOMO Loewdin population (benzoyl center)",
            "electronic.homo_loewdin_population_benzoyl_center",
            5,
        ),
        (
            "LUMO Loewdin population (benzoyl center)",
            "electronic.lumo_loewdin_population_benzoyl_center",
            5,
        ),
        ("HOMO energy / eV", "electronic.homo_energy_ev", 4),
        ("LUMO energy / eV", "electronic.lumo_energy_ev", 4),
        ("carbonyl stretch / cm⁻¹", "carbonyl_stretch.frequency_cm1", 2),
    )
    for label, metric, digits in descriptors:
        bz_item = stat(bz_summary, metric)
        sb_item = stat(sb_summary, metric)
        bz_mean = bz_item.get("weighted_mean")
        sb_mean = sb_item.get("weighted_mean")
        difference = (
            float(sb_mean) - float(bz_mean)
            if isinstance(bz_mean, (int, float))
            and isinstance(sb_mean, (int, float))
            else None
        )
        descriptor_rows.append(
            "| "
            + " | ".join(
                (
                    label,
                    mean_std(bz_summary, metric, digits),
                    mean_std(sb_summary, metric, digits),
                    number(difference, digits),
                    f"{number(bz_item.get('weight_coverage'), 3)} / "
                    f"{number(sb_item.get('weight_coverage'), 3)}",
                )
            )
            + " |"
        )

    def conformer_table(records: Sequence[Mapping[str, object]]) -> List[str]:
        rows = [
            "| conformer | ΔE / kJ mol⁻¹ | ΔG / kJ mol⁻¹ | Freq | 虚振動(<−20) | Opt由来 |",
            "| --- | ---: | ---: | --- | ---: | --- |",
        ]
        for record in sorted(
            records, key=lambda item: float(item["relative_energy_kj_mol"])
        ):
            frequency = bool(record.get("frequency_calculated"))
            imaginary = record.get("imaginary_modes")
            rows.append(
                f"| {record['conformer_id']} | "
                f"{number(record.get('relative_energy_kj_mol'), 3)} | "
                f"{number(record.get('relative_gibbs_energy_kj_mol'), 3)} | "
                f"{'完了' if frequency else '未実施'} | "
                f"{number(imaginary, 0) if frequency else 'null'} | "
                f"{record.get('optimization_provenance') or '不明'} |"
            )
        return rows

    reciprocal = cross_rmsd.get("reciprocal_nearest_pairs", [])
    reciprocal_rows = [
        "| Bz | SB | 共通骨格RMSD / Å | 反応中心RMSD / Å |",
        "| --- | --- | ---: | ---: |",
    ]
    for item in reciprocal if isinstance(reciprocal, list) else []:
        reciprocal_rows.append(
            f"| {item['first_conformer_id']} | {item['second_conformer_id']} | "
            f"{number(item.get('common_scaffold_rmsd_angstrom'), 3)} | "
            f"{number(item.get('reaction_center_rmsd_angstrom'), 3)} |"
        )

    screen = sb_pool_review.get("cheap_reranking", {})
    screen_performed = (
        isinstance(screen, Mapping)
        and screen.get("independent_electronic_reranking_performed") is True
    )
    screen_rows = [
        "| pool index | 安価な電子構造screen ΔE / kJ mol⁻¹ | 判定 |",
        "| ---: | ---: | --- |",
    ]
    if screen_performed:
        recommendations = sb_pool_review.get(
            "additional_dft_recommendations_after_screen", []
        )
        if not recommendations:
            recommendations = sb_pool_review.get(
                "additional_dft_recommendations_after_xtb", []
            )
        for item in recommendations if isinstance(recommendations, list) else []:
            screen_rows.append(
                f"| {item.get('pool_index')} | "
                f"{number(item.get('relative_to_best_retained_initial_geometry_kj_mol', item.get('relative_energy_kj_mol')), 3)} | "
                f"{item.get('additional_opt_priority', item.get('recommendation'))} |"
            )
    else:
        screen_rows.append("| — | — | 単一点順位付け未完了 |")

    bz_freq = int(bz_summary["frequency_completed_count"])
    sb_freq = int(sb_summary["frequency_completed_count"])
    bz_count = int(bz_summary["dft_representative_count"])
    sb_count = int(sb_summary["dft_representative_count"])
    bz_unique = int(
        bz_summary.get("common_scaffold_unique", {}).get(
            "representative_count", bz_count
        )
    )
    sb_unique = int(
        sb_summary.get("common_scaffold_unique", {}).get(
            "representative_count", sb_count
        )
    )
    bz_frequency = bz_summary.get("frequency_bookkeeping", {})
    sb_frequency = sb_summary.get("frequency_bookkeeping", {})
    bz_frequency = bz_frequency if isinstance(bz_frequency, Mapping) else {}
    sb_frequency = sb_frequency if isinstance(sb_frequency, Mapping) else {}
    targeted = sb_pool_review.get("targeted_followup", {})
    targeted = targeted if isinstance(targeted, Mapping) else {}
    pool_unrepresented = sb_pool_review.get("unrepresented_candidates", [])
    pool_unrepresented_count = (
        len(pool_unrepresented) if isinstance(pool_unrepresented, list) else 0
    )
    si_distance = stat(
        sb_summary, "steric_access.silicon_distance_from_carbonyl_c_angstrom"
    ).get("weighted_mean")
    lines = [
        "# 1Bz-LSD_RR / 1SB-LSD_RR 予備比較レポート",
        "",
        "## 結論の範囲",
        "",
        "この比較はB3LYP/def2-SVP・気相で得た計算候補構造の比較であり、実験構造、加水分解速度、酵素反応障壁を直接示すものではない。",
        "現時点では、差が配座内分布または計算法の不確かさ以下なら『この条件では検出できない』とし、『効果が存在しない』とは扱わない。",
        "",
        f"- Bz: DFT代表 {bz_count} / 共通骨格unique {bz_unique}、元のFreq選択 "
        f"{bz_frequency.get('original_selection_completed')}/{bz_frequency.get('original_selection_candidates')}完了、"
        f"補足Freq {bz_frequency.get('supplemental_completed')}、利用可能合計 {bz_freq}",
        f"- SB感度解析集合: DFT代表 {sb_count}（元の10 + targeted follow-up 1） / "
        f"共通骨格unique {sb_unique}、元のFreq選択 "
        f"{sb_frequency.get('original_selection_completed')}/{sb_frequency.get('original_selection_candidates')}完了、"
        f"補足Freq {sb_frequency.get('supplemental_completed')}、利用可能合計 {sb_freq}",
        f"- SB力場候補: 全重原子22 cluster → 共通骨格 {sb_pool_review.get('common_scaffold_cluster_count')} cluster、DFT未カバー {pool_unrepresented_count} cluster",
        f"- MBIS/CHELPG property対象: Bz {bz_property_scope.get('matched_conformer_count', 0)}、SB {sb_property_scope.get('matched_conformer_count', 0)}",
        "- 過去のlauncher/startup失敗はoperational_historyとして保存し、最終Opt/Freq失敗とは別集計",
        "",
        "主表は共通骨格clusterごとの最低電子エネルギー代表を1構造だけ採り、その集合を電子エネルギーで重み付けする。TMS回転だけの重複計上は避けるが、回転子の縮退や配置エントロピーを厳密に扱うものではなく、完全な分子ensembleのBoltzmann分布ではない。Freq未完了構造が残るため、Gibbs重みはFreq完了部分集合に対する条件付き値である。",
        "",
        "## conformer別結果 — Bz",
        "",
        *conformer_table(bz_records),
        "",
        "## conformer別結果 — SB",
        "",
        *conformer_table(sb_records),
        "",
        "## Bz–SB共通骨格の相互最近傍",
        "",
        *reciprocal_rows,
        "",
        "原子対応はTMS末端枝を除去したグラフ同型写像で固定した。対称原子の全置換によるRMSD最小化は行っていないため、値は再現可能だが対称性を許した最小RMSDの上限になり得る。",
        "",
        "## 指標比較",
        "",
        "値は電子エネルギー重み付き平均 ± 配座間標準偏差。coverageは当該指標を持つ構造が全電子重みの何割を占めるかをBz / SBで示す。",
        "",
        "| 指標 | Bz | SB | SB−Bz | weight coverage Bz/SB |",
        "| --- | ---: | ---: | ---: | ---: |",
        *descriptor_rows,
        "",
        "## 仮説1（Si–C結合・配置による立体障害緩和）",
        "",
        "### 支持と整合する点",
        "",
        "- 固定Bürgi–Dunitz円錐でのaccessible fractionはBz/SBとも0.1528、best clearanceのSB−Bzは−0.0009 Åで、配座内標準偏差（約0.01 Å）よりはるかに小さい。少なくとも静的気相構造では、SBのpara-TMSがcarbonyl Cへの幾何学的アクセスを大きく閉じる差は検出されなかった。",
        f"- SBのcarbonyl C–Si距離は電子エネルギー重み付きで約{number(si_distance, 2)} Åで、TMS重原子はcarbonyl Cから5 Å以内に0個だった。TMSは見かけの体積ほど反応中心を直接覆っていないという考えとは整合する。",
        "- 共通骨格RMSDにより、TMSメチル回転だけを別のLSD–benzoyl配座として数える問題を除いた。",
        "",
        "### 反する可能性のある点",
        "",
        "- SBでclearanceが増えた証拠はなく、平均はごくわずかに低い。したがって『Si–C結合が長いこと自体がアクセスを改善する』までは支持できない。現状が示すのは、para-TMSによる直接遮蔽をこの指標では検出できない、という範囲である。",
        "- accessible fractionは全構造で同一になり、角度刻みと閾値に量子化されている。best clearance、p90、probe-radius感度を併記しても、これは予備的な幾何学指標である。",
        "- この幾何学probeには水、酵素ポケット、基質誘導適合、遷移状態が含まれない。静的アクセスが同じでも加水分解障壁が同じとは限らない。",
        "",
        "## 仮説2（TMSの電子効果）",
        "",
        "### 支持と整合する点",
        "",
        "- carbonyl-CのMBIS電荷はSB−Bz = −0.00371 eで、各分子内の配座標準偏差（約0.0004 e）より大きく、低エネルギーproperty部分集合では一貫した差として検出された。TMS-Bz側でcarbonyl Cがわずかに低正電荷になる方向である。",
        "- benzoyl centerのLUMO Loewdin populationはSB−Bz = −0.02154で、配座標準偏差より大きい。C=Oが+0.00010 Å長く、Mayer C=Oが−0.00032、carbonyl stretchが−0.45 cm⁻¹という小さな変化も、弱い電子供与を想定した方向とは概ね整合する。",
        "",
        "### 反する可能性のある点",
        "",
        "- CHELPG carbonyl-C差（−0.00327 e）は配座標準偏差（約0.013–0.015 e）より小さく、結合長、Mayer bond order、dipole、stretchの差も多くは配座ばらつき以下である。MBIS/LUMO局在化ほど全指標が強く一致しているわけではない。",
        "- したがって『TMSに電子効果がある』という予備的証拠は得たが、その大きさが溶液中でも維持されるか、加水分解を速めるか遅くするかは判断できない。電子差の検出と反応速度の説明は分ける。",
        "",
        "## SB未DFTクラスタの再評価",
        "",
        *screen_rows,
        "",
        f"pool 109は事前記録後に1構造だけtargeted DFT Optし、既存current best conf002より{number(targeted.get('relative_to_existing_current_best_kj_mol'), 3)} kJ mol⁻¹高かった。共通骨格RMSD 0.75 Åでは既存構造と重複せず、RRを保持したが、Freq未実施のためimaginary_modesはnullである。",
        "残るpool 46はMMFF構造上のDFT単一点で+6.481 kJ mol⁻¹、GFN2-xTB Optで+6.882 kJ mol⁻¹だった。未緩和単一点/別モデルの順位なので、追加Opt候補の優先順位付けにだけ用いる。今回は大量投入せず保留した。",
        "",
        "## 現時点で判断できないこと",
        "",
        "- Bz→SBの置換では立体効果と電子効果が同時に変わるため、観測差の因果分離はこの二分子だけでは完結しない。",
        "- Freqで虚振動が0でも、それは当該構造が局所極小候補であることを示すだけでglobal minimumを保証しない。",
        "- 低周波モードを含む調和近似のΔGはsub-kJ/mol順位に敏感であり、Freq部分集合だけの平均を完全ensemble平均とは呼べない。",
        "",
        "## 探索範囲の非対称性",
        "",
        "Bzは最大10枠に達する前に候補が収まった一方、SBは元の全重原子cluster数が22で最大10枠に切られた。共通骨格では9 clusterに縮約し、pool 109を追加Optした結果、未カバーはpool 46の1 clusterまで減った。ただしBz/SBの初期候補選択履歴は完全対称ではないため、SB側current bestの確度はなおわずかに低い。",
        "",
        "## 次に価値が高い計算",
        "",
        "1. 直ちに追加の長時間DFTは不要。まず固定円錐より連続的な求核攻撃trajectory/SASA指標を実装し、相互最近傍のBz/SB配座でpaired差を確認する（新規量子化学計算なし）。",
        "2. 仮説2のMBIS/LUMO差を検証するなら、相互対応する低エネルギー構造だけに高い基底・dispersion・暗黙溶媒を用いたproperty単一点を行う。OptやTSより安価で、気相/基底依存性を直接判定できる。",
        "3. pool 46 Optは探索非対称性を完全に閉じたい場合だけ実行する。現screenではcurrent bestを更新する優先度は低く、見積り0.4–0.8時間、上位化した場合のみFreq約2時間を追加する。",
        "4. 上記で電子差または立体差がモデル変更後も維持されて初めて、溶媒中反応物複合体、加水分解TS、置換基対照へ進む価値を再判定する。",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _uff_single_point_energy(model: MoleculeModel, smiles: str) -> Optional[float]:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol.GetNumAtoms() != len(model.atoms) or not AllChem.UFFHasAllMoleculeParams(mol):
        return None
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for atom in model.atoms:
        conformer.SetAtomPosition(atom.index, atom.position)
    mol.RemoveAllConformers()
    mol.AddConformer(conformer, assignId=True)
    forcefield = AllChem.UFFGetMoleculeForceField(mol, confId=0)
    return float(forcefield.CalcEnergy()) if forcefield is not None else None


def review_forcefield_pool(
    ensemble_dir: Path,
    output_path: Path,
) -> Dict[str, object]:
    """Regenerate the deterministic pool and test common-core coverage cheaply."""

    payload = _load_json(Path(ensemble_dir) / "ensemble.json")
    molecule = payload["molecule"]
    search = payload["conformer_search"]
    initial = search["initial_screening"]
    smiles = str(molecule.get("original_smiles") or molecule["smiles"])
    pool_size = int(search["pool_size"])
    candidates, inspection = generate_forcefield_conformers(
        smiles,
        pool_size,
        name=str(payload["name"]),
        random_seed=int(search["random_seed"]),
        prune_rms_threshold=float(search["embedding_prune_rmsd_angstrom"]),
    )
    records = {
        int(record["pool_index"]): record for record in initial.get("records", [])
    }
    energy_matches = []
    for candidate in candidates:
        recorded = records.get(candidate.conformer_index, {}).get(
            "energy_kcal_mol"
        )
        if recorded is None or candidate.energy_kcal_mol is None:
            continue
        energy_matches.append(abs(float(recorded) - candidate.energy_kcal_mol))
    exact_reproduction = bool(energy_matches) and max(energy_matches) < 1.0e-8
    minimum = min(
        float(candidate.energy_kcal_mol)
        for candidate in candidates
        if candidate.converged and candidate.energy_kcal_mol is not None
    )
    window = float(initial["energy_window_kj_mol"])
    eligible = sorted(
        [
            candidate
            for candidate in candidates
            if candidate.converged
            and candidate.energy_kcal_mol is not None
            and (float(candidate.energy_kcal_mol) - minimum) * 4.184
            <= window + 1.0e-9
        ],
        key=lambda candidate: (
            float(candidate.energy_kcal_mol),
            candidate.conformer_index,
        ),
    )
    subsets = rmsd_atom_subsets(eligible[0].model)

    def clusters(atom_indices: Sequence[int], threshold: float):
        representatives = []
        members: List[List[int]] = []
        for candidate in eligible:
            selected = None
            for index, representative in enumerate(representatives):
                if (
                    aligned_atom_subset_rmsd(
                        candidate.model, representative.model, atom_indices
                    )
                    < threshold
                ):
                    selected = index
                    break
            if selected is None:
                selected = len(representatives)
                representatives.append(candidate)
                members.append([])
            members[selected].append(candidate.conformer_index)
        return representatives, members

    common_representatives, common_members = clusters(
        subsets.common_scaffold, 0.75
    )
    reaction_representatives, _reaction_members = clusters(
        subsets.reaction_center, 0.25
    )
    selected_pool = {
        int(index)
        for index, record in records.items()
        if record.get("selected_for_dft") is True
    }
    _current_payload, optimized_conformers = load_existing_conformers(ensemble_dir)
    uff = {
        candidate.conformer_index: _uff_single_point_energy(candidate.model, smiles)
        for candidate in eligible
    }
    available_uff = [value for value in uff.values() if value is not None]
    uff_minimum = min(available_uff) if available_uff else None
    unrepresented = []
    cluster_coverage = []
    for cluster_index, (representative, members) in enumerate(
        zip(common_representatives, common_members), 1
    ):
        optimized_rmsds = [
            (
                str(conformer.entry["conformer_id"]),
                aligned_atom_subset_rmsd(
                    representative.model,
                    conformer.model,
                    subsets.common_scaffold,
                ),
            )
            for conformer in optimized_conformers
        ]
        closest_id, closest_rmsd = min(optimized_rmsds, key=lambda item: item[1])
        represented_by_selected_pool = bool(selected_pool.intersection(members))
        represented_by_optimized_dft = closest_rmsd < 0.75
        cluster_coverage.append(
            {
                "common_scaffold_cluster_id": f"ff_core_cluster_{cluster_index:03d}",
                "representative_pool_index": representative.conformer_index,
                "represented_by_original_selected_pool": represented_by_selected_pool,
                "closest_optimized_dft_conformer_id": closest_id,
                "closest_optimized_dft_common_scaffold_rmsd_angstrom": closest_rmsd,
                "represented_by_optimized_dft_within_0_75_angstrom": represented_by_optimized_dft,
            }
        )
        if represented_by_optimized_dft:
            continue
        uff_energy = uff[representative.conformer_index]
        unrepresented.append(
            {
                "common_scaffold_cluster_id": f"ff_core_cluster_{cluster_index:03d}",
                "representative_pool_index": representative.conformer_index,
                "member_pool_indices": members,
                "mmff_energy_kcal_mol": representative.energy_kcal_mol,
                "mmff_relative_energy_kj_mol": (
                    (float(representative.energy_kcal_mol) - minimum) * 4.184
                ),
                "uff_single_point_energy_kcal_mol": uff_energy,
                "uff_relative_energy_kj_mol": (
                    (float(uff_energy) - float(uff_minimum)) * 4.184
                    if uff_energy is not None and uff_minimum is not None
                    else None
                ),
            }
        )
    unrepresented.sort(
        key=lambda item: (
            float("inf")
            if item["uff_relative_energy_kj_mol"] is None
            else float(item["uff_relative_energy_kj_mol"]),
            float(item["mmff_relative_energy_kj_mol"]),
        )
    )
    recommended = []
    for item in unrepresented[:3]:
        recommended.append(
            {
                **item,
                "recommendation": "single_point_reranking_before_any_additional_opt",
                "reason": (
                    "Common-scaffold cluster has no selected DFT representative. UFF is "
                    "only a force-field sensitivity check, not an electronic reranking."
                ),
                "estimated_opt_wall_time_if_later_selected_hours": "0.4-0.8",
            }
        )
    pool_records = []
    for candidate in candidates:
        pool_records.append(
            {
                "pool_index": candidate.conformer_index,
                "forcefield": candidate.forcefield,
                "converged": candidate.converged,
                "energy_kcal_mol": candidate.energy_kcal_mol,
                "selected_for_original_dft": candidate.conformer_index
                in selected_pool,
                "coordinates_angstrom": [
                    [atom.x, atom.y, atom.z] for atom in candidate.model.atoms
                ],
            }
        )
    report = {
        "generated_at_utc": _utc_now(),
        "name": payload["name"],
        "canonical_isomeric_smiles": inspection.canonical_isomeric_smiles,
        "pool_size": pool_size,
        "random_seed": search["random_seed"],
        "rdkit_version": search.get("rdkit_version"),
        "exact_energy_reproduction_of_retained_metadata": exact_reproduction,
        "maximum_energy_difference_kcal_mol": (
            max(energy_matches) if energy_matches else None
        ),
        "energy_window_kj_mol": window,
        "eligible_count": len(eligible),
        "original_all_heavy_cluster_count": initial.get("cluster_count"),
        "common_scaffold_cluster_count": len(common_representatives),
        "reaction_center_cluster_count": len(reaction_representatives),
        "original_selected_pool_indices": sorted(selected_pool),
        "common_scaffold_cluster_coverage": cluster_coverage,
        "unrepresented_candidates": unrepresented,
        "recommended_candidates": recommended,
        "cheap_reranking": {
            "method": "UFF single-point energy on MMFF94s-optimized geometry",
            "role": "force-field sensitivity check only",
            "independent_electronic_reranking_performed": False,
            "reason": (
                "No standalone xTB backend is installed. Launching a batch of new ORCA "
                "single points was deferred until the common-scaffold redundancy and "
                "existing-property analysis were known."
            ),
        },
        "additional_dft_opt_launched": False,
        "additional_dft_policy": (
            "Record candidates and expected impact before launching any new Opt jobs."
        ),
        "pool_records": pool_records,
    }
    _atomic_write_json(output_path, report)
    return report


def run_preliminary_comparison(
    bz_ensemble_dir: Path,
    sb_ensemble_dir: Path,
    output_dir: Path,
    *,
    update_ensemble_json: bool = True,
) -> Dict[str, object]:
    """Repair both ensembles and write reproducible Bz/SB comparison artifacts."""

    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    refresh_existing_ensemble_metadata(
        bz_ensemble_dir, write=update_ensemble_json
    )
    refresh_existing_ensemble_metadata(
        sb_ensemble_dir, write=update_ensemble_json
    )
    bz_payload, bz_conformers = load_existing_conformers(bz_ensemble_dir)
    sb_payload, sb_conformers = load_existing_conformers(sb_ensemble_dir)
    targeted_followup = load_targeted_followup_conformer(sb_payload, output_dir)
    if targeted_followup is not None:
        sb_conformers.append(targeted_followup)
        followup_rmsd = secondary_rmsd_analysis(sb_conformers)
        followup_rmsd_by_index = {
            int(item["conformer_index"]): item
            for item in followup_rmsd.get("records", [])
        }
        for conformer in sb_conformers:
            record = followup_rmsd_by_index[int(conformer.entry["conformer_index"])]
            conformer.entry.update(
                {
                    "all_heavy_rmsd_cluster_id": record[
                        "all_heavy_rmsd_cluster_id"
                    ],
                    "common_scaffold_rmsd_cluster_id": record[
                        "common_scaffold_rmsd_cluster_id"
                    ],
                    "reaction_center_rmsd_cluster_id": record[
                        "reaction_center_rmsd_cluster_id"
                    ],
                }
            )
    bz_records = [analyze_existing_conformer(item) for item in bz_conformers]
    sb_records = [analyze_existing_conformer(item) for item in sb_conformers]
    cross_rmsd = cross_ensemble_rmsd_analysis(bz_conformers, sb_conformers)
    property_summary = output_dir / "electronic_property_results.json"
    bz_property_scope = merge_property_single_points(
        "1Bz-LSD_RR", bz_records, property_summary
    )
    sb_property_scope = merge_property_single_points(
        "1SB-LSD_RR", sb_records, property_summary
    )
    bz_summary = ensemble_descriptor_summary(bz_payload, bz_records)
    sb_summary_payload = dict(sb_payload)
    if targeted_followup is not None:
        retained_entries = list(sb_payload.get("final_conformer_ensemble", []))
        retained_entries.append(targeted_followup.entry)
        sb_summary_payload["final_conformer_ensemble"] = retained_entries
    sb_summary = ensemble_descriptor_summary(sb_summary_payload, sb_records)
    sb_pool = review_forcefield_pool(
        sb_ensemble_dir, output_dir / "1SB-LSD_RR_regenerated_pool.json"
    )
    if targeted_followup is not None:
        for coverage in sb_pool.get("common_scaffold_cluster_coverage", []):
            if coverage.get("representative_pool_index") == 109:
                coverage[
                    "evaluated_by_targeted_followup_dft_opt"
                ] = "conf011_pool109"
        sb_pool["unrepresented_candidates"] = [
            item
            for item in sb_pool.get("unrepresented_candidates", [])
            if item.get("representative_pool_index") != 109
        ]
        sb_pool["recommended_candidates"] = [
            item
            for item in sb_pool.get("recommended_candidates", [])
            if item.get("representative_pool_index") != 109
        ]
        sb_pool["additional_dft_opt_launched"] = True
        sb_pool["targeted_followup"] = _load_json(
            output_dir / "sb_additional_dft_opt_result.json"
        )
    xtb_path = output_dir / "sb_pool_gfn2_xtb_screen.json"
    if xtb_path.is_file():
        xtb = _load_json(xtb_path)
        sb_pool["cheap_reranking"] = {
            "method": xtb.get("method"),
            "software": xtb.get("software"),
            "version": xtb.get("version"),
            "role": "independent relaxed screening before any additional DFT Opt",
            "independent_electronic_reranking_performed": True,
            "results": xtb.get("results", []),
            "interpretation": xtb.get("interpretation"),
        }
        sb_pool["additional_dft_recommendations_after_xtb"] = xtb.get(
            "additional_dft_opt_recommendation", []
        )
    screen_path = output_dir / "sb_pool_singlepoint_screen.json"
    if screen_path.is_file():
        screen = _load_json(screen_path)
        screen_results = screen.get("omitted_cluster_results", [])
        sb_pool["cheap_reranking"] = {
            "method": screen.get("calculation"),
            "role": "unrelaxed electronic-energy screen before any additional Opt",
            "independent_electronic_reranking_performed": True,
            "results": screen_results,
            "interpretation": screen.get("interpretation"),
        }
        sb_pool["additional_dft_recommendations_after_screen"] = [
            {
                **item,
                "additional_opt_priority": (
                    "review_for_one_additional_opt"
                    if float(
                        item.get(
                            "relative_to_best_retained_initial_geometry_kj_mol",
                            float("inf"),
                        )
                    )
                    <= 10.0
                    else "low_priority_at_current_screening_level"
                ),
            }
            for item in screen_results
            if isinstance(item, dict)
        ]
    if targeted_followup is not None:
        completed_pool_indices = {109}
        for key in (
            "additional_dft_recommendations_after_xtb",
            "additional_dft_recommendations_after_screen",
        ):
            sb_pool[key] = [
                item
                for item in sb_pool.get(key, [])
                if isinstance(item, dict)
                and item.get("pool_index") not in completed_pool_indices
            ]
        sb_pool["completed_targeted_followup_pool_indices"] = sorted(
            completed_pool_indices
        )
    report = {
        "schema_version": 1,
        "generated_at_utc": _utc_now(),
        "purpose": "Bz versus TMS-Bz preliminary hypothesis comparison",
        "calculation_scope": (
            "Retained B3LYP/def2-SVP gas-phase Opt/Freq data, two supplemental "
            "Freq-only jobs, deterministic force-field pool regeneration, cheap "
            "xTB/DFT single-point screening, one targeted SB DFT Opt (pool 109), "
            "and property-only single points; no bulk DFT launch."
        ),
        "ensembles": {
            "1Bz-LSD_RR": {
                "summary": bz_summary,
                "conformers": bz_records,
                "property_single_point_scope": bz_property_scope,
            },
            "1SB-LSD_RR": {
                "summary": sb_summary,
                "conformers": sb_records,
                "property_single_point_scope": sb_property_scope,
            },
        },
        "sb_forcefield_pool_review": {
            key: value for key, value in sb_pool.items() if key != "pool_records"
        },
        "targeted_followup_included_in_sb_electronic_sensitivity": (
            targeted_followup is not None
        ),
        "cross_ensemble_common_scaffold_rmsd": cross_rmsd,
        "interpretation_rules": {
            "detected_vs_absent": (
                "A difference smaller than conformer spread or model uncertainty is "
                "reported as not detected at this level, never as proof of absence."
            ),
            "global_minimum": (
                "Lowest-energy structure found is not called the global minimum."
            ),
            "gibbs_scope": (
                "Gibbs weights are conditional whenever Freq is incomplete."
            ),
        },
    }
    operational_history = output_dir / "operational_history.json"
    if operational_history.is_file():
        report["operational_history"] = _load_json(operational_history)
    _atomic_write_json(output_dir / "analysis.json", report)
    _write_metric_csv(output_dir / "1Bz-LSD_RR_conformer_metrics.csv", bz_records)
    _write_metric_csv(output_dir / "1SB-LSD_RR_conformer_metrics.csv", sb_records)
    cross_records = cross_rmsd.get("records", [])
    if isinstance(cross_records, list) and cross_records:
        with (output_dir / "Bz_SB_common_scaffold_rmsd_matrix.csv").open(
            "w", newline="", encoding="utf-8-sig"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(cross_records[0]))
            writer.writeheader()
            writer.writerows(cross_records)
    write_preliminary_comparison_report(
        output_dir / "PRELIMINARY_COMPARISON.md",
        bz_summary,
        sb_summary,
        bz_records,
        sb_records,
        sb_pool_review=sb_pool,
    )
    write_japanese_preliminary_report(
        output_dir / "PRELIMINARY_COMPARISON_JA.md",
        bz_summary,
        sb_summary,
        bz_records,
        sb_records,
        sb_pool_review=sb_pool,
        cross_rmsd=cross_rmsd,
        bz_property_scope=bz_property_scope,
        sb_property_scope=sb_property_scope,
    )
    return report

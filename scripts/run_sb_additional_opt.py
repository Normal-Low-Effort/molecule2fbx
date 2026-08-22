"""Run exactly one recorded additional SB DFT Opt for pool index 109."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

from molecule2fbx.comparison import load_existing_conformers, structural_descriptors
from molecule2fbx.model import Atom, Bond, MoleculeModel, normalize_bond_order
from molecule2fbx.quantum.base import QuantumSettings
from molecule2fbx.quantum.orca import OrcaBackend
from molecule2fbx.structures import (
    aligned_atom_subset_rmsd,
    rmsd_atom_subsets,
    validate_model_stereochemistry,
)


WORKSPACE = Path(__file__).resolve().parents[1]
ENSEMBLE_DIR = WORKSPACE / "outputs" / "1SB-LSD_RR"
COMPARISON_DIR = WORKSPACE / "outputs" / "Bz_vs_SB_preliminary_comparison"
POOL_PATH = COMPARISON_DIR / "1SB-LSD_RR_regenerated_pool.json"
WORK_DIR = COMPARISON_DIR / "sb_additional_dft_opt" / "conformer_011_pool_109"
STATUS_PATH = COMPARISON_DIR / "sb_additional_dft_opt_status.json"
RESULT_PATH = COMPARISON_DIR / "sb_additional_dft_opt_result.json"
ORCA = Path(r"C:\ORCA_6.1.1\orca.exe")
HARTREE_TO_KJ_MOL = 2625.4996394799


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def pool_model() -> tuple[MoleculeModel, str]:
    ensemble = json.loads((ENSEMBLE_DIR / "ensemble.json").read_text(encoding="utf-8"))
    pool = json.loads(POOL_PATH.read_text(encoding="utf-8"))
    smiles = ensemble["molecule"].get("original_smiles") or ensemble["molecule"]["smiles"]
    record = next(item for item in pool["pool_records"] if item["pool_index"] == 109)
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    coordinates = record["coordinates_angstrom"]
    if len(coordinates) != molecule.GetNumAtoms():
        raise RuntimeError("Pool-109 coordinate atom count does not match SMILES")
    atoms = tuple(
        Atom(index, atom.GetSymbol(), *map(float, coordinates[index]))
        for index, atom in enumerate(molecule.GetAtoms())
    )
    bonds = tuple(
        Bond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            normalize_bond_order(
                float(bond.GetBondTypeAsDouble()), bond.GetIsAromatic()
            ),
        )
        for bond in molecule.GetBonds()
    )
    return (
        MoleculeModel(
            None,
            "1SB-LSD_RR",
            atoms,
            bonds,
            {
                "original_smiles": smiles,
                "canonical_isomeric_smiles": Chem.MolToSmiles(
                    Chem.MolFromSmiles(smiles), isomericSmiles=True
                ),
                "conformer_index": 10,
                "conformer_pool_index": 109,
                "conformer_pool_index_provenance": "regenerated_exact_seed_pool",
                "optimization_provenance": "computed_targeted_followup",
            },
        ),
        smiles,
    )


def main() -> int:
    if not ORCA.is_file():
        raise SystemExit(f"ORCA executable not found: {ORCA}")
    model, smiles = pool_model()
    initial_stereo = validate_model_stereochemistry(smiles, model)
    if not initial_stereo.matches:
        raise RuntimeError("Pool-109 initial geometry does not retain expected RR stereo")
    if WORK_DIR.exists() and any(WORK_DIR.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty Opt directory: {WORK_DIR}")
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    os.environ["PATH"] = (
        rf"C:\Program Files\Microsoft MPI\Bin;{ORCA.parent};"
        + os.environ.get("PATH", "")
    )
    started = datetime.now(timezone.utc)
    write_json(
        STATUS_PATH,
        {
            "status": "RUNNING",
            "pool_index": 109,
            "started_at_utc": started.isoformat().replace("+00:00", "Z"),
            "finished_at_utc": None,
            "work_directory": str(WORK_DIR.resolve()),
        },
    )
    settings = QuantumSettings(
        method="dft",
        functional="B3LYP",
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
        frequency=False,
        timeout=7200,
        max_opt_steps=200,
        max_scf_iterations=200,
        nprocs=16,
        maxcore_mb=1000,
    )
    result = OrcaBackend(str(ORCA)).optimize(model, settings, WORK_DIR, 10)
    stereo = validate_model_stereochemistry(smiles, result.model)
    if not stereo.matches:
        raise RuntimeError("Additional DFT Opt changed expected RR stereochemistry")
    ensemble, existing = load_existing_conformers(ENSEMBLE_DIR)
    subsets = rmsd_atom_subsets(result.model)
    rmsds = [
        {
            "conformer_id": item.entry["conformer_id"],
            "common_scaffold_rmsd_angstrom": aligned_atom_subset_rmsd(
                result.model, item.model, subsets.common_scaffold
            ),
            "reaction_center_rmsd_angstrom": aligned_atom_subset_rmsd(
                result.model, item.model, subsets.reaction_center
            ),
        }
        for item in existing
    ]
    closest = min(rmsds, key=lambda item: item["common_scaffold_rmsd_angstrom"])
    current_best = min(
        float(item["dft_energy_hartree"])
        for item in ensemble["final_conformer_ensemble"]
    )
    payload = {
        "generated_at_utc": now(),
        "pool_index": 109,
        "conformer_id": "conf011_pool109",
        "software": "ORCA",
        "version": result.backend_version,
        "functional": "B3LYP",
        "basis": "def2-SVP",
        "charge": 0,
        "multiplicity": 1,
        "geometry_converged": result.geometry_converged,
        "scf_converged": result.scf_converged,
        "frequency_calculated": False,
        "imaginary_modes": None,
        "energy_hartree": result.final_energy_hartree,
        "relative_to_existing_current_best_kj_mol": (
            result.final_energy_hartree - current_best
        )
        * HARTREE_TO_KJ_MOL,
        "stereochemistry_validation": stereo.to_metadata(),
        "closest_existing_structure": closest,
        "unique_common_scaffold_at_0_75_angstrom": (
            float(closest["common_scaffold_rmsd_angstrom"]) >= 0.75
        ),
        "all_existing_rmsds": rmsds,
        "structural_descriptors": structural_descriptors(result.model),
        "calculation_directory": str(WORK_DIR.resolve()),
        "optimized_xyz": str((WORK_DIR / "conformer_011.xyz").resolve()),
        "interpretation": (
            "Targeted follow-up candidate; not automatically part of the original "
            "ensemble and not a frequency-confirmed local minimum."
        ),
    }
    write_json(RESULT_PATH, payload)
    write_json(
        STATUS_PATH,
        {
            "status": "SUCCESS",
            "pool_index": 109,
            "started_at_utc": started.isoformat().replace("+00:00", "Z"),
            "finished_at_utc": now(),
            "elapsed_minutes": round(
                (datetime.now(timezone.utc) - started).total_seconds() / 60.0, 2
            ),
            "work_directory": str(WORK_DIR.resolve()),
        },
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

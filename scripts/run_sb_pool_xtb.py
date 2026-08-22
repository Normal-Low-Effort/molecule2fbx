"""GFN2-xTB reoptimization/reranking of the nine SB common-core representatives."""

from __future__ import annotations

import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

from molecule2fbx.frequency import read_xyz_model
from molecule2fbx.structures import (
    aligned_atom_subset_rmsd,
    rmsd_atom_subsets,
    validate_model_stereochemistry,
)


WORKSPACE = Path(__file__).resolve().parents[1]
ENSEMBLE_DIR = WORKSPACE / "outputs" / "1SB-LSD_RR"
COMPARISON_DIR = WORKSPACE / "outputs" / "Bz_vs_SB_preliminary_comparison"
POOL_PATH = COMPARISON_DIR / "1SB-LSD_RR_regenerated_pool.json"
OUTPUT_ROOT = COMPARISON_DIR / "sb_pool_gfn2_xtb"
STATUS_PATH = COMPARISON_DIR / "sb_pool_gfn2_xtb_status.json"
SUMMARY_PATH = COMPARISON_DIR / "sb_pool_gfn2_xtb_screen.json"
XTB = Path(
    os.environ.get(
        "XTB_EXECUTABLE", r"C:\ORCA_6.1.1\xtb-6.7.1pre\xtb.exe"
    )
)
ENERGY_RE = re.compile(
    r"TOTAL ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh",
    re.IGNORECASE,
)
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


def write_xyz(path: Path, elements: list[str], coordinates: list[list[float]]) -> None:
    lines = [str(len(elements)), "retained ETKDG/MMFF94s pool geometry"]
    lines.extend(
        f"{element:<3} {xyz[0]: .12f} {xyz[1]: .12f} {xyz[2]: .12f}"
        for element, xyz in zip(elements, coordinates)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_completed(job_dir: Path) -> tuple[float, Path]:
    output_path = job_dir / "xtb.out"
    optimized_xyz = job_dir / "xtbopt.xyz"
    text = output_path.read_text(encoding="utf-8", errors="replace")
    energies = [float(value) for value in ENERGY_RE.findall(text)]
    if "normal termination of xtb" not in text or not energies:
        raise RuntimeError(f"Incomplete xTB output: {output_path}")
    if not optimized_xyz.is_file():
        raise RuntimeError(f"xTB optimized XYZ is missing: {optimized_xyz}")
    return energies[-1], optimized_xyz


def main() -> int:
    if not XTB.is_file():
        raise SystemExit(f"xTB executable not found: {XTB}")
    pool = json.loads(POOL_PATH.read_text(encoding="utf-8"))
    ensemble = json.loads((ENSEMBLE_DIR / "ensemble.json").read_text(encoding="utf-8"))
    smiles = ensemble["molecule"].get("original_smiles") or ensemble["molecule"]["smiles"]
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    records = {int(item["pool_index"]): item for item in pool["pool_records"]}
    representatives = [
        int(item["representative_pool_index"])
        for item in pool["common_scaffold_cluster_coverage"]
    ]
    omitted = {
        int(item["representative_pool_index"])
        for item in pool["unrepresented_candidates"]
    }
    state = {
        "status": "RUNNING",
        "started_at_utc": now(),
        "finished_at_utc": None,
        "current_pool_index": None,
        "jobs": [],
    }
    write_json(STATUS_PATH, state)
    results = []
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "16"
    for pool_index in representatives:
        record = records[pool_index]
        coordinates = record["coordinates_angstrom"]
        if len(coordinates) != len(elements):
            raise RuntimeError(f"Atom-count mismatch for pool {pool_index}")
        job_dir = OUTPUT_ROOT / f"pool_{pool_index:03d}"
        job_dir.mkdir(parents=True, exist_ok=True)
        input_xyz = job_dir / "input.xyz"
        output_path = job_dir / "xtb.out"
        if output_path.is_file():
            energy, optimized_xyz = parse_completed(job_dir)
            provenance = "reused_completed_xtb_job"
        else:
            unexpected = [path for path in job_dir.iterdir() if path != input_xyz]
            if unexpected:
                raise RuntimeError(
                    f"Refusing to overwrite incomplete/non-empty xTB directory: {job_dir}"
                )
            write_xyz(input_xyz, elements, coordinates)
            state["current_pool_index"] = pool_index
            write_json(STATUS_PATH, state)
            started = datetime.now(timezone.utc)
            process = subprocess.run(
                [
                    str(XTB),
                    input_xyz.name,
                    "--gfn",
                    "2",
                    "--opt",
                    "tight",
                    "--chrg",
                    "0",
                    "--uhf",
                    "0",
                    "--parallel",
                    "16",
                ],
                cwd=job_dir,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                env=env,
                check=False,
                timeout=3600,
            )
            output_path.write_text(
                process.stdout + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
                encoding="utf-8",
            )
            energy, optimized_xyz = parse_completed(job_dir)
            if process.returncode != 0:
                raise RuntimeError(f"xTB exited with code {process.returncode}")
            provenance = "computed_this_run"
            state["jobs"].append(
                {
                    "pool_index": pool_index,
                    "started_at_utc": started.isoformat().replace("+00:00", "Z"),
                    "finished_at_utc": now(),
                    "provenance": provenance,
                }
            )
            write_json(STATUS_PATH, state)
        model = read_xyz_model(optimized_xyz).with_metadata(original_smiles=smiles)
        stereo = validate_model_stereochemistry(smiles, model)
        if not stereo.matches:
            raise RuntimeError(f"xTB structure changed expected stereochemistry: pool {pool_index}")
        results.append(
            {
                "pool_index": pool_index,
                "originally_unrepresented": pool_index in omitted,
                "energy_hartree": energy,
                "stereochemistry_validation": stereo.to_metadata(),
                "optimized_xyz": str(optimized_xyz.resolve()),
                "output": str(output_path.resolve()),
                "provenance": provenance,
                "model": model,
            }
        )
    minimum = min(float(item["energy_hartree"]) for item in results)
    ranked = sorted(results, key=lambda item: float(item["energy_hartree"]))
    representatives_models = []
    subsets = rmsd_atom_subsets(ranked[0]["model"])
    for item in ranked:
        assigned = None
        assigned_rmsd = 0.0
        for cluster_index, retained in enumerate(representatives_models):
            rmsd = aligned_atom_subset_rmsd(
                item["model"], retained["model"], subsets.common_scaffold
            )
            if rmsd < 0.75:
                assigned = cluster_index
                assigned_rmsd = rmsd
                break
        if assigned is None:
            assigned = len(representatives_models)
            representatives_models.append(item)
        item["relative_energy_kj_mol"] = (
            float(item["energy_hartree"]) - minimum
        ) * HARTREE_TO_KJ_MOL
        item["xtb_common_scaffold_cluster_id"] = f"xtb_core_cluster_{assigned + 1:03d}"
        item["xtb_common_scaffold_rmsd_to_representative_angstrom"] = assigned_rmsd
    serializable = [{key: value for key, value in item.items() if key != "model"} for item in ranked]
    omitted_results = [item for item in serializable if item["originally_unrepresented"]]
    summary = {
        "generated_at_utc": now(),
        "software": "xTB",
        "version": "6.7.1pre",
        "method": "GFN2-xTB tight geometry optimization",
        "charge": 0,
        "unpaired_electrons": 0,
        "nprocs": 16,
        "purpose": "Independent cheap reranking before any additional DFT Opt",
        "original_common_scaffold_representative_count": len(representatives),
        "xtb_common_scaffold_cluster_count": len(representatives_models),
        "results": serializable,
        "originally_unrepresented_results": omitted_results,
        "additional_dft_opt_recommendation": [
            {
                "pool_index": item["pool_index"],
                "relative_energy_kj_mol": item["relative_energy_kj_mol"],
                "recommendation": (
                    "review_for_at_most_one_additional_dft_opt"
                    if float(item["relative_energy_kj_mol"]) <= 10.0
                    else "low_priority_no_additional_dft_opt_now"
                ),
            }
            for item in omitted_results
        ],
        "interpretation": (
            "GFN2-xTB rankings are a screening model and are not combined with DFT "
            "energies or thermochemistry. No additional DFT Opt was launched."
        ),
    }
    write_json(SUMMARY_PATH, summary)
    state.update(status="SUCCESS", finished_at_utc=now(), current_pool_index=None)
    write_json(STATUS_PATH, state)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

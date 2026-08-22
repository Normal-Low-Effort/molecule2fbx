"""Non-destructive B3LYP/def2-SVP screening of omitted SB pool clusters."""

from __future__ import annotations

import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

from molecule2fbx.quantum.orca import parse_orca_output


WORKSPACE = Path(__file__).resolve().parents[1]
ENSEMBLE_DIR = WORKSPACE / "outputs" / "1SB-LSD_RR"
COMPARISON_DIR = WORKSPACE / "outputs" / "Bz_vs_SB_preliminary_comparison"
POOL_PATH = COMPARISON_DIR / "1SB-LSD_RR_regenerated_pool.json"
OUTPUT_ROOT = COMPARISON_DIR / "sb_pool_singlepoints"
STATUS_PATH = COMPARISON_DIR / "sb_pool_singlepoint_status.json"
SUMMARY_PATH = COMPARISON_DIR / "sb_pool_singlepoint_screen.json"
ORCA = Path(os.environ.get("ORCA_EXECUTABLE", r"C:\ORCA_6.1.1\orca.exe"))
ENERGY_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+(?:\.\d+)?)")
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


def render_input(elements, coordinates) -> str:
    lines = [
        "! B3LYP def2-SVP TightSCF",
        "",
        "%scf",
        "  MaxIter 200",
        "end",
        "%pal",
        "  nprocs 16",
        "end",
        "%maxcore 1000",
        "",
        "* xyz 0 1",
    ]
    lines.extend(
        f"  {element:<3} {xyz[0]: .12f} {xyz[1]: .12f} {xyz[2]: .12f}"
        for element, xyz in zip(elements, coordinates)
    )
    lines.extend(["*", ""])
    return "\n".join(lines)


def retained_initial_energies(ensemble: dict) -> list[dict]:
    records = []
    for entry in ensemble["final_conformer_ensemble"]:
        index = int(entry["conformer_index"])
        stem = f"conformer_{index + 1:03d}"
        output = Path(entry["calculation_directory"]) / f"{stem}.out"
        energies = [
            float(value)
            for value in ENERGY_RE.findall(
                output.read_text(encoding="utf-8", errors="replace")
            )
        ]
        if energies:
            records.append(
                {
                    "conformer_id": entry["conformer_id"],
                    "initial_geometry_energy_hartree": energies[0],
                    "optimized_energy_hartree": energies[-1],
                }
            )
    return records


def main() -> int:
    if not ORCA.is_file():
        raise SystemExit(f"ORCA executable not found: {ORCA}")
    pool = json.loads(POOL_PATH.read_text(encoding="utf-8"))
    ensemble = json.loads((ENSEMBLE_DIR / "ensemble.json").read_text(encoding="utf-8"))
    smiles = ensemble["molecule"].get("original_smiles") or ensemble["molecule"]["smiles"]
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    records = {int(item["pool_index"]): item for item in pool["pool_records"]}
    recommended = {
        int(item["representative_pool_index"]): item
        for item in pool["recommended_candidates"]
    }
    requested = [
        int(item["representative_pool_index"]) for item in pool["recommended_candidates"]
    ]
    retained = retained_initial_energies(ensemble)
    reference = min(item["initial_geometry_energy_hartree"] for item in retained)
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
    env["PATH"] = rf"C:\Program Files\Microsoft MPI\Bin;{ORCA.parent};" + env.get("PATH", "")
    for pool_index in requested:
        record = records[pool_index]
        coordinates = record["coordinates_angstrom"]
        if len(coordinates) != len(elements):
            raise RuntimeError(f"Atom-count mismatch for pool {pool_index}")
        job_dir = OUTPUT_ROOT / f"pool_{pool_index:03d}"
        job_dir.mkdir(parents=True, exist_ok=True)
        input_path = job_dir / f"pool_{pool_index:03d}_sp.inp"
        output_path = job_dir / f"pool_{pool_index:03d}_sp.out"
        if output_path.is_file():
            parsed = parse_orca_output(
                output_path.read_text(encoding="utf-8", errors="replace"),
                require_geometry_convergence=False,
            )
            provenance = "reused_completed_property_job"
        else:
            unexpected = [path for path in job_dir.iterdir() if path != input_path]
            if unexpected:
                raise RuntimeError(
                    f"Refusing to overwrite incomplete/non-empty job directory: {job_dir}"
                )
            input_path.write_text(render_input(elements, coordinates), encoding="utf-8")
            state["current_pool_index"] = pool_index
            write_json(STATUS_PATH, state)
            started = datetime.now(timezone.utc)
            process = subprocess.run(
                [str(ORCA), input_path.name],
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
            parsed = parse_orca_output(
                output_path.read_text(encoding="utf-8", errors="replace"),
                require_geometry_convergence=False,
            )
            if process.returncode != 0:
                raise RuntimeError(f"ORCA exited with code {process.returncode}")
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
        results.append(
            {
                "pool_index": pool_index,
                "energy_hartree": parsed.energy_hartree,
                "relative_to_best_retained_initial_geometry_kj_mol": (
                    parsed.energy_hartree - reference
                )
                * HARTREE_TO_KJ_MOL,
                "mmff_relative_energy_kj_mol": record.get(
                    "mmff_relative_energy_kj_mol",
                    recommended[pool_index].get("mmff_relative_energy_kj_mol"),
                ),
                "provenance": provenance,
                "output": str(output_path.resolve()),
            }
        )
    summary = {
        "generated_at_utc": now(),
        "calculation": "B3LYP/def2-SVP single point on retained MMFF94s geometry",
        "charge": 0,
        "multiplicity": 1,
        "geometry_optimization_performed": False,
        "interpretation": (
            "Screening of omitted basins only; unrelaxed energies are not conformer free energies."
        ),
        "retained_initial_geometry_reference_hartree": reference,
        "retained_initial_geometry_energies": retained,
        "omitted_cluster_results": results,
    }
    write_json(SUMMARY_PATH, summary)
    state.update(status="SUCCESS", finished_at_utc=now(), current_pool_index=None)
    write_json(STATUS_PATH, state)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

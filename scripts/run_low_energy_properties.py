"""Run non-destructive property-only ORCA jobs on low-energy Bz/SB structures."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from molecule2fbx.comparison import parse_existing_orca_electronic
from molecule2fbx.frequency import read_xyz_model
from molecule2fbx.quantum.orca import parse_orca_output


WORKSPACE = Path(__file__).resolve().parents[1]
COMPARISON_DIR = WORKSPACE / "outputs" / "Bz_vs_SB_preliminary_comparison"
OUTPUT_ROOT = COMPARISON_DIR / "electronic_properties"
STATUS_PATH = COMPARISON_DIR / "electronic_property_status.json"
SUMMARY_PATH = COMPARISON_DIR / "electronic_property_results.json"
ORCA = Path(os.environ.get("ORCA_EXECUTABLE", r"C:\ORCA_6.1.1\orca.exe"))
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


def selected_conformers(name: str) -> list[dict]:
    payload = json.loads(
        (WORKSPACE / "outputs" / name / "ensemble.json").read_text(encoding="utf-8")
    )
    entries = payload["final_conformer_ensemble"]
    minimum = min(float(entry["dft_energy_hartree"]) for entry in entries)
    selected = []
    for entry in sorted(entries, key=lambda item: float(item["dft_energy_hartree"])):
        delta = (float(entry["dft_energy_hartree"]) - minimum) * HARTREE_TO_KJ_MOL
        if delta <= 5.0 + 1.0e-9 and len(selected) < 5:
            selected.append({**entry, "relative_energy_kj_mol": delta})
    if name == "1SB-LSD_RR":
        followup_path = COMPARISON_DIR / "sb_additional_dft_opt_result.json"
        if followup_path.is_file():
            followup = json.loads(followup_path.read_text(encoding="utf-8"))
            delta = float(followup["relative_to_existing_current_best_kj_mol"])
            if delta <= 5.0 + 1.0e-9:
                selected.append(
                    {
                        "conformer_id": "conf011_pool109",
                        "conformer_index": 10,
                        "dft_energy_hartree": float(followup["energy_hartree"]),
                        "relative_energy_kj_mol": delta,
                        "calculation_directory": followup["calculation_directory"],
                    }
                )
    return selected


def render_input(model) -> str:
    lines = [
        "! B3LYP def2-SVP TightSCF MOREAD MBIS CHELPG FMOPop",
        "",
        '%moinp "start.gbw"',
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
        f"  {atom.element:<3} {atom.x: .12f} {atom.y: .12f} {atom.z: .12f}"
        for atom in model.atoms
    )
    lines.extend(["*", ""])
    return "\n".join(lines)


def compact_properties(parsed: dict) -> dict:
    return {
        "mbis_charges": parsed.get("mbis_charges"),
        "chelpg_charges": parsed.get("chelpg_charges"),
        "mayer_bond_orders": {
            f"{first}-{second}": value
            for (first, second), value in parsed.get("mayer_bond_orders", {}).items()
        },
        "dipole_magnitude_debye": parsed.get("dipole_magnitude_debye"),
        "homo": parsed.get("homo"),
        "lumo": parsed.get("lumo"),
        "homo_lumo_gap_ev": parsed.get("homo_lumo_gap_ev"),
        "frontier_orbital_local_contributions": parsed.get(
            "frontier_orbital_local_contributions"
        ),
        "frontier_orbital_local_contributions_reason": parsed.get(
            "frontier_orbital_local_contributions_reason"
        ),
    }


def main() -> int:
    if not ORCA.is_file():
        raise SystemExit(f"ORCA executable not found: {ORCA}")
    targets = []
    for name in ("1Bz-LSD_RR", "1SB-LSD_RR"):
        targets.extend((name, entry) for entry in selected_conformers(name))
    state = {
        "status": "RUNNING",
        "started_at_utc": now(),
        "finished_at_utc": None,
        "current_job": None,
        "jobs": [],
    }
    write_json(STATUS_PATH, state)
    env = os.environ.copy()
    env["PATH"] = rf"C:\Program Files\Microsoft MPI\Bin;{ORCA.parent};" + env.get("PATH", "")
    results = []
    for name, entry in targets:
        index = int(entry["conformer_index"])
        stem = f"conformer_{index + 1:03d}"
        source_dir = Path(entry["calculation_directory"])
        source_xyz = source_dir / f"{stem}.xyz"
        source_gbw = source_dir / f"{stem}.gbw"
        if not source_xyz.is_file() or not source_gbw.is_file():
            raise RuntimeError(f"Missing retained XYZ/GBW for {name} {entry['conformer_id']}")
        model = read_xyz_model(source_xyz)
        job_dir = OUTPUT_ROOT / name / str(entry["conformer_id"])
        job_dir.mkdir(parents=True, exist_ok=True)
        input_path = job_dir / "property.inp"
        output_path = job_dir / "property.out"
        start_gbw = job_dir / "start.gbw"
        if output_path.is_file():
            parse_orca_output(
                output_path.read_text(encoding="utf-8", errors="replace"),
                require_geometry_convergence=False,
            )
            provenance = "reused_completed_property_job"
        else:
            unexpected = [
                path
                for path in job_dir.iterdir()
                if path.name not in {"property.inp", "start.gbw"}
            ]
            if unexpected:
                raise RuntimeError(
                    f"Refusing to overwrite incomplete/non-empty job directory: {job_dir}"
                )
            if not start_gbw.is_file():
                shutil.copy2(source_gbw, start_gbw)
            input_path.write_text(render_input(model), encoding="utf-8")
            state["current_job"] = f"{name} {entry['conformer_id']}"
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
            parse_orca_output(
                output_path.read_text(encoding="utf-8", errors="replace"),
                require_geometry_convergence=False,
            )
            if process.returncode != 0:
                raise RuntimeError(f"ORCA exited with code {process.returncode}")
            provenance = "computed_this_run"
            state["jobs"].append(
                {
                    "molecule": name,
                    "conformer_id": entry["conformer_id"],
                    "started_at_utc": started.isoformat().replace("+00:00", "Z"),
                    "finished_at_utc": now(),
                    "provenance": provenance,
                }
            )
            write_json(STATUS_PATH, state)
        parsed = parse_existing_orca_electronic(output_path)
        results.append(
            {
                "molecule": name,
                "conformer_id": entry["conformer_id"],
                "conformer_index": index,
                "relative_energy_kj_mol": entry["relative_energy_kj_mol"],
                "geometry_source": str(source_xyz.resolve()),
                "wavefunction_start_source": str(source_gbw.resolve()),
                "output": str(output_path.resolve()),
                "provenance": provenance,
                "properties": compact_properties(parsed),
            }
        )
    summary = {
        "generated_at_utc": now(),
        "calculation": "B3LYP/def2-SVP property-only single point",
        "simple_input": "B3LYP def2-SVP TightSCF MOREAD MBIS CHELPG FMOPop",
        "selection": (
            "Electronic delta-E <= 5 kJ/mol, maximum 5 original-ensemble "
            "structures per molecule, plus the recorded SB pool-109 follow-up"
        ),
        "geometry_optimization_performed": False,
        "frequency_calculation_performed": False,
        "charge": 0,
        "multiplicity": 1,
        "results": results,
    }
    write_json(SUMMARY_PATH, summary)
    state.update(status="SUCCESS", finished_at_utc=now(), current_job=None)
    write_json(STATUS_PATH, state)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

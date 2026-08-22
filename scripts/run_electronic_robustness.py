"""Non-destructive paired electronic-property robustness calculations.

The script never performs Opt or Freq.  It first fills B3LYP/def2-SVP gas-phase
property data for the six matched Bz/SB conformer pairs, then evaluates the two
lowest-energy matched pairs with B3LYP-D3BJ/def2-TZVP and CPCM(water).
Completed ORCA outputs are validated and reused; incomplete directories are not
overwritten.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Optional

from molecule2fbx.comparison import (
    ExistingConformer,
    load_existing_conformers,
    load_targeted_followup_conformer,
    parse_existing_orca_electronic,
    structural_descriptors,
)
from molecule2fbx.quantum.orca import parse_orca_output


WORKSPACE = Path(__file__).resolve().parents[1]
COMPARISON_DIR = WORKSPACE / "outputs" / "Bz_vs_SB_preliminary_comparison"
OUTPUT_ROOT = COMPARISON_DIR / "electronic_robustness"
STATUS_PATH = OUTPUT_ROOT / "status.json"
RESULT_PATH = OUTPUT_ROOT / "results.json"
ORCA = Path(os.environ.get("ORCA_EXECUTABLE", r"C:\ORCA_6.1.1\orca.exe"))
NPROCS = 16
MAXCORE_MB = 1000

CONDITIONS = {
    "baseline_b3lyp_def2_svp_gas": {
        "simple_input": (
            "B3LYP def2-SVP TightSCF MOREAD MBIS CHELPG FMOPop"
        ),
        "selection": "all_six_matched_pairs",
    },
    "b3lyp_d3bj_def2_tzvp_cpcm_water": {
        "simple_input": (
            "B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF MOREAD "
            "MBIS CHELPG FMOPop CPCM(Water)"
        ),
        "selection": "two_low_energy_matched_pairs",
    },
}


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


def load_pairs() -> list[dict]:
    payload = json.loads(
        (COMPARISON_DIR / "paired_continuous_steric_access.json").read_text(
            encoding="utf-8"
        )
    )
    pairs = payload["pair_selection"]["selected_pairs"]
    if len(pairs) != 6:
        raise RuntimeError(f"Expected six matched pairs, found {len(pairs)}")
    return pairs


def conformers() -> tuple[dict[str, ExistingConformer], dict[str, ExistingConformer]]:
    _, bz = load_existing_conformers(WORKSPACE / "outputs" / "1Bz-LSD_RR")
    sb_payload, sb = load_existing_conformers(
        WORKSPACE / "outputs" / "1SB-LSD_RR"
    )
    followup = load_targeted_followup_conformer(sb_payload, COMPARISON_DIR)
    if followup is not None:
        sb.append(followup)
    return (
        {str(item.entry["conformer_id"]): item for item in bz},
        {str(item.entry["conformer_id"]): item for item in sb},
    )


def render_input(model, simple_input: str) -> str:
    lines = [
        f"! {simple_input}",
        "",
        '%moinp "start.gbw"',
        "%scf",
        "  MaxIter 300",
        "end",
        "%pal",
        f"  nprocs {NPROCS}",
        "end",
        f"%maxcore {MAXCORE_MB}",
        "",
        "* xyz 0 1",
    ]
    lines.extend(
        f"  {atom.element:<3} {atom.x: .12f} {atom.y: .12f} {atom.z: .12f}"
        for atom in model.atoms
    )
    lines.extend(["*", ""])
    return "\n".join(lines)


def existing_baseline_output(molecule: str, conformer_id: str) -> Optional[Path]:
    candidate = (
        COMPARISON_DIR
        / "electronic_properties"
        / molecule
        / conformer_id
        / "property.out"
    )
    if not candidate.is_file():
        return None
    parse_orca_output(
        candidate.read_text(encoding="utf-8", errors="replace"),
        require_geometry_convergence=False,
    )
    return candidate.resolve()


def source_gbw(conformer: ExistingConformer) -> Path:
    stem = f"conformer_{int(conformer.entry['conformer_index']) + 1:03d}"
    path = conformer.xyz_path.parent / f"{stem}.gbw"
    if not path.is_file():
        raise RuntimeError(f"Missing retained wavefunction: {path}")
    return path


def run_job(
    condition_name: str,
    condition: Mapping[str, str],
    molecule: str,
    conformer: ExistingConformer,
    state: dict,
    env: Mapping[str, str],
) -> tuple[Path, str]:
    conformer_id = str(conformer.entry["conformer_id"])
    if condition_name == "baseline_b3lyp_def2_svp_gas":
        existing = existing_baseline_output(molecule, conformer_id)
        if existing is not None:
            return existing, "reused_existing_baseline_property_job"

    job_dir = OUTPUT_ROOT / condition_name / molecule / conformer_id
    input_path = job_dir / "property.inp"
    output_path = job_dir / "property.out"
    start_gbw = job_dir / "start.gbw"
    if output_path.is_file():
        parse_orca_output(
            output_path.read_text(encoding="utf-8", errors="replace"),
            require_geometry_convergence=False,
        )
        return output_path.resolve(), "reused_completed_robustness_job"

    if job_dir.exists():
        unexpected = [
            path
            for path in job_dir.iterdir()
            if path.name not in {"property.inp", "start.gbw"}
        ]
        if unexpected:
            raise RuntimeError(
                f"Refusing to overwrite incomplete job directory: {job_dir}"
            )
    job_dir.mkdir(parents=True, exist_ok=True)
    if not start_gbw.is_file():
        shutil.copy2(source_gbw(conformer), start_gbw)
    input_path.write_text(
        render_input(conformer.model, str(condition["simple_input"])),
        encoding="utf-8",
    )
    state["current_job"] = {
        "condition": condition_name,
        "molecule": molecule,
        "conformer_id": conformer_id,
        "started_at_utc": now(),
    }
    write_json(STATUS_PATH, state)
    process = subprocess.run(
        [str(ORCA), input_path.name],
        cwd=job_dir,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=dict(env),
        check=False,
        timeout=7200,
    )
    output_path.write_text(
        process.stdout
        + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
        encoding="utf-8",
    )
    parse_orca_output(
        output_path.read_text(encoding="utf-8", errors="replace"),
        require_geometry_convergence=False,
    )
    if process.returncode != 0:
        raise RuntimeError(
            f"ORCA exited with {process.returncode}: {condition_name} "
            f"{molecule} {conformer_id}"
        )
    state["completed_jobs"].append(
        {
            **state["current_job"],
            "finished_at_utc": now(),
            "provenance": "computed_this_run",
        }
    )
    state["current_job"] = None
    write_json(STATUS_PATH, state)
    return output_path.resolve(), "computed_this_run"


def compact_result(
    molecule: str,
    conformer: ExistingConformer,
    condition_name: str,
    output_path: Path,
    provenance: str,
) -> dict:
    parsed = parse_existing_orca_electronic(output_path)
    structure = structural_descriptors(conformer.model)
    indices = structure["atom_indices"]
    c = int(indices["benzoyl_carbonyl_c"])
    o = int(indices["benzoyl_carbonyl_o"])
    n = int(indices["benzoyl_amide_n"])
    ipso = int(indices["benzoyl_ipso_c"])

    def charge(kind: str, index: int):
        values = parsed.get(kind)
        return values.get(index) if isinstance(values, dict) else None

    mayer = parsed.get("mayer_bond_orders", {})
    frontier = parsed.get("frontier_orbital_local_contributions")
    populations = (
        frontier.get("atom_populations", {})
        if isinstance(frontier, dict)
        else {}
    )

    def local_population(kind: str):
        if not isinstance(populations, dict):
            return None
        values = [
            populations[index].get(kind)
            for index in (c, o, n, ipso)
            if index in populations and isinstance(populations[index], dict)
        ]
        if len(values) != 4 or any(value is None for value in values):
            return None
        return sum(float(value) for value in values)

    return {
        "condition": condition_name,
        "molecule": molecule,
        "conformer_id": conformer.entry["conformer_id"],
        "relative_energy_kj_mol": conformer.entry.get("relative_energy_kj_mol"),
        "geometry_source": str(conformer.xyz_path),
        "output": str(output_path),
        "provenance": provenance,
        "properties": {
            "mbis_charge_carbonyl_c": charge("mbis_charges", c),
            "mbis_charge_carbonyl_o": charge("mbis_charges", o),
            "mbis_charge_amide_n": charge("mbis_charges", n),
            "chelpg_charge_carbonyl_c": charge("chelpg_charges", c),
            "chelpg_charge_carbonyl_o": charge("chelpg_charges", o),
            "chelpg_charge_amide_n": charge("chelpg_charges", n),
            "loewdin_charge_carbonyl_c": charge("loewdin_charges", c),
            "loewdin_charge_carbonyl_o": charge("loewdin_charges", o),
            "loewdin_charge_amide_n": charge("loewdin_charges", n),
            "mulliken_charge_carbonyl_c": charge("mulliken_charges", c),
            "mulliken_charge_carbonyl_o": charge("mulliken_charges", o),
            "mulliken_charge_amide_n": charge("mulliken_charges", n),
            "mayer_bond_order_c_o": mayer.get(tuple(sorted((c, o)))),
            "mayer_bond_order_c_n": mayer.get(tuple(sorted((c, n)))),
            "homo_energy_ev": (
                parsed["homo"]["energy_ev"] if parsed.get("homo") else None
            ),
            "lumo_energy_ev": (
                parsed["lumo"]["energy_ev"] if parsed.get("lumo") else None
            ),
            "homo_lumo_gap_ev": parsed.get("homo_lumo_gap_ev"),
            "homo_loewdin_population_benzoyl_center": local_population(
                "homo_loewdin"
            ),
            "lumo_loewdin_population_benzoyl_center": local_population(
                "lumo_loewdin"
            ),
            "dipole_magnitude_debye": parsed.get("dipole_magnitude_debye"),
        },
    }


def main() -> int:
    if not ORCA.is_file():
        raise SystemExit(f"ORCA executable not found: {ORCA}")
    pairs = load_pairs()
    bz, sb = conformers()
    low_energy_pairs = pairs[:2]
    targets = []
    for condition_name, condition in CONDITIONS.items():
        selected = (
            pairs
            if condition["selection"] == "all_six_matched_pairs"
            else low_energy_pairs
        )
        for pair in selected:
            targets.extend(
                (
                    (condition_name, condition, "1Bz-LSD_RR", bz[str(pair["first_conformer_id"])]),
                    (condition_name, condition, "1SB-LSD_RR", sb[str(pair["second_conformer_id"])]),
                )
            )
    state = {
        "status": "RUNNING",
        "started_at_utc": now(),
        "finished_at_utc": None,
        "current_job": None,
        "completed_jobs": [],
        "planned_job_count_including_reuse": len(targets),
        "estimated_wall_time_hours": "1-2",
        "reason": (
            "Test whether the observed MBIS and local-LUMO shifts survive "
            "conformer pairing and a larger basis/dispersion/aqueous continuum."
        ),
        "effect_on_hypothesis": (
            "Persistence with one sign supports a model-robust TMS electronic "
            "effect; sign reversal or method-scale variability makes it unresolved."
        ),
    }
    write_json(STATUS_PATH, state)
    env = os.environ.copy()
    env["PATH"] = (
        rf"C:\Program Files\Microsoft MPI\Bin;{ORCA.parent};"
        + env.get("PATH", "")
    )
    results = []
    try:
        for condition_name, condition, molecule, conformer in targets:
            output, provenance = run_job(
                condition_name,
                condition,
                molecule,
                conformer,
                state,
                env,
            )
            results.append(
                compact_result(
                    molecule,
                    conformer,
                    condition_name,
                    output,
                    provenance,
                )
            )
            write_json(
                RESULT_PATH,
                {
                    "generated_at_utc": now(),
                    "status": "RUNNING",
                    "conditions": CONDITIONS,
                    "geometry_optimization_performed": False,
                    "frequency_calculation_performed": False,
                    "charge": 0,
                    "multiplicity": 1,
                    "pairs": pairs,
                    "results": results,
                },
            )
    except Exception as exc:
        state.update(status="FAILED", finished_at_utc=now(), error=str(exc))
        write_json(STATUS_PATH, state)
        raise
    summary = {
        "generated_at_utc": now(),
        "status": "SUCCESS",
        "conditions": CONDITIONS,
        "geometry_optimization_performed": False,
        "frequency_calculation_performed": False,
        "charge": 0,
        "multiplicity": 1,
        "pairs": pairs,
        "results": results,
    }
    write_json(RESULT_PATH, summary)
    state.update(status="SUCCESS", finished_at_utc=now(), current_job=None)
    write_json(STATUS_PATH, state)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

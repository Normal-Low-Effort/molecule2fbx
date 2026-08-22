"""Minimal fixed-geometry validation of the proposed TMS electronic effect.

This script performs property-only ORCA single points for two low-energy matched
Bz/SB conformer pairs.  Existing B3LYP-D3BJ/def2-TZVP/CPCM(water) results are
reused.  New calculations add an independent PBE0 condition and a fixed-scaffold
counterfactual in which SiMe3 is replaced by H along the original C--Si vector.
No geometry optimization or frequency calculation is performed, and completed
outputs are never overwritten.
"""

from __future__ import annotations

import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

from molecule2fbx.comparison import (
    ExistingConformer,
    load_existing_conformers,
    load_targeted_followup_conformer,
    parse_existing_orca_electronic,
    replace_terminal_tms_with_hydrogen,
    structural_descriptors,
)
from molecule2fbx.model import MoleculeModel
from molecule2fbx.quantum.orca import parse_orca_output


WORKSPACE = Path(__file__).resolve().parents[1]
COMPARISON_DIR = WORKSPACE / "outputs" / "Bz_vs_SB_preliminary_comparison"
PRIOR_ROOT = COMPARISON_DIR / "electronic_robustness"
OUTPUT_ROOT = COMPARISON_DIR / "hypothesis2_validation"
STATUS_PATH = OUTPUT_ROOT / "status.json"
RESULT_PATH = OUTPUT_ROOT / "results.json"
ORCA = Path(os.environ.get("ORCA_EXECUTABLE", r"C:\ORCA_6.1.1\orca.exe"))
NPROCS = 16
MAXCORE_MB = 1000

CONDITIONS = {
    "b3lyp_d3bj_def2_tzvp_cpcm_water": {
        "simple_input": (
            "B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF "
            "MBIS CHELPG FMOPop CPCM(Water)"
        ),
        "functional": "B3LYP",
        "basis": "def2-TZVP",
        "solvent": "CPCM(Water)",
    },
    "pbe0_d3bj_def2_tzvp_cpcm_water": {
        "simple_input": (
            "PBE0 D3BJ def2-TZVP def2/J RIJCOSX TightSCF "
            "MBIS CHELPG FMOPop CPCM(Water)"
        ),
        "functional": "PBE0",
        "basis": "def2-TZVP",
        "solvent": "CPCM(Water)",
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
    if len(pairs) < 2:
        raise RuntimeError("At least two matched conformer pairs are required")
    return pairs[:2]


def conformers() -> tuple[dict[str, ExistingConformer], dict[str, ExistingConformer]]:
    _, bz = load_existing_conformers(WORKSPACE / "outputs" / "1Bz-LSD_RR")
    sb_payload, sb = load_existing_conformers(WORKSPACE / "outputs" / "1SB-LSD_RR")
    followup = load_targeted_followup_conformer(sb_payload, COMPARISON_DIR)
    if followup is not None:
        sb.append(followup)
    return (
        {str(item.entry["conformer_id"]): item for item in bz},
        {str(item.entry["conformer_id"]): item for item in sb},
    )


def render_input(model: MoleculeModel, simple_input: str) -> str:
    lines = [
        f"! {simple_input}",
        "",
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


def prior_b3lyp_output(molecule: str, conformer_id: str) -> Path:
    path = (
        PRIOR_ROOT
        / "b3lyp_d3bj_def2_tzvp_cpcm_water"
        / molecule
        / conformer_id
        / "property.out"
    )
    if not path.is_file():
        raise RuntimeError(f"Missing completed B3LYP property output: {path}")
    parse_orca_output(
        path.read_text(encoding="utf-8", errors="replace"),
        require_geometry_convergence=False,
    )
    return path.resolve()


def run_job(
    *,
    condition_name: str,
    condition: Mapping[str, str],
    pair_id: str,
    species: str,
    model: MoleculeModel,
    state: dict,
    env: Mapping[str, str],
) -> tuple[Path, str]:
    job_dir = OUTPUT_ROOT / condition_name / pair_id / species
    input_path = job_dir / "property.inp"
    output_path = job_dir / "property.out"
    if output_path.is_file():
        parse_orca_output(
            output_path.read_text(encoding="utf-8", errors="replace"),
            require_geometry_convergence=False,
        )
        return output_path.resolve(), "reused_completed_hypothesis2_job"
    if job_dir.exists() and any(job_dir.iterdir()):
        raise RuntimeError(f"Refusing to overwrite incomplete job directory: {job_dir}")

    job_dir.mkdir(parents=True, exist_ok=True)
    input_path.write_text(
        render_input(model, str(condition["simple_input"])), encoding="utf-8"
    )
    state["current_job"] = {
        "condition": condition_name,
        "pair_id": pair_id,
        "species": species,
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
            f"ORCA exited with {process.returncode}: {condition_name} {pair_id} {species}"
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
    *,
    condition_name: str,
    pair_id: str,
    species: str,
    conformer_id: str,
    model: MoleculeModel,
    geometry_source: Path,
    output_path: Path,
    provenance: str,
) -> dict:
    parsed = parse_existing_orca_electronic(output_path)
    structure = structural_descriptors(model)
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
    populations = frontier.get("atom_populations", {}) if isinstance(frontier, dict) else {}

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
        "pair_id": pair_id,
        "species": species,
        "conformer_id": conformer_id,
        "geometry_source": str(geometry_source.resolve()),
        "geometry_optimized_for_this_species": species != "sb_tms_to_h_fixed",
        "output": str(output_path),
        "provenance": provenance,
        "atom_count": len(model.atoms),
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
            "mayer_bond_order_c_o": mayer.get(tuple(sorted((c, o)))),
            "mayer_bond_order_c_n": mayer.get(tuple(sorted((c, n)))),
            "homo_energy_ev": parsed["homo"]["energy_ev"] if parsed.get("homo") else None,
            "lumo_energy_ev": parsed["lumo"]["energy_ev"] if parsed.get("lumo") else None,
            "homo_lumo_gap_ev": parsed.get("homo_lumo_gap_ev"),
            "homo_loewdin_population_benzoyl_center": local_population("homo_loewdin"),
            "lumo_loewdin_population_benzoyl_center": local_population("lumo_loewdin"),
            "dipole_magnitude_debye": parsed.get("dipole_magnitude_debye"),
        },
    }


def main() -> int:
    if not ORCA.is_file():
        raise SystemExit(f"ORCA executable not found: {ORCA}")
    pairs = load_pairs()
    bz, sb = conformers()
    state = {
        "status": "RUNNING",
        "started_at_utc": now(),
        "finished_at_utc": None,
        "current_job": None,
        "completed_jobs": [],
        "planned_records_including_reuse": 12,
        "planned_new_orca_jobs": 8,
        "estimated_wall_time_hours": "3-4 based on observed def2-TZVP property jobs",
        "reason": (
            "Test whether the small carbonyl electronic shift is functional-independent "
            "and directly attributable to TMS at fixed common-scaffold geometry."
        ),
        "effect_on_hypothesis": (
            "A shared direction across B3LYP/PBE0, two matched conformers, and the "
            "fixed-geometry TMS-to-H control supports a small model-robust electronic effect."
        ),
        "geometry_optimization_performed": False,
        "frequency_calculation_performed": False,
    }
    write_json(STATUS_PATH, state)
    env = os.environ.copy()
    env["PATH"] = (
        rf"C:\Program Files\Microsoft MPI\Bin;{ORCA.parent};" + env.get("PATH", "")
    )
    results = []
    try:
        for pair_number, pair in enumerate(pairs, start=1):
            pair_id = f"pair{pair_number:03d}"
            bz_conformer = bz[str(pair["first_conformer_id"])]
            sb_conformer = sb[str(pair["second_conformer_id"])]
            counterfactual = replace_terminal_tms_with_hydrogen(sb_conformer.model)
            species_models = {
                "bz_optimized": (
                    bz_conformer.model,
                    bz_conformer.xyz_path,
                    str(bz_conformer.entry["conformer_id"]),
                    "1Bz-LSD_RR",
                ),
                "sb_optimized": (
                    sb_conformer.model,
                    sb_conformer.xyz_path,
                    str(sb_conformer.entry["conformer_id"]),
                    "1SB-LSD_RR",
                ),
                "sb_tms_to_h_fixed": (
                    counterfactual,
                    sb_conformer.xyz_path,
                    str(sb_conformer.entry["conformer_id"]),
                    None,
                ),
            }
            for condition_name, condition in CONDITIONS.items():
                for species, (model, geometry_source, conformer_id, molecule) in species_models.items():
                    if condition_name.startswith("b3lyp") and molecule is not None:
                        output_path = prior_b3lyp_output(molecule, conformer_id)
                        provenance = "reused_prior_b3lyp_property_job"
                    else:
                        output_path, provenance = run_job(
                            condition_name=condition_name,
                            condition=condition,
                            pair_id=pair_id,
                            species=species,
                            model=model,
                            state=state,
                            env=env,
                        )
                    results.append(
                        compact_result(
                            condition_name=condition_name,
                            pair_id=pair_id,
                            species=species,
                            conformer_id=conformer_id,
                            model=model,
                            geometry_source=geometry_source,
                            output_path=output_path,
                            provenance=provenance,
                        )
                    )
                    write_json(
                        RESULT_PATH,
                        {
                            "status": "RUNNING",
                            "generated_at_utc": now(),
                            "conditions": CONDITIONS,
                            "pairs": pairs,
                            "charge": 0,
                            "multiplicity": 1,
                            "geometry_optimization_performed": False,
                            "frequency_calculation_performed": False,
                            "results": results,
                        },
                    )
    except Exception as exc:
        state.update(status="FAILED", finished_at_utc=now(), error=str(exc))
        write_json(STATUS_PATH, state)
        raise

    payload = {
        "status": "SUCCESS",
        "generated_at_utc": now(),
        "conditions": CONDITIONS,
        "pairs": pairs,
        "charge": 0,
        "multiplicity": 1,
        "geometry_optimization_performed": False,
        "frequency_calculation_performed": False,
        "counterfactual": {
            "transformation": "SiMe3 replaced by H along the source aryl-C--Si vector",
            "common_atom_coordinates_retained": True,
            "replacement_c_h_length_angstrom": 1.085,
            "energy_comparison_allowed": False,
            "purpose": "local electronic descriptors only",
        },
        "results": results,
    }
    write_json(RESULT_PATH, payload)
    state.update(status="SUCCESS", finished_at_utc=now(), current_job=None)
    write_json(STATUS_PATH, state)
    print(
        json.dumps(
            {
                "status": payload["status"],
                "record_count": len(results),
                "new_job_count": len(state["completed_jobs"]),
                "results": str(RESULT_PATH),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

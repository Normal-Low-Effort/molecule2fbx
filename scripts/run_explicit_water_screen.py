"""One-water microsolvation screen for matched Bz/SB conformers.

Multiple water placements around the N1 carbonyl are optimized with
GFN2-xTB/ALPB(water).  This intentionally small screen tests local first-shell
geometry only; it is not explicit-solvent MD, a hydration free energy, or a
reaction-path calculation.
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "work"))

from analyze_hypotheses4_5 import atom_groups  # noqa: E402
from molecule2fbx.comparison import _rdkit_model_from_xyz  # noqa: E402
from molecule2fbx.frequency import read_xyz_model  # noqa: E402
from molecule2fbx.model import Atom, MoleculeModel, validate_model  # noqa: E402


XTB = Path(os.environ.get("XTB_EXECUTABLE", r"C:\ORCA_6.1.1\xtb-6.7.1pre\xtb.exe"))
PROTONATION = ROOT / "outputs" / "Bz_vs_SB_preliminary_comparison" / "protonation_screen"
OUTPUT = ROOT / "outputs" / "Bz_vs_SB_preliminary_comparison" / "explicit_water_screen"
STATUS = OUTPUT / "status.json"
RESULTS = OUTPUT / "results.json"
FINDINGS = OUTPUT / "FINDINGS_JA.md"
ENERGY_RE = re.compile(
    r"TOTAL ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh", re.IGNORECASE
)
HARTREE_TO_KJ_MOL = 2625.4996394799

PAIRS = (
    ("pair_low_energy", "conf008", "conf002"),
    ("pair_close_geometry", "conf005", "conf008"),
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def write_xyz(path: Path, elements, coordinates, comment: str) -> None:
    lines = [str(len(elements)), comment]
    lines.extend(
        f"{element:<3} {xyz[0]: .12f} {xyz[1]: .12f} {xyz[2]: .12f}"
        for element, xyz in zip(elements, coordinates)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_completed(job_dir: Path):
    output = job_dir / "xtb.out"
    optimized = job_dir / "xtbopt.xyz"
    text = output.read_text(encoding="utf-8", errors="replace")
    energies = [float(value) for value in ENERGY_RE.findall(text)]
    if "normal termination of xtb" not in text or not energies or not optimized.is_file():
        raise RuntimeError(f"Incomplete xTB job: {job_dir}")
    return energies[-1], optimized


def run_job(job_dir: Path, elements, coordinates, charge: int):
    job_dir.mkdir(parents=True, exist_ok=True)
    output = job_dir / "xtb.out"
    if output.is_file():
        energy, optimized = parse_completed(job_dir)
        return energy, optimized, "reused_completed_xtb_job"
    if any(job_dir.iterdir()):
        raise RuntimeError(f"Refusing to overwrite incomplete directory: {job_dir}")
    input_path = job_dir / "input.xyz"
    write_xyz(input_path, elements, coordinates, f"charge={charge}; explicit H2O + ALPB(water)")
    process = subprocess.run(
        [
            str(XTB), input_path.name, "--gfn", "2", "--opt", "tight",
            "--alpb", "water", "--chrg", str(charge), "--uhf", "0",
            "--parallel", "16",
        ],
        cwd=job_dir,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env={**os.environ, "OMP_NUM_THREADS": "16"},
        timeout=3600,
        check=False,
    )
    output.write_text(
        process.stdout + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
        encoding="utf-8",
    )
    energy, optimized = parse_completed(job_dir)
    if process.returncode != 0:
        raise RuntimeError(f"xTB returned {process.returncode}: {job_dir}")
    return energy, optimized, "computed_this_run"


def run_single_point(job_dir: Path, elements, coordinates, charge: int):
    job_dir.mkdir(parents=True, exist_ok=True)
    output = job_dir / "sp.out"
    if output.is_file():
        text = output.read_text(encoding="utf-8", errors="replace")
        energies = [float(value) for value in ENERGY_RE.findall(text)]
        if "normal termination of xtb" not in text or not energies:
            raise RuntimeError(f"Incomplete xTB single point: {job_dir}")
        return energies[-1], "reused_completed_xtb_job"
    if any(job_dir.iterdir()):
        raise RuntimeError(f"Refusing to overwrite incomplete SP directory: {job_dir}")
    input_path = job_dir / "input.xyz"
    write_xyz(input_path, elements, coordinates, f"fragment SP charge={charge}; ALPB(water)")
    process = subprocess.run(
        [
            str(XTB), input_path.name, "--gfn", "2", "--sp", "--alpb", "water",
            "--chrg", str(charge), "--uhf", "0", "--parallel", "16",
        ],
        cwd=job_dir,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env={**os.environ, "OMP_NUM_THREADS": "16"},
        timeout=3600,
        check=False,
    )
    output.write_text(
        process.stdout + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
        encoding="utf-8",
    )
    text = output.read_text(encoding="utf-8", errors="replace")
    energies = [float(value) for value in ENERGY_RE.findall(text)]
    if process.returncode != 0 or "normal termination of xtb" not in text or not energies:
        raise RuntimeError(f"xTB fragment SP failed: {job_dir}")
    return energies[-1], "computed_this_run"


def orthogonal(vector: np.ndarray) -> np.ndarray:
    trial = np.asarray((1.0, 0.0, 0.0))
    if abs(float(np.dot(vector, trial))) > 0.85:
        trial = np.asarray((0.0, 1.0, 0.0))
    result = trial - np.dot(trial, vector) * vector
    return result / np.linalg.norm(result)


def water_geometry(oxygen: np.ndarray, bisector: np.ndarray, plane_axis: np.ndarray):
    bisector = bisector / np.linalg.norm(bisector)
    plane_axis = plane_axis - np.dot(plane_axis, bisector) * bisector
    if np.linalg.norm(plane_axis) < 1.0e-8:
        plane_axis = orthogonal(bisector)
    plane_axis /= np.linalg.norm(plane_axis)
    half_angle = math.radians(104.5 / 2.0)
    first = math.cos(half_angle) * bisector + math.sin(half_angle) * plane_axis
    second = math.cos(half_angle) * bisector - math.sin(half_angle) * plane_axis
    return [oxygen, oxygen + 0.96 * first, oxygen + 0.96 * second]


def placements(model):
    groups = atom_groups(model)
    c = np.asarray(model.atoms[groups["carbonyl_c"]].position, dtype=float)
    o = np.asarray(model.atoms[groups["carbonyl_o"]].position, dtype=float)
    n = np.asarray(model.atoms[groups["amide_n"]].position, dtype=float)
    axis = (o - c) / np.linalg.norm(o - c)
    in_plane = n - c - np.dot(n - c, axis) * axis
    in_plane /= np.linalg.norm(in_plane)
    normal = np.cross(axis, in_plane)
    normal /= np.linalg.norm(normal)
    result = []

    water_o = o + 2.85 * axis
    # One O-H points approximately toward the carbonyl oxygen.
    result.append(("carbonyl_o_hbond", water_geometry(water_o, -axis, in_plane)))

    theta = math.radians(107.0)
    for sample in range(8):
        phi = 2.0 * math.pi * sample / 8.0
        radial = math.cos(phi) * in_plane + math.sin(phi) * normal
        direction = math.cos(theta) * axis + math.sin(theta) * radial
        water_o = c + 3.00 * direction
        # Hydrogens point away from C so the oxygen lone-pair side faces C.
        result.append(
            (
                f"attack_{sample:02d}",
                water_geometry(water_o, direction, np.cross(direction, axis)),
            )
        )
    return result


def optimized_solute_model(source, optimized_xyz: Path):
    raw = read_xyz_model(optimized_xyz)
    if len(raw.atoms) != len(source.atoms) + 3:
        raise ValueError("Complex XYZ atom count is not solute + H2O")
    atoms = tuple(
        Atom(index, source.atoms[index].element, *raw.atoms[index].position)
        for index in range(len(source.atoms))
    )
    model = MoleculeModel(
        source.cid, source.name, atoms, source.bonds, dict(source.metadata)
    )
    validate_model(model)
    return model, raw


def angle(first, center, second):
    a = first - center
    b = second - center
    cosine = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return math.degrees(math.acos(min(1.0, max(-1.0, cosine))))


def analyze_complex(source, optimized_xyz: Path):
    solute, raw = optimized_solute_model(source, optimized_xyz)
    groups = atom_groups(solute)
    c = np.asarray(solute.atoms[groups["carbonyl_c"]].position, dtype=float)
    o = np.asarray(solute.atoms[groups["carbonyl_o"]].position, dtype=float)
    n6 = np.asarray(solute.atoms[9].position, dtype=float)
    offset = len(solute.atoms)
    water_o = np.asarray(raw.atoms[offset].position, dtype=float)
    water_h = [np.asarray(raw.atoms[offset + index].position, dtype=float) for index in (1, 2)]
    water_o_h_distances = [
        float(np.linalg.norm(water_o - hydrogen)) for hydrogen in water_h
    ]
    intact_water = max(water_o_h_distances) <= 1.30
    c_ow = float(np.linalg.norm(c - water_o))
    o_ow = float(np.linalg.norm(o - water_o))
    o_hw = min(float(np.linalg.norm(o - hydrogen)) for hydrogen in water_h)
    n6_ow = float(np.linalg.norm(n6 - water_o))
    attack_angle = angle(o, c, water_o)
    tms = groups["tms"]
    water_tms_min = min(
        (float(np.linalg.norm(water_o - np.asarray(solute.atoms[index].position))) for index in tms),
        default=None,
    )
    if not intact_water:
        location = "proton_transfer_or_dissociated_water"
    elif o_hw <= 2.50:
        location = "carbonyl_o_hbond"
    elif c_ow <= 3.50 and 90.0 <= attack_angle <= 125.0:
        location = "carbonyl_attack_region"
    elif n6_ow <= 3.50:
        location = "n6_region"
    else:
        location = "other_or_dispersed"
    return {
        "final_location": location,
        "intact_water": intact_water,
        "water_o_h_distances_angstrom": water_o_h_distances,
        "carbonyl_c_water_o_distance_angstrom": c_ow,
        "carbonyl_o_water_o_distance_angstrom": o_ow,
        "carbonyl_o_nearest_water_h_distance_angstrom": o_hw,
        "o_c_water_o_angle_deg": attack_angle,
        "n6_water_o_distance_angstrom": n6_ow,
        "water_o_minimum_tms_atom_distance_angstrom": water_tms_min,
    }


def source_record(results, molecule, state, conformer_id):
    dataset = results["ensembles"][molecule][state]
    record = next(
        item for item in dataset["records"] if item["source_conformer_id"] == conformer_id
    )
    model = _rdkit_model_from_xyz(
        Path(record["optimized_xyz"]), dataset["smiles"], name=f"{molecule}_{state}"
    )
    return dataset, record, model


def aggregate(records):
    valid = [item for item in records if item["geometry"]["intact_water"]]
    local = [
        item for item in valid
        if item["geometry"]["final_location"] in ("carbonyl_o_hbond", "carbonyl_attack_region")
    ]
    formation = [item["apparent_cluster_formation_energy_kj_mol"] for item in local]
    interaction = [item["fragment_interaction_energy_kj_mol"] for item in local]
    return {
        "placement_count": len(records),
        "intact_water_count": len(valid),
        "proton_transfer_or_dissociation_count": len(records) - len(valid),
        "carbonyl_local_count": len(local),
        "carbonyl_local_fraction_of_intact": len(local) / len(valid) if valid else None,
        "location_counts": {
            location: sum(item["geometry"]["final_location"] == location for item in records)
            for location in (
                "carbonyl_o_hbond", "carbonyl_attack_region", "n6_region",
                "other_or_dispersed", "proton_transfer_or_dissociated_water",
            )
        },
        "best_apparent_carbonyl_local_cluster_formation_kj_mol": min(formation) if formation else None,
        "median_apparent_carbonyl_local_cluster_formation_kj_mol": (
            float(np.median(formation)) if formation else None
        ),
        "best_carbonyl_local_fragment_interaction_kj_mol": (
            min(interaction) if interaction else None
        ),
        "median_carbonyl_local_fragment_interaction_kj_mol": (
            float(np.median(interaction)) if interaction else None
        ),
    }


def write_findings(payload):
    rows = []
    for pair in payload["pairs"]:
        for state in ("neutral", "n6_protonated"):
            bz = pair["states"][state]["1Bz-LSD_RR"]["summary"]
            sb = pair["states"][state]["1SB-LSD_RR"]["summary"]
            first = bz["best_carbonyl_local_fragment_interaction_kj_mol"]
            second = sb["best_carbonyl_local_fragment_interaction_kj_mol"]
            difference = None if first is None or second is None else second - first
            rows.append(
                f"| {pair['pair_id']} | {state} | {bz['carbonyl_local_count']}/{bz['intact_water_count']} | "
                f"{sb['carbonyl_local_count']}/{sb['intact_water_count']} | "
                f"{'—' if first is None else f'{first:.2f}'} | "
                f"{'—' if second is None else f'{second:.2f}'} | "
                f"{'—' if difference is None else f'{difference:+.2f}'} |"
            )
    lines = [
        "# 仮説6 明示的1水分子microsolvation screen",
        "",
        "対応する2組のBz/SB配座について、neutral/N6H+それぞれのN1-carbonyl周辺へ水1分子を9方向から配置し、GFN2-xTB/ALPB(水)で最適化した。これは明示的溶媒MDでも加水分解反応計算でもない。O–Hが1.30 Åを超えた配置は水のproton transfer/dissociationとして局所水和集計から除外した。",
        "",
        "| pair | state | Bz carbonyl-local/intact | SB carbonyl-local/intact | Bz best Eint | SB best Eint | SB−Bz / kJ mol⁻¹ |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
        *rows,
        "",
        "Eintは最適化複合体座標で E(complex)−E(solute fragment)−E(water fragment) としたxTB相互作用energy proxy。溶質変形は別項に分離したが、標準状態、エントロピー、濃度、BSSE、多水分子効果を含まないため結合自由エネルギーとは呼ばない。",
        "",
        "このscreenでcarbonyl-local保持率と対応pairの差が一貫しない場合、TMSによる局所第一水和殻効果はこの解像度では検出不能とする。差が一貫しても、十分な水和ensembleと反応障壁で再検証するまで速度差へ結び付けない。",
    ]
    FINDINGS.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if not XTB.is_file():
        raise SystemExit(f"xTB not found: {XTB}")
    protonation = json.loads((PROTONATION / "results.json").read_text(encoding="utf-8"))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    # One water-reference Opt plus, per placement, one complex Opt and two
    # fragment single points at the optimized complex geometry.
    maximum_jobs = 1 + len(PAIRS) * 2 * 2 * 9 * 3
    status = {
        "status": "RUNNING", "started_at_utc": now(), "finished_at_utc": None,
        "maximum_jobs_if_all_waters_remain_intact": maximum_jobs,
        "complex_optimization_jobs_expected": 1 + len(PAIRS) * 2 * 2 * 9,
        "processed_jobs": 0, "computed_jobs": 0,
        "reused_jobs": 0, "current_job": None,
    }
    write_json(STATUS, status)
    try:
        water_dir = OUTPUT / "water_reference"
        water_elements = ["O", "H", "H"]
        water_coordinates = water_geometry(
            np.zeros(3), np.asarray((1.0, 0.0, 0.0)), np.asarray((0.0, 1.0, 0.0))
        )
        water_energy, water_xyz, provenance = run_job(
            water_dir, water_elements, water_coordinates, 0
        )
        status["processed_jobs"] += 1
        status["computed_jobs" if provenance == "computed_this_run" else "reused_jobs"] += 1
        pair_results = []
        for pair_id, bz_id, sb_id in PAIRS:
            pair = {"pair_id": pair_id, "bz_conformer": bz_id, "sb_conformer": sb_id, "states": {}}
            for state, charge in (("neutral", 0), ("n6_protonated", 1)):
                state_result = {}
                for molecule, conformer_id in (("1Bz-LSD_RR", bz_id), ("1SB-LSD_RR", sb_id)):
                    dataset, solute_record, model = source_record(
                        protonation, molecule, state, conformer_id
                    )
                    elements = [atom.element for atom in model.atoms]
                    coordinates = [np.asarray(atom.position, dtype=float) for atom in model.atoms]
                    records = []
                    for placement_id, water in placements(model):
                        status["current_job"] = f"{pair_id}/{state}/{molecule}/{placement_id}"
                        write_json(STATUS, status)
                        job_dir = OUTPUT / pair_id / state / molecule / placement_id
                        energy, optimized, provenance = run_job(
                            job_dir,
                            elements + water_elements,
                            coordinates + water,
                            charge,
                        )
                        status["processed_jobs"] += 1
                        status["computed_jobs" if provenance == "computed_this_run" else "reused_jobs"] += 1
                        status["current_job"] = None
                        write_json(STATUS, status)
                        geometry = analyze_complex(model, optimized)
                        raw = read_xyz_model(optimized)
                        solute_count = len(model.atoms)
                        raw_elements = [atom.element for atom in raw.atoms]
                        raw_coordinates = [np.asarray(atom.position, dtype=float) for atom in raw.atoms]
                        solute_sp = None
                        water_sp = None
                        if geometry["intact_water"]:
                            status["current_job"] = f"{pair_id}/{state}/{molecule}/{placement_id}/solute_sp"
                            write_json(STATUS, status)
                            solute_sp, solute_provenance = run_single_point(
                                job_dir / "fragment_solute_sp",
                                raw_elements[:solute_count],
                                raw_coordinates[:solute_count],
                                charge,
                            )
                            status["processed_jobs"] += 1
                            status["computed_jobs" if solute_provenance == "computed_this_run" else "reused_jobs"] += 1
                            status["current_job"] = f"{pair_id}/{state}/{molecule}/{placement_id}/water_sp"
                            write_json(STATUS, status)
                            water_sp, water_provenance = run_single_point(
                                job_dir / "fragment_water_sp",
                                raw_elements[solute_count:],
                                raw_coordinates[solute_count:],
                                0,
                            )
                            status["processed_jobs"] += 1
                            status["computed_jobs" if water_provenance == "computed_this_run" else "reused_jobs"] += 1
                            status["current_job"] = None
                            write_json(STATUS, status)
                        records.append(
                            {
                                "placement_id": placement_id,
                                "energy_hartree": energy,
                                "apparent_cluster_formation_energy_kj_mol": (
                                    energy - float(solute_record["energy_hartree"]) - water_energy
                                ) * HARTREE_TO_KJ_MOL,
                                "fragment_interaction_energy_kj_mol": (
                                    (energy - solute_sp - water_sp) * HARTREE_TO_KJ_MOL
                                    if solute_sp is not None and water_sp is not None else None
                                ),
                                "solute_deformation_energy_kj_mol": (
                                    (solute_sp - float(solute_record["energy_hartree"])) * HARTREE_TO_KJ_MOL
                                    if solute_sp is not None else None
                                ),
                                "water_deformation_energy_kj_mol": (
                                    (water_sp - water_energy) * HARTREE_TO_KJ_MOL
                                    if water_sp is not None else None
                                ),
                                "optimized_xyz": str(optimized.resolve()),
                                "provenance": provenance,
                                "geometry": geometry,
                            }
                        )
                    state_result[molecule] = {
                        "source_conformer_id": conformer_id,
                        "charge": charge,
                        "summary": aggregate(records),
                        "records": records,
                    }
                pair["states"][state] = state_result
            pair_results.append(pair)
        payload = {
            "schema_version": 1,
            "generated_at_utc": now(),
            "software": "xTB 6.7.1pre",
            "method": "one explicit H2O, GFN2-xTB tight Opt, ALPB(water)",
            "water_reference_energy_hartree": water_energy,
            "water_reference_xyz": str(water_xyz.resolve()),
            "scope_limit": "microsolvation geometry screen; not MD, hydration free energy, or reaction barrier",
            "pairs": pair_results,
        }
        write_json(RESULTS, payload)
        write_findings(payload)
        status.update(
            status="SUCCESS",
            finished_at_utc=now(),
            current_job=None,
            skipped_fragment_jobs_due_nonintact_water=(
                maximum_jobs - status["processed_jobs"]
            ),
        )
        write_json(STATUS, status)
        return 0
    except Exception as exc:
        status.update(status="FAILED", finished_at_utc=now(), error=f"{type(exc).__name__}: {exc}")
        write_json(STATUS, status)
        raise


if __name__ == "__main__":
    raise SystemExit(main())

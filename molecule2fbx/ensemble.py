"""Conformer screening, post-DFT deduplication, and ensemble reporting."""

from __future__ import annotations

import json
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .config import ConversionRequest
from .errors import ForceFieldError
from .model import MoleculeModel
from .quantum.base import QuantumResult
from .structures import (
    ConformerCandidate,
    aligned_atom_subset_rmsd,
    aligned_heavy_atom_rmsd,
    rmsd_atom_subsets,
)


HARTREE_TO_KJ_MOL = 2625.49963948
KCAL_TO_KJ_MOL = 4.184


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _optimization_provenance(result: QuantumResult) -> str:
    recorded = result.model.metadata.get("optimization_provenance")
    if isinstance(recorded, str) and recorded:
        return recorded
    return (
        "reused_legacy"
        if result.model.metadata.get("reused_quantum_result") is True
        else "computed_this_run"
    )


def _frequency_provenance(result: QuantumResult) -> Optional[str]:
    if not (result.frequency_requested and result.frequencies_cm1):
        return None
    recorded = result.model.metadata.get("frequency_provenance")
    if isinstance(recorded, str) and recorded:
        return recorded
    return "reused_or_preexisting_unspecified"


def _forcefield_used_summary(
    metadata: Dict[str, object], initial_screening: Dict[str, object]
) -> Optional[str]:
    """Recover the force field from screening records when the best Opt was reused.

    Externally reused quantum structures do not necessarily retain RDKit
    force-field metadata.  The force-field stage belongs to the current
    screening run, so its records are the authoritative fallback.
    """

    recorded = metadata.get("forcefield")
    if isinstance(recorded, str) and recorded:
        return recorded
    records = initial_screening.get("records", [])
    if not isinstance(records, list):
        return None
    values = sorted(
        {
            str(item["forcefield"])
            for item in records
            if isinstance(item, dict)
            and isinstance(item.get("forcefield"), str)
            and item["forcefield"]
        }
    )
    return values[0] if len(values) == 1 else "/".join(values) or None


def _cluster_by_subset(
    results: Sequence[QuantumResult],
    atom_indices: Sequence[int],
    *,
    rmsd_threshold: float,
    prefix: str,
) -> Tuple[Dict[int, int], Dict[int, float], List[QuantumResult]]:
    representatives: List[QuantumResult] = []
    assignment: Dict[int, int] = {}
    rmsd_to_representative: Dict[int, float] = {}
    for result in results:
        cluster_index: Optional[int] = None
        selected_rmsd = 0.0
        for index, representative in enumerate(representatives):
            value = aligned_atom_subset_rmsd(
                result.model, representative.model, atom_indices
            )
            if value < rmsd_threshold:
                cluster_index = index
                selected_rmsd = value
                break
        if cluster_index is None:
            cluster_index = len(representatives)
            representatives.append(result)
        assignment[result.conformer_index] = cluster_index
        rmsd_to_representative[result.conformer_index] = selected_rmsd
    return assignment, rmsd_to_representative, representatives


def screen_forcefield_candidates(
    candidates: Sequence[ConformerCandidate],
    *,
    energy_window_kj: float,
    rmsd_threshold: float,
    maximum: int,
) -> Tuple[List[ConformerCandidate], Dict[str, object]]:
    """Energy-filter and greedily cluster converged force-field structures."""

    usable = [
        candidate
        for candidate in candidates
        if candidate.converged and candidate.energy_kcal_mol is not None
    ]
    if not usable:
        raise ForceFieldError("No converged force-field conformer has a usable energy")
    minimum = min(float(candidate.energy_kcal_mol) for candidate in usable)
    ranked = sorted(
        usable,
        key=lambda candidate: (float(candidate.energy_kcal_mol), candidate.conformer_index),
    )
    eligible = [
        candidate
        for candidate in ranked
        if (float(candidate.energy_kcal_mol) - minimum) * KCAL_TO_KJ_MOL
        <= energy_window_kj + 1.0e-9
    ]

    representatives: List[ConformerCandidate] = []
    assignment: Dict[int, int] = {}
    for candidate in eligible:
        cluster_index: Optional[int] = None
        for index, representative in enumerate(representatives):
            if aligned_heavy_atom_rmsd(candidate.model, representative.model) < rmsd_threshold:
                cluster_index = index
                break
        if cluster_index is None:
            cluster_index = len(representatives)
            representatives.append(candidate)
        assignment[candidate.conformer_index] = cluster_index

    chosen = representatives[:maximum]
    selected_pool_indices = {candidate.conformer_index for candidate in chosen}
    reindexed: List[ConformerCandidate] = []
    for job_index, candidate in enumerate(chosen):
        pool_index = candidate.conformer_index
        model = candidate.model.with_metadata(
            conformer_index=job_index,
            conformer_pool_index=pool_index,
            conformer_pool_index_provenance="recorded_from_current_etkdg_pool",
            conformer_pool_size=len(candidates),
            initial_rmsd_cluster_id=f"ff_cluster_{assignment[pool_index] + 1:03d}",
            conformer_selection_rmsd_threshold_angstrom=rmsd_threshold,
            conformer_selection_threshold_relaxed=False,
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

    records: List[Dict[str, object]] = []
    for candidate in sorted(candidates, key=lambda item: item.conformer_index):
        energy = candidate.energy_kcal_mol
        delta = None if energy is None else (float(energy) - minimum) * KCAL_TO_KJ_MOL
        cluster_index = assignment.get(candidate.conformer_index)
        excluded_reason = None
        if not candidate.converged:
            excluded_reason = candidate.failure_status or "forcefield_optimization_failure"
        elif energy is None:
            excluded_reason = "forcefield_energy_unavailable"
        elif delta is not None and delta > energy_window_kj + 1.0e-9:
            excluded_reason = "forcefield_energy_window"
        elif candidate.conformer_index not in selected_pool_indices:
            representative = (
                representatives[cluster_index] if cluster_index is not None else None
            )
            excluded_reason = (
                "rmsd_duplicate"
                if representative is not None
                and representative.conformer_index != candidate.conformer_index
                else "dft_candidate_limit"
            )
        records.append(
            {
                "pool_index": candidate.conformer_index,
                "forcefield": candidate.forcefield,
                "optimization_converged": candidate.converged,
                "optimization_status": candidate.model.metadata.get(
                    "forcefield_optimization_status"
                ),
                "optimization_iterations": candidate.optimization_iterations,
                "energy_kcal_mol": energy,
                "relative_energy_kj_mol": delta,
                "rmsd_cluster_id": (
                    f"ff_cluster_{cluster_index + 1:03d}"
                    if cluster_index is not None
                    else None
                ),
                "selected_for_dft": candidate.conformer_index in selected_pool_indices,
                "excluded_reason": excluded_reason,
            }
        )
    report: Dict[str, object] = {
        "generated": len(candidates),
        "converged": len(usable),
        "energy_filtered": len(eligible),
        "energy_reference_kcal_mol": minimum,
        "energy_window_kj_mol": energy_window_kj,
        "rmsd_threshold_angstrom": rmsd_threshold,
        "cluster_count": len(representatives),
        "selected_for_dft": len(reindexed),
        "threshold_relaxed": False,
        "records": records,
    }
    return reindexed, report


def cluster_optimized_results(
    results: Sequence[QuantumResult],
    *,
    rmsd_threshold: float,
    energy_window_kj: Optional[float] = None,
) -> Tuple[List[QuantumResult], Dict[str, object]]:
    """Remove optimized structures that converged to the same heavy-atom geometry."""

    if not results:
        return [], {"records": [], "cluster_count": 0}
    minimum = min(result.final_energy_hartree for result in results)
    ranked = sorted(results, key=lambda result: result.final_energy_hartree)
    eligible = [
        result
        for result in ranked
        if energy_window_kj is None
        or (result.final_energy_hartree - minimum) * HARTREE_TO_KJ_MOL
        <= energy_window_kj + 1.0e-9
    ]
    representatives: List[QuantumResult] = []
    assignment: Dict[int, int] = {}
    for result in eligible:
        cluster_index: Optional[int] = None
        for index, representative in enumerate(representatives):
            if aligned_heavy_atom_rmsd(result.model, representative.model) < rmsd_threshold:
                cluster_index = index
                break
        if cluster_index is None:
            cluster_index = len(representatives)
            representatives.append(result)
        assignment[result.conformer_index] = cluster_index

    atom_subsets = rmsd_atom_subsets(eligible[0].model)
    secondary = {}
    for key, indices, prefix in (
        ("common_scaffold", atom_subsets.common_scaffold, "dft_core_cluster"),
        ("reaction_center", atom_subsets.reaction_center, "dft_site_cluster"),
    ):
        if len(indices) < 3:
            continue
        subset_assignment, subset_rmsd, subset_representatives = _cluster_by_subset(
            eligible,
            indices,
            rmsd_threshold=rmsd_threshold,
            prefix=prefix,
        )
        secondary[key] = {
            "assignment": subset_assignment,
            "rmsd_to_representative": subset_rmsd,
            "representatives": subset_representatives,
            "prefix": prefix,
            "atom_indices": tuple(indices),
        }

    updated_representatives: List[QuantumResult] = []
    representative_indices = {result.conformer_index for result in representatives}
    for result in representatives:
        cluster_index = assignment[result.conformer_index]
        model = result.model.with_metadata(
            dft_rmsd_cluster_id=f"dft_cluster_{cluster_index + 1:03d}",
            dft_all_heavy_rmsd_cluster_id=f"dft_cluster_{cluster_index + 1:03d}",
            dft_common_scaffold_rmsd_cluster_id=(
                f"dft_core_cluster_{secondary['common_scaffold']['assignment'][result.conformer_index] + 1:03d}"
                if "common_scaffold" in secondary
                else None
            ),
            dft_reaction_center_rmsd_cluster_id=(
                f"dft_site_cluster_{secondary['reaction_center']['assignment'][result.conformer_index] + 1:03d}"
                if "reaction_center" in secondary
                else None
            ),
            dft_cluster_representative=True,
            relative_electronic_energy_kj_mol=(
                result.final_energy_hartree - minimum
            )
            * HARTREE_TO_KJ_MOL,
        )
        updated_representatives.append(replace(result, model=model))

    records: List[Dict[str, object]] = []
    for result in sorted(results, key=lambda item: item.conformer_index):
        delta = (result.final_energy_hartree - minimum) * HARTREE_TO_KJ_MOL
        cluster_index = assignment.get(result.conformer_index)
        excluded_reason = None
        duplicate_of = None
        if energy_window_kj is not None and delta > energy_window_kj + 1.0e-9:
            excluded_reason = "dft_energy_window"
        elif result.conformer_index not in representative_indices:
            excluded_reason = "dft_rmsd_duplicate"
            if cluster_index is not None:
                duplicate_of = representatives[cluster_index].conformer_index
        records.append(
            {
                "conformer_index": result.conformer_index,
                "electronic_energy_hartree": result.final_energy_hartree,
                "relative_electronic_energy_kj_mol": delta,
                "rmsd_cluster_id": (
                    f"dft_cluster_{cluster_index + 1:03d}"
                    if cluster_index is not None
                    else None
                ),
                "all_heavy_rmsd_cluster_id": (
                    f"dft_cluster_{cluster_index + 1:03d}"
                    if cluster_index is not None
                    else None
                ),
                "common_scaffold_rmsd_cluster_id": (
                    f"dft_core_cluster_{secondary['common_scaffold']['assignment'][result.conformer_index] + 1:03d}"
                    if "common_scaffold" in secondary
                    and result.conformer_index
                    in secondary["common_scaffold"]["assignment"]
                    else None
                ),
                "common_scaffold_rmsd_to_representative_angstrom": (
                    secondary["common_scaffold"]["rmsd_to_representative"].get(
                        result.conformer_index
                    )
                    if "common_scaffold" in secondary
                    else None
                ),
                "reaction_center_rmsd_cluster_id": (
                    f"dft_site_cluster_{secondary['reaction_center']['assignment'][result.conformer_index] + 1:03d}"
                    if "reaction_center" in secondary
                    and result.conformer_index
                    in secondary["reaction_center"]["assignment"]
                    else None
                ),
                "reaction_center_rmsd_to_representative_angstrom": (
                    secondary["reaction_center"]["rmsd_to_representative"].get(
                        result.conformer_index
                    )
                    if "reaction_center" in secondary
                    else None
                ),
                "selected_as_representative": result.conformer_index
                in representative_indices,
                "duplicate_of_conformer_index": duplicate_of,
                "excluded_reason": excluded_reason,
                "source": (
                    "reused"
                    if result.model.metadata.get("reused_quantum_result") is True
                    else "generated"
                ),
                "calculation_directory": (
                    str(result.calculation_directory)
                    if result.calculation_directory is not None
                    else None
                ),
            }
        )
    report: Dict[str, object] = {
        "energy_reference_hartree": minimum,
        "energy_window_kj_mol": energy_window_kj,
        "rmsd_threshold_angstrom": rmsd_threshold,
        "cluster_count": len(representatives),
        "rmsd_analyses": {
            "all_heavy_fixed_order": {
                "cluster_count": len(representatives),
                "threshold_angstrom": rmsd_threshold,
                "atom_indices": list(atom_subsets.all_heavy),
                "symmetry_permutations": False,
            },
            **{
                key: {
                    "cluster_count": len(value["representatives"]),
                    "threshold_angstrom": rmsd_threshold,
                    "atom_indices": list(value["atom_indices"]),
                    "symmetry_permutations": False,
                }
                for key, value in secondary.items()
            },
            "atom_subset_definition": atom_subsets.to_metadata(),
        },
        "records": records,
    }
    return updated_representatives, report


def _local_minimum_status(result: QuantumResult) -> str:
    if not result.frequency_requested or not result.frequencies_cm1:
        return "not_evaluated"
    if result.imaginary_frequencies_cm1:
        return "not_a_confirmed_local_minimum"
    return "local_minimum_candidate"


def build_ensemble_report(
    request: ConversionRequest,
    results: Sequence[QuantumResult],
    *,
    initial_screening: Dict[str, object],
    post_dft_screening: Dict[str, object],
    frequency_selection: Optional[Dict[str, object]],
    dft_failures: Sequence[Dict[str, object]] = (),
) -> Dict[str, object]:
    """Create the UI-neutral JSON payload for one complete ensemble run."""

    electronic_minimum = min(result.final_energy_hartree for result in results)
    gibbs_values = [
        result.thermochemistry.gibbs_free_energy_hartree
        for result in results
        if result.thermochemistry is not None
        and result.thermochemistry.gibbs_free_energy_hartree is not None
    ]
    gibbs_minimum = min(gibbs_values) if gibbs_values else None
    entries = []
    for result in sorted(results, key=lambda item: item.final_energy_hartree):
        thermo = result.thermochemistry
        gibbs = thermo.gibbs_free_energy_hartree if thermo is not None else None
        frequency_calculated = bool(
            result.frequency_requested and result.frequencies_cm1
        )
        entries.append(
            {
                "conformer_id": f"conf{result.conformer_index + 1:03d}",
                "conformer_index": result.conformer_index,
                "stereochemistry": result.model.metadata.get("stereochemistry"),
                "dft_energy_hartree": result.final_energy_hartree,
                "relative_energy_kj_mol": (
                    result.final_energy_hartree - electronic_minimum
                )
                * HARTREE_TO_KJ_MOL,
                "gibbs_energy_hartree": gibbs,
                "relative_gibbs_energy_kj_mol": (
                    (gibbs - gibbs_minimum) * HARTREE_TO_KJ_MOL
                    if gibbs is not None and gibbs_minimum is not None
                    else None
                ),
                "frequency_calculated": frequency_calculated,
                "frequency_provenance": _frequency_provenance(result),
                "number_of_modes": (
                    len(result.frequencies_cm1) if frequency_calculated else None
                ),
                "lowest_frequencies_cm1": (
                    sorted(result.frequencies_cm1)[:10] if frequency_calculated else []
                ),
                "imaginary_modes": (
                    len(result.imaginary_frequencies_cm1)
                    if frequency_calculated
                    else None
                ),
                "imaginary_frequencies_cm1": list(result.imaginary_frequencies_cm1),
                "translation_rotation_near_zero_modes_cm1": [
                    value for value in result.frequencies_cm1 if abs(value) < 1.0
                ],
                "low_frequency_modes_cm1": [
                    value
                    for value in result.frequencies_cm1
                    if abs(value) >= 1.0
                    and value < request.low_frequency_threshold_cm1
                ],
                "local_minimum_assessment": _local_minimum_status(result),
                "rmsd_cluster_id": result.model.metadata.get("dft_rmsd_cluster_id"),
                "all_heavy_rmsd_cluster_id": result.model.metadata.get(
                    "dft_all_heavy_rmsd_cluster_id",
                    result.model.metadata.get("dft_rmsd_cluster_id"),
                ),
                "common_scaffold_rmsd_cluster_id": result.model.metadata.get(
                    "dft_common_scaffold_rmsd_cluster_id"
                ),
                "reaction_center_rmsd_cluster_id": result.model.metadata.get(
                    "dft_reaction_center_rmsd_cluster_id"
                ),
                "source": (
                    "reused"
                    if result.model.metadata.get("reused_quantum_result") is True
                    else "generated"
                ),
                "optimization_provenance": _optimization_provenance(result),
                "conformer_pool_index": result.model.metadata.get(
                    "conformer_pool_index"
                ),
                "conformer_pool_index_provenance": result.model.metadata.get(
                    "conformer_pool_index_provenance",
                    (
                        "recorded_from_current_etkdg_pool"
                        if result.model.metadata.get("conformer_pool_index") is not None
                        else "unavailable"
                    ),
                ),
                "calculation_directory": (
                    str(result.calculation_directory)
                    if result.calculation_directory is not None
                    else None
                ),
                "thermochemistry": thermo.to_metadata() if thermo is not None else None,
                "post_dft_stereochemistry_validation": result.model.metadata.get(
                    "post_dft_stereochemistry_validation"
                ),
            }
        )
    best = min(results, key=lambda result: result.final_energy_hartree)
    metadata = best.model.metadata
    etkdg = metadata.get("etkdg", {})
    if not isinstance(etkdg, dict):
        etkdg = {}
    frequency_selection = frequency_selection or {}
    frequency_failures = frequency_selection.get("failures", [])
    if not isinstance(frequency_failures, list):
        frequency_failures = []
    current_dft_failures = [
        failure
        for failure in dft_failures
        if failure.get("phase") in {"dft_optimization", "post_dft_validation"}
    ]
    recovery_history = [
        failure for failure in dft_failures if failure.get("phase") == "resume_recovery"
    ]
    status = "SUCCESS"
    if current_dft_failures or frequency_failures:
        status = "PARTIAL"
    selected_frequency_indices = set(
        frequency_selection.get("selected_conformer_indices", [])
    )
    completed_frequency_indices = {
        entry["conformer_index"]
        for entry in entries
        if entry["frequency_calculated"]
    }
    completed_selected_indices = (
        selected_frequency_indices & completed_frequency_indices
    )
    summary = {
        "etkdg_candidates": initial_screening.get("generated", 0),
        "forcefield_valid": initial_screening.get("converged", 0),
        "energy_filtered": initial_screening.get("energy_filtered", 0),
        "rmsd_clusters": initial_screening.get("cluster_count", 0),
        "dft_candidates": len(post_dft_screening.get("records", []))
        + len(current_dft_failures),
        "dft_completed": len(post_dft_screening.get("records", [])),
        "unique_dft_structures": len(results),
        "frequency_candidates": len(selected_frequency_indices),
        "frequency_completed": len(completed_selected_indices),
        "frequency_selected_this_run": len(selected_frequency_indices),
        "frequency_completed_for_this_run_selection": len(
            completed_selected_indices
        ),
        "frequency_available_total": len(completed_frequency_indices),
        "frequency_preexisting_or_reused_total": sum(
            1
            for entry in entries
            if entry["frequency_calculated"]
            and entry["frequency_provenance"]
            != "computed_this_run"
        ),
        "frequency_computed_this_run": sum(
            1
            for entry in entries
            if entry["frequency_calculated"]
            and entry["frequency_provenance"] == "computed_this_run"
        ),
    }
    report = {
        "schema_version": 2,
        "run_type": "conformer_ensemble_dft_screening",
        "calculation_status": status,
        "summary": summary,
        "generated_at_utc": _utc_now(),
        "interpretation": {
            "structure_provenance": "computed_quantum",
            "best_structure_claim": (
                "Lowest electronic energy among structures evaluated in this run; "
                "not a proven global minimum and not an experimental structure."
            ),
            "frequency_scope": (
                "Frequency confirms stationary-point character only for structures on which "
                "Freq was completed; it does not establish a global minimum."
            ),
            "gibbs_ensemble_scope": (
                "Relative Gibbs energies and any Gibbs-weighted statistics are conditional "
                "on conformers with completed frequency calculations unless every final "
                "conformer has thermochemistry."
            ),
        },
        "molecule": {
            "name": best.model.name,
            "cid": best.model.cid,
            "smiles": metadata.get("canonical_isomeric_smiles"),
            "original_smiles": metadata.get("original_smiles"),
            "charge": request.effective_charge,
            "multiplicity": request.effective_multiplicity,
            "stereochemistry": metadata.get("stereochemistry"),
            "stereochemistry_specified": True,
            "configuration": request.stereochemistry_label,
        },
        "conformer_search": {
            "method": "ETKDG",
            "pool_size": request.effective_conformer_pool,
            "random_seed": request.random_seed,
            "embedding_prune_rmsd_angstrom": request.embedding_prune_rmsd,
            "embedding_failures": etkdg.get("embedding_failures"),
            "rdkit_version": etkdg.get("rdkit_version"),
            "forcefield_preference": ["MMFF94s", "UFF"],
            "forcefield_used": _forcefield_used_summary(
                metadata, initial_screening
            ),
            "initial_screening": initial_screening,
        },
        "dft": {
            "software": "ORCA",
            "version": best.backend_version,
            "method": request.method,
            "functional": request.effective_functional if request.method == "dft" else None,
            "basis": request.effective_basis,
            "solvent_model": "gas_phase",
            "dispersion_correction": "none",
            "charge": request.effective_charge,
            "multiplicity": request.effective_multiplicity,
            "nprocs": request.nprocs,
            "maxcore_mb_per_process": request.maxcore_mb,
            "post_optimization_screening": post_dft_screening,
            "failures": current_dft_failures,
            "recovery_history": recovery_history,
        },
        "frequency": {
            "window_kj_mol": request.frequency_window_kj,
            "max_structures": request.frequency_max,
            "imaginary_threshold_cm1": request.imaginary_threshold_cm1,
            "low_frequency_threshold_cm1": request.low_frequency_threshold_cm1,
            "translation_rotation_near_zero_threshold_abs_cm1": 1.0,
            "low_frequency_modes_exclude_near_zero_translation_rotation": True,
            "selection": frequency_selection,
        },
        "best_conformer_id": f"conf{best.conformer_index + 1:03d}",
        "final_conformer_ensemble": entries,
    }
    report.update(
        {
            "name": best.model.name,
            "charge": request.effective_charge,
            "multiplicity": request.effective_multiplicity,
            "stereochemistry": {
                "specified": True,
                "configuration": request.stereochemistry_label,
            },
            "conformer_pool": request.effective_conformer_pool,
            "forcefield": "MMFF94s/UFF",
            "energy_window_kj": request.effective_forcefield_energy_window_kj,
            "rmsd_threshold_angstrom": request.conformer_rmsd_threshold,
            "dft_max_conformers": request.conformers,
            "frequency_window_kj": request.frequency_window_kj,
            "frequency_max": request.frequency_max,
            "conformers": entries,
        }
    )
    return report


def write_ensemble_report(
    payload: Dict[str, object],
    output_dir: Path,
    stem: str,
    *,
    filename: Optional[str] = None,
) -> Path:
    path = output_dir / (filename or f"{stem}_ensemble.json")
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)
    return path


def _write_run_summary_legacy(payload: Dict[str, object], output_dir: Path) -> Path:
    """Write a compact morning-review report next to ensemble.json."""

    summary = payload["summary"]
    conformers = payload["final_conformer_ensemble"]
    best_id = payload["best_conformer_id"]
    best = next(item for item in conformers if item["conformer_id"] == best_id)
    low_frequency_count = sum(
        len(item["low_frequency_modes_cm1"])
        for item in conformers
        if item["frequency_calculated"]
    )
    imaginary = best["imaginary_modes"] if best["frequency_calculated"] else "not evaluated"
    lines = [
        f"# {payload['name']} Ensemble Night Run",
        "",
        "## Calculation status",
        "",
        str(payload["calculation_status"]),
        "",
        "## Summary",
        "",
        "```text",
        f"ETKDG candidates:       {summary['etkdg_candidates']}",
        f"MMFF/UFF valid:         {summary['forcefield_valid']}",
        f"Energy-filtered:        {summary['energy_filtered']}",
        f"RMSD clusters:          {summary['rmsd_clusters']}",
        f"DFT candidates:         {summary['dft_candidates']}",
        f"DFT completed:          {summary['dft_completed']}",
        f"Unique DFT structures:  {summary['unique_dft_structures']}",
        f"Freq candidates:        {summary['frequency_candidates']}",
        f"Freq completed:         {summary['frequency_completed']}",
        "```",
        "",
        "## Lowest-energy structure found",
        "",
        "```text",
        f"Conformer:       {best['conformer_id']}",
        f"Energy:          {best['dft_energy_hartree']:.12f} Eh",
        f"Relative ΔE:     {best['relative_energy_kj_mol']:.3f} kJ/mol",
        f"Gibbs:           {best['gibbs_energy_hartree']}",
        f"Relative ΔG:     {best['relative_gibbs_energy_kj_mol']}",
        f"Imaginary modes: {imaginary}",
        "```",
        "",
        "## Morning checklist",
        "",
        f"- [{'x' if payload['calculation_status'] == 'SUCCESS' else ' '}] ジョブ全体がSUCCESSか",
        f"- [{'x' if payload['molecule'].get('configuration') else ' '}] stereochemistryがRRのままか",
        f"- [x] DFT Opt成功数: {summary['dft_completed']}",
        f"- [x] 最低エネルギー構造: {best['conformer_id']}",
        f"- [{'x' if imaginary == 0 else ' '}] 最低エネルギー構造の虚振動が0か",
        "",
        "## Interpretation",
        "",
        "This is the lowest-energy structure found among the conformers evaluated in this run. ",
        "It is not a proven global minimum and is not an experimentally determined structure. ",
        "Zero imaginary modes only supports local-minimum character for that calculated structure.",
    ]
    if low_frequency_count:
        lines.extend(
            [
                "",
                "## Low-frequency caution",
                "",
                f"Freq-completed structures contain {low_frequency_count} modes below the configured ",
                "50 cm^-1 reporting threshold. Their Gibbs energies can be sensitive to the ",
                "harmonic/quasi-RRHO treatment and should be interpreted cautiously.",
            ]
        )
    path = output_dir / "RUN_SUMMARY.md"
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path


def write_run_summary(payload: Dict[str, object], output_dir: Path) -> Path:
    """Write an unambiguous morning-review report next to ensemble.json."""

    summary = payload["summary"]
    conformers = payload["final_conformer_ensemble"]
    best_id = payload["best_conformer_id"]
    best = next(item for item in conformers if item["conformer_id"] == best_id)
    low_frequency_count = sum(
        len(item["low_frequency_modes_cm1"])
        for item in conformers
        if item["frequency_calculated"]
    )
    imaginary = (
        best["imaginary_modes"]
        if best["frequency_calculated"]
        else "not evaluated"
    )
    lines = [
        f"# {payload['name']} Ensemble Night Run",
        "",
        "## Calculation status",
        "",
        str(payload["calculation_status"]),
        "",
        "## Summary",
        "",
        "```text",
        f"ETKDG candidates:       {summary['etkdg_candidates']}",
        f"MMFF/UFF valid:         {summary['forcefield_valid']}",
        f"Energy-filtered:        {summary['energy_filtered']}",
        f"All-heavy RMSD clusters:{summary['rmsd_clusters']:>9}",
        f"DFT candidates:         {summary['dft_candidates']}",
        f"DFT completed:          {summary['dft_completed']}",
        f"Unique DFT structures:  {summary['unique_dft_structures']}",
        f"Freq selected this run: {summary['frequency_selected_this_run']}",
        "Freq completed for selection: "
        f"{summary['frequency_completed_for_this_run_selection']}",
        f"Freq available total:   {summary['frequency_available_total']}",
        f"Freq new supplemental:  {summary['frequency_computed_this_run']}",
        "Freq preexisting local:  "
        f"{summary.get('frequency_preexisting_local_original_ensemble', 0)}",
        "Freq reused external:    "
        f"{summary.get('frequency_reused_external_read_only', 0)}",
        "Freq preexisting total:  "
        f"{summary['frequency_preexisting_or_reused_total']}",
        "```",
        "",
        "## Lowest-energy structure found",
        "",
        "```text",
        f"Conformer:       {best['conformer_id']}",
        f"Energy:          {best['dft_energy_hartree']:.12f} Eh",
        f"Relative Delta E:{best['relative_energy_kj_mol']:>9.3f} kJ/mol",
        f"Gibbs:           {best['gibbs_energy_hartree']}",
        f"Relative Delta G:{best['relative_gibbs_energy_kj_mol']}",
        f"Imaginary modes: {imaginary}",
        "```",
        "",
        "## Morning checklist",
        "",
        f"- [{'x' if payload['calculation_status'] == 'SUCCESS' else ' '}] Job status is SUCCESS",
        f"- [{'x' if payload['molecule'].get('configuration') else ' '}] Explicit stereochemistry was checked",
        f"- [x] DFT Opt completed: {summary['dft_completed']}",
        f"- [x] Lowest electronic-energy structure found: {best['conformer_id']}",
        f"- [{'x' if imaginary == 0 else ' '}] Lowest structure has zero significant imaginary modes",
        "",
        "## Interpretation",
        "",
        "This is the lowest-energy structure found among the conformers evaluated in this run. ",
        "It is not a proven global minimum and is not an experimentally determined structure. ",
        "Zero imaginary modes only supports local-minimum character for that calculated structure.",
    ]
    if low_frequency_count:
        lines.extend(
            [
                "",
                "## Low-frequency caution",
                "",
                f"Freq-completed structures contain {low_frequency_count} modes below the configured ",
                "50 cm^-1 reporting threshold. Their Gibbs energies can be sensitive to the ",
                "harmonic/quasi-RRHO treatment and should be interpreted cautiously.",
            ]
        )
    path = output_dir / "RUN_SUMMARY.md"
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path

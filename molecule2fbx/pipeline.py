"""UI-neutral conversion orchestration."""

from __future__ import annotations

import json
import tempfile
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Type

from .blender_export import export_model, find_blender
from .config import ConversionRequest, XYZConversionRequest
from .electronic import validate_electronic_state, validate_metal_policy
from .ensemble import (
    build_ensemble_report,
    cluster_optimized_results,
    screen_forcefield_candidates,
    write_ensemble_report,
    write_run_summary,
)
from .errors import (
    APIError,
    ConfigurationError,
    ForceFieldError,
    No3DConformerError,
    QuantumCalculationError,
)
from .model import MoleculeModel, safe_filename
from .pubchem import download_and_parse, fetch_compound_properties
from .quantum import OrcaBackend, QuantumResult, QuantumSettings, find_orca
from .quantum.reuse import (
    archive_incomplete_calculation,
    discover_reusable_orca_results,
    ensure_quantum_result_frequency,
    load_reusable_orca_result,
)
from .structures import (
    ConformerCandidate,
    aligned_heavy_atom_rmsd,
    generate_forcefield_conformers,
    inspect_smiles,
    select_diverse_conformers,
    validate_model_stereochemistry,
)
from .xyz import load_xyz_for_export


Logger = Callable[[str], None]
HARTREE_TO_KJ_MOL = 2625.49963948


@dataclass(frozen=True)
class ExportedArtifact:
    fbx_path: Path
    metadata_path: Path
    structure_origin: str
    conformer_index: Optional[int]


@dataclass(frozen=True)
class ConversionOutcome:
    selected_model: MoleculeModel
    artifacts: Tuple[ExportedArtifact, ...]
    method: str
    ensemble_report_path: Optional[Path] = None
    run_summary_path: Optional[Path] = None


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _rename_model(model: MoleculeModel, name: str) -> MoleculeModel:
    return MoleculeModel(model.cid, name, model.atoms, model.bonds, model.metadata)


def _write_metadata(
    model: MoleculeModel,
    fbx_path: Path,
    metadata_path: Optional[Path] = None,
) -> Path:
    metadata_path = metadata_path or fbx_path.with_suffix(".metadata.json")
    payload = model.to_dict()
    payload.update(
        {
            "schema_version": 1,
            "metadata_generated_at_utc": _utc_now(),
            "fbx_file": fbx_path.name,
        }
    )
    metadata_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return metadata_path


def _export_artifact(
    model: MoleculeModel,
    *,
    output_dir: Path,
    stem: str,
    blender_executable: str,
    blender_timeout: float,
    exporter: Callable[..., Path],
    metadata_path: Optional[Path] = None,
) -> ExportedArtifact:
    fbx_path = exporter(
        model,
        output_path=output_dir / f"{stem}.fbx",
        blender_executable=blender_executable,
        timeout=blender_timeout,
    )
    metadata_path = _write_metadata(model, fbx_path, metadata_path)
    return ExportedArtifact(
        fbx_path=fbx_path,
        metadata_path=metadata_path,
        structure_origin=str(model.metadata.get("structure_origin", "unknown")),
        conformer_index=(
            int(model.metadata["conformer_index"])
            if "conformer_index" in model.metadata
            else None
        ),
    )


def run_xyz_conversion(
    request: XYZConversionRequest,
    *,
    log: Optional[Logger] = None,
    blender_finder: Callable[[Optional[str]], str] = find_blender,
    exporter: Callable[..., Path] = export_model,
) -> ConversionOutcome:
    """Export an existing XYZ geometry without force-field or quantum work."""

    request.validate()
    emit = log or (lambda _message: None)
    xyz_path = Path(request.xyz_path).expanduser().resolve()
    if not xyz_path.is_file():
        raise ConfigurationError(f"XYZ file not found: {xyz_path}")
    model = load_xyz_for_export(
        xyz_path,
        name=request.name,
        charge=request.charge,
    )
    output_dir = (
        Path(request.output_dir).expanduser().resolve()
        if request.output_dir is not None
        else xyz_path.parent
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    blender = blender_finder(request.blender_executable)
    base = safe_filename(model.name, xyz_path.stem)
    emit("Structure origin: imported_xyz")
    inference = model.metadata.get("xyz_bond_inference")
    if isinstance(inference, dict) and inference.get("warning"):
        emit(f"Warning: {inference['warning']}")
    artifact = _export_artifact(
        model,
        output_dir=output_dir,
        stem=base,
        blender_executable=blender,
        blender_timeout=request.blender_timeout,
        exporter=exporter,
        metadata_path=output_dir / f"{base}.fbx.metadata.json",
    )
    return ConversionOutcome(model, (artifact,), "xyz")


def _properties_for_cid(cid: int, timeout: float) -> Tuple[str, str]:
    properties = fetch_compound_properties(cid, timeout=timeout)
    smiles = properties.get("smiles")
    if not isinstance(smiles, str) or not smiles:
        raise APIError("PubChem did not return a stereochemistry-preserving SMILES")
    return str(properties.get("title") or f"CID_{cid}"), smiles


def _forcefield_candidates(
    request: ConversionRequest,
    *,
    count: Optional[int] = None,
    index_offset: int = 0,
) -> Tuple[List[ConformerCandidate], str]:
    if request.smiles is not None:
        smiles = request.smiles
        name = request.name
    else:
        name, smiles = _properties_for_cid(request.cid, request.api_timeout)  # type: ignore[arg-type]
        if request.name:
            name = request.name
    candidates, _ = generate_forcefield_conformers(
        smiles,
        count or request.conformers,
        cid=request.cid,
        name=name,
        index_offset=index_offset,
        random_seed=request.random_seed,
        prune_rms_threshold=request.embedding_prune_rmsd,
    )
    return candidates, smiles


def _select_forcefield(candidates: Sequence[ConformerCandidate]) -> ConformerCandidate:
    converged = [candidate for candidate in candidates if candidate.converged]
    if not converged:
        raise ForceFieldError("No force-field conformer converged")
    return min(
        converged,
        key=lambda candidate: (
            float("inf") if candidate.energy_kcal_mol is None else candidate.energy_kcal_mol
        ),
    )


def _reusable_conformer_indices(request: ConversionRequest) -> Tuple[int, ...]:
    if request.reuse_calculations is None:
        return ()
    root = Path(request.reuse_calculations).expanduser().resolve()
    indices = []
    for index in range(request.conformers):
        if not (root / f"conformer_{index + 1:03d}").exists():
            break
        indices.append(index)
    return tuple(indices)


def _select_quantum_candidates(
    request: ConversionRequest,
    candidates: Sequence[ConformerCandidate],
    log: Logger,
    *,
    mandatory_indices: Sequence[int] = (),
) -> List[ConformerCandidate]:
    converged = [candidate for candidate in candidates if candidate.converged]
    if request.effective_conformer_pool == request.conformers:
        return converged
    required = tuple(
        sorted(set((*mandatory_indices, *_reusable_conformer_indices(request))))
    )
    selected = select_diverse_conformers(
        candidates,
        request.conformers,
        request.conformer_rmsd_threshold,
        required_indices=required,
    )
    pool_indices = [
        int(candidate.model.metadata["conformer_pool_index"])
        for candidate in selected
    ]
    log(
        f"Generated {len(candidates)} force-field conformers and selected "
        f"{len(selected)} for quantum optimization; pool indices={pool_indices}"
    )
    return selected


def _quantum_initial_candidates(
    request: ConversionRequest,
    log: Logger,
) -> tuple:
    if request.smiles is not None:
        inspection = inspect_smiles(request.smiles)
        if inspection.has_unspecified_stereo:
            detail = (
                f"unresolved atom centers={list(inspection.unspecified_stereocenters)}, "
                f"double bonds={list(inspection.unspecified_double_bonds)}"
            )
            if request.strict_stereochemistry:
                raise ConfigurationError(
                    "SMILES contains unresolved stereochemistry; conformer and stereoisomer "
                    f"searches must remain separate ({detail})"
                )
            log(
                "Warning: SMILES contains unresolved stereochemistry; generated coordinates "
                f"do not establish a unique stereoisomer ({detail})"
            )
        if (
            request.expected_stereocenters is not None
            and len(inspection.specified_stereocenters)
            != request.expected_stereocenters
        ):
            raise ConfigurationError(
                "SMILES stereocenter count does not match the explicit expectation: "
                f"observed {len(inspection.specified_stereocenters)}, expected "
                f"{request.expected_stereocenters}"
            )
        candidates, smiles = _forcefield_candidates(
            request, count=request.effective_conformer_pool
        )
        if request.ensemble:
            selected, screening = screen_forcefield_candidates(
                candidates,
                energy_window_kj=request.effective_forcefield_energy_window_kj,
                rmsd_threshold=request.conformer_rmsd_threshold,
                maximum=request.conformers,
            )
            log(
                f"Ensemble pre-screening retained {len(selected)} DFT candidates from "
                f"{len(candidates)} embedded conformers"
            )
            return selected, smiles, screening
        return _select_quantum_candidates(request, candidates, log), smiles

    assert request.cid is not None
    pubchem_model = None
    try:
        pubchem_model = download_and_parse(request.cid, timeout=request.api_timeout)
        if request.name:
            pubchem_model = _rename_model(pubchem_model, request.name)
        log("Using PubChem 3D as quantum-optimization conformer 1")
    except No3DConformerError:
        log("PubChem has no 3D conformer; generating quantum starting geometries from SMILES")

    title = pubchem_model.name if pubchem_model else f"CID_{request.cid}"
    smiles: Optional[str] = None
    try:
        property_title, smiles = _properties_for_cid(request.cid, request.api_timeout)
        if not pubchem_model:
            title = request.name or property_title
    except APIError:
        if pubchem_model is None or request.effective_conformer_pool > 1:
            raise

    candidates: List[ConformerCandidate] = []
    if pubchem_model is not None:
        if smiles:
            pubchem_model = pubchem_model.with_metadata(**inspect_smiles(smiles).to_metadata())
        pubchem_model = pubchem_model.with_metadata(
            initial_structure_origin="pubchem_3d",
            conformer_index=0,
        )
        candidates.append(ConformerCandidate(pubchem_model, 0, None, "PubChem 3D", True))

    remaining = request.effective_conformer_pool - len(candidates)
    if remaining:
        if not smiles:
            raise APIError("SMILES is required to generate additional conformers")
        generated, _ = generate_forcefield_conformers(
            smiles,
            remaining,
            cid=request.cid,
            name=request.name or title,
            index_offset=len(candidates),
            random_seed=request.random_seed,
            prune_rms_threshold=request.embedding_prune_rmsd,
        )
        candidates.extend(candidate for candidate in generated if candidate.converged)
    if not candidates:
        raise ForceFieldError("No valid initial conformer is available for quantum optimization")
    mandatory = (0,) if pubchem_model is not None else ()
    return (
        _select_quantum_candidates(
            request,
            candidates,
            log,
            mandatory_indices=mandatory,
        ),
        smiles,
    )


def _validate_quantum_request(request: ConversionRequest, model: MoleculeModel, log: Logger) -> None:
    charge = request.effective_charge
    multiplicity = request.effective_multiplicity
    formal_charge = model.metadata.get("formal_charge")
    if request.charge is None and isinstance(formal_charge, int) and formal_charge != 0:
        raise ConfigurationError(
            f"The input structure has formal charge {formal_charge:+d}; pass --charge {formal_charge} explicitly"
        )
    if request.charge is not None and isinstance(formal_charge, int) and charge != formal_charge:
        log(
            f"Warning: requested charge {charge:+d} differs from input formal charge {formal_charge:+d}"
        )
    state = validate_electronic_state(
        (atom.element for atom in model.atoms),
        charge,
        multiplicity,
    )
    validate_metal_policy(
        state,
        method=request.method,
        allow_metals=request.allow_metals,
        basis_explicit=request.basis is not None,
        functional_explicit=request.functional is not None,
        charge_explicit=request.charge is not None,
        multiplicity_explicit=request.multiplicity is not None,
    )


def _match_external_reuse_to_selected_candidates(
    imported: Sequence[QuantumResult],
    candidates: Sequence[ConformerCandidate],
    *,
    rmsd_threshold: float,
    log: Logger,
    reuse_description: str = "external",
) -> List[QuantumResult]:
    """Retain only reusable Opt results represented in the current pre-screen.

    Legacy calculation numbering is not a conformer identity: for example,
    ``conformer_001`` in an older four-structure run may correspond to a pool
    member that is outside the current force-field energy window.  Match by
    optimized/current-candidate heavy-atom RMSD instead, without renaming or
    modifying the retained ORCA files.
    """

    edges = []
    for result_position, result in enumerate(imported):
        for candidate_position, candidate in enumerate(candidates):
            rmsd = aligned_heavy_atom_rmsd(result.model, candidate.model)
            if rmsd < rmsd_threshold:
                edges.append(
                    (rmsd, result_position, candidate_position)
                )
    used_results = set()
    used_candidates = set()
    matches: Dict[int, int] = {}
    for _rmsd, result_position, candidate_position in sorted(edges):
        if result_position in used_results or candidate_position in used_candidates:
            continue
        used_results.add(result_position)
        used_candidates.add(candidate_position)
        matches[result_position] = candidate_position

    retained = []
    copied_metadata = (
        "conformer_pool_index",
        "conformer_pool_index_provenance",
        "conformer_pool_size",
        "initial_rmsd_cluster_id",
        "conformer_selection_rmsd_threshold_angstrom",
        "conformer_selection_threshold_relaxed",
    )
    for result_position, result in enumerate(imported):
        candidate_position = matches.get(result_position)
        if candidate_position is None:
            log(
                f"Ignoring compatible {reuse_description} calculation "
                f"conformer_{result.conformer_index + 1:03d} for this ensemble: "
                "it does not match a candidate retained by the current "
                "force-field energy/RMSD pre-screen. The source files remain unchanged."
            )
            continue
        candidate = candidates[candidate_position]
        updates = {
            key: candidate.model.metadata[key]
            for key in copied_metadata
            if key in candidate.model.metadata
        }
        updates.update(
            {
                "matched_current_candidate_index": candidate.conformer_index,
                "reuse_selection_provenance": (
                    "matched_current_forcefield_prescreen_by_heavy_atom_rmsd"
                ),
            }
        )
        retained.append(
            replace(result, model=result.model.with_metadata(**updates))
        )
    return retained


def _run_quantum(
    request: ConversionRequest,
    candidates: Sequence[ConformerCandidate],
    backend,
    work_root: Path,
    log: Logger,
) -> Tuple[List[QuantumResult], List[dict]]:
    settings = QuantumSettings(
        method=request.method,
        functional=request.effective_functional,
        basis=request.effective_basis,
        charge=request.effective_charge,
        multiplicity=request.effective_multiplicity,
        frequency=request.frequency and not request.selective_frequency,
        timeout=request.quantum_timeout,
        max_opt_steps=request.max_opt_steps,
        max_scf_iterations=request.max_scf_iterations,
        nprocs=request.nprocs,
        maxcore_mb=request.maxcore_mb,
        imaginary_threshold_cm1=request.imaginary_threshold_cm1,
    )
    results: List[QuantumResult] = []
    failures: List[dict] = []
    work_root = Path(work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    run_candidates = list(candidates)

    def validate_stereo(result: QuantumResult) -> Optional[QuantumResult]:
        if request.smiles is None:
            return result
        validation = validate_model_stereochemistry(request.smiles, result.model)
        if not validation.matches:
            failures.append(
                {
                    "conformer_index": result.conformer_index,
                    "phase": "post_dft_validation",
                    "error_type": "StereochemistryMismatch",
                    "message": (
                        "Optimized geometry does not retain the explicit input "
                        "stereochemistry or atom ordering"
                    ),
                    "validation": validation.to_metadata(),
                }
            )
            log(
                f"Warning: conformer {result.conformer_index + 1} excluded because "
                "post-DFT stereochemistry validation failed"
            )
            return None
        return replace(
            result,
            model=result.model.with_metadata(
                post_dft_stereochemistry_validation=validation.to_metadata()
            ),
        )

    if request.ensemble:
        retained: List[QuantumResult] = []
        if request.reuse_calculations is not None:
            reuse_root = Path(request.reuse_calculations).expanduser().resolve()
            imported, _ = discover_reusable_orca_results(
                candidates[0], settings, backend, reuse_root, log, strict=True
            )
            imported = _match_external_reuse_to_selected_candidates(
                imported,
                candidates,
                rmsd_threshold=request.conformer_rmsd_threshold,
                log=log,
            )
            imported = [
                replace(
                    result,
                    model=result.model.with_metadata(
                        optimization_provenance="reused_external_read_only"
                    ),
                )
                for result in imported
            ]
            retained.extend(imported)
            log(
                f"Imported {len(imported)} completed calculations read-only from "
                f"{reuse_root}"
            )
        resumed, incomplete = discover_reusable_orca_results(
            candidates[0], settings, backend, work_root, log, strict=False
        )
        resumed = [
            replace(
                result,
                model=result.model.with_metadata(
                    optimization_provenance="reused_local_resume"
                ),
            )
            for result in resumed
        ]
        externally_matched_pool_indices = {
            result.model.metadata.get("conformer_pool_index")
            for result in retained
            if result.model.metadata.get("conformer_pool_index") is not None
        }
        locally_available_candidates = [
            candidate
            for candidate in candidates
            if candidate.model.metadata.get("conformer_pool_index")
            not in externally_matched_pool_indices
        ]
        resumed = _match_external_reuse_to_selected_candidates(
            resumed,
            locally_available_candidates,
            rmsd_threshold=request.conformer_rmsd_threshold,
            log=log,
            reuse_description="local resumed",
        )
        retained.extend(resumed)
        for record in incomplete:
            source = Path(str(record["directory"]))
            archived = archive_incomplete_calculation(source, work_root)
            record = dict(record)
            record["phase"] = "resume_recovery"
            record["archived_to"] = str(archived)
            record["recovery_action"] = "incomplete_attempt_archived_before_retry"
            failures.append(record)
            log(f"Archived incomplete attempt without deletion: {archived}")

        by_index: Dict[int, QuantumResult] = {}
        for retained_result in retained:
            if retained_result.conformer_index in by_index:
                raise ConfigurationError(
                    "Two reusable calculation roots contain the same conformer ID: "
                    f"{retained_result.conformer_index + 1}"
                )
            checked = validate_stereo(retained_result)
            if checked is not None:
                by_index[checked.conformer_index] = checked
        results.extend(by_index[index] for index in sorted(by_index))
        if len(results) > request.conformers:
            raise ConfigurationError(
                f"Found {len(results)} reusable conformers but --conformers is "
                f"{request.conformers}; refusing to discard completed calculations implicitly"
            )

        remaining = request.conformers - len(results)
        next_index = max(by_index, default=-1) + 1
        selected_new: List[ConformerCandidate] = []
        references = [result.model for result in results]
        for candidate in candidates:
            if len(selected_new) == remaining:
                break
            if any(
                aligned_heavy_atom_rmsd(candidate.model, model)
                < request.conformer_rmsd_threshold
                for model in references
            ):
                log(
                    "Skipping force-field candidate pool index "
                    f"{candidate.model.metadata.get('conformer_pool_index')} because it "
                    "duplicates a reusable optimized structure by heavy-atom RMSD"
                )
                continue
            index = next_index + len(selected_new)
            model = candidate.model.with_metadata(
                conformer_index=index,
                resumed_run_candidate=True,
            )
            selected_new.append(
                ConformerCandidate(
                    model,
                    index,
                    candidate.energy_kcal_mol,
                    candidate.forcefield,
                    candidate.converged,
                    candidate.optimization_iterations,
                    candidate.failure_status,
                )
            )
        run_candidates = selected_new
        log(
            f"Ensemble DFT plan: {len(results)} completed/reused + "
            f"{len(run_candidates)} new = {len(results) + len(run_candidates)} total"
        )

    for candidate in run_candidates:
        index = candidate.conformer_index
        if request.reuse_calculations is not None and not request.ensemble:
            reused = load_reusable_orca_result(
                candidate,
                settings,
                backend,
                work_root,
                log,
            )
            if reused is not None:
                reused = replace(
                    reused,
                    model=reused.model.with_metadata(
                        optimization_provenance="reused_requested_root"
                    ),
                )
                checked = validate_stereo(reused)
                if checked is not None:
                    results.append(checked)
                log(
                    f"Reused quantum optimization {index + 1}/{len(candidates)}; "
                    f"energy={reused.final_energy_hartree:.12f} Eh"
                )
                continue
        log(f"Quantum optimization {index + 1}/{len(candidates)} started")
        try:
            result = backend.optimize(
                candidate.model,
                settings,
                work_root / f"conformer_{index + 1:03d}",
                index,
            )
        except QuantumCalculationError as exc:
            failures.append(
                {
                    "conformer_index": index,
                    "phase": "dft_optimization",
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                }
            )
            log(f"Warning: conformer {index + 1} failed: {exc}")
            if len(candidates) == 1:
                raise
            continue
        checked = validate_stereo(result)
        if checked is None:
            continue
        checked = replace(
            checked,
            model=checked.model.with_metadata(
                optimization_provenance="computed_this_run"
            ),
        )
        results.append(checked)
        log(
            f"Quantum optimization {index + 1} converged; "
            f"energy={result.final_energy_hartree:.12f} Eh"
        )
    if not results:
        raise QuantumCalculationError("No quantum geometry optimization converged")
    return results, failures


def _run_selective_frequencies(
    request: ConversionRequest,
    results: Sequence[QuantumResult],
    backend,
    calculation_root: Path,
    log: Logger,
) -> Tuple[List[QuantumResult], Optional[dict]]:
    if not request.selective_frequency:
        return list(results), None

    lowest_energy = min(result.final_energy_hartree for result in results)
    window = request.frequency_window_kj
    ranked = sorted(results, key=lambda result: result.final_energy_hartree)
    eligible = [
        result
        for result in ranked
        if window is None
        or (result.final_energy_hartree - lowest_energy) * HARTREE_TO_KJ_MOL
        <= window + 1.0e-9
    ]
    limit = request.frequency_max or len(eligible)
    selected = eligible[:limit]
    selected_by_index = {result.conformer_index: result for result in selected}
    ranked_by_index = {result.conformer_index: result for result in ranked}
    for index in request.frequency_include:
        retained = ranked_by_index.get(index)
        if retained is None:
            log(
                f"Warning: conformer {index + 1} requested by --frequency-include is not "
                "a final post-DFT representative"
            )
        elif index not in selected_by_index:
            selected.append(retained)
            selected_by_index[index] = retained
            log(
                f"Conformer {index + 1} retained for Freq by --frequency-include "
                "despite the automatic energy/count selection"
            )
    selected_indices = {result.conformer_index for result in selected}
    frequency_settings = QuantumSettings(
        method=request.method,
        functional=request.effective_functional,
        basis=request.effective_basis,
        charge=request.effective_charge,
        multiplicity=request.effective_multiplicity,
        frequency=True,
        timeout=request.quantum_timeout,
        max_opt_steps=request.max_opt_steps,
        max_scf_iterations=request.max_scf_iterations,
        nprocs=request.nprocs,
        maxcore_mb=request.maxcore_mb,
        imaginary_threshold_cm1=request.imaginary_threshold_cm1,
    )
    completed = {}
    computed_this_run = set()
    reused_or_preexisting = set()
    frequency_failures = []
    for position, result in enumerate(selected, 1):
        delta_kj = (
            result.final_energy_hartree - lowest_energy
        ) * HARTREE_TO_KJ_MOL
        log(
            f"Selective frequency {position}/{len(selected)} for conformer "
            f"{result.conformer_index + 1}; deltaE={delta_kj:.3f} kJ/mol"
        )
        try:
            completed_result = ensure_quantum_result_frequency(
                result,
                frequency_settings,
                backend,
                calculation_root,
                log,
            )
            completed[result.conformer_index] = completed_result
            if (
                completed_result.model.metadata.get("frequency_provenance")
                == "computed_this_run"
            ):
                computed_this_run.add(result.conformer_index)
            else:
                reused_or_preexisting.add(result.conformer_index)
        except QuantumCalculationError as exc:
            frequency_failures.append(
                {
                    "conformer_index": result.conformer_index,
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                }
            )
            log(f"Warning: conformer {result.conformer_index + 1} Freq failed: {exc}")
    updated = [completed.get(result.conformer_index, result) for result in results]
    selection = {
        "mode": "electronic_energy_window",
        "window_kj_mol": window,
        "maximum_jobs": request.frequency_max,
        "selected_conformer_indices": sorted(selected_indices),
        "completed_conformer_indices": sorted(completed),
        "computed_this_run_conformer_indices": sorted(computed_this_run),
        "reused_or_preexisting_conformer_indices": sorted(reused_or_preexisting),
        "skipped_conformer_indices": sorted(
            result.conformer_index
            for result in results
            if result.conformer_index not in selected_indices
        ),
        "energy_reference_hartree": lowest_energy,
        "explicitly_included_conformer_indices": list(request.frequency_include),
        "failures": frequency_failures,
    }
    return updated, selection


def _execute_quantum(
    request: ConversionRequest,
    candidates: Sequence[ConformerCandidate],
    backend,
    calculation_root: Path,
    log: Logger,
) -> Tuple[List[QuantumResult], List[dict], Optional[dict], Optional[dict]]:
    results, failures = _run_quantum(
        request, candidates, backend, calculation_root, log
    )
    post_dft_screening = None
    if request.ensemble:
        optimized_count = len(results)
        results, post_dft_screening = cluster_optimized_results(
            results,
            rmsd_threshold=request.effective_dft_rmsd_threshold,
            energy_window_kj=request.dft_energy_window_kj,
        )
        log(
            f"Post-DFT screening retained {len(results)} representatives from "
            f"{optimized_count} converged optimizations"
        )
    results, frequency_selection = _run_selective_frequencies(
        request, results, backend, calculation_root, log
    )
    return results, failures, frequency_selection, post_dft_screening


def run_conversion(
    request: ConversionRequest,
    *,
    log: Optional[Logger] = None,
    blender_finder: Callable[[Optional[str]], str] = find_blender,
    orca_finder: Callable[[Optional[str]], str] = find_orca,
    orca_backend_class: Type[OrcaBackend] = OrcaBackend,
    exporter: Callable[..., Path] = export_model,
) -> ConversionOutcome:
    """Execute one conversion request without depending on argparse or a GUI."""

    request.validate()
    emit = log or (lambda _message: None)
    output_dir = Path(request.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    blender = blender_finder(request.blender_executable)

    if request.method in {"auto", "pubchem"} and request.cid is not None:
        try:
            model = download_and_parse(request.cid, timeout=request.api_timeout)
            if request.name:
                model = _rename_model(model, request.name)
            emit("Structure origin: pubchem_3d")
            base = safe_filename(model.name, f"CID_{request.cid}")
            artifact = _export_artifact(
                model,
                output_dir=output_dir,
                stem=base,
                blender_executable=blender,
                blender_timeout=request.blender_timeout,
                exporter=exporter,
            )
            return ConversionOutcome(model, (artifact,), "pubchem")
        except No3DConformerError:
            if request.method == "pubchem":
                raise
            emit("PubChem has no 3D conformer; auto mode is falling back to force field")

    if request.method in {"auto", "forcefield"}:
        candidates, _ = _forcefield_candidates(request)
        selected = _select_forcefield(candidates)
        converged = [candidate for candidate in candidates if candidate.converged]
        search = {
            "requested_conformers": request.conformers,
            "converged_conformers": len(converged),
            "selected_conformer_index": selected.conformer_index,
            "energies_kcal_mol": [
                {
                    "conformer_index": candidate.conformer_index,
                    "energy": candidate.energy_kcal_mol,
                }
                for candidate in converged
            ],
        }
        selected_model = selected.model.with_metadata(conformer_search=search)
        emit("Structure origin: computed_forcefield")
        base = safe_filename(selected_model.name, f"CID_{request.cid}" if request.cid else "molecule")
        artifacts = [
            _export_artifact(
                selected_model,
                output_dir=output_dir,
                stem=f"{base}_forcefield",
                blender_executable=blender,
                blender_timeout=request.blender_timeout,
                exporter=exporter,
            )
        ]
        if request.save_all_conformers and len(converged) > 1:
            for candidate in converged:
                artifacts.append(
                    _export_artifact(
                        candidate.model.with_metadata(conformer_search=search),
                        output_dir=output_dir,
                        stem=f"{base}_forcefield_conf{candidate.conformer_index + 1:03d}",
                        blender_executable=blender,
                        blender_timeout=request.blender_timeout,
                        exporter=exporter,
                    )
                )
        return ConversionOutcome(selected_model, tuple(artifacts), "forcefield")

    # Explicit HF/DFT is the only path that may launch a long quantum calculation.
    orca_executable = orca_finder(request.quantum_executable)
    backend = orca_backend_class(orca_executable)
    initial_payload = _quantum_initial_candidates(request, emit)
    candidates, _smiles = initial_payload[:2]
    initial_screening: Optional[Dict[str, object]] = (
        initial_payload[2] if len(initial_payload) > 2 else None
    )
    _validate_quantum_request(request, candidates[0].model, emit)
    if len(candidates) > 1 and request.reuse_calculations is None:
        emit(
            f"Warning: {len(candidates)} independent quantum optimizations will run; "
            "runtime is approximately multiplied by the conformer count"
        )
    elif len(candidates) > 1:
        emit(
            f"Requested {len(candidates)} optimized conformers; compatible completed "
            "calculations will be reused and only missing conformers will run"
        )

    base = safe_filename(
        candidates[0].model.name,
        f"CID_{request.cid}" if request.cid else "molecule",
    )
    persistent_calculations = (
        request.keep_calculation_files or request.reuse_calculations is not None
    )
    if request.ensemble:
        calculation_root = output_dir / "conformers"
        calculation_root.mkdir(parents=True, exist_ok=True)
        if request.reuse_calculations is not None:
            emit(
                "Reading retained ORCA calculations without writing to them: "
                f"{Path(request.reuse_calculations).expanduser().resolve()}"
            )
        emit(f"New and resumed calculations directory: {calculation_root}")
        results, failures, frequency_selection, post_dft_screening = _execute_quantum(
            request, candidates, backend, calculation_root, emit
        )
    elif request.reuse_calculations is not None:
        calculation_root = Path(request.reuse_calculations).expanduser().resolve()
        calculation_root.mkdir(parents=True, exist_ok=True)
        emit(f"Reusing ORCA calculation directory: {calculation_root}")
        results, failures, frequency_selection, post_dft_screening = _execute_quantum(
            request, candidates, backend, calculation_root, emit
        )
    elif request.keep_calculation_files:
        calculation_root = output_dir / f"{base}_{request.method}_calculations"
        calculation_root.mkdir(parents=True, exist_ok=True)
        results, failures, frequency_selection, post_dft_screening = _execute_quantum(
            request, candidates, backend, calculation_root, emit
        )
    else:
        with tempfile.TemporaryDirectory(prefix="molecule2fbx-quantum-") as temporary:
            results, failures, frequency_selection, post_dft_screening = _execute_quantum(
                request, candidates, backend, Path(temporary), emit
            )

    selected_result = min(results, key=lambda result: result.final_energy_hartree)
    search = {
        "requested_conformers": request.conformers,
        "converged_conformers": len(results),
        "selected_conformer_index": selected_result.conformer_index,
        "energies_hartree": [
            {
                "conformer_index": result.conformer_index,
                "energy": result.final_energy_hartree,
            }
            for result in results
        ],
        "failures": failures,
        "reused_conformers": [
            result.conformer_index
            for result in results
            if result.model.metadata.get("reused_quantum_result") is True
        ],
    }
    if request.conformer_pool is not None:
        search["conformer_pool"] = {
            "generated": request.effective_conformer_pool,
            "selected_for_quantum": request.conformers,
            "rmsd_threshold_angstrom": request.conformer_rmsd_threshold,
            "selected_pool_indices": [
                result.model.metadata.get("conformer_pool_index")
                for result in results
            ],
        }
    if frequency_selection is not None:
        search["frequency_selection"] = frequency_selection
    ensemble_report_path = None
    run_summary_path = None
    ensemble_payload = None
    if request.ensemble:
        if initial_screening is None or post_dft_screening is None:
            raise QuantumCalculationError("Ensemble screening metadata is incomplete")
        ensemble_payload = build_ensemble_report(
            request,
            results,
            initial_screening=initial_screening,
            post_dft_screening=post_dft_screening,
            frequency_selection=frequency_selection,
            dft_failures=failures,
        )
        ensemble_report_path = write_ensemble_report(
            ensemble_payload,
            output_dir,
            f"{base}_{request.method}",
            filename="ensemble.json",
        )
        run_summary_path = write_run_summary(ensemble_payload, output_dir)
        emit(f"Ensemble report: {ensemble_report_path}")
        emit(f"Run summary: {run_summary_path}")
    selected_model = selected_result.model.with_metadata(
        conformer_search=search,
        ensemble_summary=(
            {
                "report": str(ensemble_report_path),
                "best_conformer_id": ensemble_payload["best_conformer_id"],
                "best_structure_claim": ensemble_payload["interpretation"][
                    "best_structure_claim"
                ],
            }
            if ensemble_payload is not None and ensemble_report_path is not None
            else None
        ),
    )
    if persistent_calculations:
        selected_model = selected_model.with_metadata(
            calculation_files_directory=str(calculation_root)
        )
    emit("Structure origin: computed_quantum")
    artifact_output_dir = output_dir / "structures" if request.ensemble else output_dir
    artifact_output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = [
        _export_artifact(
            selected_model,
            output_dir=artifact_output_dir,
            stem=f"{base}_{request.method}",
            blender_executable=blender,
            blender_timeout=request.blender_timeout,
            exporter=exporter,
        )
    ]
    if request.save_all_conformers and len(results) > 1:
        for result in results:
            conformer_model = result.model.with_metadata(conformer_search=search)
            if persistent_calculations:
                conformer_model = conformer_model.with_metadata(
                    calculation_files_directory=str(calculation_root)
                )
            artifacts.append(
                _export_artifact(
                    conformer_model,
                    output_dir=artifact_output_dir,
                    stem=f"{base}_{request.method}_conf{result.conformer_index + 1:03d}",
                    blender_executable=blender,
                    blender_timeout=request.blender_timeout,
                    exporter=exporter,
                )
            )
    return ConversionOutcome(
        selected_model,
        tuple(artifacts),
        request.method,
        ensemble_report_path,
        run_summary_path,
    )

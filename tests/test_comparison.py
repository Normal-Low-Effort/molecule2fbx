import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from molecule2fbx.comparison import (
    _recorded_pool_index,
    _successful_selection_repair_pool_indices,
    _synchronize_ensemble_conformer_aliases,
    carbonyl_continuous_steric_access,
    carbonyl_steric_access,
    ensemble_descriptor_summary,
    parse_existing_orca_electronic,
    structural_descriptors,
    replace_terminal_tms_with_hydrogen,
)
from molecule2fbx.config import ConversionRequest
from molecule2fbx.ensemble import build_ensemble_report, cluster_optimized_results
from molecule2fbx.quantum.base import QuantumResult
from molecule2fbx.structures import (
    aligned_atom_subset_rmsd,
    aligned_mapped_atom_rmsd,
    common_scaffold_atom_mapping,
    generate_forcefield_conformers,
    rmsd_atom_subsets,
)


SB_SMILES = (
    "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn("
    "C(=O)c4ccc(cc4)[Si](C)(C)C)c5cccc(C2=C1)c35"
)
BZ_SMILES = (
    "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn("
    "C(=O)c4ccccc4)c5cccc(C2=C1)c35"
)


def _quantum_result(model, index, energy, *, frequency=False, reused=False):
    return QuantumResult(
        model=model.with_metadata(
            conformer_index=index,
            reused_quantum_result=reused,
            optimization_provenance=(
                "reused_external_read_only" if reused else "computed_this_run"
            ),
        ),
        backend="ORCA",
        backend_version="6.1.1",
        method="dft",
        functional="B3LYP",
        basis="def2-SVP",
        charge=0,
        multiplicity=1,
        final_energy_hartree=energy,
        geometry_converged=True,
        scf_converged=True,
        conformer_index=index,
        frequency_requested=frequency,
        frequencies_cm1=(25.0, 100.0) if frequency else (),
        imaginary_frequencies_cm1=(),
        calculation_directory=Path(f"conformer_{index + 1:03d}"),
    )


def test_tms_common_scaffold_excludes_terminal_si_and_three_carbons():
    candidate = generate_forcefield_conformers(SB_SMILES, 1)[0][0]
    subsets = rmsd_atom_subsets(candidate.model)
    excluded_elements = [
        candidate.model.atoms[index].element
        for index in subsets.excluded_terminal_substituent
    ]
    assert excluded_elements.count("Si") == 1
    assert excluded_elements.count("C") == 3
    assert len(subsets.common_scaffold) == len(subsets.all_heavy) - 4
    assert subsets.benzoyl_carbonyl_c in subsets.reaction_center
    assert subsets.benzoyl_carbonyl_o in subsets.reaction_center


def test_fixed_geometry_tms_to_h_counterfactual_preserves_common_atoms():
    sb = generate_forcefield_conformers(SB_SMILES, 1)[0][0].model
    bz = generate_forcefield_conformers(BZ_SMILES, 1)[0][0].model
    counterfactual = replace_terminal_tms_with_hydrogen(sb)

    assert all(atom.element != "Si" for atom in counterfactual.atoms)
    assert len(counterfactual.atoms) == len(bz.atoms)
    assert counterfactual.metadata["counterfactual_geometry_optimized"] is False
    removed = set(counterfactual.metadata["removed_source_atom_indices"])
    index_mapping = counterfactual.metadata[
        "source_to_counterfactual_atom_indices"
    ]
    retained_source_atoms = [atom for atom in sb.atoms if atom.index not in removed]
    for source in retained_source_atoms:
        retained = counterfactual.atoms[index_mapping[source.index]]
        assert retained.element == source.element
        assert retained.position == source.position

    anchor_source_index = counterfactual.metadata["source_aryl_anchor_atom_index"]
    anchor_new_index = index_mapping[anchor_source_index]
    anchor = counterfactual.atoms[anchor_new_index]
    hydrogen = counterfactual.atoms[
        counterfactual.metadata["replacement_hydrogen_atom_index"]
    ]
    distance = sum(
        (first - second) ** 2
        for first, second in zip(anchor.position, hydrogen.position)
    ) ** 0.5
    assert distance == pytest.approx(1.085)
    assert structural_descriptors(counterfactual)[
        "carbonyl_c_o_length_angstrom"
    ] == pytest.approx(
        structural_descriptors(sb)["carbonyl_c_o_length_angstrom"]
    )


def test_tms_only_motion_does_not_change_common_scaffold_rmsd():
    candidate = generate_forcefield_conformers(SB_SMILES, 1)[0][0]
    subsets = rmsd_atom_subsets(candidate.model)
    moved_atoms = []
    excluded = set(subsets.excluded_terminal_substituent)
    for atom in candidate.model.atoms:
        if atom.index in excluded:
            moved_atoms.append(
                type(atom)(atom.index, atom.element, atom.x + 2.0, atom.y, atom.z)
            )
        else:
            moved_atoms.append(atom)
    moved = type(candidate.model)(
        candidate.model.cid,
        candidate.model.name,
        tuple(moved_atoms),
        candidate.model.bonds,
        candidate.model.metadata,
    )
    assert aligned_atom_subset_rmsd(
        candidate.model, moved, subsets.all_heavy
    ) > 0.1
    assert aligned_atom_subset_rmsd(
        candidate.model, moved, subsets.common_scaffold
    ) == pytest.approx(0.0, abs=1.0e-8)


def test_bz_sb_common_scaffold_mapping_preserves_reaction_center():
    bz = generate_forcefield_conformers(BZ_SMILES, 1)[0][0].model
    sb = generate_forcefield_conformers(SB_SMILES, 1)[0][0].model
    mapping = common_scaffold_atom_mapping(bz, sb)
    bz_subsets = rmsd_atom_subsets(bz)
    sb_subsets = rmsd_atom_subsets(sb)
    assert len(mapping.first_indices) == len(bz_subsets.all_heavy)
    assert len(mapping.second_indices) == len(sb_subsets.common_scaffold)
    assert set(mapping.first_reaction_indices) == set(bz_subsets.reaction_center)
    assert set(mapping.second_reaction_indices) == set(sb_subsets.reaction_center)
    assert aligned_mapped_atom_rmsd(
        bz, sb, mapping.first_indices, mapping.second_indices
    ) >= 0.0


def test_structural_and_steric_analysis_identifies_benzoyl_site():
    candidate = generate_forcefield_conformers(SB_SMILES, 1)[0][0]
    structure = structural_descriptors(candidate.model)
    steric = carbonyl_steric_access(candidate.model, azimuth_samples=12)
    assert 1.1 < structure["carbonyl_c_o_length_angstrom"] < 1.4
    assert 1.2 < structure["amide_c_n_length_angstrom"] < 1.6
    assert 1.7 < structure["aryl_c_si_length_angstrom"] < 2.2
    assert 0.0 <= steric["overall_accessible_fraction"] <= 1.0
    assert steric["clearance_percentiles_angstrom"]["p90"] <= steric[
        "best_clearance_angstrom"
    ]
    assert steric["silicon_distance_from_carbonyl_c_angstrom"] > 5.0
    assert steric["terminal_tms_heavy_atoms_within_5_angstrom"] == 0
    assert steric["face_assignment"]["re_si_assigned"] is False


def test_continuous_steric_access_separates_local_and_remote_shielding():
    candidate = generate_forcefield_conformers(SB_SMILES, 1)[0][0]
    steric = carbonyl_continuous_steric_access(
        candidate.model,
        azimuth_samples=72,
        probe_radii_angstrom=(0.0, 0.5, 1.0),
        include_trajectory=True,
    )
    assert steric["face_assignment"]["re_si_assigned"] is True
    assert steric["face_assignment"]["coplanar_sample_count"] <= 2
    assert len(steric["trajectory"]) == 72
    assert len(steric["atom_indices"]["terminal_tms_heavy"]) == 4
    total = steric["scopes"]["total"]["probe_radius_sensitivity"]
    remote = steric["scopes"]["nonlocal_environment"][
        "probe_radius_sensitivity"
    ]
    assert remote["probe_radius_0.50_angstrom"][
        "best_clearance_angstrom"
    ] >= total["probe_radius_0.50_angstrom"]["best_clearance_angstrom"]
    assert total["probe_radius_0.00_angstrom"][
        "positive_clearance_fraction"
    ] >= total["probe_radius_1.00_angstrom"][
        "positive_clearance_fraction"
    ]
    assert steric["terminal_tms_as_total_limiter_fraction"] is not None
    direct = steric["direct_terminal_tms_counterfactual"][
        "probe_radius_sensitivity"
    ]["probe_radius_0.50_angstrom"]["total"]
    assert direct[
        "positive_clearance_integral_with_minus_without_tms_angstrom"
    ] <= 0.0


def test_continuous_steric_access_has_no_tms_scope_for_bz_control():
    candidate = generate_forcefield_conformers(BZ_SMILES, 1)[0][0]
    steric = carbonyl_continuous_steric_access(
        candidate.model,
        azimuth_samples=24,
    )
    assert steric["atom_indices"]["terminal_tms_heavy"] == []
    assert steric["terminal_tms_as_total_limiter_fraction"] is None
    assert steric["scopes"]["terminal_tms"]["clearance_without_probe"][
        "sample_count"
    ] == 0
    assert steric["direct_terminal_tms_counterfactual"][
        "probe_radius_sensitivity"
    ] == {}


def test_existing_orca_property_parser_uses_final_population_sections(tmp_path):
    output = tmp_path / "job.out"
    output.write_text(
        """
MULLIKEN ATOMIC CHARGES
  0 C : 0.100
  1 O : -0.100

LOEWDIN ATOMIC CHARGES
  0 C : 0.200
  1 O : -0.200

ORBITAL ENERGIES
 NO OCC E(Eh) E(eV)
  0 2.0000 -0.3000 -8.16
  1 0.0000 -0.1000 -2.72

MAYER POPULATION ANALYSIS
B(  0-C ,  1-O ) : 2.0500

DIPOLE MOMENT
Magnitude (Debye) : 4.2500

FRONTIER MOLECULAR ORBITAL POPULATION ANALYSIS
ANALYZING ORBITALS: HOMO= 10 LUMO= 11
  0-C  0.10 0.20 0.30 0.40
  1-O  0.50 0.60 0.70 0.80

MBIS ANALYSIS
  ATOM CHARGE POPULATION SPIN
  0 C 0.321 5.679 0.0
  1 O -0.321 8.321 0.0

CHELPG Charges
  0 C : 0.456
  1 O : -0.456
""",
        encoding="utf-8",
    )
    parsed = parse_existing_orca_electronic(output)
    assert parsed["mulliken_charges"][0] == pytest.approx(0.1)
    assert parsed["loewdin_charges"][0] == pytest.approx(0.2)
    assert parsed["mayer_bond_orders"][(0, 1)] == pytest.approx(2.05)
    assert parsed["dipole_magnitude_debye"] == pytest.approx(4.25)
    assert parsed["homo_lumo_gap_ev"] == pytest.approx(5.44)
    assert parsed["mbis_charges"][0] == pytest.approx(0.321)
    assert parsed["chelpg_charges"][0] == pytest.approx(0.456)
    frontier = parsed["frontier_orbital_local_contributions"]
    assert frontier["homo_index"] == 10
    assert frontier["atom_populations"][1]["lumo_loewdin"] == pytest.approx(0.8)


def test_ensemble_report_uses_null_for_uncomputed_imaginary_modes():
    candidates, _inspection = generate_forcefield_conformers("CCCC", 2)
    results = [
        _quantum_result(candidates[0].model, 0, -10.0, frequency=True),
        _quantum_result(candidates[1].model, 1, -9.999, frequency=False),
    ]
    post = {
        "records": [
            {"conformer_index": 0},
            {"conformer_index": 1},
        ],
        "cluster_count": 2,
    }
    request = ConversionRequest(
        smiles="CCCC",
        method="dft",
        ensemble=True,
        conformers=2,
        conformer_pool=100,
        frequency=True,
        frequency_window_kj=5.0,
        frequency_max=1,
    )
    payload = build_ensemble_report(
        request,
        results,
        initial_screening={"records": []},
        post_dft_screening=post,
        frequency_selection={
            "selected_conformer_indices": [0],
            "completed_conformer_indices": [0],
        },
        dft_failures=[
            {"phase": "resume_recovery", "message": "historical"}
        ],
    )
    entries = {item["conformer_index"]: item for item in payload["conformers"]}
    assert entries[1]["imaginary_modes"] is None
    assert payload["summary"]["frequency_available_total"] == 1
    assert payload["summary"]["frequency_completed_for_this_run_selection"] == 1
    assert payload["dft"]["failures"] == []
    assert len(payload["dft"]["recovery_history"]) == 1


def test_metadata_repair_synchronizes_legacy_conformer_alias():
    payload = {
        "final_conformer_ensemble": [
            {
                "conformer_id": "conf001",
                "frequency_calculated": False,
                "imaginary_modes": None,
            }
        ],
        "conformers": [
            {
                "conformer_id": "conf001",
                "frequency_calculated": False,
                "imaginary_modes": 0,
            }
        ],
    }

    _synchronize_ensemble_conformer_aliases(payload)

    assert payload["conformers"] == payload["final_conformer_ensemble"]
    assert payload["conformers"][0]["imaginary_modes"] is None
    assert payload["conformers"] is not payload["final_conformer_ensemble"]


def test_pool_index_is_recovered_from_external_reuse_metadata(tmp_path):
    calculation_root = tmp_path / "calculations"
    calculation_dir = calculation_root / "conformer_002"
    calculation_dir.mkdir(parents=True)
    (tmp_path / "conf002.metadata.json").write_text(
        json.dumps(
            {
                "metadata": {
                    "calculation_files_directory": str(calculation_root.resolve()),
                    "conformer_index": 1,
                    "conformer_pool_index": 78,
                }
            }
        ),
        encoding="utf-8",
    )

    assert _recorded_pool_index(
        tmp_path / "ensemble",
        "conf002",
        calculation_dir,
    ) == 78


def test_successful_selection_repair_pool_is_preserved(tmp_path):
    (tmp_path / "strict_selection_repair_plan.json").write_text(
        json.dumps(
            {
                "status": "SUCCESS",
                "target": {"pool_index": 151},
                "result": {"geometry_optimization_converged": True},
            }
        ),
        encoding="utf-8",
    )

    assert _successful_selection_repair_pool_indices(tmp_path) == {151}


def test_descriptor_summary_deduplicates_common_scaffold_clusters():
    records = [
        {
            "conformer_id": "conf001",
            "common_scaffold_rmsd_cluster_id": "core001",
            "dft_energy_hartree": -10.0,
            "gibbs_energy_hartree": None,
            "frequency_calculated": False,
            "structure": {"carbonyl_c_o_length_angstrom": 1.20},
            "electronic": {"mbis_charge_carbonyl_o": -0.50},
        },
        {
            "conformer_id": "conf002",
            "common_scaffold_rmsd_cluster_id": "core001",
            "dft_energy_hartree": -9.9999,
            "gibbs_energy_hartree": None,
            "frequency_calculated": False,
            "structure": {"carbonyl_c_o_length_angstrom": 1.30},
            "electronic": {"mbis_charge_carbonyl_o": -0.40},
        },
        {
            "conformer_id": "conf003",
            "common_scaffold_rmsd_cluster_id": "core002",
            "dft_energy_hartree": -9.999,
            "gibbs_energy_hartree": None,
            "frequency_calculated": False,
            "structure": {"carbonyl_c_o_length_angstrom": 1.25},
            "electronic": {"mbis_charge_carbonyl_o": -0.45},
        },
    ]
    summary = ensemble_descriptor_summary(
        {"name": "test", "final_conformer_ensemble": records}, records
    )
    unique = summary["common_scaffold_unique"]
    assert unique["representative_count"] == 2
    assert unique["representative_conformer_ids"] == ["conf001", "conf003"]
    oxygen = unique["metrics"]["electronic.mbis_charge_carbonyl_o"][
        "electronic_energy_weighted"
    ]
    assert oxygen["conformer_count"] == 2
    assert -0.50 <= oxygen["weighted_mean"] <= -0.45
    assert records[1]["common_scaffold_unique_electronic_weight"] is None

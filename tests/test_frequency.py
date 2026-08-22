import json

from molecule2fbx.config import FrequencyOnlyRequest
from molecule2fbx.frequency import read_orca_input_settings, run_frequency_only
from molecule2fbx.quantum.base import FrequencyResult


def test_reads_original_charge_multiplicity_and_method(tmp_path):
    source = tmp_path / "conformer_001.inp"
    source.write_text(
        """! PBE0 def2-TZVP Opt TightSCF
%pal
  nprocs 12
end
%maxcore 2000
* xyz -1 2
H 0 0 0
*
""",
        encoding="utf-8",
    )
    parsed = read_orca_input_settings(source)
    assert parsed.method == "dft"
    assert parsed.functional == "PBE0"
    assert parsed.basis == "def2-TZVP"
    assert parsed.charge == -1
    assert parsed.multiplicity == 2
    assert parsed.nprocs == 12
    assert parsed.maxcore_mb == 2000


def test_frequency_only_updates_related_metadata_without_opt(tmp_path):
    output_dir = tmp_path / "output"
    calculation_root = output_dir / "Hydrogen_dft_calculations"
    conformer_dir = calculation_root / "conformer_001"
    conformer_dir.mkdir(parents=True)
    xyz_path = conformer_dir / "conformer_001.xyz"
    xyz_path.write_text("2\noptimized\nH 0 0 0\nH 0.74 0 0\n", encoding="utf-8")
    xyz_path.with_suffix(".inp").write_text(
        """! B3LYP def2-SVP Opt TightSCF
* xyz 0 1
H 0 0 0
H 0.74 0 0
*
""",
        encoding="utf-8",
    )
    metadata_path = output_dir / "Hydrogen_dft.metadata.json"
    metadata_path.write_text(
        json.dumps(
            {
                "name": "Hydrogen",
                "cid": None,
                "metadata": {
                    "calculation_files_directory": str(calculation_root.resolve()),
                    "quantum_method": "dft",
                    "functional": "B3LYP",
                    "basis": "def2-SVP",
                    "charge": 0,
                    "multiplicity": 1,
                },
            }
        ),
        encoding="utf-8",
    )
    captured = {}

    class FakeBackend:
        def __init__(self, executable):
            captured["executable"] = executable

        def frequency(self, model, settings, work_dir, stem):
            captured["model"] = model
            captured["settings"] = settings
            input_path = work_dir / f"{stem}.inp"
            output_path = work_dir / f"{stem}.out"
            input_path.write_text("! B3LYP def2-SVP Freq TightSCF\n", encoding="utf-8")
            output_path.write_text("mock", encoding="utf-8")
            return FrequencyResult(
                backend="ORCA",
                backend_version="6.1.0",
                method="dft",
                functional="B3LYP",
                basis="def2-SVP",
                charge=0,
                multiplicity=1,
                final_energy_hartree=-1.125,
                frequencies_cm1=(-25.0, -5.0, 100.0),
                imaginary_frequencies_cm1=(-25.0,),
                warnings=("imaginary frequency",),
                input_path=input_path,
                output_path=output_path,
            )

    outcome = run_frequency_only(
        FrequencyOnlyRequest(xyz_path=xyz_path, nprocs=8),
        orca_finder=lambda explicit: "orca",
        orca_backend_class=FakeBackend,
    )
    assert captured["settings"].charge == 0
    assert captured["settings"].multiplicity == 1
    assert captured["settings"].nprocs == 8
    assert [atom.element for atom in captured["model"].atoms] == ["H", "H"]
    assert captured["model"].atoms[1].x == 0.74
    assert outcome.updated_metadata_path == metadata_path
    updated = json.loads(metadata_path.read_text(encoding="utf-8"))
    analysis = updated["metadata"]["frequency_analyses"][-1]
    assert analysis["geometry_reoptimized"] is False
    assert analysis["imaginary_frequencies_cm1"] == [-25.0]
    assert outcome.result_metadata_path.is_file()

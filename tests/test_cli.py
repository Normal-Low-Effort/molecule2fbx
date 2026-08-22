from pathlib import Path

import pytest

from molecule2fbx.cli import _frequency_request_from_args, _request_from_args, build_parser
from molecule2fbx.config import ConversionRequest


def test_nprocs_defaults_to_logical_processor_count(monkeypatch):
    monkeypatch.setattr("molecule2fbx.cli.os.cpu_count", lambda: 12)
    args = build_parser().parse_args(["--smiles", "O", "--method", "dft"])
    request = _request_from_args(args)
    assert request.nprocs == 12


@pytest.mark.parametrize("count", ["8", "12"])
def test_explicit_nprocs_overrides_auto_detection(monkeypatch, count):
    monkeypatch.setattr("molecule2fbx.cli.os.cpu_count", lambda: 12)
    args = build_parser().parse_args(
        ["--smiles", "O", "--method", "dft", "--nprocs", count]
    )
    request = _request_from_args(args)
    assert request.nprocs == int(count)


@pytest.mark.parametrize("value", ["0", "-1", "1.5", "abc"])
def test_invalid_nprocs_is_rejected_by_cli(value):
    with pytest.raises(SystemExit) as exc_info:
        build_parser().parse_args(["--smiles", "O", "--method", "dft", "--nprocs", value])
    assert exc_info.value.code == 2


def test_frequency_only_cli_does_not_require_cid_or_smiles(tmp_path):
    xyz = tmp_path / "conformer_001.xyz"
    args = build_parser().parse_args(
        ["--frequency-only", str(xyz), "--nprocs", "8"]
    )
    request = _frequency_request_from_args(args)
    assert request.xyz_path == xyz
    assert request.nprocs == 8
    assert request.method is None


def test_reuse_calculations_cli_option_is_preserved():
    args = build_parser().parse_args(
        [
            "--smiles",
            "O",
            "--method",
            "dft",
            "--conformers",
            "4",
            "--reuse-calculations",
            "saved-calculations",
        ]
    )
    request = _request_from_args(args)
    assert request.reuse_calculations == Path("saved-calculations")
    assert request.conformers == 4


def test_staged_conformer_and_frequency_options_are_parsed():
    args = build_parser().parse_args(
        [
            "--smiles",
            "O",
            "--method",
            "dft",
            "--conformer-pool",
            "200",
            "--conformers",
            "10",
            "--conformer-rmsd-threshold",
            "0.8",
            "--frequency",
            "--frequency-window-kj",
            "5",
            "--frequency-max",
            "3",
        ]
    )
    request = _request_from_args(args)
    request.validate()
    assert request.conformer_pool == 200
    assert request.conformers == 10
    assert request.conformer_rmsd_threshold == pytest.approx(0.8)
    assert request.frequency_window_kj == pytest.approx(5.0)
    assert request.frequency_max == 3
    assert request.selective_frequency is True


def test_conformer_pool_cannot_be_smaller_than_quantum_count():
    request = ConversionRequest(
        smiles="O", method="dft", conformer_pool=4, conformers=10
    )
    with pytest.raises(ValueError, match="conformer-pool"):
        request.validate()


def test_selective_frequency_requires_frequency_flag():
    request = ConversionRequest(
        smiles="O", method="dft", frequency_window_kj=5.0, frequency_max=3
    )
    with pytest.raises(ValueError, match="require --frequency"):
        request.validate()

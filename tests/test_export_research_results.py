from pathlib import Path

from scripts.export_research_results import (
    ENSEMBLE_NAMES,
    SNAPSHOT_FILES,
    ensure_within,
    portable_snapshot_bytes,
)


def test_public_snapshot_includes_all_comparison_ensembles():
    sources = {source for source, _ in SNAPSHOT_FILES}

    for name in ENSEMBLE_NAMES:
        assert f"{name}/ensemble.json" in sources
        assert f"{name}/RUN_SUMMARY.md" in sources


def test_snapshot_destinations_are_unique_and_relative():
    destinations = [destination for _, destination in SNAPSHOT_FILES]
    assert len(destinations) == len(set(destinations))
    assert all(not Path(destination).is_absolute() for destination in destinations)
    assert all(".." not in Path(destination).parts for destination in destinations)


def test_portable_snapshot_redacts_normal_and_json_escaped_workspace(tmp_path):
    source = tmp_path / "result.json"
    normal = str(tmp_path)
    escaped = normal.replace("\\", "\\\\")
    source.write_text(
        '{"normal": "' + normal + '", "escaped": "' + escaped + '"}',
        encoding="utf-8",
    )

    data, replacements = portable_snapshot_bytes(source, tmp_path)
    text = data.decode("utf-8")

    assert replacements == 2
    assert normal not in text
    assert escaped not in text
    assert text.count("${WORKSPACE}") == 2


def test_portable_snapshot_removes_leading_utf8_bom(tmp_path):
    source = tmp_path / "powershell-result.json"
    source.write_text('\ufeff{"status": "SUCCESS"}', encoding="utf-8")

    data, replacements = portable_snapshot_bytes(source, tmp_path)

    assert replacements == 0
    assert data == b'{"status": "SUCCESS"}'


def test_ensure_within_rejects_parent_escape(tmp_path):
    inside = tmp_path / "inside" / "result.json"
    ensure_within(inside, tmp_path)

    outside = tmp_path.parent / "outside.json"
    try:
        ensure_within(outside, tmp_path)
    except ValueError as exc:
        assert "escapes expected root" in str(exc)
    else:
        raise AssertionError("Expected a parent-directory escape to be rejected")

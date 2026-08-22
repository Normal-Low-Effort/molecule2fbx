# Curated research results

This directory contains small, Git-friendly snapshots copied from the local
calculation tree under `outputs/`.

- `outputs/` remains the canonical source of ORCA, xTB, runEDDB, Freq, and FBX
  calculation records and is intentionally ignored by Git.
- Files here are summaries, parsed tables, or selected optimized XYZ files.
- They are computational results, not experimental structures or measurements.
- Local workspace prefixes in text files are replaced with `${WORKSPACE}` so
  the Git snapshot does not publish a user-specific absolute path.
- `SNAPSHOT_MANIFEST.json` records each source path, byte size, source SHA-256,
  snapshot SHA-256, and the number of path replacements.

Refresh the snapshot from the repository root with:

```powershell
python scripts\export_research_results.py
```

Validate the source set without copying with:

```powershell
python scripts\export_research_results.py --check
```

No ORCA/xTB calculation is started by this export operation.

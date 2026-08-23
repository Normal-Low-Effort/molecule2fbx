# Research scripts

These scripts operate on the local, Git-ignored `outputs/` calculation tree.
Run them from the repository root so the installed `molecule2fbx` package and
external ORCA/xTB/runEDDB programs can be located consistently.

## Groups

- `night_run_*.ps1`, `run_missing_freqs.ps1`: ensemble and supplemental Freq runs.
- `run_para_alkyl_control_ensembles.ps1`: sequential, resumable p-Me, p-iPr,
  and p-tBu control ensembles using the same B3LYP/def2-SVP settings as the
  existing 1Bz-LSD/1SB-LSD comparison.
- `run_*.py`: non-destructive ORCA, xTB, and runEDDB follow-up calculations.
- `analyze_*.py`: parsing and comparison of existing calculations.
- `analyze_para_substituent_series.py`: read-only five-series common-scaffold
  and carbonyl-access analysis; it also refreshes derived ensemble metadata.
- `export_research_results.py`: copy the approved lightweight result set into
  the Git-managed `research_results/` snapshot.

The calculation scripts continue to write to `outputs/`. They do not use
`research_results/` as a calculation source.

PowerShell run scripts resolve the repository from `$PSScriptRoot`. Set
`ORCA_EXECUTABLE` and `BLENDER_EXECUTABLE` before ensemble runs. Set
`MOLECULE2FBX_PYTHON` only when a Python executable other than the local
`work/test-venv` or `python` on PATH should be used.

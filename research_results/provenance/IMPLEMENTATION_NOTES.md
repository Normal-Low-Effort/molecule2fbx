# Preliminary-comparison implementation notes

## Changed source files

- `molecule2fbx/structures.py`
  - adds all-heavy/common-scaffold/reaction-centre atom subsets;
  - adds deterministic cross-molecule Bz/TMS-Bz graph mapping;
  - keeps fixed input atom identity and does not silently permute symmetric atoms.
- `molecule2fbx/ensemble.py`
  - writes `imaginary_modes: null` when Freq is absent;
  - separates selected/completed/available/reused Freq counts;
  - separates final DFT failures from resume-recovery history;
  - records Opt/Freq and conformer-pool provenance;
  - writes secondary common-scaffold and reaction-centre RMSD metadata.
- `molecule2fbx/pipeline.py` and `molecule2fbx/quantum/reuse.py`
  - distinguish newly computed, locally resumed, and externally reused results;
  - avoid assigning current-pool force-field metadata to legacy reused results.
- `molecule2fbx/comparison.py`
  - repairs retained ensemble metadata non-destructively (with one-time backup);
  - parses retained ORCA population, Mayer, dipole, orbital, Freq, Hessian, and
    thermochemistry data;
  - calculates carbonyl-centred directional steric access;
  - calculates within-ensemble and cross-ensemble RMSDs;
  - distinguishes all-representative weights, common-scaffold-unique weights,
    and conditional Freq/Gibbs weights;
  - writes JSON, CSV, English summary, and Japanese assessment artifacts.
- `tests/test_comparison.py`
  - tests TMS exclusion, TMS-only motion, Bz/SB graph mapping, steric descriptors,
    ORCA property parsing, null Freq semantics, provenance/failure separation,
    and common-scaffold weight de-duplication.
- `README.md`
  - documents RMSD atom correspondence, symmetry assumptions, weight scope,
    and analysis regeneration.

## Analysis runners

```powershell
work\test-venv\Scripts\python.exe work\analyze_bz_sb_preliminary.py
work\test-venv\Scripts\python.exe work\run_sb_pool_xtb.py
work\test-venv\Scripts\python.exe work\run_low_energy_properties.py
```

`work/run_sb_pool_singlepoints.py` is a conditional fallback for the two
omitted SB basins. It is run only if the GFN2-xTB screen justifies it.

## Test status

```text
69 passed
```

No existing Opt result is overwritten or deleted by these changes. New Freq,
xTB, and property results are written to separate directories. Historical
incomplete launch attempts remain archived and recoverable.

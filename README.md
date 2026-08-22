# molecule2fbx 0.4.0

[Japanese version](READMEja.md)

`molecule2fbx` is a command-line tool that generates ball-and-stick molecular
models from a PubChem CID or SMILES and exports them to FBX through Blender.
It distinguishes PubChem 3D structures, RDKit force-field structures, and ORCA
quantum-chemistry structures, and records their provenance and calculation
settings in JSON metadata.

Version 0.4.0 adds an `--ensemble` mode that runs ETKDG conformer generation,
DFT geometry optimization, post-DFT duplicate removal, selective frequency
calculations, and relative energy and Gibbs-energy analysis as one workflow.

Computed structures are not experimental structures. The lowest-energy
structure in a run is only the current best among the candidates examined; it
is not proof of the true global minimum.

## Existing-ensemble comparison analysis

`molecule2fbx.comparison` adds analysis views without replacing the historical
all-heavy-atom clustering:

- `all_heavy_rmsd_cluster_id`: fixed input atom order, all heavy atoms.
- `common_scaffold_rmsd_cluster_id`: fixed input atom order with a terminal
  aryl-Si(CH3)3 branch, consisting of Si and its three methyl carbons, excluded.
- `reaction_center_rmsd_cluster_id`: fixed input atom order for heavy atoms no
  more than two graph bonds from the aromatic N-benzoyl carbonyl carbon.

No symmetry permutation is applied. ORCA XYZ files retain the input atom order,
while automatic permutation can silently change atom identity. The
common-scaffold metric removes the whole terminal TMS branch so that rotation
or permutation of its equivalent methyl groups does not create a new
LSD-benzoyl scaffold cluster. This additional descriptor does not rewrite
coordinates or alter the original all-heavy cluster IDs.

For direct Bz-versus-TMS-Bz RMSD, the analysis removes the terminal TMS branch
from the RDKit graph, finds a chirality-preserving graph isomorphism, anchors
the mapping at the benzoyl C/O/N/ipso atoms, and selects the mapping with the
smallest displacement in retained input-atom order. Kabsch alignment is applied
to that deterministic mapping. Symmetry-equivalent benzene or ethyl-branch
permutations are not searched to minimize RMSD. Both the selected atom map and
this assumption are recorded in `analysis.json`. The cross-molecule RMSD is
therefore reproducible and preserves atom identity, but it can be an upper
bound when an unsearched symmetry-equivalent mapping would be lower.

Descriptor reports retain the historical all-representative weights and also
provide `common_scaffold_unique`, in which the lowest-electronic-energy
structure from each common-scaffold cluster is counted once. The latter is used
for the main Bz/SB sensitivity comparison so that TMS-only rotational
duplicates are not counted as independent scaffold conformers. It is not a
rigorous treatment of rotor degeneracy or configurational entropy and is
labelled accordingly.

Regenerate the preliminary Bz/SB analysis with:

```powershell
work\test-venv\Scripts\python.exe scripts\analyze_bz_sb_preliminary.py
```

The analysis uses `imaginary_modes: null` when no frequency calculation exists,
separates final failures from resume-recovery history, records Opt/Freq
provenance, regenerates the deterministic force-field pool, extracts retained
ORCA properties, and reports directional carbonyl access on a sampled
Burgi-Dunitz cone. The steric probe is a reactant-geometry descriptor only; it
is not an enzyme, solvent, or activation-barrier model. Gibbs-weighted values
are labelled conditional unless frequency results exist for every final DFT
representative.

## Repository and research data

The repository separates Git-managed files from large local calculation data.

```text
molecule2fbx/      CLI implementation
tests/             automated tests
scripts/           ensemble, follow-up calculation, and analysis scripts
research_results/  lightweight, Git-managed analysis snapshots
outputs/           ORCA, xTB, runEDDB, and FBX records excluded from Git
work/              virtual environments, tests, and packaging scratch space
tools/             locally installed external tools excluded from Git
```

`outputs/` is the canonical calculation record. `research_results/` contains
only explicitly selected small summaries, tables, and current-best XYZ files
copied by `scripts/export_research_results.py`; it is not a replacement for the
ORCA source data.

```powershell
python scripts\export_research_results.py --check
python scripts\export_research_results.py
```

Local absolute workspace paths in the snapshot are replaced with
`${WORKSPACE}`. `research_results/SNAPSHOT_MANIFEST.json` records the source,
size, source SHA-256, and snapshot SHA-256 for every copied file.

## Requirements

- Python 3.9 or later
- Blender 3.x or 4.x
- RDKit, NumPy, and requests
- A separately installed ORCA 6.x for DFT or HF calculations

```powershell
python -m pip install .
```

For development:

```powershell
python -m pip install -e ".[dev]"
python -m pytest
```

ORCA is not bundled because of its licensing and distribution conditions.
Specify it with `--orca`, `ORCA_EXECUTABLE`, `ORCADIR`, or `PATH`.

```powershell
$env:ORCA_EXECUTABLE = "C:\Orca_6.1.1\orca.exe"
$env:BLENDER_EXECUTABLE = "C:\Program Files\Blender Foundation\Blender 4.0\blender.exe"
```

## Existing CLI

The original positional CID interface remains supported.

```powershell
molecule2fbx 5360697
molecule2fbx --cid 5360697 --method auto
```

`auto` prefers PubChem 3D data and falls back to ETKDG plus MMFF94s or UFF only
when PubChem 3D is unavailable. It never starts a quantum-chemistry calculation.

```powershell
molecule2fbx --smiles "CCO" --name Ethanol --method forcefield
molecule2fbx --smiles "O" --name Water --method dft
```

The default ORCA electronic-structure model is `B3LYP/def2-SVP`, charge 0, and
multiplicity 1.

```powershell
molecule2fbx --smiles "[O]" --method dft `
  --functional B3LYP --basis "6-31G(d)" `
  --charge 0 --multiplicity 3
```

When `--nprocs` is omitted, `os.cpu_count()` is used. If Windows exposes a
Ryzen 5 5600 as 12 logical processors, the generated ORCA input contains
`%pal nprocs 12 end`.

## Conformer ensemble mode

Provide a SMILES string with complete stereochemistry. In ensemble mode,
omitting `--method` selects DFT because `--ensemble` itself is explicit
authorization for a potentially long calculation.

```powershell
molecule2fbx --ensemble `
  --smiles "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn(C(=O)c4ccc(cc4)[Si](C)(C)C)c5cccc(C2=C1)c35" `
  --name 1SB-LSD_RR `
  --nprocs 16 --maxcore 1000
```

Default ensemble workflow:

```text
Validate SMILES, charge, multiplicity, and stereochemistry
  -> generate 200 ETKDG conformers
  -> optimize with MMFF94s, falling back to UFF
  -> retain structures within 10 kJ/mol of the lowest force-field energy
  -> cluster by heavy-atom RMSD at 0.75 angstrom
  -> optimize up to 10 representatives with ORCA B3LYP/def2-SVP
  -> recluster DFT structures at 0.75 angstrom
  -> run Freq on up to 3 structures within 5 kJ/mol of the lowest electronic energy
  -> write frequencies, ZPE, thermal corrections, H, G, relative E, and relative G
```

Important options:

| Option | Ensemble default | Description |
|---|---:|---|
| `--conformer-pool` | 200 | Number of ETKDG structures; range is 100 to 500 |
| `--conformers` | 10 | Maximum cluster representatives sent to DFT |
| `--random-seed` | 61453 | ETKDG random seed |
| `--embedding-prune-rmsd` | -1 | Embedding pruning; -1 disables it |
| `--forcefield-energy-window-kj` | 10 | Force-field energy window before DFT |
| `--conformer-rmsd-threshold` | 0.75 angstrom | Pre-DFT heavy-atom RMSD threshold |
| `--dft-energy-window-kj` | unlimited | Optional post-DFT energy window |
| `--dft-rmsd-threshold` | 0.75 angstrom | Post-DFT heavy-atom RMSD threshold |
| `--frequency-window-kj` | 5 | Electronic-energy window for Freq selection |
| `--frequency-max` | 3 | Maximum automatically selected Freq structures |
| `--frequency-include N` | none | Retain conformer N for Freq even at higher energy |
| `--imaginary-threshold-cm1` | -20 | Frequencies below this count as imaginary modes |
| `--low-frequency-threshold-cm1` | 50 | Frequencies below this are reported as low modes |

RMSD thresholds are never relaxed automatically. If fewer than 10 independent
clusters exist, only the available representatives are sent to DFT.

The RDKit bulk MMFF/UFF optimizer does not expose the actual iteration count
through its Python API. Metadata therefore records
`forcefield_optimization_iterations: null` and
`forcefield_iteration_count_available: false`, together with the status and
iteration limit. No iteration count is guessed.

## Stereochemistry

Normal mode warns about unspecified stereocenters. `--strict-stereochemistry`
and `--ensemble` stop before ORCA if an unspecified stereocenter or E/Z
candidate is detected.

Conformer search and stereoisomer search are separate operations. The program
retains the `@` and `@@` information in the input SMILES and does not enumerate
R/S isomers from unspecified centers. Canonical isomeric SMILES, CIP labels,
and unspecified centers are stored in metadata.

## Reusing existing ORCA calculations

```powershell
molecule2fbx --ensemble --smiles "..." --name Molecule `
  --reuse-calculations output\Molecule_dft_calculations `
  --nprocs 16 --maxcore 1000
```

Reuse validation checks the ORCA input, output, and XYZ; normal termination;
Opt convergence; atom order; canonical isomeric SMILES; method; functional;
basis; charge; and multiplicity. Compatible Opt results are not rerun, and only
missing conformers are calculated. Incomplete or incompatible existing
directories are never overwritten.

The current ensemble implementation uses gas-phase optimization without an
additional dispersion keyword. These settings are recorded in ensemble JSON,
and calculations with different conditions are not silently compared as one
ensemble.

Run only Freq on an existing optimized XYZ with:

```powershell
molecule2fbx --frequency-only path\to\conformer_001.xyz `
  --nprocs 16 --maxcore 1000
```

This generates a Freq input without Opt or `%geom`. Charge, multiplicity,
functional, and basis are recovered from metadata or the matching `.inp` file.
Freq-only parsing does not require an Opt-convergence marker.

## Output

Example ensemble output:

```text
output/
|-- Molecule_dft.fbx
|-- Molecule_dft.metadata.json
|-- Molecule_dft_conf001.fbx
|-- Molecule_dft_conf001.metadata.json
|-- Molecule_dft_ensemble.json
`-- Molecule_dft_calculations/
    |-- conformer_001/
    `-- frequency_additions/
```

Ensemble JSON records force-field exclusions, initial clusters, relative
electronic energies for all converged DFT structures, post-DFT clusters and
duplicates, final representatives, Freq selection, all frequencies,
low-frequency and imaginary modes, thermochemistry, relative Gibbs energies,
reuse status, calculation directories, software settings, and timestamps.

A structure without Freq is never labelled as a confirmed local minimum. The
assessment has three states:

- `local_minimum_candidate`: Freq completed with no mode below -20 cm-1.
- `not_a_confirmed_local_minimum`: one or more imaginary modes were found.
- `not_evaluated`: Freq was not run or could not be parsed.

Low-frequency modes can strongly affect Gibbs energies under the harmonic
approximation. Metadata records whether ORCA used Quasi-RRHO, but this does not
guarantee thermochemical accuracy or chemical correctness.

## Metals and special molecules

Organic-molecule defaults are not applied automatically when a metal is
detected. The calculation requires `--allow-metals` and explicit functional,
basis, charge, and multiplicity settings.

```powershell
molecule2fbx --smiles "..." --method dft --allow-metals `
  --functional B3LYP --basis def2-TZVP --charge 0 --multiplicity 1
```

This guard does not establish that the basis, ECP, oxidation state, or
open-shell state is appropriate. Metal complexes require specialist
electronic-structure assessment.

## Architecture

- `config.py`: shared request and validation model for the CLI and a future GUI.
- `cli.py`: command-line argument parsing.
- `structures.py`: SMILES validation, ETKDG, MMFF94s/UFF, and heavy-atom RMSD.
- `ensemble.py`: energy filtering, clustering, relative energies, and ensemble JSON.
- `pipeline.py`: user-interface-independent workflow integration.
- `quantum/base.py`: backend-independent settings, results, and thermochemistry.
- `quantum/orca.py`: ORCA input generation, execution, and output parsing.
- `quantum/reuse.py`: non-destructive Opt and Freq reuse.
- `frequency.py`: Freq-only workflow for an optimized XYZ.
- `blender_export.py` and `blender_worker.py`: FBX generation.

ORCA was selected because it is practical as an external executable on Windows
and supports DFT geometry optimization, open-shell systems, frequency
calculations, and broad basis/ECP choices without embedding a large quantum
chemistry runtime in the Python environment.

`outputs/1SB-LSD_RR_redo` remains local validation data. Its conf002 is the
current best among the original four structures, not a global minimum. A new
structure becomes the current-best candidate only if a broader ensemble search
finds a lower-energy result under comparable conditions.

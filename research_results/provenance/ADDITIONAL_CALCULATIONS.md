# Additional calculation record

## Missing frequency calculations

Before execution, two frequency-only jobs were selected:

| Molecule | Conformer | Relative electronic energy | Reason |
| --- | ---: | ---: | --- |
| 1Bz-LSD_RR | conf008 | 4.785 kJ/mol | Inside the configured 5 kJ/mol window but omitted by the maximum-three limit |
| 1SB-LSD_RR | conf006 | 1.095 kJ/mol | Low-energy structure lacking thermochemistry and local-minimum assessment |

Conditions are unchanged from the ensemble runs: ORCA 6.1.1,
B3LYP/def2-SVP, charge 0, multiplicity 1, `nprocs 16`, and
`maxcore 1000`. The optimized XYZ coordinates are used directly; Opt is not
rerun. New files are written only below each run's
`conformers/frequency_additions/<conformer>` directory.

Estimated wall time from the completed jobs on this host is about 80 minutes
for 1Bz-LSD_RR and 115 minutes for 1SB-LSD_RR, or about 3--4 hours when run
sequentially. These jobs can change the conditional Gibbs weighting and the
local-minimum status of two already-low-energy structures. They cannot resolve
the incomplete SB conformer search or prove a global minimum.

No additional geometry optimization is authorized by this record.

Two earlier Bz conf008 launch attempts are retained non-destructively below
`conformers/frequency_additions/incomplete_attempts`. One was a deliberately
terminated startup check and one failed in the PowerShell wrapper before a
meaningful ORCA run because CLI progress on stderr was treated as a terminating
error. They are operational history, not final quantum-calculation failures;
the distinction is recorded in `operational_history.json`.

## SB force-field pool screening

Common-scaffold re-clustering reduces the 22 all-heavy-atom force-field
clusters to 9 scaffold clusters. Seven are represented within 0.75 Angstrom
by an optimized DFT structure; the representatives of two remaining clusters
are ETKDG/MMFF pool indices 46 and 109.

Before any additional geometry optimization, run one B3LYP/def2-SVP
single-point calculation on each retained MMFF94s geometry. Compare these
unrelaxed energies against the first-cycle single-point energies already
present in the retained Opt outputs. This is a screening comparison only:
an unoptimized single-point energy is not a conformer free energy and cannot
replace an Opt/Freq result.

- New jobs: 2 single points; no Opt and no Freq.
- Charge/multiplicity: 0/1.
- Resources: nprocs 16, maxcore 1000 MB.
- Estimated wall time after the current Freq queue: 5-20 minutes total.
- Decision impact: an omitted candidate comparable to or below the best
  retained initial-geometry energy justifies one additional Opt; a clearly
  higher result argues against spending hours on it, but does not prove the
  basin is absent.

After the initial environment check, `xTB 6.7.1pre` was found in the installed
ORCA 6.1.1 directory. Therefore the first independent reranking step is updated
to GFN2-xTB tight geometry optimization for all 9 common-scaffold cluster
representatives. This compares relaxed structures at one consistent cheap
electronic-structure level and is preferable to comparing MMFF rankings alone.

- New jobs: 9 GFN2-xTB Opt jobs; these are not DFT Opt jobs.
- Estimated wall time: 10-45 minutes total on this host.
- The original `(6aR,9R)` assignment is validated after every xTB Opt.
- The two B3LYP/def2-SVP single points (pool 46 and 109) are conditional: run
  them only if xTB leaves an omitted basin within 10 kJ/mol or otherwise makes
  it competitive. This avoids unnecessary work.

### Additional DFT Opt decision after screening

Both independent screens retain pool 109 as the strongest omitted candidate:

| pool | GFN2-xTB relative energy | B3LYP/def2-SVP initial-geometry SP relative to retained best initial geometry |
| ---: | ---: | ---: |
| 109 | 5.04 kJ/mol | 4.62 kJ/mol |
| 46 | 6.88 kJ/mol | 6.48 kJ/mol |

Existing DFT Opt jobs on this host required about 24-46 minutes and relaxed
the initial structures by amounts that vary by several kJ/mol. Pool 109 can
therefore still enter the low-energy DFT set. Run exactly one additional
B3LYP/def2-SVP Opt from the retained pool-109 MMFF94s geometry, with charge 0,
multiplicity 1, nprocs 16, and maxcore 1000 MB.

- Estimated wall time: 25-50 minutes.
- Output is written to the comparison directory, not into or over an existing
  ensemble conformer directory.
- Pool 46 is not optimized in this pass.
- Decision impact: if pool 109 converges to a unique structure near the current
  best, the SB search asymmetry materially affects the hypothesis comparison;
  if it converges to an existing minimum or remains clearly higher, additional
  DFT expansion has lower immediate value.
- Freq is not pre-authorized by this decision. If the new structure is unique
  and low enough to affect the Gibbs comparison, its Freq need and estimated
  cost are reported separately before execution.

### Pool-109 Opt outcome and Freq decision

The targeted Opt converged in 23.5 minutes at -1767.603665704513 Eh. It retained
the RR stereochemistry, lies 4.03 kJ/mol above the existing current best, and
is a new common-scaffold structure (closest existing core RMSD 1.141 Angstrom).

Five existing SB structures lie between 0 and 1.10 kJ/mol, so pool 109 is not a
top-three electronic-energy Freq target even though it is inside the 5 kJ/mol
window. A new Freq would be expected to take about 118 minutes on this host.
It is not run in this preliminary pass. The structure is included in the
electronic-energy sensitivity analysis with `frequency_calculated: false` and
`imaginary_modes: null`; it is excluded from Gibbs weighting and is not called
a confirmed local minimum.

## Robust charge/frontier-property calculations

The retained Opt outputs do not contain MBIS or CHELPG charges and do not
contain atom-resolved frontier-orbital populations. After the missing Freq
jobs, property-only single points may be run on the low-electronic-energy
subsets (within 5 kJ/mol, capped at five structures per molecule) using the
same B3LYP/def2-SVP model and retained optimized coordinates.

- Requested properties: MBIS, CHELPG, and ORCA frontier-MO population output.
- Geometry changes: none; no Opt and no Freq.
- Estimated wall time: 30-90 minutes total for the matched low-energy subset.
- Decision impact: this can make the preliminary electronic comparison more
  robust than Mulliken/Loewdin alone. It still does not calculate a hydrolysis
  barrier or separate substituent effects from solvation/enzyme effects.

### Property-only outcome

The property-only jobs completed for four Bz and six SB structures, including
the targeted pool-109 follow-up. All jobs reused the retained optimized XYZ
coordinates and did not perform Opt or Freq. The resulting MBIS/CHELPG values
cover about 93% of the common-scaffold-unique electronic weight in each
molecule; uncovered structures are retained as null rather than imputed.

The preliminary analysis detects a consistent lower carbonyl-C MBIS charge and
lower benzoyl-center LUMO population in SB. CHELPG, bond-length, Mayer-order,
dipole, and carbonyl-stretch shifts are smaller relative to conformer spread.
This is recorded as an electronic-structure signal at the present gas-phase
B3LYP/def2-SVP level, not as a hydrolysis-rate prediction.

### Remaining pool-46 decision

Pool 46 is now the only common-scaffold cluster in the regenerated SB
force-field window without a DFT-optimized representative. It remained 6.88
kJ/mol above the xTB-screen minimum and 6.48 kJ/mol above the best retained
initial geometry in the unrelaxed B3LYP/def2-SVP single-point comparison.
No further Opt was started. If complete coverage is later prioritized, one Opt
is estimated at 0.4--0.8 hours; Freq (about 2 hours) is warranted only if that
Opt reaches the final low-energy/Freq selection region.

# 1SB-LSD_RR Ensemble Night Run

## Calculation status

SUCCESS

## Summary

```text
ETKDG candidates:       200
MMFF/UFF valid:         200
Energy-filtered:        40
All-heavy RMSD clusters:       22
DFT candidates:         10
DFT completed:          10
Unique DFT structures:  10
Freq selected this run: 3
Freq completed for selection: 3
Freq available total:   7
Freq new supplemental:  1
Freq preexisting local:  2
Freq reused external:    4
Freq preexisting total:  6
```

## Lowest-energy structure found

```text
Conformer:       conf002
Energy:          -1767.605199644220 Eh
Relative Delta E:    0.000 kJ/mol
Gibbs:           -1767.05772569
Relative Delta G:0.0
Imaginary modes: 0
```

## Morning checklist

- [x] Job status is SUCCESS
- [x] Explicit stereochemistry was checked
- [x] DFT Opt completed: 10
- [x] Lowest electronic-energy structure found: conf002
- [x] Lowest structure has zero significant imaginary modes

## Interpretation

This is the lowest-energy structure found among the conformers evaluated in this run. 
It is not a proven global minimum and is not an experimentally determined structure. 
Zero imaginary modes only supports local-minimum character for that calculated structure.

## Low-frequency caution

Freq-completed structures contain 47 modes below the configured 
50 cm^-1 reporting threshold. Their Gibbs energies can be sensitive to the 
harmonic/quasi-RRHO treatment and should be interpreted cautiously.

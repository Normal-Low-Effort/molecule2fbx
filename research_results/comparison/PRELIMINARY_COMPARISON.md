# 1Bz-LSD_RR vs 1SB-LSD_RR preliminary comparison

Generated from retained B3LYP/def2-SVP gas-phase structures. These are computed candidate structures, not experimental structures.

## Scope

- Bz: 9 final DFT representatives; original Freq selection 3/3 completed, 1 supplemental, 4 total available.
- SB sensitivity set: 11 DFT representatives including one targeted follow-up; original Freq selection 3/3 completed, 1 supplemental, 7 total available.
- SB force-field pool review: 9 common-scaffold clusters within the force-field window; 1 clusters lacked a selected DFT representative.
- Values in the main table use one lowest-energy representative per common-scaffold cluster, then electronic-energy weighting. This is not claimed to be a complete Boltzmann ensemble or a rotor-degeneracy treatment.
- Gibbs-weighted values in analysis.json are conditional on structures with Freq unless every final representative has Freq.

## Electronically weighted descriptor comparison

| Descriptor | 1Bz-LSD_RR | 1SB-LSD_RR | SB - Bz |
| --- | ---: | ---: | ---: |
| C=O length (A) | 1.21559 | 1.21568 | 0.00010 |
| Amide C-N length (A) | 1.40226 | 1.40263 | 0.00037 |
| Benzoyl-carbonyl torsion (deg) | 40.141 | 40.536 | 0.395 |
| BD-cone accessible fraction | 0.1528 | 0.1528 | 0.0000 |
| Best BD-cone clearance (A) | 0.2860 | 0.2851 | -0.0009 |
| MBIS carbonyl-C charge | 0.65068 | 0.64697 | -0.00371 |
| MBIS carbonyl-O charge | -0.51157 | -0.51134 | 0.00023 |
| CHELPG carbonyl-C charge | 0.53219 | 0.52892 | -0.00327 |
| CHELPG carbonyl-O charge | -0.48241 | -0.48465 | -0.00224 |
| Mayer C=O bond order | 2.08854 | 2.08821 | -0.00032 |
| Mayer C-N bond order | 1.03927 | 1.03792 | -0.00135 |
| Dipole magnitude (D) | 5.2785 | 5.3726 | 0.0941 |
| LUMO Loewdin population at benzoyl center | 0.46325 | 0.44171 | -0.02154 |
| Carbonyl stretch (cm-1) | 1777.46 | 1777.01 | -0.45 |

## Evidence supporting hypothesis 1 (steric relief from the longer aryl-Si bond)

- SB retains essentially the same sampled directional carbonyl access as Bz. The terminal para-TMS atoms remain remote from the carbonyl center. This is consistent with absence of direct static occlusion at the carbonyl, but it does not prove that the longer Si-C bond is the cause.

## Evidence against hypothesis 1

- SB does not show an increased best clearance; the mean difference is much smaller than the conformer spread. The fixed-cone accessible fraction is also quantized and identical for every retained structure, so it is not sufficiently discriminating on its own.
- The probe omits enzyme-pocket reorganization, water, and an actual nucleophile trajectory. Comparable static access therefore cannot be converted into a hydrolysis-rate claim.

## Evidence supporting hypothesis 2 (electronic effect)

- Property-only single points provide MBIS/CHELPG charges and atom-resolved frontier populations for the low-energy subset. SB shows a coherent lower carbonyl-C MBIS charge and lower benzoyl-center LUMO population; smaller shifts in C=O length, Mayer C=O order, and carbonyl stretch point in a compatible direction.

## Evidence against hypothesis 2

- CHELPG, bond-length, Mayer-order, dipole, and stretch differences are small relative to conformer spread or expected model sensitivity. The property subset covers about 93% of the electronic weight, not every structure. The present calculation therefore detects an electronic-structure signal but does not establish its magnitude in solution or its kinetic consequence.

## Not decidable from the current data

- Two molecules change steric and electronic factors simultaneously; causality cannot be separated without a control or reaction-energy calculation.
- A static reactant structure does not determine hydrolysis rate or an enzymatic activation barrier.
- Differences comparable to conformer spread or method uncertainty must be reported as not detected at this level, not as absent.
- Re/Si face labels are intentionally not assigned by the geometric probe implementation; signed faces are reported instead.

## Method limitations

- B3LYP/def2-SVP, gas phase, no explicit dispersion correction and neutral singlet only.
- Common-scaffold RMSD uses fixed input atom order. It does not permute symmetric atoms; the entire terminal aryl-Si(CH3)3 branch is excluded to prevent methyl rotation from defining scaffold clusters.
- Low-frequency modes make sub-kJ/mol Gibbs rankings sensitive to the thermochemistry treatment.

## Highest-value next calculation

No further long calculation is required for this preliminary comparison. Pool 109 was already optimized as the single targeted follow-up and remained 4.03 kJ/mol above SB conf002; pool 46 is the only uncovered common-scaffold cluster and was 6.48 kJ/mol above the best retained initial geometry in the DFT single-point screen. The highest-value next checks are a finer continuous steric trajectory descriptor and matched higher-basis/solvent property single points. Pool 46 Opt is optional if closing the residual search asymmetry becomes more important than model validation. Reaction TS and substituted controls remain later-stage work.

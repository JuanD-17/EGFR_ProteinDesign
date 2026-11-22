# Step 12 - Rotamer feasibility inside the re-predicted variant

Generated 2026-11-22T14:58:55+00:00 by `scripts/12_rotamer_in_variant.py`.

Step 11 showed the engineered histidines are buried without a counter-charge and that their pKa falls. This asks whether the contact was reachable at all in the predicted variant, or whether AlphaFold2 simply packed the side chain away from it.

The probe is step 07's, unchanged: 200 orientations on spheres of 2.53 and 3.67 A around CB, clashes discarded, a contact counted between 2.6 and 4.0 A. Only the structure differs: the variant complex as predicted, with the histidine's own side chain removed and the target left in place.

## Result

**12 of 16** histidines can reach their intended handle by some rotamer inside the predicted complex.

| Design | Variant | His | Handle | Best possible | AF2 built | Reaches |
|---|---|---:|---|---:|---:|---|
| l68 | H7 M10H H11 | 10 | H370 | 2.91 | 7.81 | yes |
| l68 | M10H H11 | 10 | H370 | 2.91 | 7.81 | yes |
| l74 | M25H | 25 | D460 | 2.99 | 5.01 | yes |
| l66 | Q23H | 23 | H433 | 3.04 | 8.62 | yes |
| l68 | H7 H11 | 7 | D368 | 3.05 | 2.99 | yes |
| l55 | M55H | 55 | D460 | 3.08 | 6.14 | yes |
| l68 | H7 M10H H11 | 7 | D368 | 3.16 | 3.18 | yes |
| l74 | W63H | 63 | H433 | 3.27 | 10.47 | yes |
| l74 | W24H M25H | 24 | H433 | 3.27 | 11.15 | yes |
| l73 | F59H | 59 | H433 | 3.38 | 11.19 | yes |
| l74 | W24H | 24 | H433 | 3.47 | 10.5 | yes |
| l74 | W24H M25H | 25 | D460 | 3.58 | 5.54 | yes |
| l68 | H7 H11 | 11 | H433 | 4.41 | 9.84 | no |
| l68 | H7 M10H H11 | 11 | H433 | 4.89 | 9.86 | no |
| l68 | M10H H11 | 11 | H433 | 4.89 | 9.86 | no |
| l63 | R49H | 49 | H370 | 6.06 | 10.67 | no |

## Reading this

*Best possible* is the closest approach any unclashed rotamer could make. *AF2 built* is the distance in the predicted structure. Where the two differ widely the predictor had a better option and did not take it; where both are large, no rotamer was available.

## Comparison with step 07

Step 07 ran this probe on the parent backbone and 18 of 23 substitutions passed. Step 12 runs it on the variant complex and 12 of 16 pass.

The gap between those two numbers is the finding. A position that can accept a reaching histidine before the substitution is not the same as a position that still can once the sequence has changed and the complex has been predicted again.

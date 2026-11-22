# Step 12 - Rotamer feasibility inside the re-predicted variant

Generated 2026-11-22T16:34:37+00:00 by `scripts/12_rotamer_in_variant.py`.

Step 11 showed the engineered histidines are buried without a counter-charge and that their pKa falls. This asks whether the contact was reachable at all in the predicted variant, or whether AlphaFold2 simply packed the side chain away from it.

The probe is step 07's, unchanged: 200 orientations on spheres of 2.53 and 3.67 A around CB, clashes discarded, a contact counted between 2.6 and 4.0 A. Only the structure differs: the variant complex as predicted, with the histidine's own side chain removed and the target left in place.

## Result

**12 of 16** histidines can reach their intended handle by some rotamer inside the predicted complex.

| Design | Variant | His | Handle | Best possible | AF2 built | Reaches |
|---|---|---:|---|---:|---:|---|
| l68 | M10H H11 | 10 | H370 | 2.91 | 8.05 | yes |
| l74 | W24H M25H | 24 | H433 | 2.93 | 10.2 | yes |
| l68 | H7 M10H H11 | 7 | D368 | 2.98 | 3.11 | yes |
| l66 | Q23H | 23 | H433 | 3.06 | 8.53 | yes |
| l68 | H7 M10H H11 | 11 | H433 | 3.08 | 10.77 | yes |
| l68 | H7 M10H H11 | 10 | H370 | 3.11 | 7.27 | yes |
| l68 | H7 H11 | 7 | D368 | 3.14 | 3.99 | yes |
| l63 | R49H | 49 | H370 | 3.16 | 7.21 | yes |
| l68 | H7 H11 | 11 | H433 | 3.16 | 9.58 | yes |
| l55 | M55H | 55 | D460 | 3.18 | 5.03 | yes |
| l74 | W63H | 63 | H433 | 3.31 | 11.7 | yes |
| l73 | F59H | 59 | H433 | 3.34 | 11.34 | yes |
| l68 | M10H H11 | 11 | H370 | 5.08 | 9.97 | no |
| l74 | W24H M25H | 25 | H433 | 11.64 | 9.24 | no |
| l74 | W24H | 24 | H433 | 13.21 | 4.56 | no |
| l74 | M25H | 25 | D460 | 16.48 | 8.51 | no |

## Reading this

*Best possible* is the closest approach any unclashed rotamer could make. *AF2 built* is the distance in the predicted structure. Where the two differ widely the predictor had a better option and did not take it; where both are large, no rotamer was available.

## Comparison with step 07

Step 07 ran this probe on the parent backbone and 18 of 23 substitutions passed. Step 12 runs it on the variant complex and 12 of 16 pass.

The gap between those two numbers is the finding. A position that can accept a reaching histidine before the substitution is not the same as a position that still can once the sequence has changed and the complex has been predicted again.

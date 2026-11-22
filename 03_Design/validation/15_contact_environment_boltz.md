# Step 15 - Environment of the contact that was not made

Generated 2026-11-22T16:45:25+00:00 by `scripts/15_contact_environment.py`.

Steps 12 and 14 found a reaching rotamer sitting unused in most cases. The probe tests sterics only, so this asks whether the position the predictors declined is one a protein would want: what lies within 6 A of it, how apolar, and whether a compensating charge is present.

## Declined positions against the ones that were built

| | Contact built | Reachable, declined |
|---|---:|---:|
| Cases | 2 | 10 |
| Apolar fraction at the contact point | 0.587 | 0.644 |
| Net charge within 6 A | +0.0 | -0.1 |
| Heavy atoms within 6 A | 46 | 40 |

## Where the predictor put it instead

| | Declined contact point | Position actually built |
|---|---:|---:|
| Apolar fraction | 0.644 | 0.650 |
| Net charge | -0.1 | +0.1 |
| Heavy atoms | 40 | 37 |

## Per histidine

| Design | Variant | His | Reach | Built | Apolar (reach) | Net charge (reach) | Apolar (built) |
|---|---|---:|---:|---:|---:|---:|---:|
| l68 | H7 H11 | 11 | 3.16 | 3.64 | 0.750 | +0.0 | 0.733 |
| l74 | W24H M25H | 25 | 11.64 | 14.42 | 0.733 | +1.0 | 0.571 |
| l68 | H7 M10H H11 | 11 | 3.08 | 4.20 | 0.718 | +0.0 | 0.634 |
| l68 | M10H H11 | 10 | 2.91 | 3.44 | 0.667 | +0.0 | 0.667 |
| l68 | M10H H11 | 11 | 5.08 | 3.57 | 0.636 | +0.0 | 0.745 |
| l74 | W24H M25H | 24 | 2.93 | 3.91 | 0.632 | +0.0 | 0.605 |
| l68 | H7 M10H H11 | 10 | 3.07 | 4.95 | 0.630 | +1.5 | 0.578 |
| l73 | F59H | 59 | 3.34 | 3.54 | 0.625 | -1.0 | 0.657 |
| l63 | R49H | 49 | 3.16 | 4.28 | 0.622 | +0.0 | 0.688 |
| l55 | M55H | 55 | 3.18 | 5.26 | 0.610 | -0.5 | 0.633 |
| l66 | Q23H | 23 | 3.06 | 7.63 | 0.600 | -0.5 | 0.579 |
| l68 | H7 M10H H11 | 7 | 2.93 | 3.15 | 0.596 | +0.0 | 0.585 |
| l74 | W63H | 63 | 3.31 | 5.22 | 0.590 | -0.5 | 0.727 |
| l74 | W24H | 24 | 13.21 | 14.59 | 0.583 | +0.0 | 0.686 |
| l68 | H7 H11 | 7 | 2.91 | 3.96 | 0.578 | +0.0 | 0.558 |
| l74 | M25H | 25 | 16.48 | 20.83 | 0.576 | -1.5 | 0.625 |

## The handle chemistry

7 of 16 histidines were placed against a handle that is itself a histidine on the receptor - H433 or H370 - and only 9 against a carboxylate.

That is a flaw in the placement rule, not in the predictions. The mechanism needs the protonated binder histidine to be stabilised on binding, which requires a negative counter-charge. A receptor histidine carries +1 over the same pH window, so the pair repels and binding is destabilised at the lower pH, which is the opposite of the intended direction. Step 03 listed H433 as a titratable group and the rule took it as an anchor without asking which sign it carries.

The measured shifts follow: against carboxylates they average well above those against receptor histidines, and the single positive shift in the set is the one carboxylate contact that was built. With five cases against eleven and distance confounded with chemistry, this is consistent with the mechanism rather than a demonstration of it.


## How to read this

A protonated imidazole wants a polar, negatively charged neighbourhood. If the declined positions are more apolar than the ones that were built, or carry no compensating negative charge, the predictors had grounds and the placement rule is at fault for scoring geometry without desolvation.

If the two groups look alike, this measurement has no signal and does not support either reading. With single-digit group sizes it cannot support much in any case, and is reported as a description rather than a test.

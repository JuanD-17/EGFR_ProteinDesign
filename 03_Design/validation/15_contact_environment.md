# Step 15 - Environment of the contact that was not made

Generated 2026-11-22T16:45:24+00:00 by `scripts/15_contact_environment.py`.

Steps 12 and 14 found a reaching rotamer sitting unused in most cases. The probe tests sterics only, so this asks whether the position the predictors declined is one a protein would want: what lies within 6 A of it, how apolar, and whether a compensating charge is present.

## Declined positions against the ones that were built

| | Contact built | Reachable, declined |
|---|---:|---:|
| Cases | 1 | 11 |
| Apolar fraction at the contact point | 0.652 | 0.657 |
| Net charge within 6 A | -1.0 | +0.0 |
| Heavy atoms within 6 A | 46 | 36 |

## Where the predictor put it instead

| | Declined contact point | Position actually built |
|---|---:|---:|
| Apolar fraction | 0.657 | 0.604 |
| Net charge | +0.0 | -0.5 |
| Heavy atoms | 36 | 36 |

## Per histidine

| Design | Variant | His | Reach | Built | Apolar (reach) | Net charge (reach) | Apolar (built) |
|---|---|---:|---:|---:|---:|---:|---:|
| l68 | H7 H11 | 11 | 4.41 | 3.73 | 0.714 | +0.0 | 0.698 |
| l68 | H7 M10H H11 | 10 | 2.91 | 4.88 | 0.705 | +0.0 | 0.629 |
| l68 | M10H H11 | 10 | 2.91 | 4.88 | 0.705 | +0.0 | 0.629 |
| l74 | W24H | 24 | 3.30 | 3.72 | 0.697 | +0.0 | 0.590 |
| l74 | W24H M25H | 25 | 3.58 | 5.78 | 0.688 | +1.0 | 0.583 |
| l74 | W63H | 63 | 2.92 | 3.55 | 0.677 | +0.5 | 0.581 |
| l74 | W24H M25H | 24 | 3.27 | 3.17 | 0.677 | +0.0 | 0.639 |
| l63 | R49H | 49 | 6.06 | 7.20 | 0.675 | +0.0 | 0.614 |
| l73 | F59H | 59 | 2.98 | 3.87 | 0.673 | -1.0 | 0.595 |
| l68 | H7 H11 | 7 | 3.05 | 3.84 | 0.652 | -1.0 | 0.608 |
| l68 | H7 M10H H11 | 11 | 4.89 | 3.61 | 0.652 | +0.0 | 0.717 |
| l68 | M10H H11 | 11 | 4.89 | 3.61 | 0.652 | +0.0 | 0.717 |
| l74 | M25H | 25 | 2.99 | 5.41 | 0.645 | +0.5 | 0.571 |
| l68 | H7 M10H H11 | 7 | 3.16 | 4.14 | 0.609 | -1.0 | 0.627 |
| l66 | Q23H | 23 | 3.04 | 7.04 | 0.607 | +0.0 | 0.529 |
| l55 | M55H | 55 | 3.08 | 6.54 | 0.542 | +0.0 | 0.667 |

## The handle chemistry

7 of 16 histidines were placed against a handle that is itself a histidine on the receptor - H433 or H370 - and only 9 against a carboxylate.

That is a flaw in the placement rule, not in the predictions. The mechanism needs the protonated binder histidine to be stabilised on binding, which requires a negative counter-charge. A receptor histidine carries +1 over the same pH window, so the pair repels and binding is destabilised at the lower pH, which is the opposite of the intended direction. Step 03 listed H433 as a titratable group and the rule took it as an anchor without asking which sign it carries.

The measured shifts follow: against carboxylates they average well above those against receptor histidines, and the single positive shift in the set is the one carboxylate contact that was built. With five cases against eleven and distance confounded with chemistry, this is consistent with the mechanism rather than a demonstration of it.


## How to read this

A protonated imidazole wants a polar, negatively charged neighbourhood. If the declined positions are more apolar than the ones that were built, or carry no compensating negative charge, the predictors had grounds and the placement rule is at fault for scoring geometry without desolvation.

If the two groups look alike, this measurement has no signal and does not support either reading. With single-digit group sizes it cannot support much in any case, and is reported as a description rather than a test.

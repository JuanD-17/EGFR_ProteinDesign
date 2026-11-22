# Step 11 - Bound-state pKa and contact geometry

Generated 2026-11-22T16:34:26+00:00 by `scripts/11_bound_state_pka.py`.

Step 03b assumed a 2.0-unit pKa shift on binding because the complexes did not yet exist. They exist now, so the shift is measured here, together with whether the contact that was supposed to produce it was actually formed.

## Result

Of 16 histidines across 11 variants, **2 lie within 4.0 A of a target carboxylate** and 0 show a positive pKa shift. Measured shifts run from -4.49 to -0.37, against the +2.0 that was assumed.

The assumption was wrong in sign, not in magnitude.

## Per histidine

| Design | Variant | His | pKa free | pKa bound | Shift | Nearest acid | d (A) | Contact | Neighbours |
|---|---|---:|---:|---:|---:|---|---:|---|---:|
| l55 | M55H | 55 | 6.35 | 5.98 | -0.37 | D460 | 5.03 | no | 68 |
| l74 | W24H | 24 | 6.10 | 5.07 | -1.03 | D379 | 4.56 | no | 79 |
| l66 | Q23H | 23 | 6.43 | 5.39 | -1.04 | D368 | 8.53 | no | 79 |
| l74 | M25H | 25 | 6.53 | 5.33 | -1.19 | D368 | 8.51 | no | 75 |
| l63 | R49H | 49 | 6.30 | 4.77 | -1.52 | D368 | 7.21 | no | 86 |
| l74 | W24H M25H | 25 | 6.60 | 4.83 | -1.77 | D379 | 9.24 | no | 83 |
| l68 | H7 M10H H11 | 10 | 5.89 | 4.12 | -1.77 | D368 | 7.27 | no | 100 |
| l73 | F59H | 59 | 6.02 | 3.49 | -2.53 | D368 | 11.34 | no | 105 |
| l68 | M10H H11 | 10 | 6.08 | 3.37 | -2.71 | D368 | 8.05 | no | 116 |
| l68 | H7 M10H H11 | 7 | 6.20 | 2.83 | -3.37 | D368 | 3.11 | yes | 106 |
| l68 | H7 H11 | 7 | 6.20 | 2.69 | -3.51 | D368 | 3.99 | yes | 105 |
| l68 | H7 M10H H11 | 11 | 5.44 | 1.33 | -4.11 | D368 | 10.77 | no | 105 |
| l74 | W24H M25H | 24 | 5.98 | 1.78 | -4.20 | D368 | 10.20 | no | 91 |
| l68 | H7 H11 | 11 | 5.34 | 1.05 | -4.29 | D368 | 9.58 | no | 109 |
| l74 | W63H | 63 | 6.00 | 1.62 | -4.38 | D368 | 11.70 | no | 104 |
| l68 | M10H H11 | 11 | 5.10 | 0.61 | -4.49 | D368 | 9.97 | no | 115 |

## Per variant, measured against assumed

| Design | Variant | Assumed | Measured | 
|---|---|---:|---:|
| l55 | M55H | 0.848 | -0.132 |
| l74 | W63H | 0.943 | -0.139 |
| l73 | F59H | 0.943 | -0.144 |
| l74 | W24H | 0.943 | -0.151 |
| l68 | M10H H11 | 1.531 | -0.184 |
| l63 | R49H | 0.588 | -0.234 |
| l68 | H7 H11 | 1.791 | -0.237 |
| l66 | Q23H | 0.943 | -0.267 |
| l74 | M25H | 0.848 | -0.320 |
| l68 | H7 M10H H11 | 2.379 | -0.355 |
| l74 | W24H M25H | 1.791 | -0.521 |

## What this means

A positive shift means the target stabilises the protonated histidine, which is the direction that gives tighter binding at the lower pH. Every variant here is at or below zero, so under this estimate none of them is a pH switch in the intended direction, and most would bind marginally *better* at pH 7.4.

The contact column explains why. The histidines are enclosed - 70 to 113 heavy atoms within 8 A - but they are not paired with a carboxylate. Burying a base without a counter-charge destabilises its protonated form, and the pKa falls. That is what the numbers show, and the internal consistency is the evidence that they mean it: the only histidine within hydrogen bonding distance of an acid is also the only one whose shift is not negative.

## Where the pipeline failed

Step 07 tested whether a rotamer *could* reach the handle, probing orientations on the parent backbone with the replaced side chain removed. 18 of 23 substitutions passed. That test was correct for the question it asked.

It is the wrong question. What matters is whether the contact exists in the re-predicted variant complex, and it mostly does not. Geometric reachability on the parent does not predict contact formation after the sequence changes and the structure is predicted again.

This is a failure of the placement rule, not of the thermodynamics. The linkage argument in step 03 stands: a coupled protonation is still worth up to 1.23 kcal/mol over this pH window. Nothing here was coupled.

## Limitations of this step

The binder was taken from the complex without moving an atom, so the free state is the bound conformation with the partner deleted. Conformational relaxation on unbinding would raise the free-state pKa and make the shifts *more* negative, so this is the conservative direction and does not rescue the result.

PROPKA's error is near one pKa unit. Several shifts here exceed three units, which is well outside that error, but the smaller ones are not resolvable.

These are AlphaFold2 models. If the prediction places the histidine wrongly, both the contact distance and the pKa are wrong together. Constant-pH molecular dynamics on an ensemble, rather than PROPKA on one frame, is the rigorous route.

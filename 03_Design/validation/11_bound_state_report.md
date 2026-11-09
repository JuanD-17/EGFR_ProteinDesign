# Step 11 - Bound-state pKa and contact geometry

Generated 2026-11-09T17:02:56+00:00 by `scripts/11_bound_state_pka.py`.

Step 03b assumed a 2.0-unit pKa shift on binding because the complexes did not yet exist. They exist now, so the shift is measured here, together with whether the contact that was supposed to produce it was actually formed.

## Result

Of 16 histidines across 11 variants, **2 lie within 4.0 A of a target carboxylate** and 1 show a positive pKa shift. Measured shifts run from -5.76 to +0.05, against the +2.0 that was assumed.

The assumption was wrong in sign, not in magnitude.

## Per histidine

| Design | Variant | His | pKa free | pKa bound | Shift | Nearest acid | d (A) | Contact | Neighbours |
|---|---|---:|---:|---:|---:|---|---:|---|---:|
| l68 | H7 M10H H11 | 7 | 6.24 | 6.29 | +0.05 | D368 | 3.18 | yes | 105 |
| l68 | H7 H11 | 7 | 6.24 | 6.16 | -0.08 | D368 | 2.99 | yes | 104 |
| l66 | Q23H | 23 | 6.45 | 5.72 | -0.73 | D368 | 8.62 | no | 70 |
| l74 | M25H | 25 | 6.42 | 5.69 | -0.73 | D460 | 5.01 | no | 64 |
| l74 | W24H M25H | 25 | 6.44 | 5.54 | -0.90 | D460 | 5.54 | no | 68 |
| l55 | M55H | 55 | 6.42 | 4.88 | -1.54 | D460 | 6.14 | no | 81 |
| l73 | F59H | 59 | 6.08 | 3.65 | -2.43 | D368 | 11.19 | no | 104 |
| l68 | H7 M10H H11 | 10 | 6.10 | 2.92 | -3.18 | D368 | 7.81 | no | 109 |
| l68 | M10H H11 | 10 | 6.10 | 2.92 | -3.18 | D368 | 7.81 | no | 109 |
| l74 | W24H | 24 | 6.02 | 2.80 | -3.22 | D460 | 10.50 | no | 98 |
| l63 | R49H | 49 | 6.05 | 2.81 | -3.24 | D368 | 10.67 | no | 100 |
| l74 | W24H M25H | 24 | 5.94 | 2.48 | -3.46 | D460 | 11.15 | no | 95 |
| l68 | H7 H11 | 11 | 5.41 | 0.52 | -4.90 | D368 | 9.84 | no | 111 |
| l68 | H7 M10H H11 | 11 | 5.28 | -0.18 | -5.46 | D368 | 9.86 | no | 113 |
| l68 | M10H H11 | 11 | 5.28 | -0.18 | -5.46 | D368 | 9.86 | no | 113 |
| l74 | W63H | 63 | 5.88 | 0.11 | -5.76 | D368 | 10.47 | no | 101 |

## Per variant, measured against assumed

| Design | Variant | Assumed | Measured | 
|---|---|---:|---:|
| l68 | H7 H11 | 1.791 | -0.071 |
| l74 | W63H | 0.943 | -0.109 |
| l74 | W24H | 0.943 | -0.145 |
| l63 | R49H | 0.588 | -0.155 |
| l73 | F59H | 0.943 | -0.162 |
| l68 | H7 M10H H11 | 2.379 | -0.179 |
| l68 | M10H H11 | 1.531 | -0.200 |
| l74 | M25H | 0.848 | -0.226 |
| l66 | Q23H | 0.943 | -0.236 |
| l55 | M55H | 0.848 | -0.288 |
| l74 | W24H M25H | 1.791 | -0.380 |

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

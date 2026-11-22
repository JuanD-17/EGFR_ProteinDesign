# Step 14 - AlphaFold2 against Boltz-2

Generated 2026-11-22T16:35:33+00:00 by `scripts/14_compare_predictors.py`.

## Agreement

| | |
|---|---:|
| Histidines compared | 16 |
| Both place it beyond 4 Å of the handle | **14** |
| Both place it in contact | 2 |
| Disagree | 0 |
| Both give a pKa shift at or below zero | **15** |
| A rotamer reached, and neither model used it | **10** |

Boltz-2 interface confidence across 19 complexes: median 0.914, range 0.435 to 0.950. The complexes are not being rejected; the models are confident and still do not build the contacts.

## Per histidine

| Design | Variant | His | Handle | AF2 d | Boltz d | AF2 shift | Boltz shift |
|---|---|---:|---|---:|---:|---:|---:|
| l68 | H7 H11 | 7 | D368 | 2.99 | 3.99 | -0.08 | -3.51 |
| l68 | H7 M10H H11 | 7 | D368 | 3.18 | 3.11 | +0.05 | -3.37 |
| l74 | M25H | 25 | D460 | 5.01 | 8.51 | -0.73 | -1.19 |
| l74 | W24H M25H | 25 | D460 | 5.54 | 9.24 | -0.90 | -1.77 |
| l55 | M55H | 55 | D460 | 6.14 | 5.03 | -1.54 | -0.37 |
| l68 | H7 M10H H11 | 10 | H370 | 7.81 | 7.27 | -3.18 | -1.77 |
| l68 | M10H H11 | 10 | H370 | 7.81 | 8.05 | -3.18 | -2.71 |
| l66 | Q23H | 23 | H433 | 8.62 | 8.53 | -0.73 | -1.04 |
| l68 | H7 H11 | 11 | H433 | 9.84 | 9.58 | -4.90 | -4.29 |
| l68 | H7 M10H H11 | 11 | H433 | 9.86 | 10.77 | -5.46 | -4.11 |
| l68 | M10H H11 | 11 | H433 | 9.86 | 9.97 | -5.46 | -4.49 |
| l74 | W63H | 63 | H433 | 10.47 | 11.70 | -5.76 | -4.38 |
| l74 | W24H | 24 | H433 | 10.50 | 4.56 | -3.22 | -1.03 |
| l63 | R49H | 49 | H370 | 10.67 | 7.21 | -3.24 | -1.52 |
| l74 | W24H M25H | 24 | H433 | 11.15 | 10.20 | -3.46 | -4.20 |
| l73 | F59H | 59 | H433 | 11.19 | 11.34 | -2.43 | -2.53 |

## What this settles, and what it does not

Two predictors of different architectures, given the same sequences under the same conditions, place the engineered histidines away from their handles in 14 of 16 cases and agree on the direction of the resulting pKa shift in 15 of 16. The behaviour is not an artefact of one implementation.

It does not establish that the contacts cannot form. Both models are trained on the Protein Data Bank and inherit its biases together, and neither is an experiment. What the agreement removes is the hypothesis that a single predictor packed the interface carelessly.

The practical consequence stands either way. Validating a designed polar contact by re-prediction, with either model, will report the contact absent - and any electrostatic quantity computed on that structure will inherit the error.

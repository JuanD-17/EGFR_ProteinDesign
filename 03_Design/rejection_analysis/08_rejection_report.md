# Step 08 - Where the campaign loses candidates

Generated 2026-09-30T15:11:25+00:00 by `scripts/08_analyze_rejections.py` over 2 run directories.

## Per run

| Run | Started | Successful | Accepted | Success % | Accept % | Median s |
|-----|--------:|-----------:|---------:|----------:|---------:|---------:|
| 03_Design/run1/gpu0 | 15 | 12 | 3 | 80.0 | 20.0 | 1800 |
| 03_Design/run1/gpu1 | 14 | 10 | 2 | 71.4 | 14.29 | 2029 |

Overall, 5 of 29 trajectories yielded an accepted design, 17.2 percent.

## What rejects them

273 of 322 rejections (85%) are on AlphaFold2's confidence in the re-predicted complex rather than on composition or geometry.

| Filter | Category | Rejections | % |
|--------|----------|-----------:|--:|
| i_pAE | confidence | 108 | 33.5 |
| i_pTM | confidence | 89 | 27.6 |
| pLDDT | confidence | 65 | 20.2 |
| n_InterfaceUnsatHbonds | other | 21 | 6.5 |
| Surface_Hydrophobicity | other | 9 | 2.8 |
| ShapeComplementarity | other | 7 | 2.2 |
| n_InterfaceHbonds | other | 5 | 1.6 |
| Binder_RMSD | other | 5 | 1.6 |
| Binder_pLDDT | confidence | 4 | 1.2 |
| Trajectory_final_pLDDT | confidence | 3 | 0.9 |
| Trajectory_Clashes | other | 2 | 0.6 |
| Trajectory_one-hot_pLDDT | confidence | 2 | 0.6 |
| Trajectory_logits_pLDDT | confidence | 1 | 0.3 |
| Trajectory_softmax_pLDDT | confidence | 1 | 0.3 |

## Reading this

The split between confidence filters and the rest is what decides whether the acceptance rate can be improved cheaply.

Rejections concentrated on peripheral criteria, such as surface hydrophobicity or unsatisfied hydrogen bonds, would suggest a threshold worth revisiting: those measure properties of a design that is otherwise believed to bind.

Rejections concentrated on i_pAE, i_pTM and pLDDT do not. Those are the model's confidence that the complex forms at all, and relaxing them means accepting designs AlphaFold2 does not believe in, which moves the failure from the filter to the assay. The honest response to a low rate of this kind is to accept it and run more trajectories.

## Between runs

A rejection profile that holds across runs indicates the acceptance rate is a property of the target and protocol rather than of chance in any one campaign.

| Filter | gpu0 | gpu1 |
|--------|---:|---:|
| i_pAE | 62 | 46 |
| i_pTM | 56 | 33 |
| pLDDT | 28 | 37 |
| n_InterfaceUnsatHbonds | 13 | 8 |
| Surface_Hydrophobicity | 9 | 0 |
| ShapeComplementarity | 2 | 5 |
| n_InterfaceHbonds | 4 | 1 |
| Binder_RMSD | 3 | 2 |
| Binder_pLDDT | 1 | 3 |
| Trajectory_final_pLDDT | 1 | 2 |

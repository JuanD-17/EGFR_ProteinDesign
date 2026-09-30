# Step 08 - Where the campaign loses candidates

Generated 2026-09-30T23:32:48+00:00 by `scripts/08_analyze_rejections.py` over 4 run directories.

## Per run

| Run | Started | Successful | Accepted | Success % | Accept % | Median s |
|-----|--------:|-----------:|---------:|----------:|---------:|---------:|
| 03_Design/run1/gpu0 | 15 | 12 | 3 | 80.0 | 20.0 | 1800 |
| 03_Design/run1/gpu1 | 14 | 10 | 2 | 71.4 | 14.29 | 2029 |
| 03_Design/run2/gpu0 | 15 | 10 | 1 | 66.7 | 6.67 | 1687 |
| 03_Design/run2/gpu1 | 15 | 12 | 2 | 80.0 | 13.33 | 1645 |

Overall, 8 of 59 trajectories yielded an accepted design, 13.6 percent.

## What rejects them

483 of 628 rejections (77%) are on AlphaFold2's confidence in the re-predicted complex rather than on composition or geometry.

| Filter | Category | Rejections | % |
|--------|----------|-----------:|--:|
| i_pAE | confidence | 196 | 31.2 |
| i_pTM | confidence | 157 | 25.0 |
| pLDDT | confidence | 106 | 16.9 |
| n_InterfaceUnsatHbonds | other | 75 | 11.9 |
| Surface_Hydrophobicity | other | 30 | 4.8 |
| ShapeComplementarity | other | 20 | 3.2 |
| Trajectory_final_pLDDT | confidence | 9 | 1.4 |
| n_InterfaceHbonds | other | 9 | 1.4 |
| Binder_RMSD | other | 8 | 1.3 |
| Trajectory_one-hot_pLDDT | confidence | 7 | 1.1 |
| Binder_pLDDT | confidence | 5 | 0.8 |
| Trajectory_Clashes | other | 3 | 0.5 |
| Trajectory_logits_pLDDT | confidence | 2 | 0.3 |
| Trajectory_softmax_pLDDT | confidence | 1 | 0.2 |

## Reading this

The split between confidence filters and the rest is what decides whether the acceptance rate can be improved cheaply.

Rejections concentrated on peripheral criteria, such as surface hydrophobicity or unsatisfied hydrogen bonds, would suggest a threshold worth revisiting: those measure properties of a design that is otherwise believed to bind.

Rejections concentrated on i_pAE, i_pTM and pLDDT do not. Those are the model's confidence that the complex forms at all, and relaxing them means accepting designs AlphaFold2 does not believe in, which moves the failure from the filter to the assay. The honest response to a low rate of this kind is to accept it and run more trajectories.

## Between runs

A rejection profile that holds across runs indicates the acceptance rate is a property of the target and protocol rather than of chance in any one campaign.

| Filter | run1 | run2 |
|--------|---:|---:|
| i_pAE | 108 | 88 |
| i_pTM | 89 | 68 |
| pLDDT | 65 | 41 |
| n_InterfaceUnsatHbonds | 21 | 54 |
| Surface_Hydrophobicity | 9 | 21 |
| ShapeComplementarity | 7 | 13 |
| Trajectory_final_pLDDT | 3 | 6 |
| n_InterfaceHbonds | 5 | 4 |
| Binder_RMSD | 5 | 3 |
| Trajectory_one-hot_pLDDT | 2 | 5 |

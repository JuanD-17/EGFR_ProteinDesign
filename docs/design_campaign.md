# Design campaign

Binder generation on rented GPU, and where the campaign lost its candidates.

---

### `scripts/08_analyze_rejections.py` — where candidates are lost

A campaign accepting 5 designs has discarded the rest, and where they died
decides whether more compute helps and what to change. Candidates are lost
in three places, and they mean different things: a trajectory that never
converges is cheap, a backbone whose redesigned sequences all fail
re-prediction is expensive, and a scored complex rejected by a threshold is
informative, because it says which property limits the campaign.

```bash
python scripts/08_analyze_rejections.py --runs 03_Design/run1/gpu0 \
                                               03_Design/run1/gpu1
```

## Findings from step 08

**The acceptance rate is 17%, not the 9% first reported.** Counting
trajectories with `grep "Starting trajectory"` double-counts, because
BindCraft writes both `Starting trajectory: <name>` when one begins and
`Starting trajectory took: …` when its hallucination stage ends. The first
run started 29 trajectories, not 58; 22 reached the end of the four-stage
protocol, a 76% completion rate, and 5 were accepted.

This is the error the script exists to prevent, and it is worth recording
that the ad-hoc command produced a figure nearly two-fold wrong while the
parser with a precise pattern did not.

| Run | Started | Completed | Accepted | Median trajectory |
|-----|--------:|----------:|---------:|------------------:|
| gpu0 | 15 | 12 | 3 | 1800 s |
| gpu1 | 14 | 10 | 2 | 2029 s |

**85% of rejections are on AlphaFold2 confidence.**

| Filter | Category | Rejections |
|--------|----------|-----------:|
| i_pAE | confidence | 108 |
| i_pTM | confidence | 89 |
| pLDDT | confidence | 65 |
| n_InterfaceUnsatHbonds | other | 21 |
| Surface_Hydrophobicity | other | 9 |
| ShapeComplementarity | other | 7 |

The split is what decides whether the rate can be improved cheaply.
Rejections concentrated on peripheral criteria would point to a threshold
worth revisiting, since those measure properties of a design otherwise
believed to bind. Rejections concentrated on i_pAE, i_pTM and pLDDT do not:
those are the model's confidence that the complex forms at all, and
relaxing them accepts designs AlphaFold2 does not believe in, which moves
the failure from the filter to the assay. The rate is therefore accepted as
a property of this target and protocol.


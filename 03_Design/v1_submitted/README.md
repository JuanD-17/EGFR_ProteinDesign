# The submitted campaign, frozen

What was entered into the Anthropic × Adaptyv 2026 competition on 3 October
2026, and the analysis that followed it. **Nothing in this directory is
regenerated.** It is the record of what was actually submitted, including
the error it contained.

The live pipeline writes to `03_Design/ph_variants/`,
`03_Design/negative_control/` and `03_Design/validation/`. Those hold the
corrected campaign.

## What is here

| | |
|---|---|
| `ph_variants/` | the 16 variants step 05 proposed, under the original rule |
| `negative_control/` | step 06 against the matched null |
| `validation/predictions/` | AlphaFold2 re-predictions, 19 complexes |
| `validation/boltz_predictions/` | Boltz-2 re-predictions of the same 19 |
| `validation/11_*` | bound-state pKa, both predictors |
| `validation/12_*` | rotamer feasibility inside the variant, both predictors |
| `validation/14_*` | AlphaFold2 against Boltz-2 |
| `validation/15_*` | environment of the declined contacts |

## The error this campaign contains

Step 03b labels each target handle by role. A carboxylate is an *anchor for
an engineered histidine*; H370 and H433 are *receptor-side switches*, a
different mechanism. Step 05 read that table and ignored the column,
placing a histidine against every handle.

Seventeen of twenty-three placements went to H370 or H433. A protonated
binder histidine next to a protonated receptor histidine is two cations in
contact: they repel, and binding is destabilised at the lower pH, which is
the opposite of the objective.

The measured consequences are in `validation/11_*` and are what led to the
correction. They are kept because the campaign with the error is the
evidence that the error existed and what it cost.

## Why it is frozen rather than deleted

Three reasons.

It is what was submitted. A record of an entry that no longer matches the
entry is worse than no record.

The comparison between this campaign and the corrected one is the result.
Same binders, same epitope, same geometry, one chemistry changed — that is
as close to a controlled comparison as this project will get without a
laboratory.

And the failure modes found here are transferable whether or not the
correction works: that re-prediction does not build designed polar
contacts, that two independent predictors agree on declining them, and
that a matched negative control validated a placement rule whose premise
was wrong. The control drew both arms from the same handle list, so the
null shared the error and could not see it.

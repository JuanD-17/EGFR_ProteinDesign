# Step 07 - Negative control

Generated 2026-09-30T14:23:11+00:00 by `scripts/07_negative_control.py`, seed 20260930, 3000 draws per design and substitution count.

## The question

Step 05 reports a coupling estimate for each designed variant. That figure is uninterpretable alone. What matters is whether placing histidines by the mechanistic rule beats placing the same number at random among the positions that were available.

## The null model

The null draws from exactly the pool step 05 could have chosen from: binder positions at the interface, exposed in the unbound binder, and chemically substitutable. Only the targeting is removed.

This is deliberately the hard version of the control. A null that scattered histidines anywhere on the binder would be easy to beat, and beating it would show only that interface positions differ from the rest, which nobody doubts.

| Design | Pool size | Handles in reach |
|--------|-----------|------------------|
| EGFRd3_p2g0_l66_s408238_mpnn4_model2 | 29 | 7 |
| EGFRd3_p2g0_l73_s455990_mpnn4_model1 | 29 | 7 |
| EGFRd3_p2g0_l74_s203049_mpnn7_model1 | 25 | 7 |
| EGFRd3_p2g1_l63_s42253_mpnn4_model2 | 24 | 7 |
| EGFRd3_p2g1_l68_s39881_mpnn3_model1 | 27 | 7 |

## Results

| Design | Substitutions | Observed | Null mean ± SD | Null max | p |
|--------|---------------|---------:|----------------|--------:|--:|
| EGFRd3_p2g0_l66_s408238_mpnn | Q23H V12H | 1.791 | 0.251 ± 0.383 | 1.791 | 0.003 |
| EGFRd3_p2g0_l66_s408238_mpnn | Q23H V12H | 1.436 | 0.251 ± 0.383 | 1.791 | 0.01899 |
| EGFRd3_p2g0_l66_s408238_mpnn | Q23H | 0.943 | 0.118 ± 0.271 | 0.943 | 0.03299 |
| EGFRd3_p2g0_l73_s455990_mpnn | Y63H F59H M9H | 2.379 | 0.39 ± 0.493 | 2.379 | 0.00267 |
| EGFRd3_p2g0_l73_s455990_mpnn | F59H M9H | 1.791 | 0.259 ± 0.406 | 1.791 | 0.005 |
| EGFRd3_p2g0_l73_s455990_mpnn | Y63H F59H | 1.531 | 0.259 ± 0.406 | 1.791 | 0.01466 |
| EGFRd3_p2g0_l74_s203049_mpnn | H59 V2H W70H | 2.379 | 0.556 ± 0.536 | 2.379 | 0.003 |
| EGFRd3_p2g0_l74_s203049_mpnn | V2H W70H | 1.791 | 0.381 ± 0.481 | 1.791 | 0.015 |
| EGFRd3_p2g0_l74_s203049_mpnn | H59 V2H | 1.531 | 0.381 ± 0.481 | 1.791 | 0.02932 |
| EGFRd3_p2g1_l63_s42253_mpnn4 | R49H Q57H W60H | 2.639 | 0.638 ± 0.549 | 2.639 | 0.00233 |
| EGFRd3_p2g1_l63_s42253_mpnn4 | R49H T52H Q57H | 2.379 | 0.638 ± 0.549 | 2.639 | 0.00633 |
| EGFRd3_p2g1_l63_s42253_mpnn4 | T52H Q57H W60H | 2.379 | 0.638 ± 0.549 | 2.639 | 0.00633 |
| EGFRd3_p2g1_l68_s39881_mpnn3 | M10H H11 | 1.791 | 0.492 ± 0.522 | 1.791 | 0.03632 |
| EGFRd3_p2g1_l68_s39881_mpnn3 | M10H H11 | 1.531 | 0.492 ± 0.522 | 1.791 | 0.04499 |
| EGFRd3_p2g1_l68_s39881_mpnn3 | H11 | 0.943 | 0.258 ± 0.396 | 0.943 | 0.12163 |

14 of 15 variants reach p < 0.05. 0 exceed every random draw.

## Reading this

`p` is the fraction of random draws achieving coupling at least as high as the designed variant, with the usual +1 correction so that zero never appears. `null_fraction_zero` is the proportion of random placements that engage no handle at all, which measures how much of the pool is useless for this purpose.

A small p means the mechanistic rule found something chance rarely finds. A large p means the interface is small enough that almost any exposed substitution lands near a handle, in which case the targeting contributes little and the honest conclusion is that the pool, not the rule, is doing the work.

## What this does not show

Both arms are scored by the same model, with the same assumed pKa shift. The control tests whether the placement rule beats chance *under that model*; it does not validate the model. A designed variant that beats the null is a better candidate, not a demonstrated pH switch.

Nor does this test structural plausibility. A variant can win here and still fail to fold or to bind, which is what step 06 exists to check.

Two statistical caveats. No designed variant exceeds the maximum of its null distribution: chance does occasionally find placements as good, it just does so rarely, so the claim is about frequency and not about reaching something unreachable. And the p-values are uncorrected. Fifteen tests were run, and variants within a design share positions so they are not independent, which makes a clean correction awkward; under a conservative Bonferroni threshold the strongest variants survive and the marginal ones do not.

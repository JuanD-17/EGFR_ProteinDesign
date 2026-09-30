# Step 05 - Histidine engineering

Generated 2026-09-30T14:05:30+00:00 by `scripts/05_engineer_ph_switch.py`.

5 designs examined, 15 variants proposed.

## Why substitutions come in sets

Wyman linkage caps the pH dependence of binding at about 1.23 kcal/mol per coupled protonation across the 0.9 pH units between 7.4 and 6.5, which is roughly eight-fold in affinity. The challenge asks for binding at 6.5 and none detectable at 7.4, so a single histidine cannot deliver it. Variants engaging two or three handles are therefore ranked above single substitutions regardless of individual geometry.

## Target handles

| Handle | Role | pKa free | kcal/mol |
|--------|------|---------:|---------:|
| D368 | anchor for an engineered histidine | 4.67 | 0.848 |
| H370 | receptor-side switch | 4.93 | 0.588 |
| H433 | receptor-side switch | 6.22 | 0.943 |
| E455 | anchor for an engineered histidine | 4.31 | 0.848 |
| D458 | anchor for an engineered histidine | 1.95 | 0.848 |
| D460 | anchor for an engineered histidine | 3.7 | 0.848 |
| E496 | anchor for an engineered histidine | 4.88 | 0.848 |

## Proposed variants

| Design | Substitutions | Engages | Mean distance Å | kcal/mol | Fold |
|--------|---------------|---------|----------------|---------:|-----:|
| EGFRd3_p2g0_l66_s408238_mpnn4_model2 | Q23H V12H | H433 D460 | 7.77 | 1.791 | 20.6× |
| EGFRd3_p2g0_l66_s408238_mpnn4_model2 | Q23H V12H | H370 D460 | 8.56 | 1.436 | 11.3× |
| EGFRd3_p2g0_l66_s408238_mpnn4_model2 | Q23H | H433 | 6.26 | 0.943 | 4.9× |
| EGFRd3_p2g0_l73_s455990_mpnn4_model1 | Y63H F59H M9H | H370 H433 D460 | 8.0 | 2.379 | 55.4× |
| EGFRd3_p2g0_l73_s455990_mpnn4_model1 | F59H M9H | H433 D460 | 7.32 | 1.791 | 20.6× |
| EGFRd3_p2g0_l73_s455990_mpnn4_model1 | Y63H F59H | H370 H433 | 7.07 | 1.531 | 13.3× |
| EGFRd3_p2g0_l74_s203049_mpnn7_model1 | H59 V2H W70H | H370 H433 D460 | 8.67 | 2.379 | 55.4× |
| EGFRd3_p2g0_l74_s203049_mpnn7_model1 | V2H W70H | H433 D460 | 8.7 | 1.791 | 20.6× |
| EGFRd3_p2g0_l74_s203049_mpnn7_model1 | H59 V2H | H370 H433 | 8.9 | 1.531 | 13.3× |
| EGFRd3_p2g1_l63_s42253_mpnn4_model2 | R49H Q57H W60H | D368 H433 D460 | 7.33 | 2.639 | 86.0× |
| EGFRd3_p2g1_l63_s42253_mpnn4_model2 | R49H T52H Q57H | D368 H370 H433 | 7.24 | 2.379 | 55.4× |
| EGFRd3_p2g1_l63_s42253_mpnn4_model2 | T52H Q57H W60H | H370 H433 D460 | 6.75 | 2.379 | 55.4× |
| EGFRd3_p2g1_l68_s39881_mpnn3_model1 | M10H H11 | D368 H433 | 7.79 | 1.791 | 20.6× |
| EGFRd3_p2g1_l68_s39881_mpnn3_model1 | M10H H11 | H370 H433 | 5.89 | 1.531 | 13.3× |
| EGFRd3_p2g1_l68_s39881_mpnn3_model1 | H11 | H433 | 5.8 | 0.943 | 4.9× |

## Selection rules

- A binder position qualifies when its CB lies 4.0–10.0 Å from the handle's titratable atom, which is the range a histidine rotamer can bridge, and its side chain points towards the handle.
- Cysteine, proline and glycine are not substituted: C (may be disulfide-bonded), P (backbone conformation depends on proline), G (may occupy a position requiring positive phi).
- Positions below 0.2 relative SASA in the unbound binder are excluded as core: a charge introduced there destabilises the fold instead of forming a switch. SASA is measured on the binder alone, so a position that is exposed when free and contacts the target when bound still qualifies.
- A position that is already histidine and faces a handle is kept as an existing switch rather than discarded, since it supplies a protonation event at no cost in substitutions.
- One binder residue cannot serve two handles, so combinations reusing a position are discarded.

## Limitations

The coupling figures are carried over from step 03b's estimates for the target handles and assume the engineered histidine achieves the same pKa shift on binding. That shift is the quantity that actually decides the outcome and it cannot be computed before the variant complex exists. Nothing here establishes that a proposed substitution preserves binding at all: introducing a histidine changes the interface and may simply break it.

**Every variant must be refolded and rescored before it is believed.** This ranking is a prioritisation for that validation, not a result.

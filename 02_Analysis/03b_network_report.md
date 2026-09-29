# Step 03b - pKa and protonation networks

Generated 2026-09-29T17:34:32+00:00 by `scripts/03b_protonation_networks.py`.

## Selectivity comes from the pKa shift, not the pKa

The linkage coefficient is the number of protons taken up on binding, `Δn(pH) = f_bound(pH) − f_free(pH)`. A group whose pKa does not change on binding contributes exactly nothing, however favourable its interaction. The design requirement is therefore to place each histidine so that **binding stabilises its protonated form** — burying it against a carboxylate is the mechanism that achieves this.

Estimates below assume a shift of 2.0 pKa units on binding, which is modest for a buried salt bridge. Coupling is integrated numerically over the 6.5–7.4 window rather than taken from the two-state ceiling.

## Free-state pKa of the reachable groups

| Group | Role | pKa free | Titrating pKa | Bound | kcal/mol | Fold | Verdict |
|-------|------|---------:|--------------:|------:|---------:|-----:|---------|
| D368 | anchor for an engineered histidine | 4.67 | 6.5 | 8.5 | 0.848 | 4.19× | good anchor |
| H370 | receptor-side switch | 4.93 | 4.93 | 6.93 | 0.588 | 2.7× | usable |
| H433 | receptor-side switch | 6.22 | 6.22 | 8.22 | 0.943 | 4.91× | usable |
| E455 | anchor for an engineered histidine | 4.31 | 6.5 | 8.5 | 0.848 | 4.19× | excellent anchor, fully ionised |
| D458 | anchor for an engineered histidine | 1.95 | 6.5 | 8.5 | 0.848 | 4.19× | excellent anchor, fully ionised |
| D460 | anchor for an engineered histidine | 3.7 | 6.5 | 8.5 | 0.848 | 4.19× | excellent anchor, fully ionised |
| E496 | anchor for an engineered histidine | 4.88 | 6.5 | 8.5 | 0.848 | 4.19× | good anchor |

## Geometrically viable networks

A network is viable when its groups lie between 5.0 and 22.0 Å apart — close enough for one binder face to reach all of them, far enough that they titrate against the binder rather than each other — and their outward normals fall within 80.0°.

| Groups | Receptor His | Engineered His | Span Å | Angle ° | kcal/mol | Fold | Patches |
|--------|-------------:|---------------:|-------:|--------:|---------:|-----:|---------|
| D368 H433 D458 | 1 | 2 | 15.9 | 20.0 | 2.639 | 86.0× | 1 2 3 |
| D368 H433 D460 | 1 | 2 | 15.1 | 16.5 | 2.639 | 86.0× | 1 2 3 |
| H433 D458 D460 | 1 | 2 | 15.9 | 20.0 | 2.639 | 86.0× | 1 2 3 |
| D368 E455 D458 | 0 | 3 | 17.1 | 30.0 | 2.544 | 73.2× | 5 |
| D368 E455 D460 | 0 | 3 | 17.1 | 30.0 | 2.544 | 73.2× | 5 |
| D368 D458 D460 | 0 | 3 | 15.1 | 18.2 | 2.544 | 73.2× | 1 2 3 5 |
| E455 D458 D460 | 0 | 3 | 16.2 | 14.2 | 2.544 | 73.2× | 5 |
| D368 H370 H433 | 2 | 1 | 11.3 | 8.5 | 2.379 | 55.4× | 1 3 |
| H370 H433 D458 | 2 | 1 | 16.9 | 25.8 | 2.379 | 55.4× | 1 3 |
| H370 H433 D460 | 2 | 1 | 17.3 | 21.8 | 2.379 | 55.4× | 1 3 |
| D368 H370 E455 | 1 | 2 | 21.3 | 35.9 | 2.284 | 47.2× | none |
| D368 H370 D458 | 1 | 2 | 16.9 | 25.8 | 2.284 | 47.2× | 1 3 |
| D368 H370 D460 | 1 | 2 | 17.3 | 21.8 | 2.284 | 47.2× | 1 3 |
| H370 E455 D458 | 1 | 2 | 21.3 | 35.9 | 2.284 | 47.2× | none |
| H370 E455 D460 | 1 | 2 | 21.3 | 35.9 | 2.284 | 47.2× | none |
| H370 D458 D460 | 1 | 2 | 17.3 | 25.8 | 2.284 | 47.2× | 1 3 |
| D368 H433 | 1 | 1 | 11.3 | 3.5 | 1.791 | 20.6× | 1 2 3 |
| H433 D458 | 1 | 1 | 15.9 | 20.0 | 1.791 | 20.6× | 1 2 3 |
| H433 D460 | 1 | 1 | 12.8 | 16.5 | 1.791 | 20.6× | 1 2 3 |
| D368 E455 | 0 | 2 | 17.1 | 30.0 | 1.696 | 17.5× | 5 |
| D368 D458 | 0 | 2 | 12.5 | 18.2 | 1.696 | 17.5× | 1 2 3 5 |
| D368 D460 | 0 | 2 | 15.1 | 16.1 | 1.696 | 17.5× | 1 2 3 5 |
| E455 D458 | 0 | 2 | 11.4 | 14.0 | 1.696 | 17.5× | 5 |
| E455 D460 | 0 | 2 | 16.2 | 14.2 | 1.696 | 17.5× | 5 |
| D458 D460 | 0 | 2 | 8.1 | 7.9 | 1.696 | 17.5× | 1 2 3 5 |
| D460 E496 | 0 | 2 | 21.3 | 24.6 | 1.696 | 17.5× | none |
| H370 H433 | 2 | 0 | 9.8 | 5.8 | 1.531 | 13.3× | 1 3 |
| D368 H370 | 1 | 1 | 5.6 | 8.5 | 1.436 | 11.3× | 1 3 |
| H370 E455 | 1 | 1 | 21.3 | 35.9 | 1.436 | 11.3× | none |
| H370 D458 | 1 | 1 | 16.9 | 25.8 | 1.436 | 11.3× | 1 3 |
| H370 D460 | 1 | 1 | 17.3 | 21.8 | 1.436 | 11.3× | 1 3 |

## Limitations

pKa values in the complex cannot be computed before the complex exists, so the assumed shift is a stand-in for the quantity that actually decides the outcome. PROPKA is empirical with errors near one pKa unit, comparable to the effects reasoned about here. Rigorous coupling free energies would require constant-pH molecular dynamics, which is out of scope.

Every anchor is scored with the same engineered-histidine model, so all anchors return an identical coupling energy and the ranking among anchor-only networks is driven by how many groups they contain rather than by their quality. In this model anchor pKa governs *reliability* — whether the carboxylate is dependably ionised across the window — not the magnitude of the shift it induces. A more careful treatment would let anchor pKa and burial depth modulate the achievable shift, and would separate the networks that here score identically.

These numbers are a screen that removes groups whose free-state pKa already rules them out, and a geometric filter on which combinations a single binder could engage. They are not predictions of experimental pH selectivity.

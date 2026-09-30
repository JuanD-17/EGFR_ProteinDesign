# Step 07 - Structural validation, preparation

Generated 2026-09-30T14:49:24+00:00 by `scripts/07_validate_ph_variants.py --prepare`.

## Why the scores stay separate

Introducing charged side chains into a packed interface can make it worse. A variant with an elegant protonation network and a broken interface is worthless, so the network score is never combined with the structural scores. Three layers are reported independently: folding plausibility, interface quality, and whether the designed network survives.

## Rotamer feasibility

For each substitution, an imidazole nitrogen is placed at 200 orientations on spheres of 2.53 and 3.67 Å around CB, the CB-ND1 and CB-NE2 distances measured from the 23 histidines resolved in 6ARU. Orientations clashing with the binder or the target within 2.9 Å are discarded, and the substitution passes only if a surviving orientation lands 2.6–4.0 Å from the handle.

This is a feasibility test rather than a rotamer library search. It answers whether the contact is geometrically possible at all, which is a veto AlphaFold2 does not provide and which is worth applying before spending GPU time.

12 of 14 substitutions can reach their handle; 2 cannot.

## Prepared for prediction

7 variants and 5 unmodified parents, as complex records in `07_variants_for_af2.fasta`.

The parents matter as much as the variants. Without predicting the unmodified design under identical settings there is no baseline, and any change in interface confidence could not be attributed to the substitutions rather than to prediction noise.

## Next

Re-predict every record on GPU, then:

```bash
python scripts/07_validate_ph_variants.py --score --predictions <dir>
```

# Step 09 - Novelty of the designed binders

Generated 2026-09-30T23:33:14+00:00 by `scripts/09_check_novelty.py`.

## Why this is measured

The competition filters out designs that are not de novo and zero-shot, or that lack adequate sequence and structural diversity from known proteins. These binders come from AlphaFold2 hallucination with no starting binder, so novelty is expected, but an expectation is not evidence.

## Against the Protein Data Bank

**No design returned a hit.** The search covers every polymer entity in the PDB, so this is the strongest available statement that these sequences do not resemble a protein of known structure.

## Composition

Hallucinated sequences can drift into compositions no natural protein has, which is a different novelty problem: too novel to express or fold. Total variation distance from mean UniProtKB frequencies, where 0 is identical and 1 is disjoint.

| Design | Deviation | Largest departures |
|--------|----------:|--------------------|
| EGFRd3_p2g0_l66_s408238_mpnn4_ | 0.270 | E 15% (nat 7%), A 2% (nat 8%), G 2% (nat 7%) |
| EGFRd3_p2g0_l73_s455990_mpnn4_ | 0.320 | L 0% (nat 10%), K 15% (nat 6%), E 15% (nat 7%) |
| EGFRd3_p2g0_l74_s203049_mpnn7_ | 0.235 | E 14% (nat 7%), R 11% (nat 6%), V 11% (nat 7%) |
| EGFRd3_p2g0_l74_s258459_mpnn1_ | 0.348 | E 23% (nat 7%), K 16% (nat 6%), L 4% (nat 10%) |
| EGFRd3_p2g1_l55_s473498_mpnn1_ | 0.229 | E 16% (nat 7%), G 0% (nat 7%), A 2% (nat 8%) |
| EGFRd3_p2g1_l63_s42253_mpnn4_m | 0.247 | E 17% (nat 7%), K 14% (nat 6%), G 3% (nat 7%) |
| EGFRd3_p2g1_l66_s66715_mpnn4_m | 0.357 | E 21% (nat 7%), K 18% (nat 6%), G 0% (nat 7%) |
| EGFRd3_p2g1_l68_s39881_mpnn3_m | 0.384 | E 28% (nat 7%), A 1% (nat 8%), D 0% (nat 5%) |

## Against each other

A portfolio of near-identical binders is not a portfolio, whatever its distance from the PDB.

| Design A | Design B | Identity |
|----------|----------|---------:|
| EGFRd3_p2g0_l74_s258459_mp | EGFRd3_p2g1_l63_s42253_mpn | 27.0% |
| EGFRd3_p2g0_l66_s408238_mp | EGFRd3_p2g1_l63_s42253_mpn | 25.8% |
| EGFRd3_p2g0_l74_s258459_mp | EGFRd3_p2g1_l66_s66715_mpn | 25.7% |
| EGFRd3_p2g0_l66_s408238_mp | EGFRd3_p2g1_l55_s473498_mp | 24.2% |
| EGFRd3_p2g1_l66_s66715_mpn | EGFRd3_p2g1_l68_s39881_mpn | 23.5% |
| EGFRd3_p2g0_l73_s455990_mp | EGFRd3_p2g1_l66_s66715_mpn | 23.3% |
| EGFRd3_p2g0_l73_s455990_mp | EGFRd3_p2g1_l68_s39881_mpn | 21.9% |
| EGFRd3_p2g0_l74_s203049_mp | EGFRd3_p2g0_l74_s258459_mp | 21.6% |
| EGFRd3_p2g0_l74_s258459_mp | EGFRd3_p2g1_l68_s39881_mpn | 21.6% |
| EGFRd3_p2g0_l66_s408238_mp | EGFRd3_p2g1_l66_s66715_mpn | 21.2% |
| EGFRd3_p2g1_l63_s42253_mpn | EGFRd3_p2g1_l66_s66715_mpn | 21.2% |
| EGFRd3_p2g1_l63_s42253_mpn | EGFRd3_p2g1_l68_s39881_mpn | 20.6% |
| EGFRd3_p2g0_l73_s455990_mp | EGFRd3_p2g1_l63_s42253_mpn | 20.5% |
| EGFRd3_p2g0_l73_s455990_mp | EGFRd3_p2g0_l74_s203049_mp | 18.9% |
| EGFRd3_p2g0_l73_s455990_mp | EGFRd3_p2g0_l74_s258459_mp | 18.9% |
| EGFRd3_p2g0_l74_s203049_mp | EGFRd3_p2g1_l63_s42253_mpn | 18.9% |
| EGFRd3_p2g1_l55_s473498_mp | EGFRd3_p2g1_l66_s66715_mpn | 18.2% |
| EGFRd3_p2g0_l74_s203049_mp | EGFRd3_p2g1_l66_s66715_mpn | 17.6% |
| EGFRd3_p2g1_l55_s473498_mp | EGFRd3_p2g1_l63_s42253_mpn | 17.5% |
| EGFRd3_p2g0_l66_s408238_mp | EGFRd3_p2g0_l73_s455990_mp | 16.4% |
| EGFRd3_p2g0_l66_s408238_mp | EGFRd3_p2g0_l74_s203049_mp | 16.2% |
| EGFRd3_p2g0_l66_s408238_mp | EGFRd3_p2g1_l68_s39881_mpn | 16.2% |
| EGFRd3_p2g0_l73_s455990_mp | EGFRd3_p2g1_l55_s473498_mp | 15.1% |
| EGFRd3_p2g0_l74_s203049_mp | EGFRd3_p2g1_l55_s473498_mp | 14.9% |
| EGFRd3_p2g0_l74_s203049_mp | EGFRd3_p2g1_l68_s39881_mpn | 14.9% |
| EGFRd3_p2g1_l55_s473498_mp | EGFRd3_p2g1_l68_s39881_mpn | 14.7% |
| EGFRd3_p2g0_l66_s408238_mp | EGFRd3_p2g0_l74_s258459_mp | 13.5% |
| EGFRd3_p2g0_l74_s258459_mp | EGFRd3_p2g1_l55_s473498_mp | 9.5% |

## Limitations

Sequence identity is a weak proxy for structural novelty. Two proteins can share a fold at 15% identity, and the competition asks for structural diversity as well as sequence diversity. A complete answer would compare the designed backbones against the PDB structurally, with something like Foldseek, which is not run here.

No automatic threshold is applied. The competition publishes no numeric cutoff, so hits are reported and the judgement is left explicit rather than hidden in a filter.

# Methods

Conditional EGFR binder design for Challenge 01 of the Anthropic × Adaptyv 2026 competition, Track 3.

Generated 2026-10-01T14:54:20+00:00 from the pipeline outputs by `scripts/10_build_submission.py`. Every figure below is read from the analysis files rather than transcribed.

---

## Summary

16 single-chain de novo proteins of 55–74 residues, targeting domain III of the human EGFR extracellular region. 8 carry an engineered histidine network intended to make binding pH-dependent. 1 already presented such a network without modification and is submitted unchanged. 7 are parent binders with no network, included as a fallback should the switch fail.

That a design arrived with a usable protonation network already in place is reported rather than absorbed into the engineered count. BindCraft optimises interface confidence and has no notion of pH, so the network is incidental, and counting it as engineered would claim credit for something the pipeline did not do. It is submitted regardless, because an incidental network is as real as a designed one.

## Objectives addressed

The challenge asks for three properties, ranked with pH selectivity first.

*Cross-species reactivity* is handled upstream rather than designed for. A residue entered the designable pool only if it is **identical** between human and mouse, not merely conserved. The bar is deliberately absolute: the cetuximab epitope is 20 of 27 residues identical with none differing outright, and cetuximab still fails to recognise murine EGFR, so a handful of conservative substitutions across a footprint is enough to abolish antibody binding. Every submitted binder therefore targets a surface with no human-mouse difference at all.

*pH selectivity* is engineered, since the target cannot be mutated and no binder-design pipeline optimises for it.

*Affinity* is what BindCraft optimises, and is the property these designs are least differentiated on.

## Software

| Tool | Version | Role |
|------|---------|------|
| BindCraft | `martinpacesa/BindCraft`, cloned 2026-09-29 | binder generation |
| AlphaFold2 | `alphafold_params_2022-12-06`, multimer v3 | backbone hallucination and complex prediction |
| ProteinMPNN | as vendored by BindCraft | sequence design |
| PyRosetta | `2026.29+release.quarterly` | interface relaxation and scoring |
| PROPKA | 3.5.1 | free-state pKa estimation |
| FreeSASA | 2.2.1 | solvent accessibility |
| Biopython | 1.88 | structure and sequence handling |

PyRosetta and the AlphaFold2 parameters carry non-commercial licences; this work is academic. Full terms are recorded in the repository.

## Target definition

Human EGFR (UniProt P00533, entry version 301) and mouse EGFR (Q01279, version 253) were retrieved from UniProtKB release 2026_03, and the EGFR ectodomain in complex with a cetuximab Fab mutant (PDB 6ARU) from the RCSB PDB, all on 29 September 2026. SHA-256 checksums are recorded in the accompanying repository.

Every residue of the extracellular region (25–645, precursor numbering) was annotated on four layers: human/mouse conservation from a global BLOSUM62 alignment; solvent accessibility computed with FreeSASA on chain A of 6ARU; occlusion by N-glycans, using the sugars the crystal resolved and the glycosylated asparagine where it did not; and disulfide-derived architecture. A residue was called designable when exposed above 0.2 relative SASA, identical between species, clear of glycan, and not disulfide-bonded. 168 of 621 qualified.

### Numbering

Structure 6ARU is numbered by the mature protein, 24 below UniProt. The offset was derived by alignment and verified independently: mapping UniProt's 25 annotated disulfide pairs through it placed 24 of them at 2.03–2.05 Å SG–SG in the crystal, against the 2.05 Å of a covalent bond. That test is geometric and does not reuse the alignment it validates.

## Epitope selection

Designable residues were grouped into contiguous surface patches and ranked on pH-switch potential, conservation of the footprint including a shell beyond it, overlap with the cetuximab interface as evidence of druggability, and buriable area. Patches were then characterised geometrically, which showed that several were not single connected surfaces; residues close in a distance matrix but detached from the patch face were removed from the hotspot lists.

Patch 2 was selected. It is the only one among the top-ranked patches forming a single connected surface, it contacts the cetuximab Fab directly, and it contains H433, Q435 and K489, three positions with prior experimental support that the selection recovered without using the literature as input.

## pH-switch mechanism

Wyman linkage caps the pH dependence of binding at 2.303·R·T·Δn per pH unit, so across the 0.9 units between 7.4 and 6.5 each fully coupled protonation contributes at most 1.23 kcal/mol. A single histidine therefore buys about eightfold, which is not the requested behaviour, and the design families were organised around the number of coupled protonation events rather than around any one position.

Selectivity arises from the pKa *shift* on binding, not from the pKa: a group whose protonation state is unchanged on binding contributes nothing. Free-state pKa values were estimated with PROPKA. H433 came out at 6.22, titrating almost ideally for this window; the target carboxylates D458, D460 and E455 sit between 1.95 and 4.31 and are therefore dependable counter-charges that upshift an engineered histidine's pKa on burial.

## Binder generation

Domain III was trimmed from 6ARU chain A to 204 residues, which is what makes the job fit in the memory of the available hardware. Binders were generated with BindCraft, which hallucinates a backbone by AlphaFold2 backpropagation, designs its sequence with ProteinMPNN, re-predicts the complex and scores it with PyRosetta. Optimisation iterations were reduced from the published defaults to fit a 30 GPU-hour weekly budget on two T4 cards; validation settings were left untouched, since weakening the stage that decides whether a design is good would only move the failure downstream.

59 trajectories were run across two batches, of which 8 produced accepted designs. 483 of 628 rejections were on AlphaFold2 confidence in the re-predicted complex rather than on composition or geometry, so the filters were left as published: relaxing them would accept designs the model does not believe in.

## Histidine engineering

Histidines were placed at binder positions facing a target ionizable group, within the span a histidine can actually cover. That span was measured on the 23 histidines resolved in 6ARU, giving CB–ND1 at 2.53 Å and CB–NE2 at 3.67 Å; with a 4.0 Å charge-assisted contact the ceiling is 7.7 Å from CB. Substitutions were proposed in sets rather than singly, since the thermodynamics make the count the governing variable, and positions already carrying a histidine were kept as existing switches rather than discarded.

Each proposal was then tested for rotamer feasibility: an imidazole nitrogen was placed at 200 orientations on spheres of both measured radii around CB, orientations clashing with the binder or target discarded, and the substitution accepted only if a surviving orientation reached hydrogen-bonding distance of the handle.

## Controls and validation

**Negative control.** Designed placements were compared against histidines placed at random among the positions the method could have chosen from: interface, exposed in the unbound binder, and chemically substitutable, with only the targeting removed. Across 3000 draws per design and substitution count, designed variants average 1.82 kcal/mol against 0.40 for random placement. 5 of the 8 submitted variants reach p < 0.05. No designed variant exceeds the maximum of its null distribution, so the claim is about frequency rather than about reaching something unreachable, and the p-values are uncorrected.

**Structural revalidation.** Every variant and its unmodified parent were re-predicted from sequence under identical settings. 6 of the 9 submitted designs that carry a network retain interface confidence within 0.05 i_pTM of their parent. Damage does not track the number of substitutions: the lead variant, which engages three handles, changes i_pTM by +0.002, while one single substitution collapsed its complex from 0.819 to 0.243 and was excluded. The variants that cost least are those built on histidines the design already carried.

**Novelty.** No binder returns a hit from the RCSB sequence service across every polymer entity in the PDB, and the highest pairwise identity within the portfolio is 27%.

## Ranking

Submissions are ordered to follow the competition's own criteria. Variants whose substitutions destroyed the interface were excluded outright. The remainder are ordered by estimated coupling, then by whether the interface survived intact rather than merely tolerably, then by the control p-value, then by composition risk. Unmodified parents follow the variants: they carry no switch but are validated binders and cost nothing within the submission limit.

## Limitations

The coupling estimates are model-dependent upper bounds under an assumed 2.0-unit pKa shift on binding. That shift is the quantity that actually decides the outcome and it cannot be computed before the complex exists; a rigorous treatment would need constant-pH molecular dynamics. PROPKA's error is near one pKa unit, comparable to the effects reasoned about.

The binders were generated by AlphaFold2 and re-predicted by AlphaFold2, so absolute confidence values are not independent evidence of binding. The variant-versus-parent comparison is sound, since no variant was optimised and both are predicted identically, but only the assay can establish that any of these bind at all.

Every design is strongly acidic, with glutamate between 14 and 28 percent against a natural average near 7. This is characteristic of hallucinated sequences rather than of any one trajectory, and it is an expression and solubility liability carried by the whole portfolio. It is weighed in the ranking but not solved.

Sequence identity is a weak proxy for structural novelty, and no structural comparison against the PDB was performed.

## Reproducibility

All code, parameters, intermediate results and the provenance of every input are in the accompanying repository. Analysis scripts are numbered in execution order and run on CPU in minutes; only binder generation and the revalidation require a GPU. Raw inputs are reconstructed by a download script and pinned by SHA-256 checksums that can be re-verified with a single command.

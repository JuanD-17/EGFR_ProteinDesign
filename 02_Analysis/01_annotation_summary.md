# Step 01 - EGFR ectodomain annotation

Generated 2026-09-29T15:56:47+00:00 by `scripts/01_annotate_target.py`.

Ectodomain residues 25-645 (621 aa), UniProt precursor numbering. Structure 6ARU chain A, numbering offset UniProt - PDB = +24.

## Domain architecture, derived from the disulfide pattern

UniProt does not annotate the ectodomain subdomains. These boundaries come from the distribution of disulfide bonds: cysteine-free stretches mark the L-domain cores.

| Range | Length | Character | Assignment |
|-------|--------|-----------|------------|
| 25-58 | 34 aa | cysteine-rich | I/II (N-terminal, cysteine-rich) |
| 59-156 | 98 aa | cysteine-free | I (L-domain core) |
| 157-362 | 206 aa | cysteine-rich | II (furin-like) |
| 363-469 | 107 aa | cysteine-free | III (L-domain core) |
| 470-645 | 176 aa | cysteine-rich | IV (furin-like) |

## Coverage

- Resolved in the crystal structure: 609/621
- Exposed (relative SASA >= 0.2): 360
- Cetuximab interface: 27 residues
- **Designable** (exposed, conserved, unglycosylated, not disulfide-bonded): 168

## Human/mouse conservation across the ectodomain

- identical: 551 (88.7%)
- similar: 52 (8.4%)
- different: 18 (2.9%)
- unaligned: 0 (0.0%)

## Conservation of the cetuximab epitope

This is the control on the decision to avoid it.

- identical: 20 (74.1%)
- similar: 7 (25.9%)
- different: 0 (0.0%)
- unaligned: 0 (0.0%)

## Designable surface per domain

| Domain | Residues | Designable | Fraction |
|--------|----------|------------|----------|
| II (furin-like) | 206 | 74 | 35.9% |
| IV (furin-like) | 176 | 52 | 29.5% |
| III (L-domain core) | 107 | 20 | 18.7% |
| I/II (N-terminal, cysteine-rich) | 34 | 13 | 38.2% |
| I (L-domain core) | 98 | 9 | 9.2% |

## Glycosylation sites in the ectodomain

| Position | Domain | Description |
|----------|--------|-------------|
| N56 | I/II (N-terminal, cysteine-rich) | N-linked (GlcNAc...) (complex) asparagine; atypical; partial |
| N73 | I (L-domain core) | N-linked (GlcNAc...) asparagine; atypical |
| N128 | I (L-domain core) | N-linked (GlcNAc...) asparagine |
| N175 | II (furin-like) | N-linked (GlcNAc...) asparagine |
| N196 | II (furin-like) | N-linked (GlcNAc...) asparagine |
| N352 | II (furin-like) | N-linked (GlcNAc...) asparagine |
| N361 | II (furin-like) | N-linked (GlcNAc...) asparagine |
| N413 | III (L-domain core) | N-linked (GlcNAc...) asparagine |
| N444 | III (L-domain core) | N-linked (GlcNAc...) asparagine |
| N528 | IV (furin-like) | N-linked (GlcNAc...) asparagine |
| N568 | IV (furin-like) | N-linked (GlcNAc...) asparagine; partial |
| N603 | IV (furin-like) | N-linked (GlcNAc...) asparagine; partial |
| N623 | IV (furin-like) | N-linked (GlcNAc...) (high mannose) asparagine |

Surface within 12.0 A of a glycan is treated as occluded. Glycan trees are rarely resolved crystallographically, so this surface looks accessible in the structure but is not.

## Histidines on the designable surface

Relevant to objective 1. A target histidine near the interface can pair with a binder carboxylate, adding a second protonation-sensitive contact on top of the histidines engineered into the binder itself.

| Position | Domain | Relative SASA | Conservation | Glycan shadow |
|----------|--------|---------------|--------------|---------------|
| H145 | I (L-domain core) | 0.468 | different | yes |
| H183 | II (furin-like) | 0.864 | identical | no |
| H233 | II (furin-like) | 0.334 | identical | no |
| H304 | II (furin-like) | 0.775 | identical | no |
| H358 | II (furin-like) | 0.583 | identical | yes |
| H383 | III (L-domain core) | 0.882 | similar | yes |
| H433 | III (L-domain core) | 0.801 | identical | no |
| H507 | IV (furin-like) | 0.387 | similar | yes |
| H559 | IV (furin-like) | 0.325 | identical | no |
| H584 | IV (furin-like) | 0.47 | identical | no |
| H615 | IV (furin-like) | 0.51 | similar | no |
| H618 | IV (furin-like) | 0.463 | identical | no |
| H621 | IV (furin-like) | 0.342 | identical | yes |

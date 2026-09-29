# Step 02 - Epitope patch selection

Generated 2026-09-29T16:54:05+00:00 by `scripts/02_select_epitope.py`.

Scope: UniProt 335-538. 110 exposed residues, 45 in the designable pool, 12 distinct candidate surfaces.

A pool residue is one a binder is *permitted* to touch. A patch is a contiguous surface a binder could actually land on. The pool is not a hotspot list.

## Ranked patches

| Rank | Centre | Residues | Area Å² | His | Cetuximab overlap | Conserved % | Score |
|------|--------|----------|---------|-----|-------------------|-------------|-------|
| 1 | 430 | 10 | 785 | H433 | 3 | 94.3 | 0.8461 |
| 2 | 460 | 11 | 846 | H433 | 5 | 88.2 | 0.8001 |
| 3 | 432 | 8 | 566 | H433 | 4 | 85.7 | 0.755 |
| 4 | 490 | 9 | 638 | — | 6 | 70.8 | 0.6809 |
| 5 | 427 | 11 | 761 | — | 0 | 95.1 | 0.6504 |
| 6 | 487 | 8 | 720 | — | 4 | 86.1 | 0.6113 |
| 7 | 483 | 8 | 518 | — | 1 | 90.6 | 0.5502 |
| 8 | 473 | 9 | 646 | — | 3 | 69.6 | 0.4418 |
| 9 | 452 | 8 | 753 | — | 0 | 88.6 | 0.4241 |
| 10 | 479 | 8 | 710 | — | 0 | 85.3 | 0.4163 |
| 11 | 476 | 8 | 621 | — | 0 | 72.7 | 0.3786 |
| 12 | 455 | 10 | 698 | — | 0 | 92.1 | 0.3673 |

## Top 5 in detail

### Patch 1 — centred on residue 430

- **Score 0.8461** (pH 1.0, conservation 0.943, functional 0.375, geometry 0.883)
- 10 pool residues, 785 Å² buriable
- Histidines in patch: H433
- Histidines adjacent: H358
- Overlaps the cetuximab interface at 3 residues
- Nearest glycan 15.3 Å
- Footprint 94.3% identical in mouse, counting the surrounding shell

Residues, UniProt precursor numbering:

```
S366 R427 R429 T430 K431 Q432 H433 G434 Q435 D458
```

BindCraft `target_hotspot_residues`, PDB numbering:

```
A342,A403,A405,A406,A407,A408,A409,A410,A411,A434
```

### Patch 2 — centred on residue 460

- **Score 0.8001** (pH 0.801, conservation 0.882, functional 0.625, geometry 0.9)
- 11 pool residues, 846 Å² buriable
- Histidines in patch: H433
- Histidines adjacent: none
- Overlaps the cetuximab interface at 5 residues
- Nearest glycan 16.0 Å
- Footprint 88.2% identical in mouse, counting the surrounding shell

Residues, UniProt precursor numbering:

```
K431 Q432 H433 G434 Q435 D458 D460 T483 K487 T488 K489
```

BindCraft `target_hotspot_residues`, PDB numbering:

```
A407,A408,A409,A410,A411,A434,A436,A459,A463,A464,A465
```

### Patch 3 — centred on residue 432

- **Score 0.755** (pH 0.801, conservation 0.857, functional 0.5, geometry 0.774)
- 8 pool residues, 566 Å² buriable
- Histidines in patch: H433
- Histidines adjacent: none
- Overlaps the cetuximab interface at 4 residues
- Nearest glycan 12.1 Å
- Footprint 85.7% identical in mouse, counting the surrounding shell

Residues, UniProt precursor numbering:

```
Q408 T430 K431 Q432 H433 G434 Q435 D460
```

BindCraft `target_hotspot_residues`, PDB numbering:

```
A384,A406,A407,A408,A409,A410,A411,A436
```

### Patch 4 — centred on residue 490

- **Score 0.6809** (pH 0.594, conservation 0.708, functional 0.75, geometry 0.807)
- 9 pool residues, 638 Å² buriable
- Histidines in patch: none
- Histidines adjacent: H433, H507
- Overlaps the cetuximab interface at 6 residues
- Nearest glycan 12.3 Å
- Footprint 70.8% identical in mouse, counting the surrounding shell

Residues, UniProt precursor numbering:

```
S464 N473 T474 K487 T488 K489 I490 N493 E496
```

BindCraft `target_hotspot_residues`, PDB numbering:

```
A440,A449,A450,A463,A464,A465,A466,A469,A472
```

### Patch 5 — centred on residue 427

- **Score 0.6504** (pH 0.692, conservation 0.951, functional 0.0, geometry 0.883)
- 11 pool residues, 761 Å² buriable
- Histidines in patch: none
- Histidines adjacent: H358, H433
- Overlaps the cetuximab interface at 0 residues
- Nearest glycan 15.3 Å
- Footprint 95.1% identical in mouse, counting the surrounding shell

Residues, UniProt precursor numbering:

```
S366 E424 R427 R429 T430 K431 E455 S457 D458 G482 T483
```

BindCraft `target_hotspot_residues`, PDB numbering:

```
A342,A400,A403,A405,A406,A407,A431,A433,A434,A458,A459
```

## Scoring

Weights follow the competition's ranking order, in which pH selectivity outranks cross-reactivity, which outranks affinity.

| Term | Weight | Meaning |
|------|--------|---------|
| pH potential | 0.4 | exposure-weighted histidines in and beside the patch |
| Conservation | 0.3 | fraction of the footprint identical in mouse, shell included |
| Functional | 0.2 | overlap with the cetuximab interface, counted as positive evidence of druggability |
| Geometry | 0.1 | buriable area in the 600-900 Å² window, and clearance from glycans |

The weights are a judgement, not a measurement. They are written in one place in the source so they can be changed and the ranking re-derived.

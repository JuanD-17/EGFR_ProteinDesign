# Epitope selection

How the designable surface was narrowed to one patch, and why that patch.

Depends on [target_analysis.md](target_analysis.md); feeds [ph_mechanism.md](ph_mechanism.md).

---

### `scripts/02_select_epitope.py` — epitope patch selection

Step 01 yields a *pool* of residues a binder is structurally permitted to
touch. That pool is not a hotspot list, and the distinction drives this step:
a binder does not contact 45 scattered residues, it lands on a contiguous
surface of roughly 600–900 Å², on the order of 15–25 residues.

Every pool residue is taken as a patch centre and the patch is the pool
residues within 11 Å of it, CB to CB. Patches are scored, ranked, and greedily
deduplicated so the selection spans distinct surfaces rather than many views
of one. Four terms, weighted to follow the competition's ranking order:

| Term | Weight | Meaning |
|------|--------|---------|
| pH potential | 0.40 | exposure-weighted histidines in and beside the patch |
| Conservation | 0.30 | fraction of the footprint identical in mouse, shell included |
| Functional | 0.20 | overlap with the cetuximab interface |
| Geometry | 0.10 | buriable area in window, clearance from glycans |

Conservation is measured over the patch *and* a 13 Å shell. The patch core is
identical by construction, so it carries no information; a binder footprint
spills past it, and the shell is what actually decides cross-reactivity.

Overlap with the cetuximab interface is scored as **positive**. The
competition recommends targeting a functional epitope and names cetuximab's,
and a surface with a therapeutic antibody bound to it is demonstrably
druggable. The non-conserved residues within that interface are already
excluded by the pool definition, so the risk they pose to objective 2 is
handled upstream rather than by avoiding the region.

```bash
python scripts/02_select_epitope.py
python scripts/02_select_epitope.py --scope 157 362     # domain II instead
```

Outputs a patch table, a readable report, a PyMOL session for visual checking,
and `02_hotspots_bindcraft.json` carrying the `target_hotspot_residues` string
**in PDB numbering**, since BindCraft indexes into the structure file it is
given and 6ARU is numbered 24 lower than UniProt.

The weights are a judgement, not a measurement, and the ranking is sensitive
to them. They live in one dictionary in the source so the ranking can be
re-derived under different assumptions.

## Findings from step 02

Searching domain III (335–538) gives 45 pool residues, 23 patches of at least
8 residues, and 12 distinct surfaces after deduplication.

**Three of the six literature residues cannot be interface contacts.** They
were used as a blind test: the pipeline never reads them, so their fate is a
check on whether structural filtering recovers what experiments found.

| Residue | Rel. SASA | Exposed | Conservation | Glycan shadow | In pool |
|---------|-----------|---------|--------------|---------------|---------|
| H370 | 0.171 | no | identical | **yes** | no |
| R377 | 0.704 | yes | **similar** | **yes** | no |
| L406 | 0.053 | **no** | identical | no | no |
| H433 | 0.801 | yes | identical | no | **yes** |
| Q435 | 0.489 | yes | identical | no | **yes** |
| K489 | 0.544 | yes | identical | no | **yes** |

H370, described in the literature as a mechanistic pH-switch histidine, sits
below the exposure threshold and is glycan-shadowed in this structure. L406
sits at 5% relative accessibility, in the hydrophobic core.

*(Step 03 revisits H370 with a finer test and finds it marginal rather than
unavailable — see* Findings from step 03 *below. L406 is unaffected.)*

What follows from that is bounded, and worth stating carefully. Low
structural accessibility makes it improbable that these residues act as
direct interface contacts *in the crystallographic state analysed here*. It
does not establish what mechanism produced their alanine-scan phenotypes;
structural analysis alone cannot demonstrate that. An indirect effect such
as destabilisation is one consistent explanation among others, and H370 may
well be mechanistically relevant to binding or pH response without being an
interface hotspot.

The usable conclusion is narrower than a mechanistic claim and still
decision-relevant: these positions are not candidates for direct contact, so
designing a binder *at* them would be a mistake that the literature alone
would have invited.

**The three that survive all land in the same patch.** Patch 2 contains
H433, Q435 and K489 together. The scoring function never saw the literature;
it saw exposure, conservation, glycan clearance, histidine content and
cetuximab overlap.

This is independent convergence and it raises the interest of patch 2. It
does not make patch 2 the correct epitope. Two filters agreeing is evidence,
not proof, and both are looking at the same crystal structure.

The two leading candidates trade off against each other:

| | Patch 1 (centre 430) | Patch 2 (centre 460) |
|---|---|---|
| Score | 0.846 | 0.800 |
| Area | 785 Å² | 846 Å² |
| Histidine | H433, H358 adjacent | H433 |
| Conservation of footprint | 94.3% | 88.2% |
| Cetuximab overlap | 3 residues | 5 residues |
| Literature residues | H433, Q435 | H433, Q435, K489 |

Patch 1 ranks higher on conservation and on pH potential, helped by H358
sitting just outside it. Patch 2 is larger, overlaps the validated interface
more, and contains every literature residue that survived filtering. Which
one leads depends on the weights, so both are carried forward rather than
resolved on the score alone. Step 02b then measures the geometry that
composition scores cannot see.

### `scripts/02b_characterize_patches.py` — patch geometry

Step 02 groups residues by proximity and scores them by composition. Neither
operation can tell whether a patch is a single usable binding face, and the
residue-overlap deduplication cannot tell two genuinely different faces from
one face found twice from different seeds. Two patches can share few
residues and still sit on the same surface.

This step measures compactness (diameter, radius of gyration, RMS deviation
from the best-fit plane), continuity (connected components under a 5 Å
heavy-atom contact rule), orientation (the outward normal of each patch),
and a graded rather than binary relationship to the Fab.

```bash
python scripts/02b_characterize_patches.py --top 6
```

## Findings from step 02b

**Three of the six top patches are not continuous surfaces.** Reporting
*which* residues detach makes the defect actionable:

| Patch | Components | Core | Detached |
|-------|------------|------|----------|
| 1 | 2 | 427–435, 458 | **S366** |
| 2 | **1** | all 11 | — |
| 3 | 2 | 430–435, 460 | Q408 |
| 4 | **1** | all 9 | — |
| 5 | 2 | 424–458, 482–483 | **S366** |
| 6 | 3 | 435–490 | 478, 483 |

S366 detaches from both patches that contain it. It lies some 60 sequence
positions from the rest and was grouped by CB proximity, but its surface
does not touch theirs. This is precisely the failure mode that a distance
matrix cannot detect and a composition score cannot penalise.

Patch 2 is the only one of the top three that forms a single connected
surface, at 11 residues, 17.4 Å across and 1.89 Å from planar, in direct
contact with the Fab.

**The face count is soft, and is reported as such.** Two pairs fall within
10° of the 45° orientation cutoff:

| Pair | Centroid separation | Normal angle |
|------|--------------------:|-------------:|
| 1 vs 2 | 7.2 Å | 47.2° |
| 2 vs 3 | 5.5 Å | 48.3° |

At a 50° threshold, patches 1, 2 and 3 merge. The honest reading is that
they are not three faces but one curved region around H433, sampled from
three seeds. The script flags these pairs as unresolved rather than
reporting a verdict that a small change of threshold would reverse.

What survives that caveat:

- **The H433 region** — patches 1, 2 and 3, most likely one surface, all in
  direct Fab contact
- **The 487–490 face** — patches 2, 4 and 6 group unambiguously
- **Patch 5** — normal 135° from patch 1 and 9.1 Å from the Fab, so a
  genuinely different surface that does not overlap the cetuximab epitope

Patch 5 is the only real diversification available within domain III. For a
20-design portfolio judged partly on novelty, that matters more than its
rank.


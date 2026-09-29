# Conditional EGFR binder design

Work for **Challenge 01 (EGFR)** of the Anthropic × Adaptyv 2026 protein design
competition, Track 3 (open track).

The goal is to design, *de novo*, a protein that binds the human epidermal
growth factor receptor (EGFR) while satisfying three conditions at once:

| # | Objective | Ranking weight |
|---|-----------|----------------|
| 1 | Selective binding at pH 6.5, no detectable binding at pH 7.4 | **highest** |
| 2 | Cross-reactivity with mouse EGFR | intermediate |
| 3 | Affinity for human EGFR | lowest |

That order is deliberate and counterintuitive. The competition rules state
explicitly that a weak but clearly pH-sensitive binder ranks above a
high-affinity binder that is not pH-sensitive: demonstrated difficulty is
rewarded over the magnitude of the dissociation constant.

## Design consequences

Two constraints follow directly from that ordering, and they shape everything
downstream.

**pH sensitivity is engineered into the binder, not discovered in the target.**
EGFR cannot be mutated. The available mechanism is to place histidines in the
binder's interface: the imidazole side chain has a pKa near 6.5, so it is
protonated in the acidic tumour microenvironment and neutral at physiological
pH. That change in charge is what should break the contact. No available
binder-design pipeline optimises this property, so it has to be added as a
stage of its own.

**Mouse cross-reactivity rules out the cetuximab epitope.** Cetuximab is
human-specific and does not recognise murine EGFR, precisely because its
domain III epitope is not conserved across species. Designing against that
patch — the natural reflex, given that the reference structure is a cetuximab
complex — works against objective 2. Structure 6ARU is used to understand the
domain III surface, not to copy its epitope.

For that reason the first analysis is not the Fab contact map but a
human/mouse conservation map crossed with solvent accessibility.

## Repository layout

```
EGFR_ProteinDesign/
├── README.md                  this file
├── requirements.txt           pinned dependencies
├── .gitignore
├── scripts/                   code, numbered in execution order
│   ├── 00_download_data.py    download and verification of raw data
│   ├── 01_annotate_target.py  per-residue annotation of the ectodomain
│   ├── 01b_verify_numbering.py  independent checks on 01's derived claims
│   ├── 02_select_epitope.py   candidate epitope patches and hotspot lists
│   ├── 02b_characterize_patches.py  patch geometry and face grouping
│   └── citation.py            manuscript provenance text from the manifest
├── 00_Competition/            challenge rules, FAQ, submission requirements
├── 01_Target/                 raw target data (not version-controlled)
│   ├── MANIFEST.json          provenance and checksums  ← version-controlled
│   ├── human/                 human EGFR: FASTA and UniProt entry
│   ├── mouse/                 mouse EGFR: FASTA and UniProt entry
│   └── structure/             6ARU in PDB and mmCIF format
├── 02_Analysis/               target analysis results
├── 03_Design/                 generated designs and their validation
├── 04_Submission/             final CSV and methods report
└── docs/                      notes, decisions, figures
```

Raw data is not version-controlled: it is reconstructed by the download script
and its integrity is guaranteed by the manifest checksums. What is
version-controlled is the code and the provenance, not the megabytes.

## Getting started

```bash
git clone <repository-url>
cd EGFR_ProteinDesign

python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt

python scripts/00_download_data.py
```

To confirm later that the data has not changed:

```bash
python scripts/00_download_data.py --verify
```

Developed on macOS (Apple Silicon) with Python 3.14. The analysis stages run
on CPU; the generative design stages require a CUDA GPU and are run on rented
hardware (see *Compute*).

## Pipeline

### `scripts/00_download_data.py` — raw data acquisition

Retrieves the six files the project starts from and records where each one
came from.

Sources are declared as a list of dictionaries at the top of the module, each
carrying the URL, the destination path, the database accession and a prose
description of why the file is needed. Adding a dataset means adding one entry
and nothing else.

For each source the script:

1. **Downloads** with up to three attempts, backing off on network errors.
2. **Validates the content** before trusting it. A failing server can return an
   HTML error page under status 200, so the script rejects anything that
   starts with `<`, requires FASTA files to begin with `>`, requires PDB files
   to carry a `HEADER` or `ATOM` record, requires mmCIF files to begin with
   `data_`, and parses JSON files to confirm they are well-formed.
3. **Computes a SHA-256 checksum**, reading in 1 MiB blocks.
4. **Captures the database release** from the HTTP response headers. UniProt
   exposes `X-UniProt-Release`, which pins the exact release the sequence came
   from — the single most useful field for reproducibility, and one that is
   lost if the file is downloaded by hand.
5. **Summarises the file** so the manifest can be audited without opening the
   data: sequence length for FASTA, atom counts per chain and the structure
   title for PDB, entry version and last sequence update for UniProt JSON.

Everything is written to `01_Target/MANIFEST.json`.

Two modes support later use. `--verify` recomputes every checksum and compares
it against the manifest, reporting missing or altered files and exiting
non-zero if any fail; this is what turns the manifest from documentation into
a test. `--force` re-downloads everything, which is how the manifest is
regenerated when a database release changes.

The script depends only on the Python standard library, so it runs before the
environment is built and cannot break because of a dependency upgrade.

### `scripts/citation.py` — provenance text for the manuscript

Renders the contents of the manifest as text ready to paste into a paper, in
three formats: a Materials and Methods paragraph (`--format paragraph`), a
Markdown provenance table for supplementary material (`--format table`), and
BibTeX entries for the database records (`--format bibtex`).

The reason this is a script and not a paragraph typed once is drift. Provenance
written by hand stops matching the data the moment a database release changes
and the download is repeated, and nothing flags the discrepancy until peer
review. Generating the text from the manifest makes that failure impossible.

```bash
python scripts/citation.py
python scripts/citation.py --format bibtex > docs/data_sources.bib
```

### `scripts/01_annotate_target.py` — per-residue annotation

Builds the table every downstream decision rests on. For each of the 621
ectodomain residues it crosses four layers:

1. **Conservation** — human versus mouse, from a global BLOSUM62 alignment,
   classified as identical, similar or different.
2. **Accessibility** — solvent-accessible surface area computed with FreeSASA
   from chain A of 6ARU, absolute and relative to the residue maximum. A
   residue is called exposed above 20% relative SASA.
3. **Glycan occlusion** — distance to the nearest glycan, using the sugar
   residues the crystal resolved and falling back to the ND2 atom of the
   glycosylated asparagine where it did not. Surface within 12 Å is treated as
   buried.
4. **Architecture** — region assignment from the disulfide pattern. See
   *Domain boundaries* below: this layer locates L-domain cores, which are
   narrower than the domains themselves.

It also computes the cetuximab interface by ΔSASA between the free receptor
and the complex. That measurement is used as an exclusion and as a positive
control that the geometry code works, not as a design target.

A residue is called **designable** when it is exposed, identical between
human and mouse, free of glycan shadow, and not disulfide-bonded.

```bash
python scripts/01_annotate_target.py
python scripts/01_annotate_target.py --glycan-radius 15 --exposure 0.25
```

Outputs `02_Analysis/01_residue_annotation.csv` (the full table),
`01_annotation_summary.md` (a readable report) and
`01_annotation_metadata.json` (parameters and provenance).

### `scripts/01b_verify_numbering.py` — verification

Step 01 rests on two derived claims: a UniProt-to-PDB offset obtained by
sequence alignment, and domain boundaries inferred from the disulfide
pattern. Both are load-bearing, and an error in either would shift every
downstream residue number without raising an exception. This script tests
both by routes that do not reuse the reasoning that produced them.

- **Test 0** — is the offset uniform across all mapped residues?
- **Test 1** — mapping UniProt's annotated disulfide pairs through the offset,
  do the SG–SG distances in the crystal come out at the 2.05 Å of a real
  covalent bond? This is geometric, not sequence-based, so it is independent
  of the alignment being checked.
- **Test 2** — do residues quoted in the EGFR literature have the expected
  identity in both numbering frames?
- **Test 3** — how much disulfide content do the competing domain III ranges
  contain?
- **Test 4** — what fraction of each range's residue contacts stay inside it?
  A genuine domain is compact; a range drawn across a boundary leaks contacts
  to its neighbours.

```bash
python scripts/01b_verify_numbering.py
```

Output is kept at `02_Analysis/01b_numbering_verification.txt`.

#### Domain boundaries

UniProt does not annotate the ectodomain subdomains I–IV. Its only `Domain`
feature is the intracellular kinase, at 712–979, and the two `Repeat`
features covering the L-domains, at 75–300 and 390–600, are flagged
*Approximate*.

Step 01 therefore infers regions from the disulfide distribution. Domains II
and IV are furin-like modules packed with disulfides; domains I and III are
leucine-rich L-domains almost free of them, so long gaps between consecutive
disulfide-bonded cysteines locate the L-domain cores:

| Range | Length | Character | Assignment |
|-------|--------|-----------|------------|
| 25–58 | 34 aa | cysteine-rich | I/II N-terminal |
| 59–156 | 98 aa | cysteine-free | I (L-domain core) |
| 157–362 | 206 aa | cysteine-rich | II (furin-like) |
| 363–469 | 107 aa | cysteine-free | III (L-domain **core**) |
| 470–645 | 176 aa | cysteine-rich | IV (furin-like) |

**These are cores, not domains, and the distinction matters.** Test 4
measures compactness directly, and the derived 363–469 range loses:

| Range for domain III | Internal contacts | External | Fraction inside |
|----------------------|-------------------|----------|-----------------|
| 363–469, from the disulfide gap | 495 | 167 | 0.748 |
| **335–538** | 1027 | **32** | **0.970** |
| 390–600, UniProt approximate | 1006 | 128 | 0.887 |

A short range scores well on compactness partly by being short, so a *longer*
range scoring higher cannot be a length artefact: it means the shorter one
was cutting through a structural unit. The 363–469 core leaks a quarter of
its contacts to residues outside itself, which a domain does not do.

The reading is that cysteine content and contact topology answer different
questions. The disulfide gap locates the leucine-rich solenoid core of
domain III, which is a real substructure; the domain also carries flanking
disulfide-bonded segments that pack against that core. **Domain III is taken
as 335–538 for epitope scoping**, and 363–469 is referred to as its core.

An earlier version of this file claimed 335–538 was unsupported. It is not;
Test 4 supports it over the range derived here.

## Findings from step 01

**The structure uses mature numbering.** The alignment gives an offset of
exactly +24 across all 609 resolved residues, so 6ARU residue *n* is UniProt
residue *n*+24. Nothing in the file announces this. Any analysis that assumed
the two numberings agreed would be wrong by 24 positions throughout, and
would still run without error.

**Glycans remove more surface than expected.** Thirteen N-glycosylation sites
sit in the ectodomain, two of them (N413, N444) inside the domain III core.
The crystal resolved only 28 sugar atoms, so most of that occlusion is
invisible in the structure and has to be reconstructed from the sequence
annotation.

**Only 168 of 621 residues are designable.** Counted over the domain ranges
that survived verification, with L-domain cores shown separately:

| Region | Residues | Designable | Fraction | Designable His |
|--------|----------|------------|----------|----------------|
| II (157–362, furin-like) | 206 | 74 | 35.9% | H183, H233, H304 |
| IV (470–645, furin-like) | 176 | 52 | 29.5% | H559, H584, H618 |
| **III (335–538)** | 204 | **45** | 22.1% | H433 |
| III core (363–469) | 107 | 20 | 18.7% | H433 |
| I core (59–156) | 98 | 9 | 9.2% | — |

Domain III carries 45 designable residues, which is enough surface for a
binder epitope. Domain II carries more, and is the dimerization arm, so
blocking it would be a functional mechanism rather than an epitope of
convenience — but the margin is not large enough on its own to justify
departing from the recommended epitope.

The 20-residue figure quoted before verification came from scoping domain III
to its core, and understated the available surface by more than half.

**The cetuximab epitope is more conserved than expected.** Of its 27
residues, 20 are identical between human and mouse, 7 are conservative
substitutions and none differ outright. Cetuximab nonetheless fails to
recognise murine EGFR, which is the useful lesson: seven substitutions across
a 27-residue footprint are enough to abolish antibody binding. For objective
2 the bar is therefore not "mostly conserved" but *zero* substitutions, which
is what the `designable` criterion enforces.

**Candidate histidines.** Exposed, conserved and unshadowed histidines are
worth cataloguing for objective 1, since a target histidine facing a binder
carboxylate adds a protonation-sensitive contact on top of the histidines
engineered into the binder. Domain III contributes H433 (relative SASA 0.80);
domain II contributes H183 (0.86), H304 (0.78) and H233 (0.33).

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

H370, described in the literature as a mechanistic pH-switch histidine, is
both buried and glycan-shadowed in this structure. L406 sits at 5% relative
accessibility, in the hydrophobic core.

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

## Target data

Downloaded 29 September 2026 from **UniProt release 2026_03** and the RCSB PDB.
SHA-256 checksums for every file are in `01_Target/MANIFEST.json`.

| Data | Accession | Detail |
|------|-----------|--------|
| Human EGFR | [P00533](https://www.uniprot.org/uniprotkb/P00533) | 1210 aa precursor, entry version 301 |
| Mouse EGFR | [Q01279](https://www.uniprot.org/uniprotkb/Q01279) | 1210 aa precursor, entry version 253 |
| Structure | [6ARU](https://www.rcsb.org/structure/6ARU) | EGFR ectodomain with a cetuximab Fab mutant |

In 6ARU chain **A** is EGFR (4686 atoms), chain **B** is the Fab light chain
and chain **C** the heavy chain. The challenge specifies the extracellular
region as the target, residues 25–645 in precursor numbering (621 aa);
residues 1–24 are the signal peptide and are not part of the mature protein.

### A note on residue numbering

This is the most common source of error in EGFR work. The literature alternates
between precursor numbering and mature-protein numbering, which differ by 24
positions. The same histidine can appear as H409 or H433 depending on the
paper, and a design targeted at the wrong one is silently wrong.

In this repository **every position is given in UniProt precursor numbering**
unless explicitly stated otherwise. Correspondences with PDB numbering are
derived by sequence alignment, never by assuming a constant offset — a
crystal structure with unresolved loops does not preserve one.

## Compute

Track 3 is self-supported: no competition Claude or Modal credits.

Sequence and structure analysis runs locally on CPU. Generative binder design
(AlphaFold2-based pipelines such as BindCraft) requires a CUDA GPU and roughly
5 GB of model weights, so it runs on rented hardware. The target is trimmed to
domain III before design to keep GPU memory and runtime tractable.

## Submission requirements

- **Deadline:** Sunday 4 October 2026, 23:59 AoE
- **Format:** CSV ordered by our own ranking, best design in the first row.
  Minimum columns: `name`, `sequence`, `molecule_class`
- **Length:** 10 to 250 amino acids, single chain
- **Count:** at most 20 designs (Track 3)
- **Constraint:** de novo and zero-shot. Starting from an existing binder or
  modifying one is not permitted

In Track 3 submissions are pooled and the selection for experimental
validation is made by a model from the information submitted, weighing
predicted design quality, design novelty and method novelty. The methods
report is not supporting documentation: it is part of what gets evaluated.

## Status

- [x] Target data downloaded and verified
- [x] Domain annotation and domain III boundaries
- [x] Human/mouse conservation map
- [x] Solvent accessibility and glycan occlusion
- [x] Numbering and domain boundaries verified independently
- [x] Epitope patch selection
- [ ] Binder generation
- [ ] pH-sensitivity engineering
- [ ] Filtering, ranking and submission

## Licence and publication

Experimental data and validated sequences produced by the competition are
published on Proteinbase under an ODC-BY licence. The code in this repository
is developed with a view to accompanying a publication.

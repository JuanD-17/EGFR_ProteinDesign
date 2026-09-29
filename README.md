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
4. **Architecture** — domain boundaries derived from the disulfide pattern.

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

#### Why domain boundaries are derived rather than looked up

UniProt does not annotate the ectodomain subdomains I–IV. Its only `Domain`
feature is the intracellular kinase, at 712–979. The two `Repeat` features
that cover the L-domains, at 75–300 and 390–600, are explicitly flagged
*Approximate* and do not agree with the boundaries commonly quoted in the
literature.

Domains are therefore inferred from the disulfide distribution. Domains II
and IV are furin-like modules packed with disulfides; domains I and III are
leucine-rich L-domains almost free of them. Long gaps between consecutive
disulfide-bonded cysteines locate the L-domain cores directly from the data:

| Range | Length | Character | Assignment |
|-------|--------|-----------|------------|
| 25–58 | 34 aa | cysteine-rich | I/II N-terminal |
| 59–156 | 98 aa | cysteine-free | I (L-domain core) |
| 157–362 | 206 aa | cysteine-rich | II (furin-like) |
| 363–469 | 107 aa | cysteine-free | **III (L-domain core)** |
| 470–645 | 176 aa | cysteine-rich | IV (furin-like) |

This places domain III at **363–469**, which matches neither the 335–538
range circulating in derived material nor UniProt's approximate 390–600.

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

**Only 168 of 621 residues are designable.** Broken down by domain:

| Domain | Residues | Designable | Fraction |
|--------|----------|------------|----------|
| II (furin-like) | 206 | 74 | 35.9% |
| IV (furin-like) | 176 | 52 | 29.5% |
| III (L-domain core) | 107 | 20 | 18.7% |
| I (L-domain core) | 98 | 9 | 9.2% |

Domain III, the epitope the competition recommends, has the *lowest*
designable density of the large domains. It is small, and two glycosylation
sites land in its core. Domain II offers roughly three times as much usable
conserved surface, and it is the dimerization arm — blocking it is a
functional mechanism rather than an epitope of convenience.

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
- [ ] Epitope patch selection
- [ ] Binder generation
- [ ] pH-sensitivity engineering
- [ ] Filtering, ranking and submission

## Licence and publication

Experimental data and validated sequences produced by the competition are
published on Proteinbase under an ODC-BY licence. The code in this repository
is developed with a view to accompanying a publication.

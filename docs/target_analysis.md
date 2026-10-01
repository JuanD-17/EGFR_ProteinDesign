# Target analysis

Definition and annotation of the EGFR extracellular region, and independent verification of the numbering every later step depends on.

See also [epitope_selection.md](epitope_selection.md).

---

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


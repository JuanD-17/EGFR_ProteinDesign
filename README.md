# Conditional EGFR binder design

Work for **Challenge 01 (EGFR)** of the Anthropic × Adaptyv 2026 protein
design competition, Track 3 (open track).

The goal is to design, *de novo*, a protein that binds the human epidermal
growth factor receptor while satisfying three conditions at once:

| # | Objective | Ranking weight |
|---|-----------|----------------|
| 1 | Selective binding at pH 6.5, no detectable binding at pH 7.4 | **highest** |
| 2 | Cross-reactivity with mouse EGFR | intermediate |
| 3 | Affinity for human EGFR | lowest |

That order is deliberate and counterintuitive. The competition rules state
that a weak but clearly pH-sensitive binder ranks above a high-affinity
binder that is not pH-sensitive: demonstrated difficulty is rewarded over
the magnitude of the dissociation constant.

## What has been established

| | Result | Detail |
|---|--------|--------|
| Target | 621-residue ectodomain annotated on four layers, 168 residues designable | [target_analysis](docs/target_analysis.md) |
| Numbering | 6ARU runs 24 below UniProt; verified geometrically on 24 of 25 disulfides | [target_analysis](docs/target_analysis.md) |
| Epitope | Patch 2 of domain III, the only top-ranked patch forming one connected surface | [epitope_selection](docs/epitope_selection.md) |
| Mechanism | One histidine is worth ~8-fold; two to three coupled protonations are needed | [ph_mechanism](docs/ph_mechanism.md) |
| Handles | H433 at pKa 6.22 titrates almost ideally; D458, D460, E455 are stable anchors | [ph_mechanism](docs/ph_mechanism.md) |
| Generation | 59 trajectories, 44 completed, 8 accepted; 77% of rejections on AF2 confidence | [design_campaign](docs/design_campaign.md) |
| Engineering | 16 variants proposed, 11 geometrically feasible | [validation](docs/validation.md) |
| Control | Designed placement 1.82 kcal/mol against 0.40 for random, 5 of 8 at p < 0.05 | [validation](docs/validation.md) |
| Revalidation | 6 of 9 networked designs retain interface confidence within 0.05 i_pTM | [validation](docs/validation.md) |
| Novelty | No hit against any polymer entity in the PDB; 27% maximum internal identity | [validation](docs/validation.md) |

**Lead candidate.** A 68-residue binder engaging D368, H370 and H433, with
an estimated coupling of 2.38 kcal/mol, p = 0.0033 against random
placement, all three rotamers geometrically feasible, and interface
confidence unchanged from its parent (Δ i_pTM = +0.002).

**What is not established.** That any of these bind at all, or bind
pH-selectively. Every number above is computed. The chain this project
follows, and does not short-circuit, is:

```
computed  →  structural or mechanistic hypothesis  →  experimental validation
```

The last step belongs to the assay.

## Design rationale

Two constraints follow from the ranking order and shape everything
downstream.

**pH sensitivity is engineered into the binder, not discovered in the
target.** EGFR cannot be mutated. The available mechanism is to place
histidines in the binder's interface, where the imidazole pKa near 6.5
means protonation in the acidic tumour microenvironment and neutrality at
physiological pH. No binder-design pipeline optimises this property, so it
is added as a stage of its own.

**Cross-species reactivity is handled by an absolute criterion.** A residue
enters the designable pool only if it is *identical* between human and
mouse, not merely conserved. The bar looks arbitrarily strict until one
notes that the cetuximab epitope is 20 of 27 residues identical with none
differing outright, and cetuximab still fails to recognise murine EGFR. A
handful of conservative substitutions across a footprint is enough to
abolish antibody binding.

Structure 6ARU is used to understand the domain III surface, not to copy
the cetuximab epitope.

## Repository layout

```
EGFR_ProteinDesign/
├── README.md                  this file
├── LICENSES.md                third-party terms and outstanding actions
├── requirements.txt           pinned dependencies
├── docs/                      detailed results, one file per stage
│   ├── target_analysis.md
│   ├── epitope_selection.md
│   ├── ph_mechanism.md
│   ├── design_campaign.md
│   └── validation.md
├── scripts/                   code, numbered in execution order
├── 01_Target/                 raw target data (not tracked)
│   └── MANIFEST.json          provenance and checksums  ← tracked
├── 02_Analysis/               target analysis results
├── 03_Design/                 designs, variants, controls, validation
│   ├── target/                the trimmed domain III given to BindCraft
│   ├── accepted/              the eight generated binders
│   └── validation/predictions/  re-predicted variants and their parents
└── 04_Submission/             submission CSV and methods report
```

Raw target data is not tracked: the download script reconstructs it and
`MANIFEST.json` pins it by checksum. Of the model output, what is tracked is
the coordinates behind a reported number — the designs, the re-predicted
complexes, the run logs. The 213 MB of raw BindCraft trajectories that
produced them is not; `03_Design/rejection_analysis/` is what was learned
from it.

## Getting started

```bash
git clone https://github.com/JuanD-17/EGFR_ProteinDesign.git
cd EGFR_ProteinDesign
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python scripts/00_download_data.py
```

To confirm later that the inputs have not changed:

```bash
python scripts/00_download_data.py --verify
```

Developed on macOS (Apple Silicon) with Python 3.14. Analysis runs on CPU in
minutes; binder generation and revalidation need a CUDA GPU.

## Running the pipeline

Scripts are numbered in execution order and each requires the outputs of
those before it. Arguments that vary between runs are required rather than
defaulted, so a stale input directory cannot be analysed by accident.

**Target analysis** — CPU, a few minutes, no arguments.

```bash
python scripts/00_download_data.py
python scripts/01_annotate_target.py
python scripts/01b_verify_numbering.py
python scripts/02_select_epitope.py
python scripts/02b_characterize_patches.py
python scripts/03_ph_mechanism.py --patches 1 2 3 4 5
python scripts/03b_protonation_networks.py
python scripts/04_prepare_target.py
```

**Binder generation** — CUDA GPU, run separately. Consumes
`03_Design/target/` and returns accepted complexes, placed in
`03_Design/accepted/`.

**pH engineering, controls and submission** — CPU, seconds each.

```bash
python scripts/05_engineer_ph_switch.py --designs 03_Design/accepted/
python scripts/06_negative_control.py --samples 3000
python scripts/07_validate_ph_variants.py --prepare
# GPU: re-predict the complexes written to 03_Design/validation/
python scripts/07_validate_ph_variants.py --score
python scripts/08_analyze_rejections.py --runs 03_Design/run*/gpu*
python scripts/09_check_novelty.py
python scripts/10_build_submission.py
```

Steps 05 to 07 are rerun as a group whenever new designs arrive, since each
depends on the previous one's output. `scripts/citation.py` renders the
input provenance as manuscript text, in prose, table or BibTeX.

## Target data

Retrieved 29 September 2026 from **UniProt release 2026_03** and the RCSB
PDB. SHA-256 checksums for every file are in `01_Target/MANIFEST.json`.

| Data | Accession | Detail |
|------|-----------|--------|
| Human EGFR | [P00533](https://www.uniprot.org/uniprotkb/P00533) | 1210 aa precursor, entry version 301 |
| Mouse EGFR | [Q01279](https://www.uniprot.org/uniprotkb/Q01279) | 1210 aa precursor, entry version 253 |
| Structure | [6ARU](https://www.rcsb.org/structure/6ARU) | EGFR ectodomain with a cetuximab Fab mutant |

In 6ARU chain **A** is EGFR, **B** the Fab light chain and **C** the heavy
chain. The challenge specifies residues 25–645 in precursor numbering;
1–24 are the signal peptide.

### A note on residue numbering

This is the most common source of error in EGFR work. The literature
alternates between precursor and mature numbering, which differ by 24
positions, so the same histidine appears as H409 or H433 depending on the
paper and a design aimed at the wrong one is silently wrong.

**Every position in this repository is in UniProt precursor numbering**
unless stated otherwise. Correspondences with PDB numbering are derived by
alignment and verified geometrically, never assumed — see
[target_analysis](docs/target_analysis.md).

## Compute

Track 3 is self-supported: no competition credits. Sequence and structure
analysis runs locally on CPU. Binder generation used BindCraft on two
Tesla T4 cards within a 30 GPU-hour weekly allowance, which is why
optimisation iterations were reduced from the published defaults. The
target was trimmed to domain III, 204 residues, to fit in 15 GB of device
memory. Validation settings were left untouched: weakening the stage that
decides whether a design is good would only move the failure downstream.

## Submission

- **Deadline:** 4 October 2026, 23:59 AoE
- **Format:** CSV ordered by our own ranking, columns `name`, `sequence`,
  `molecule_class`
- **Constraints:** 10–250 residues, single chain, at most 20 designs, de
  novo and zero-shot

`04_Submission/` holds the CSV, a metadata table with every metric behind
the ranking, and `METHODS.md`. The methods report is generated from the
analysis outputs rather than transcribed, so its figures cannot drift from
the data.

In Track 3 submissions are pooled and selected by a model from the
information submitted, weighing predicted design quality, design novelty
and method novelty. The methods report is part of what is evaluated.

## Status

- [x] Target data, annotation and numbering verification
- [x] Epitope selection and patch geometry
- [x] pH mechanism, pKa estimates and protonation networks
- [x] Binder generation, two batches
- [x] Histidine engineering, negative control, rotamer feasibility
- [x] Structural revalidation against unmodified parents
- [x] Novelty check and submission build
- [ ] Submitted
- [ ] Experimental validation (competition)

## Licence and publication

Experimental data and validated sequences from the competition are
published on Proteinbase under ODC-BY. Third-party terms, including the
non-commercial licences on PyRosetta and the AlphaFold2 parameters, are
recorded in [LICENSES.md](LICENSES.md) along with the actions still
outstanding.

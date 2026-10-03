# Licences and data provenance

Record of every third-party tool, model and dataset this project depends on,
with its licence and what that licence permits.

**This file is documentation, not compliance.** Writing a licence down does
not grant it. Two dependencies below carry non-commercial terms and require
action, marked ⚠ — see *Outstanding actions*.

The code and documentation written for this project are released under the
MIT licence, in [LICENSE](LICENSE). That covers what is ours; everything
below belongs to someone else and keeps its own terms.

Prepared by Juan David Hoyos Trejos, Universidad Icesi, for the Anthropic ×
Adaptyv 2026 protein design competition, Challenge 01. Last updated
29 September 2026.

---

## Non-commercial dependencies requiring attention

### ⚠ PyRosetta

| | |
|---|---|
| Licence | Free for academic and non-commercial use; commercial use requires purchase |
| Licensor | Rosetta Commons / University of Washington (UW CoMotion) |
| Contact | `license@uw.edu` |
| Build in use | `PyRosetta4.Release.python310.ubuntu.cxx11thread.serialization` |
| Version | `2026.29+release.quarterly.80a0635615099e1b918474a63acba7b1de6fd107` |
| Build date | 2026-07-14 |
| Installed | 29 September 2026, via the BindCraft installer, in a Kaggle notebook |
| Purpose | Interface relaxation and scoring inside the BindCraft pipeline |

The banner printed on import states: *"USE OF PyRosetta FOR COMMERCIAL
PURPOSES REQUIRES PURCHASE OF A LICENSE. See LICENSE.PyRosetta.md or email
license@uw.edu for details."*

Citation: Chaudhury S, Lyskov S, Gray JJ (2010). PyRosetta: a script-based
interface for implementing molecular modeling algorithms using Rosetta.
*Bioinformatics* 26(5):689–691.

### ⚠ AlphaFold2 model parameters

| | |
|---|---|
| Licence | CC BY-NC 4.0 — attribution required, **non-commercial only** |
| Source | `storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar` |
| Size | 5,587,968,000 bytes |
| Downloaded | 29 September 2026 |
| Files | 15 `.npz`, including `params_model_{1..5}_multimer_v3.npz` |
| Purpose | Backbone hallucination and complex validation inside BindCraft |

The *parameters* and the *code* are licensed differently: the AlphaFold2
source is Apache 2.0, but the trained weights are CC BY-NC 4.0. It is the
weights this project uses, so the non-commercial term applies.

Citation: Jumper J, Evans R, Pritzel A, et al. (2021). Highly accurate
protein structure prediction with AlphaFold. *Nature* 596:583–589.
For the multimer weights: Evans R, O'Neill M, Pritzel A, et al. (2021).
Protein complex prediction with AlphaFold-Multimer. *bioRxiv*.

---

## Permissive dependencies

| Tool | Version | Licence | Role |
|------|---------|---------|------|
| BindCraft | commit from `martinpacesa/BindCraft`, cloned 29 Sep 2026 | MIT | Binder design pipeline |
| ProteinMPNN | as vendored by BindCraft | MIT | Sequence design on generated backbones |
| ColabDesign | as vendored by BindCraft | Apache 2.0 | AlphaFold2 design interface |
| JAX / jaxlib | 0.6.0 | Apache 2.0 | Numerical backend, CUDA |
| Biopython | 1.88 | Biopython Licence (BSD-style) | Structure and sequence parsing |
| NumPy | 2.5.3 | BSD 3-Clause | Numerics |
| pandas | 3.0.6 | BSD 3-Clause | Tabular analysis |
| matplotlib | 3.11.2 | Python Software Foundation Licence | Figures |
| FreeSASA | 2.2.1 | MIT | Solvent-accessible surface area |
| PROPKA | 3.5.1 | **LGPL v2.1** | pKa estimation |

Licences for the Python packages were read from installed package metadata
rather than from memory. BindCraft, ProteinMPNN and ColabDesign run on the
GPU host and their licences are recorded from their repositories; confirm
against the `LICENSE` file of the exact commit cloned.

PROPKA is LGPL v2.1, which is weaker copyleft than GPL: using it as a library
without modifying it, as this project does, does not impose LGPL terms on
this project's own code. The compatibility shim in
`scripts/03b_protonation_networks.py` patches PROPKA's behaviour at runtime
without modifying or redistributing its source.

Citation: Olsson MHM, Søndergaard CR, Rostkowski M, Jensen JH (2011).
PROPKA3: consistent treatment of internal and surface residues in empirical
pKa predictions. *J Chem Theory Comput* 7(2):525–537.

---

## Data

| Source | Identifier | Licence | Retrieved |
|--------|-----------|---------|-----------|
| UniProtKB | P00533 (human EGFR), entry version 301 | CC BY 4.0 | 29 Sep 2026, release 2026_03 |
| UniProtKB | Q01279 (mouse EGFR), entry version 253 | CC BY 4.0 | 29 Sep 2026, release 2026_03 |
| RCSB PDB | 6ARU | Public domain (CC0) | 29 Sep 2026 |

SHA-256 checksums for every input file are in `01_Target/MANIFEST.json` and
can be re-verified with `python scripts/00_download_data.py --verify`.

---

## Competition terms

Designs submitted to the competition, and any experimental data generated
from them, are published on Proteinbase under **ODC-BY**. Submitting is
therefore a decision to release those sequences under an open licence.

---

## Outstanding actions

Documentation is not a licence. Two things remain to be done, and both are
cheap now and awkward later.

**1. Confirm the PyRosetta licence.** Being able to download PyRosetta does
not establish that a licence has been granted. Check whether Universidad
Icesi already holds an institutional Rosetta/PyRosetta licence — many
universities do — and if not, register as an academic user. Academic
licences are free to degree-granting institutions.

**2. Resolve the non-commercial question for this specific context.** Both
PyRosetta and the AlphaFold2 weights are restricted to non-commercial use.
This work is academic and its outputs are published openly, which reads as
non-commercial. But the competition is co-sponsored by companies, and its
stated aim includes advancing commercially relevant design capability. That
combination is ambiguous enough to be worth a written answer rather than an
interpretation.

Emailing `license@uw.edu` with the context — academic researcher, university
affiliation, open publication of results, competition sponsored by
industry — obtains that answer in writing. The same question applies to the
AlphaFold2 weights under CC BY-NC 4.0 and should be raised with Google
DeepMind if certainty is needed there too.

Neither action blocks the current work. Both should be settled before
publication.

*This file records licence terms as understood by the project. It is not
legal advice.*

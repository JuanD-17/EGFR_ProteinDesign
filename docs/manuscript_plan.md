# Manuscript plan

Working plan for a publication built on this repository. Written 3 October
2026, the day the designs were submitted and before any experimental result
exists.

---

## 1. What the results actually support

A paper can claim only what its data carries. Sorted by how well supported
each statement is.

### Solid

**The numbering verification.** The 24-residue offset between 6ARU and
UniProt was confirmed by a route independent of the alignment that produced
it: 24 of 25 annotated disulfide pairs land at 2.03–2.05 Å SG–SG. This is
geometry, not inference, and it is reproducible from the deposited structure.

**The reach correction.** Across the 23 histidines resolved in 6ARU, CB–ND1
is 2.53 Å and CB–NE2 is 3.67 Å, giving a ceiling near 7.7 Å for a
charge-assisted contact. An earlier 10 Å bound vetoed 33 of 33 proposals,
every one for failing to reach. This is a measurement with a clear
consequence, and it is the most transferable result in the project: any
pipeline placing titratable residues by distance needs this bound.

**The negative control design.** Drawing the null from the same constrained
pool — interface, exposed, substitutable, with only the targeting removed —
is the right control and is rarely run in this literature. 9 of 16 variants
reach p < 0.05; none exceeds the maximum of its own null.

**Composition pathology.** Glutamate between 14% and 28% against a natural
7%, worst in the lead design (27.9%, zero aspartate). Total variation
distance from UniProtKB frequencies 0.23–0.38. This is a measured property
of hallucinated sequences, reported rather than hidden.

### Supported but model-dependent

**Coupling estimates.** 2.38 kcal/mol for the lead, 55-fold. These follow
correctly from Wyman linkage *given* an assumed 2.0-unit pKa shift on
binding. The thermodynamics are right; the input is a stipulation.

**Interface survival.** Δ i_pTM of +0.002 for the lead against its parent,
and the finding that damage does not track substitution count (R49H alone
collapsed a complex from 0.819 to 0.243). The parent comparison is the right
design, but both arms come from the same predictor.

### Not supported

That any design binds EGFR. That any binds pH-selectively. That any folds.
That any expresses. Every number in this project is computed, and the
composition data suggest expression is a live risk.

---

## 2. Journal: an honest assessment of PEDS

**Protein Engineering, Design and Selection** publishes protein engineering
with a strong experimental tradition. A manuscript reporting 15 designs with
no wet-lab data, whose central claim rests on an assumed pKa shift, is a hard
sell there as it stands. The likely referee objection writes itself: *the
authors assume the quantity that determines their result.*

That is not a reason to abandon PEDS. It is a reason to choose the moment.

### Three paths

**A. Preprint now, journal after the competition.** Post to bioRxiv in the
next weeks with explicit, falsifiable predictions. Adaptyv results arrive
around December 2026. If any design binds, the paper becomes a different
object and PEDS becomes straightforward. If none binds, the preprint is still
a published prediction that was tested — which is worth more than silence.

**B. Replace the assumption with a computation.** Constant-pH molecular
dynamics on the lead complex would turn the assumed Δ pKa into a computed
one. This is the single change that would make the paper publishable at PEDS
*without* experimental data, because it removes the objection above. It is
substantial work — weeks, and GPU time this project has not had — but it is
the scientifically correct next step regardless of the competition.

**C. Submit now to a computational venue.** The work would fit *Journal of
Chemical Information and Modeling*, *Bioinformatics Advances*, or *Frontiers
in Bioinformatics* as a methods contribution. Lower profile, faster, and it
forecloses the stronger paper that A or B would produce.

**Recommendation: A, with B started in parallel.** The preprint costs nothing
and establishes priority on the reach correction and the control design,
which are the parts someone else could independently publish. The CpHMD makes
the paper stand on its own if the competition yields nothing.

### What makes this publishable at all

Not the designs. Fifteen unvalidated binders are not a contribution; there
are thousands. The contribution is methodological, and it is real:

1. A geometric bound on titratable-residue placement, derived by measurement,
   that invalidated the project's own earlier result (86-fold → 55-fold).
2. A matched null model for switch placement, with the negative result
   reported (no design exceeds its null maximum; p-values uncorrected).
3. A worked demonstration that interface damage from a substitution does not
   scale with substitution count.

A paper framed around those three, using the EGFR campaign as the vehicle, is
defensible. A paper framed around "we designed pH-selective EGFR binders" is
not, until something binds.

---

## 3. Title

Framed on the method, which is what the data support:

1. **Measuring the reach of a histidine: a geometric bound and a matched null
   model for placing pH switches on de novo binder interfaces**
2. **pH-conditional binders against EGFR domain III: histidine network
   placement tested against its own null distribution**
3. **A null-model test for engineered protonation switches in computationally
   designed protein binders**

Option 1 leads with the transferable result and names both contributions.
Option 2 is the conventional target-first framing and will attract the EGFR
readership, at the cost of promising a binder the data do not demonstrate.
Option 3 is the most honest and the least likely to be read.

Prefer 1. If the competition returns a validated binder, retitle around it.

---

## 4. Materials and Methods

Most of this already exists in `04_Submission/METHODS.md`, generated from the
analysis outputs. Structure for the manuscript:

**2.1 Target definition and annotation.** UniProt P00533 (v301) and Q01279
(v253), release 2026_03; PDB 6ARU. Four-layer annotation: human/mouse
identity from global BLOSUM62 alignment, FreeSASA relative accessibility,
N-glycan occlusion at 12 Å, disulfide architecture. Designability criteria
and the 168/621 outcome. State the absolute species criterion and the
cetuximab argument for it.

**2.2 Numbering verification.** The disulfide-distance test. Report all 25
pairs including the one outlier.

**2.3 Epitope selection.** Patch construction at 11 Å, the four-term weighted
score, geometric characterisation, and the removal of residues detached from
the patch face. State that H433, Q435 and K489 were recovered without using
the literature as input — this is a useful internal validation.

**2.4 Protonation analysis.** PROPKA 3.5.1 free-state pKa. Wyman linkage:
∂ΔG/∂pH = 2.303·R·T·Δn, giving 1.23 kcal/mol per coupled protonation across
pH 7.4 → 6.5. Be explicit that selectivity arises from the pKa *shift* on
binding, not the pKa, and that the 2.0-unit shift is assumed.

**2.5 Binder generation.** BindCraft against the 204-residue domain III trim;
AF2 multimer v3 backpropagation, ProteinMPNN, PyRosetta scoring. Reduced
iterations and why. Report the campaign honestly: 59 trajectories, 8 accepted
(13.6%), 483 of 628 rejections on AF2 confidence.

**2.6 Histidine placement.** The 7.7 Å bound and its derivation from 23
histidines in 6ARU. Exposure threshold 0.20 relative SASA. Protected
residues. Substitutions proposed in sets because count is the governing
variable.

**2.7 Rotamer feasibility.** 200 orientations on spheres of both measured
radii; clash rejection; acceptance at 2.6–4.0 Å from the handle. State that
the environment excludes the side chain being replaced, and why an earlier
version that did not vetoed everything.

**2.8 Negative control.** 3000 draws per design and substitution count from
the matched pool. Define p as the fraction of draws reaching the observed
coupling, with +1 correction. State the pool sizes (22–29 positions, 7
handles).

**2.9 Structural revalidation.** Variant and parent re-predicted under
identical settings. Three layers reported separately — folding plausibility,
interface quality, network survival — and never merged.

**2.10 Novelty.** RCSB sequence service over all polymer entities.
Composition deviation as total variation distance from UniProtKB means.
Pairwise identity within the portfolio. State that no structural comparison
was run, and report that Proteinbase's own check scored one design 2/4.

**2.11 Data and code availability.** Repository URL, SHA-256 manifest,
licence. Note the non-commercial terms on PyRosetta and the AF2 parameters.

---

## 5. Results and Discussion

Order the results so each section answers the objection raised by the last.

**3.1 A designable surface on domain III.** 168 of 621 residues qualify;
patch 2 is the only top-ranked patch forming a single connected surface. The
recovery of H433, Q435 and K489 without literature input. → Figure 1.

**3.2 The numbering is verified, not assumed.** 24 of 25 disulfides at
2.03–2.05 Å. Short, but it belongs in Results because it is a measurement and
because the field gets this wrong. → Figure S1.

**3.3 One histidine is not enough.** The linkage ceiling, 1.23 kcal/mol per
protonation, ~8-fold for one. This reframes the design problem from "find a
good position" to "find several". → Figure 2.

**3.4 A histidine reaches 7.7 Å, not 10.** The measurement, the 33/33 veto
that exposed the error, and the consequence: the best variant fell from
86-fold to 55-fold. Discuss this as a general hazard — a distance criterion
chosen for convenience silently admits contacts no rotamer can form. → Figure 3.

**3.5 Targeted placement beats chance, with limits.** Null distributions and
observed values. 9 of 16 at p < 0.05, mean 1.82 vs 0.40. Then the two
honest caveats: no variant exceeds its null maximum, and the p-values are
uncorrected because variants within a design share positions. → Figure 4,
the central figure of the paper.

**3.6 Damage does not track substitution count.** The lead engages three
handles at Δ i_pTM +0.002; R49H alone collapses its complex to 0.243.
Discuss: what matters is which residue is replaced, and the cheapest
variants are those built on histidines the design already carried. That is a
design principle, and it generalises. → Figure 5.

**3.7 One network was not engineered.** H7/H11 arrived from BindCraft, which
has no notion of pH. Reported as incidental rather than counted as a result.
Discuss briefly what it implies about how often usable protonation geometry
appears by chance on designed interfaces — with n = 1, as a question, not a
rate.

**3.8 The portfolio is compositionally abnormal.** E 14–28% against 7%
natural, worst in the lead. Discuss as an expression liability and as a
general property of hallucinated sequences, not of this campaign. Connect to
the design that failed the platform's novelty check at 2/4 — sequence-level
novelty passed where a structural check did not.

**3.9 Limitations.** Give this its own numbered section, not a paragraph.
The assumed pKa shift. The AF2→AF2 circularity. n = 8 accepted designs.
No constant-pH MD. No structural novelty comparison. No experimental data.

---

## 6. Figures

| # | Content | Source |
|---|---------|--------|
| 1 | Domain III surface: four annotation layers, patch 2 outlined, cetuximab footprint overlaid | `02_Analysis/01_residue_annotation.csv`, `02_patches.csv`, 6ARU |
| 2 | Fold-change vs number of coupled protonations, 1–4, with the 8× and 55× points marked | Wyman linkage, analytic |
| 3 | Measured CB–ND1/CB–NE2 distances across 23 histidines; the 7.7 Å ceiling; veto rate before and after | `07_rotamer_check.csv`, 6ARU |
| 4 | **Central figure.** Null distributions per design with observed coupling marked; p annotated | `06_null_distribution.csv` |
| 5 | Δ i_pTM against substitution count, points labelled; R49H as the outlier | `07_variant_scores.csv` |
| 6 | Lead candidate: the D368·H370·H433 network, with distances | `pred_..._H7_M10H_H11.pdb` |
| 7 | Amino-acid composition per design against UniProtKB means | `09_novelty_report.md` |
| S1 | Disulfide SG–SG distances, 25 pairs, with the 2.05 Å reference | `01b_numbering_verification.txt` |
| S2 | Campaign funnel: 59 trajectories → 44 completed → 8 accepted, with rejection filters | `08_rejection_by_filter.csv` |

**Tables.** T1: the 15 submitted designs with every metric. T2: the
protonation handles on patch 2 with PROPKA values and behaviour.

---

## 7. The Adaptyv possibility

The designs were submitted 3 October 2026. Around 375 designs per track and
challenge go to the wet lab; results are published on Proteinbase including
failures.

**Do not write the paper contingent on being selected.** A manuscript whose
conclusion waits on someone else's decision cannot be finished, and selection
is not in your control.

**Do pre-register the prediction.** This is the strongest move available and
costs nothing. In the preprint, state before any result exists:

> Fifteen designs were submitted to the Anthropic × Adaptyv 2026 competition
> (Challenge 01, Track 3) on 3 October 2026. If any are selected, binding
> will be measured against human and mouse EGFR at pH 7.4 and 6.5. We predict
> that `EGFRd3_01_H7M10HH11_l68` shows the largest pH-dependent difference of
> the set, and that designs ranked below 9 show none. We further predict that
> expression failure, if it occurs, correlates with acidic fraction.

A published prediction that is later tested is a complete scientific act
whether it is confirmed or refuted. A refuted prediction with a declared
mechanism is more useful than no prediction, and it protects you from the
charge of reporting only what worked.

**Three scenarios.**

*Selected and something binds.* The paper is rewritten around the measurement
and goes to PEDS or higher. The computational work becomes the method section
of a validated result.

*Selected and nothing binds.* Still publishable, and more interesting than it
sounds: a pre-registered prediction, a measured refutation, and a pipeline
whose failure mode is identified. State which of the stated limitations the
failure is consistent with — composition, the assumed pKa shift, or the
AF2 circularity.

*Not selected.* Path B becomes the plan: add constant-pH MD and publish the
method on its own.

In all three the preprint is already out and dated.

---

## 8. Immediate next steps

1. Draft the preprint around Title 1, Results order in §5, Figures 1–7.
2. Make Figure 4 first. If the null-distribution figure is not convincing,
   the paper has no core and that should be known early.
3. Start constant-pH MD on the lead complex. It is the one addition that
   makes the work independent of the competition.
4. Resolve the PyRosetta licence with `license@uw.edu` before submission.
5. Reconcile `METHODS.md` (16 designs) with what was submitted (15).

---

## 9. Venues

The right journal depends on what exists when the manuscript is written, so
the list is organised by scenario rather than by prestige. APC is stated
because it decides more submissions than it should.

### Always, first, regardless of scenario

**bioRxiv** — free, immediate, citable, and it timestamps the reach
correction and the null-model design. Not a journal and not a substitute for
one. Post before any journal submission; every venue below permits it.

### Scenario A — a design binds

| Venue | Fit | APC |
|---|---|---|
| **Protein Science** (Wiley) | Protein engineering with computation; takes de novo design routinely | Hybrid — free via subscription route |
| **PEDS** (OUP) | The natural home once there is a measurement | Hybrid — free via subscription route |
| **ACS Synthetic Biology** | If framed as a conditional-activity module rather than an EGFR paper | Hybrid |
| *Nature Communications / Science Advances* | Only if the switch works cleanly and cross-species holds. A functioning de novo conditional binder is genuinely high-profile; the bar is a clean dose–response at both pH values | ~$6,000+ |

Start at Protein Science. PEDS is the fallback, not the ceiling.

### Scenario B — computational only, with constant-pH MD added

This is the realistic target if the competition yields nothing.

| Venue | Fit | APC |
|---|---|---|
| **J. Chemical Information and Modeling** (ACS) | Strong. A measured geometric bound plus a matched null is exactly their methods profile | Hybrid — free via subscription route |
| **Proteins: Structure, Function, Bioinformatics** (Wiley) | Classic fit; long history of pKa and electrostatics work | Hybrid — free route |
| **PEDS** | Works once the pKa shift is computed rather than assumed | Hybrid — free route |
| **PLOS Computational Biology** | Higher bar; needs the method framed as general, not EGFR-specific | ~$2,800, LMIC waivers available |

Start at JCIM.

### Scenario C — computational only, as it stands

Possible, but it spends the result for less than it is worth. Listed for
completeness.

**Scientific Reports** (~$2,790), **Frontiers in Bioinformatics** (~$2,000+),
**PeerJ** (~$1,500), **BMC Bioinformatics**. All open access with an APC, all
broad-scope. Note the inversion: the journals that are easiest to enter are
the ones that cost money, while the selective hybrids above can be published
free through the subscription route.

### On cost

Check whether Universidad Icesi holds a transformative agreement with Wiley,
ACS, OUP or Springer Nature before choosing. Many Colombian institutions do,
and it changes which venues are reachable. The library will know. Ask before
submitting, not after acceptance.

PLOS, Frontiers and some others operate income-based waivers for authors in
lower- and middle-income countries; these are applied for at submission.

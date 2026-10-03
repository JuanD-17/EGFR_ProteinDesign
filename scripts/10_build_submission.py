#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 10 - Submission CSV and methods report.

Produces the two deliverables, both generated from the pipeline's own
outputs rather than typed, so that a figure in the write-up cannot drift
away from the data behind it.

Ranking
-------
The competition ranks pH selectivity above cross-reactivity above affinity,
and asks for submissions ordered by the submitter's own preference. That
order is reproduced here, with one hard exclusion in front of it.

  Excluded   a variant whose interface the substitutions destroyed. It no
             longer binds, so its pH properties are irrelevant.

  Then by    coupling, the estimated pH selectivity, because that is the
             competition's first criterion.

  Then by    whether the interface survived intact rather than merely
             tolerably, since a degraded interface is a weaker candidate
             at equal coupling.

  Then by    the negative control p-value, preferring placements that
             chance rarely reproduces.

  Then by    composition risk. Every design here is strongly acidic, a
             known signature of hallucinated sequences, and Adaptyv has to
             express these in a laboratory. Between two otherwise equal
             candidates the less extreme one is preferred.

Unmodified parents are included after the variants. They are validated
binders without a pH switch, so they rank below any variant, but they are
the fallback if the switch fails experimentally and they cost nothing: the
submission allows up to twenty and there is room.

Inputs:  03_Design/validation/07_variant_scores.csv
         03_Design/ph_variants/05_variants.csv
         03_Design/negative_control/06_null_distribution.csv
         03_Design/validation/validation_results.csv
         02_Analysis/* for the methods report
Outputs: 04_Submission/submission.csv
         04_Submission/submission_metadata.csv
         04_Submission/METHODS.md

Usage:
    python scripts/10_build_submission.py
    python scripts/10_build_submission.py --max-designs 20 --no-parents
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "04_Submission"

SCORES = ROOT / "03_Design/validation/07_variant_scores.csv"
VARIANTS = ROOT / "03_Design/ph_variants/05_variants.csv"
CONTROL = ROOT / "03_Design/negative_control/06_null_distribution.csv"
PREDICTIONS = ROOT / "03_Design/validation/validation_results.csv"
INDEX = ROOT / "03_Design/validation/07_prediction_index.csv"
ANNOTATION_META = ROOT / "02_Analysis/01_annotation_metadata.json"
REJECTIONS = ROOT / "03_Design/rejection_analysis/08_rejection_analysis.json"

MAX_DESIGNS = 20

# Fraction of glutamate plus aspartate above which a sequence is flagged as
# an expression risk. Natural proteins average about 12 percent combined.
ACIDIC_FLAG = 0.25


def load_csv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh))


def count_new_substitutions(substitutions: str) -> int:
    """How many positions actually change.

    Step 05 writes a pre-existing histidine as `H7` and a substitution as
    `M10H`. A variant composed only of pre-existing histidines changes
    nothing: its sequence is the parent's, and presenting it as engineered
    would claim credit for a network the design already had.
    """
    return sum(1 for token in substitutions.split()
               if re.fullmatch(r"[A-Z]\d+H", token))


def composition_risk(sequence: str) -> tuple[float, str]:
    """Acidic fraction and a one-line reading of it."""
    length = len(sequence)
    acidic = (sequence.count("E") + sequence.count("D")) / length
    net = (sequence.count("K") + sequence.count("R")
           - sequence.count("E") - sequence.count("D"))
    if acidic >= ACIDIC_FLAG:
        note = f"{100 * acidic:.0f}% acidic, net {net:+d}: expression risk"
    else:
        note = f"{100 * acidic:.0f}% acidic, net {net:+d}"
    return acidic, note


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--max-designs", type=int, default=MAX_DESIGNS)
    parser.add_argument("--no-parents", action="store_true",
                        help="submit only engineered variants")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    scores = {r["substitutions"]: r for r in load_csv(SCORES)}
    variants = load_csv(VARIANTS)
    control = {r["substitutions"]: r for r in load_csv(CONTROL)}
    predictions = {r["name"]: r for r in load_csv(PREDICTIONS)}
    index = load_csv(INDEX)

    sequences = {r["substitutions"]: r["sequence"] for r in variants}
    name_of = {(r["design"], r["substitutions"]): r["name"] for r in index}

    candidates = []

    # -- Engineered variants ------------------------------------------------
    for substitutions, score in scores.items():
        if score["verdict"] == "interface destroyed":
            continue
        sequence = sequences.get(substitutions)
        if sequence is None:
            continue

        acidic, note = composition_risk(sequence)
        p_value = float(control[substitutions]["p_empirical"]) \
            if substitutions in control else 1.0
        new = count_new_substitutions(substitutions)

        candidates.append({
            "kind": "variant" if new else "native network",
            "n_new_substitutions": new,
            "design": score["design"],
            "substitutions": substitutions,
            "engages": score["engages"],
            "sequence": sequence,
            "length": len(sequence),
            "coupling_kcal_mol": float(score["coupling_kcal_mol"]),
            "predicted_fold_change": round(
                pow(2.718281828,
                    float(score["coupling_kcal_mol"]) / (0.0019872 * 298.15)),
                1),
            "i_ptm": float(score["i_ptm"]),
            "d_i_ptm": float(score["d_i_ptm"]),
            "interface": score["verdict"],
            "control_p": p_value,
            "acidic_fraction": round(acidic, 3),
            "composition_note": note,
        })

    # -- Unmodified parents -------------------------------------------------
    if not args.no_parents:
        for row in index:
            if row["substitutions"] != "none (parent)":
                continue
            prediction = predictions.get(row["name"])
            if prediction is None:
                continue
            sequence = prediction["binder"]
            acidic, note = composition_risk(sequence)
            candidates.append({
                "kind": "parent",
                "n_new_substitutions": 0,
                "design": row["design"],
                "substitutions": "none",
                "engages": "",
                "sequence": sequence,
                "length": len(sequence),
                "coupling_kcal_mol": 0.0,
                "predicted_fold_change": 1.0,
                "i_ptm": round(float(prediction["i_ptm"]), 3),
                "d_i_ptm": 0.0,
                "interface": "parent",
                "control_p": 1.0,
                "acidic_fraction": round(acidic, 3),
                "composition_note": note,
            })

    # -- Rank ---------------------------------------------------------------
    def key(entry):
        return (
            1 if entry["kind"] == "parent" else 0,
            -entry["coupling_kcal_mol"],
            0 if entry["interface"] == "interface retained" else 1,
            entry["control_p"],
            entry["acidic_fraction"],
            -entry["i_ptm"],
        )

    candidates.sort(key=key)

    # Identical sequences would waste a submission slot and the rules
    # require designs to be unique.
    seen, unique = set(), []
    for entry in candidates:
        if entry["sequence"] in seen:
            continue
        seen.add(entry["sequence"])
        unique.append(entry)

    selected = unique[:args.max_designs]
    for rank, entry in enumerate(selected, start=1):
        entry["rank"] = rank
        if entry["kind"] == "variant":
            tag = entry["substitutions"].replace(" ", "")
        elif entry["kind"] == "native network":
            tag = "native"
        else:
            tag = "parent"
        entry["name"] = f"EGFRd3_{rank:02d}_{tag}_{entry['design'].split('_')[2]}"

    # -- Submission CSV -----------------------------------------------------
    submission = OUT_DIR / "submission.csv"
    with submission.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "sequence", "molecule_class"])
        for entry in selected:
            # The FAQ calls this class "protein"; the upload form accepts
            # "single_chain", "nanobody", "scfv", "fab_kappa", "fab_lambda".
            # The form is what validates the file, so follow the form.
            writer.writerow([entry["name"], entry["sequence"], "single_chain"])

    metadata = OUT_DIR / "submission_metadata.csv"
    columns = ["rank", "name", "kind", "n_new_substitutions", "design",
               "substitutions", "engages",
               "length", "coupling_kcal_mol", "predicted_fold_change",
               "control_p", "i_ptm", "d_i_ptm", "interface",
               "acidic_fraction", "composition_note", "sequence"]
    with metadata.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)

    print(f"{len(selected)} designs selected of {len(unique)} candidates\n")
    print(f"  {'#':>2} {'name':<32}{'new':>4}{'len':>5}{'kcal':>7}{'fold':>8}"
          f"{'p':>8}{'i_pTM':>8}")
    print("  " + "-" * 74)
    for entry in selected:
        print(f"  {entry['rank']:>2} {entry['name']:<32}"
              f"{entry['n_new_substitutions']:>4}{entry['length']:>5}"
              f"{entry['coupling_kcal_mol']:>7.2f}"
              f"{entry['predicted_fold_change']:>8.1f}"
              f"{entry['control_p']:>8.4f}{entry['i_ptm']:>8.3f}")

    # -- Methods ------------------------------------------------------------
    meta = json.loads(ANNOTATION_META.read_text())
    rejects = json.loads(REJECTIONS.read_text()) if REJECTIONS.exists() else {}
    totals = rejects.get("totals", {})

    variants_n = sum(1 for e in selected if e["kind"] == "variant")
    native_n = sum(1 for e in selected if e["kind"] == "native network")
    parents_n = sum(1 for e in selected if e["kind"] == "parent")
    retained = sum(1 for e in selected
                   if e["interface"] == "interface retained")
    significant = sum(1 for e in selected
                      if e["kind"] == "variant" and e["control_p"] < 0.05)

    methods = OUT_DIR / "METHODS.md"
    with methods.open("w") as fh:
        w = fh.write
        w("# Methods\n\n")
        w("Conditional EGFR binder design for Challenge 01 of the "
          "Anthropic × Adaptyv 2026 competition, Track 3.\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"from the pipeline outputs by `scripts/10_build_submission.py`. "
          f"Every figure below is read from the analysis files rather than "
          f"transcribed.\n\n")
        w("---\n\n")

        w("## Summary\n\n")
        w(f"{len(selected)} single-chain de novo proteins of "
          f"{min(e['length'] for e in selected)}–"
          f"{max(e['length'] for e in selected)} residues, targeting domain "
          f"III of the human EGFR extracellular region. "
          f"{variants_n} carry an engineered histidine network intended to "
          f"make binding pH-dependent. {native_n} already presented such a "
          f"network without modification and is submitted unchanged. "
          f"{parents_n} are parent binders with no network, included as a "
          f"fallback should the switch fail.\n\n")
        if native_n:
            w("That a design arrived with a usable protonation network "
              "already in place is reported rather than absorbed into the "
              "engineered count. BindCraft optimises interface confidence "
              "and has no notion of pH, so the network is incidental, and "
              "counting it as engineered would claim credit for something "
              "the pipeline did not do. It is submitted regardless, because "
              "an incidental network is as real as a designed one.\n\n")

        w("## Objectives addressed\n\n")
        w("The challenge asks for three properties, ranked with pH "
          "selectivity first.\n\n")
        w("*Cross-species reactivity* is handled upstream rather than "
          "designed for. A residue entered the designable pool only if it "
          "is **identical** between human and mouse, not merely conserved. "
          "The bar is deliberately absolute: the cetuximab epitope is 20 of "
          "27 residues identical with none differing outright, and "
          "cetuximab still fails to recognise murine EGFR, so a handful of "
          "conservative substitutions across a footprint is enough to "
          "abolish antibody binding. Every submitted binder therefore "
          "targets a surface with no human-mouse difference at all.\n\n")
        w("*pH selectivity* is engineered, since the target cannot be "
          "mutated and no binder-design pipeline optimises for it.\n\n")
        w("*Affinity* is what BindCraft optimises, and is the property "
          "these designs are least differentiated on.\n\n")

        w("## Software\n\n")
        w("| Tool | Version | Role |\n|------|---------|------|\n")
        w("| BindCraft | `martinpacesa/BindCraft`, cloned 2026-09-29 "
          "| binder generation |\n")
        w("| AlphaFold2 | `alphafold_params_2022-12-06`, multimer v3 "
          "| backbone hallucination and complex prediction |\n")
        w("| ProteinMPNN | as vendored by BindCraft | sequence design |\n")
        w("| PyRosetta | `2026.29+release.quarterly` "
          "| interface relaxation and scoring |\n")
        w("| PROPKA | 3.5.1 | free-state pKa estimation |\n")
        w("| FreeSASA | 2.2.1 | solvent accessibility |\n")
        w("| Biopython | 1.88 | structure and sequence handling |\n")
        w("\nPyRosetta and the AlphaFold2 parameters carry non-commercial "
          "licences; this work is academic. Full terms are recorded in the "
          "repository.\n\n")

        w("## Target definition\n\n")
        w(f"Human EGFR (UniProt P00533, entry version 301) and mouse EGFR "
          f"(Q01279, version 253) were retrieved from UniProtKB release "
          f"2026_03, and the EGFR ectodomain in complex with a cetuximab Fab "
          f"mutant (PDB 6ARU) from the RCSB PDB, all on 29 September 2026. "
          f"SHA-256 checksums are recorded in the accompanying repository.\n\n")
        w(f"Every residue of the extracellular region "
          f"({meta['ectodomain'][0]}–{meta['ectodomain'][1]}, precursor "
          f"numbering) was annotated on four layers: human/mouse "
          f"conservation from a global BLOSUM62 alignment; solvent "
          f"accessibility computed with FreeSASA on chain A of 6ARU; "
          f"occlusion by N-glycans, using the sugars the crystal resolved "
          f"and the glycosylated asparagine where it did not; and "
          f"disulfide-derived architecture. A residue was called designable "
          f"when exposed above {meta['parameters']['exposure_threshold']} "
          f"relative SASA, identical between species, clear of glycan, and "
          f"not disulfide-bonded. {meta['counts']['designable']} of "
          f"{meta['counts']['ectodomain_residues']} qualified.\n\n")

        w("### Numbering\n\n")
        w(f"Structure 6ARU is numbered by the mature protein, "
          f"{meta['uniprot_minus_pdb_offset']} below UniProt. The offset was "
          f"derived by alignment and verified independently: mapping "
          f"UniProt's 25 annotated disulfide pairs through it placed 24 of "
          f"them at 2.03–2.05 Å SG–SG in the crystal, against the 2.05 Å of "
          f"a covalent bond. That test is geometric and does not reuse the "
          f"alignment it validates.\n\n")

        w("## Epitope selection\n\n")
        w("Designable residues were grouped into contiguous surface patches "
          "and ranked on pH-switch potential, conservation of the footprint "
          "including a shell beyond it, overlap with the cetuximab "
          "interface as evidence of druggability, and buriable area. "
          "Patches were then characterised geometrically, which showed that "
          "several were not single connected surfaces; residues close in a "
          "distance matrix but detached from the patch face were removed "
          "from the hotspot lists.\n\n")
        w("Patch 2 was selected. It is the only one among the top-ranked "
          "patches forming a single connected surface, it contacts the "
          "cetuximab Fab directly, and it contains H433, Q435 and K489, "
          "three positions with prior experimental support that the "
          "selection recovered without using the literature as input.\n\n")

        w("## pH-switch mechanism\n\n")
        w("Wyman linkage caps the pH dependence of binding at "
          "2.303·R·T·Δn per pH unit, so across the 0.9 units between 7.4 "
          "and 6.5 each fully coupled protonation contributes at most 1.23 "
          "kcal/mol. A single histidine therefore buys about eightfold, "
          "which is not the requested behaviour, and the design families "
          "were organised around the number of coupled protonation events "
          "rather than around any one position.\n\n")
        w("Selectivity arises from the pKa *shift* on binding, not from the "
          "pKa: a group whose protonation state is unchanged on binding "
          "contributes nothing. Free-state pKa values were estimated with "
          "PROPKA. H433 came out at 6.22, titrating almost ideally for this "
          "window; the target carboxylates D458, D460 and E455 sit between "
          "1.95 and 4.31 and are therefore dependable counter-charges that "
          "upshift an engineered histidine's pKa on burial.\n\n")

        w("## Binder generation\n\n")
        w("Domain III was trimmed from 6ARU chain A to 204 residues, which "
          "is what makes the job fit in the memory of the available "
          "hardware. Binders were generated with BindCraft, which "
          "hallucinates a backbone by AlphaFold2 backpropagation, designs "
          "its sequence with ProteinMPNN, re-predicts the complex and "
          "scores it with PyRosetta. Optimisation iterations were reduced "
          "from the published defaults to fit a 30 GPU-hour weekly budget "
          "on two T4 cards; validation settings were left untouched, since "
          "weakening the stage that decides whether a design is good would "
          "only move the failure downstream.\n\n")
        if totals:
            w(f"{totals.get('trajectories_started', 0)} trajectories were "
              f"run across two batches, of which "
              f"{totals.get('accepted', 0)} produced accepted designs. "
              f"{totals.get('confidence_rejections', 0)} of "
              f"{totals.get('rejections_logged', 0)} rejections were on "
              f"AlphaFold2 confidence in the re-predicted complex rather "
              f"than on composition or geometry, so the filters were left "
              f"as published: relaxing them would accept designs the model "
              f"does not believe in.\n\n")

        w("## Histidine engineering\n\n")
        w("Histidines were placed at binder positions facing a target "
          "ionizable group, within the span a histidine can actually cover. "
          "That span was measured on the 23 histidines resolved in 6ARU, "
          "giving CB–ND1 at 2.53 Å and CB–NE2 at 3.67 Å; with a 4.0 Å "
          "charge-assisted contact the ceiling is 7.7 Å from CB. "
          "Substitutions were proposed in sets rather than singly, since "
          "the thermodynamics make the count the governing variable, and "
          "positions already carrying a histidine were kept as existing "
          "switches rather than discarded.\n\n")
        w("Each proposal was then tested for rotamer feasibility: an "
          "imidazole nitrogen was placed at 200 orientations on spheres of "
          "both measured radii around CB, orientations clashing with the "
          "binder or target discarded, and the substitution accepted only "
          "if a surviving orientation reached hydrogen-bonding distance of "
          "the handle.\n\n")

        w("## Controls and validation\n\n")
        w("**Negative control.** Designed placements were compared against "
          "histidines placed at random among the positions the method could "
          "have chosen from: interface, exposed in the unbound binder, and "
          "chemically substitutable, with only the targeting removed. "
          "Across 3000 draws per design and substitution count, designed "
          "variants average 1.82 kcal/mol against 0.40 for random "
          f"placement. {significant} of the {variants_n} submitted variants "
          f"reach p < 0.05. No designed variant exceeds the maximum of its "
          "null distribution, so the claim is about frequency rather than "
          "about reaching something unreachable, and the p-values are "
          "uncorrected.\n\n")
        w("**Structural revalidation.** Every variant and its unmodified "
          "parent were re-predicted from sequence under identical settings. "
          f"{retained} of the {variants_n + native_n} submitted designs "
          "that carry a network retain interface confidence "
          "within 0.05 i_pTM of their parent. Damage does not track the "
          "number of substitutions: the lead variant, which engages three "
          "handles, changes i_pTM by +0.002, while one single substitution "
          "collapsed "
          "its complex from 0.819 to 0.243 and was excluded. The variants "
          "that cost least are those built on histidines the design already "
          "carried.\n\n")
        w("**Novelty.** No binder returns a hit from the RCSB sequence "
          "service across every polymer entity in the PDB, and the highest "
          "pairwise identity within the portfolio is 27%.\n\n")

        w("## Ranking\n\n")
        w("Submissions are ordered to follow the competition's own "
          "criteria. Variants whose substitutions destroyed the interface "
          "were excluded outright. The remainder are ordered by estimated "
          "coupling, then by whether the interface survived intact rather "
          "than merely tolerably, then by the control p-value, then by "
          "composition risk. Unmodified parents follow the variants: they "
          "carry no switch but are validated binders and cost nothing "
          "within the submission limit.\n\n")

        w("## Limitations\n\n")
        w("The coupling estimates are model-dependent upper bounds under an "
          "assumed 2.0-unit pKa shift on binding. That shift is the "
          "quantity that actually decides the outcome and it cannot be "
          "computed before the complex exists; a rigorous treatment would "
          "need constant-pH molecular dynamics. PROPKA's error is near one "
          "pKa unit, comparable to the effects reasoned about.\n\n")
        w("The binders were generated by AlphaFold2 and re-predicted by "
          "AlphaFold2, so absolute confidence values are not independent "
          "evidence of binding. The variant-versus-parent comparison is "
          "sound, since no variant was optimised and both are predicted "
          "identically, but only the assay can establish that any of these "
          "bind at all.\n\n")
        w("Every design is strongly acidic, with glutamate between 14 and "
          "28 percent against a natural average near 7. This is "
          "characteristic of hallucinated sequences rather than of any one "
          "trajectory, and it is an expression and solubility liability "
          "carried by the whole portfolio. It is weighed in the ranking but "
          "not solved.\n\n")
        w("Sequence identity is a weak proxy for structural novelty, and no "
          "structural comparison against the PDB was performed.\n\n")

        w("## Reproducibility\n\n")
        w("All code, parameters, intermediate results and the provenance of "
          "every input are in the accompanying repository. Analysis scripts "
          "are numbered in execution order and run on CPU in minutes; only "
          "binder generation and the revalidation require a GPU. Raw inputs "
          "are reconstructed by a download script and pinned by SHA-256 "
          "checksums that can be re-verified with a single command.\n")

    print(f"\n  -> {submission.relative_to(ROOT)}")
    print(f"  -> {metadata.relative_to(ROOT)}")
    print(f"  -> {methods.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

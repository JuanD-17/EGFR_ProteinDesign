#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 09 - Novelty of the designed binders.

The competition requires designs to be de novo and zero-shot, with adequate
sequence and structural diversity from known proteins, and states that
submissions failing this are filtered out before screening. That is a hard
criterion, not a preference, so it is checked rather than assumed.

These binders come from AlphaFold2 hallucination against a target
structure, with no starting binder, so they are expected to be novel. The
point of the check is to have a measurement instead of an expectation: a
methods section that reports "the closest match in the PDB shares 23%
identity over 40 residues" is evidence, and one that says "the designs are
de novo by construction" is a claim.

What is measured
----------------
Each binder sequence is searched against the Protein Data Bank using the
RCSB sequence search service, which runs MMseqs2 over all polymer entities.
The PDB is the right reference here: it is the set of proteins whose
structures are known, which is what "structural diversity from known
proteins" refers to, and a designed backbone resembling a solved structure
is exactly what the criterion excludes.

Two further checks come free. Sequence composition is compared against
natural protein averages, since hallucinated sequences sometimes drift into
compositions no natural protein has. And the designs are compared against
each other, because a portfolio of near-identical binders is not a
portfolio, whatever its diversity from the PDB.

Inputs:  03_Design/accepted/*.pdb and the variants from step 05
Outputs: 02_Analysis/09_novelty_pdb_hits.csv
         02_Analysis/09_novelty_report.md

Usage:
    python scripts/09_check_novelty.py
    python scripts/09_check_novelty.py --identity-cutoff 0.25
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "02_Analysis"
DESIGNS = ROOT / "03_Design" / "accepted"
VARIANTS = ROOT / "03_Design" / "ph_variants" / "05_variants.csv"

_spec = importlib.util.spec_from_file_location(
    "step05", ROOT / "scripts" / "05_engineer_ph_switch.py")
step05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step05)

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
USER_AGENT = "EGFR-ProteinDesign/1.0 (novelty check)"

# Report anything above this identity. The competition does not publish a
# numeric threshold, so nothing is filtered automatically: hits are listed
# and the judgement is left explicit.
REPORT_IDENTITY = 0.20

# Mean amino acid frequencies in natural proteins, for a composition sanity
# check. Approximate values from UniProtKB, adequate for spotting a design
# that has drifted somewhere no natural protein goes.
NATURAL_FREQUENCY = {
    "A": 0.0825, "R": 0.0553, "N": 0.0406, "D": 0.0545, "C": 0.0137,
    "Q": 0.0393, "E": 0.0675, "G": 0.0707, "H": 0.0227, "I": 0.0596,
    "L": 0.0966, "K": 0.0584, "M": 0.0242, "F": 0.0386, "P": 0.0470,
    "S": 0.0656, "T": 0.0534, "W": 0.0108, "Y": 0.0292, "V": 0.0687,
}


def search_pdb(sequence: str, identity_cutoff: float,
               rows: int = 10) -> list[dict]:
    """Sequence search against every polymer entity in the PDB."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 10,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": sequence,
            },
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "paginate": {"start": 0, "rows": rows},
        },
        "return_type": "polymer_entity",
    }
    request = urllib.request.Request(
        RCSB_SEARCH,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json",
                 "User-Agent": USER_AGENT},
    )
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            raw = response.read()
            status = response.status
        # RCSB answers 204 with an empty body when nothing matches, which
        # urllib does not raise on. For a novelty check that is the desired
        # outcome, not a failure.
        if status == 204 or not raw.strip():
            return []
        body = json.loads(raw)
    except urllib.error.HTTPError as error:
        # The service answers 204 when nothing matches, which for a novelty
        # check is the desired outcome rather than a failure.
        if error.code in (204, 404):
            return []
        raise
    except urllib.error.URLError as error:
        raise RuntimeError(f"RCSB search unreachable: {error}") from None

    hits = []
    for result in body.get("result_set", []):
        context = {}
        for service in result.get("services", []):
            for node in service.get("nodes", []):
                for match in node.get("match_context", []):
                    context = match
        hits.append({
            "pdb_entity": result.get("identifier", ""),
            "score": round(result.get("score", 0.0), 4),
            "identity": round(context.get("sequence_identity", 0.0), 4),
            "evalue": context.get("evalue", ""),
            "alignment_length": context.get("alignment_length", ""),
            "query_coverage": round(
                context.get("query_length", 0)
                and context.get("alignment_length", 0)
                / context["query_length"], 3) if context.get("query_length")
                else "",
        })
    return hits


def composition_deviation(sequence: str) -> tuple[float, list[tuple]]:
    """How far the composition sits from natural protein averages."""
    length = len(sequence)
    counts = {aa: sequence.count(aa) / length for aa in NATURAL_FREQUENCY}
    deviations = [
        (aa, counts[aa], NATURAL_FREQUENCY[aa],
         counts[aa] - NATURAL_FREQUENCY[aa])
        for aa in NATURAL_FREQUENCY
    ]
    deviations.sort(key=lambda d: -abs(d[3]))
    total = sum(abs(d[3]) for d in deviations) / 2      # total variation
    return total, deviations[:5]


def identity_between(a: str, b: str) -> float:
    """Identity of two binder sequences over their global alignment."""
    aligner = step05.Align.PairwiseAligner()
    aligner.substitution_matrix = step05.substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    alignment = aligner.align(a, b)[0]
    matches = sum(
        1
        for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned)
        for offset in range(a_end - a_start)
        if a[a_start + offset] == b[b_start + offset]
    )
    return matches / max(len(a), len(b))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--identity-cutoff", type=float, default=0.0,
                        help="lowest identity the PDB search returns "
                             "(default 0.0, i.e. report everything)")
    parser.add_argument("--offline", action="store_true",
                        help="skip the PDB search and run only the local "
                             "composition and cross-design checks")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    binders: dict[str, str] = {}
    for path in sorted(DESIGNS.glob("*.pdb")):
        model = step05.PDBParser(QUIET=True).get_structure(
            path.stem, str(path))[0]
        _, binder_id = step05.identify_chains(model, 204)
        sequence, _ = step05.chain_sequence(model[binder_id])
        binders[path.stem] = sequence

    if not binders:
        raise SystemExit(f"no designs in {DESIGNS}")

    print(f"{len(binders)} binder sequences\n")

    # -- PDB search ---------------------------------------------------------
    rows = []
    if not args.offline:
        print("Searching the PDB (RCSB sequence service)")
        for name, sequence in binders.items():
            try:
                hits = search_pdb(sequence, args.identity_cutoff)
            except RuntimeError as error:
                print(f"  {name}: {error}")
                print("  network unavailable; rerun later or use --offline")
                return 1

            if hits:
                best = max(hits, key=lambda h: h["identity"])
                print(f"  {name[:40]:<42} {len(hits)} hit(s), "
                      f"best {100 * best['identity']:.1f}% identity "
                      f"to {best['pdb_entity']}")
            else:
                print(f"  {name[:40]:<42} no hits")

            for hit in hits:
                rows.append({"design": name, "binder_length": len(sequence),
                             **hit})
            time.sleep(1)       # courtesy to a public service
        print()

    # -- Composition --------------------------------------------------------
    print("Composition against natural protein averages")
    compositions = {}
    for name, sequence in binders.items():
        total, top = composition_deviation(sequence)
        compositions[name] = (total, top)
        flagged = ", ".join(
            f"{aa} {100 * observed:.0f}% vs {100 * natural:.0f}%"
            for aa, observed, natural, _ in top[:3])
        print(f"  {name[:40]:<42} deviation {total:.3f}   {flagged}")
    print()

    # -- Designs against each other -----------------------------------------
    print("Designs against each other")
    pairs = []
    for (name_a, seq_a), (name_b, seq_b) in combinations(binders.items(), 2):
        identity = identity_between(seq_a, seq_b)
        pairs.append((name_a, name_b, identity))
        if identity >= 0.30:
            print(f"  {name_a[:28]} / {name_b[:28]}: "
                  f"{100 * identity:.1f}%  <-- similar")
    if pairs:
        worst = max(pairs, key=lambda p: p[2])
        print(f"  highest pairwise identity: {100 * worst[2]:.1f}% "
              f"({worst[0][:26]} / {worst[1][:26]})")
    print()

    # -- Outputs ------------------------------------------------------------
    if rows:
        hits_path = OUT_DIR / "09_novelty_pdb_hits.csv"
        with hits_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    else:
        hits_path = None

    report = OUT_DIR / "09_novelty_report.md"
    with report.open("w") as fh:
        w = fh.write
        w("# Step 09 - Novelty of the designed binders\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/09_check_novelty.py`.\n\n")

        w("## Why this is measured\n\n")
        w("The competition filters out designs that are not de novo and "
          "zero-shot, or that lack adequate sequence and structural "
          "diversity from known proteins. These binders come from "
          "AlphaFold2 hallucination with no starting binder, so novelty is "
          "expected, but an expectation is not evidence.\n\n")

        w("## Against the Protein Data Bank\n\n")
        if args.offline:
            w("Skipped: run without `--offline` to query RCSB.\n\n")
        elif rows:
            w("Searched with the RCSB sequence service, which runs MMseqs2 "
              "over all polymer entities.\n\n")
            w("| Design | Length | Closest PDB entity | Identity | E-value "
              "| Aligned |\n")
            w("|--------|-------:|--------------------|---------:|--------:"
              "|--------:|\n")
            for name in binders:
                design_hits = [r for r in rows if r["design"] == name]
                if not design_hits:
                    w(f"| {name[:30]} | {len(binders[name])} | no hits "
                      f"| — | — | — |\n")
                    continue
                best = max(design_hits, key=lambda h: h["identity"])
                w(f"| {name[:30]} | {len(binders[name])} "
                  f"| {best['pdb_entity']} "
                  f"| {100 * best['identity']:.1f}% | {best['evalue']} "
                  f"| {best['alignment_length']} |\n")
            w("\nFull results in `09_novelty_pdb_hits.csv`.\n\n")
        else:
            w("**No design returned a hit.** The search covers every polymer "
              "entity in the PDB, so this is the strongest available "
              "statement that these sequences do not resemble a protein of "
              "known structure.\n\n")

        w("## Composition\n\n")
        w("Hallucinated sequences can drift into compositions no natural "
          "protein has, which is a different novelty problem: too novel to "
          "express or fold. Total variation distance from mean UniProtKB "
          "frequencies, where 0 is identical and 1 is disjoint.\n\n")
        w("| Design | Deviation | Largest departures |\n")
        w("|--------|----------:|--------------------|\n")
        for name, (total, top) in compositions.items():
            departures = ", ".join(
                f"{aa} {100 * observed:.0f}% (nat {100 * natural:.0f}%)"
                for aa, observed, natural, _ in top[:3])
            w(f"| {name[:30]} | {total:.3f} | {departures} |\n")

        w("\n## Against each other\n\n")
        w("A portfolio of near-identical binders is not a portfolio, "
          "whatever its distance from the PDB.\n\n")
        w("| Design A | Design B | Identity |\n")
        w("|----------|----------|---------:|\n")
        for name_a, name_b, identity in sorted(pairs, key=lambda p: -p[2]):
            w(f"| {name_a[:26]} | {name_b[:26]} | {100 * identity:.1f}% |\n")

        w("\n## Limitations\n\n")
        w("Sequence identity is a weak proxy for structural novelty. Two "
          "proteins can share a fold at 15% identity, and the competition "
          "asks for structural diversity as well as sequence diversity. A "
          "complete answer would compare the designed backbones against the "
          "PDB structurally, with something like Foldseek, which is not run "
          "here.\n\n")
        w("No automatic threshold is applied. The competition publishes no "
          "numeric cutoff, so hits are reported and the judgement is left "
          "explicit rather than hidden in a filter.\n")

    print(f"  -> {report.relative_to(ROOT)}")
    if hits_path:
        print(f"  -> {hits_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

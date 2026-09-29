#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 05 - Histidine engineering for pH-dependent binding.

BindCraft produces binders that bind. It has no notion of pH: it optimises
AlphaFold2 confidence and interface quality, and a design that binds well at
pH 7.4 is exactly what it is trying to make. Objective 1 of the challenge
asks for the opposite behaviour at that pH, so the switch has to be added
here, on interfaces that already exist.

The mechanism, from step 03b
----------------------------
Selectivity comes from the pKa *shift* on binding, not from the pKa. A
histidine whose protonation state is the same free and bound contributes
nothing. The switch requires that binding stabilise the protonated form,
which is what burying a histidine against a carboxylate does: the anion
raises the histidine's pKa, so the complex takes up a proton, and by Wyman
linkage the affinity becomes pH-dependent in the direction we want.

Step 03b identified the anchors on the target: D458 and D460 at pKa 1.95 and
3.70 are reliably ionised across the whole window, E455 and D368 nearly so.
H433, at pKa 6.22, can act as a switch from the receptor side against a
carboxylate placed on the binder.

One protonation event buys roughly 8-fold at best. Two or three are needed,
so this script looks for *sets* of substitutions rather than single ones.

What it does
------------
For each designed complex it locates the binder positions that face a target
handle at a distance a side chain could bridge, judges whether each position
tolerates a histidine, and proposes combinations that engage two or three
handles at once. Output is a ranked table of variants with their rationale,
plus the mutated sequences ready for refolding.

What it does not do
-------------------
It does not predict pKa in the designed complex, which would need the
complex to exist first and, done properly, constant-pH molecular dynamics.
Every variant it proposes must be refolded and rescored before being
believed. The ranking is a prioritisation for that validation, not a result.

Inputs:  a directory of designed complex PDBs from BindCraft
         02_Analysis/03b_pka_estimates.csv
Outputs: 03_Design/ph_variants/05_variants.csv
         03_Design/ph_variants/05_report.md
         03_Design/ph_variants/sequences/*.fasta

Usage:
    python scripts/05_engineer_ph_switch.py --designs 03_Design/accepted/
    python scripts/05_engineer_ph_switch.py --designs <dir> --max-substitutions 3
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "03_Design" / "ph_variants"
PKA_TABLE = ROOT / "02_Analysis" / "03b_pka_estimates.csv"
TARGET_PDB = ROOT / "03_Design" / "target" / "EGFR_domainIII.pdb"

UNIPROT_MINUS_PDB = 24

TITRATABLE_ATOMS = {
    "HIS": ["ND1", "NE2"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

# A histidine imidazole reaches about 6 A from CB. Allowing 4-10 A between
# the binder CB and the target's titratable atom covers rotamers that can
# form the contact without demanding one that already does.
BRIDGE_RANGE = (4.0, 10.0)

# Positions whose side chain points away from the handle cannot reach it
# however favourable the distance.
MIN_ORIENTATION = 0.0

# Residues that should not be mutated to histidine, and why.
PROTECTED = {
    "C": "may be disulfide-bonded",
    "P": "backbone conformation depends on proline",
    "G": "may occupy a position requiring positive phi",
}

# Burial above this fraction means the position is in the binder core, where
# an introduced charge is destabilising rather than switch-forming.
MAX_BURIAL = 0.75

MAX_SUBSTITUTIONS = 3


def load_handles() -> dict[int, dict]:
    """Target handles from step 03b, keyed by UniProt position."""
    handles = {}
    with PKA_TABLE.open() as fh:
        for row in csv.DictReader(fh):
            verdict = row["verdict"]
            if verdict.startswith("weak anchor") or "cannot switch" in verdict:
                continue
            if verdict.startswith("too acidic"):
                continue
            handles[int(row["uniprot_pos"])] = {
                "uniprot_pos": int(row["uniprot_pos"]),
                "pdb_resnum": int(row["pdb_resnum"]),
                "residue": row["residue"],
                "aa": row["aa"],
                "role": row["role"],
                "coupling": float(row["coupling_kcal_mol"]),
                "pka_free": float(row["pka_free"]),
            }
    return handles


def chain_sequence(chain) -> tuple[str, list]:
    letters, residues = [], []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        try:
            letters.append(protein_letters_3to1[residue.get_resname()])
        except KeyError:
            continue
        residues.append(residue)
    return "".join(letters), residues


def identify_chains(model, target_length: int) -> tuple[str, str]:
    """Work out which chain is the target and which the binder.

    The target is the longer chain: domain III is 204 residues and the
    binders are 55-75. Relying on chain identifiers would break the moment
    BindCraft names them differently.
    """
    lengths = {}
    for chain in model:
        sequence, _ = chain_sequence(chain)
        if sequence:
            lengths[chain.id] = len(sequence)
    if len(lengths) < 2:
        raise ValueError(f"expected two chains, found {list(lengths)}")
    ordered = sorted(lengths.items(), key=lambda kv: -kv[1])
    return ordered[0][0], ordered[1][0]


def side_chain_direction(residue) -> np.ndarray | None:
    """Unit vector from CA towards the side chain."""
    if "CA" not in residue:
        return None
    if "CB" in residue:
        vector = residue["CB"].coord - residue["CA"].coord
    else:
        return None
    norm = np.linalg.norm(vector)
    return vector / norm if norm > 0 else None


def relative_burial(residue, all_atoms: np.ndarray) -> float:
    """Crude burial: neighbouring heavy atoms within 10 A, normalised.

    Full SASA on every design would be slow and is not needed: the point is
    only to exclude positions packed into the binder core, and a neighbour
    count separates those from surface positions reliably enough.
    """
    if "CB" not in residue and "CA" not in residue:
        return 1.0
    centre = (residue["CB"] if "CB" in residue else residue["CA"]).coord
    distances = np.linalg.norm(all_atoms - centre, axis=1)
    neighbours = int((distances < 10.0).sum())
    return min(1.0, neighbours / 120.0)


def evaluate_design(path: Path, handles: dict[int, dict],
                    max_substitutions: int) -> dict | None:
    model = PDBParser(QUIET=True).get_structure(path.stem, str(path))[0]

    try:
        target_id, binder_id = identify_chains(model, 204)
    except ValueError as error:
        print(f"  {path.name}: {error}")
        return None

    target_chain = model[target_id]
    binder_chain = model[binder_id]
    binder_sequence, binder_residues = chain_sequence(binder_chain)

    binder_atoms = np.array([a.coord for r in binder_residues for a in r
                             if a.element != "H"])

    # Locate the handles that are actually present in this complex.
    present = {}
    for position, handle in handles.items():
        try:
            residue = target_chain[(" ", handle["pdb_resnum"], " ")]
        except KeyError:
            continue
        atoms = [residue[n] for n in TITRATABLE_ATOMS[handle["residue"]]
                 if n in residue]
        if atoms:
            present[position] = {
                **handle,
                "coord": np.mean([a.coord for a in atoms], axis=0),
            }

    if not present:
        return None

    # For each handle, which binder positions could bridge to it.
    candidates = []
    for position, handle in present.items():
        for index, residue in enumerate(binder_residues):
            if "CB" not in residue and "CA" not in residue:
                continue
            centre = (residue["CB"] if "CB" in residue
                      else residue["CA"]).coord
            separation = float(np.linalg.norm(centre - handle["coord"]))
            if not (BRIDGE_RANGE[0] <= separation <= BRIDGE_RANGE[1]):
                continue

            direction = side_chain_direction(residue)
            if direction is None:
                continue
            towards = handle["coord"] - centre
            towards = towards / np.linalg.norm(towards)
            orientation = float(np.dot(direction, towards))
            if orientation <= MIN_ORIENTATION:
                continue

            current = binder_sequence[index]
            burial = relative_burial(residue, binder_atoms)

            blocked = PROTECTED.get(current)
            if burial > MAX_BURIAL:
                blocked = f"buried in the binder core ({burial:.2f})"
            if current == "H":
                blocked = "already histidine"

            candidates.append({
                "binder_index": index,
                "binder_position": index + 1,
                "current_aa": current,
                "handle_uniprot": position,
                "handle_aa": handle["aa"],
                "handle_role": handle["role"],
                "handle_coupling": handle["coupling"],
                "distance_a": round(separation, 2),
                "orientation": round(orientation, 3),
                "burial": round(burial, 2),
                "blocked": blocked,
            })

    usable = [c for c in candidates if c["blocked"] is None]
    if not usable:
        return {
            "design": path.stem,
            "binder_chain": binder_id,
            "target_chain": target_id,
            "binder_sequence": binder_sequence,
            "handles_present": sorted(present),
            "candidates": candidates,
            "variants": [],
        }

    # Best binder position per handle, judged on geometry.
    best_per_handle: dict[int, dict] = {}
    for candidate in usable:
        handle = candidate["handle_uniprot"]
        incumbent = best_per_handle.get(handle)
        if incumbent is None or candidate["orientation"] > incumbent["orientation"]:
            best_per_handle[handle] = candidate

    # Combinations engaging several handles. One protonation is not enough,
    # so single substitutions are recorded but ranked below pairs and triples.
    variants = []
    handle_list = sorted(best_per_handle)
    for size in range(1, min(max_substitutions, len(handle_list)) + 1):
        for combination in itertools.combinations(handle_list, size):
            picks = [best_per_handle[h] for h in combination]
            positions = [p["binder_position"] for p in picks]
            if len(set(positions)) != len(positions):
                continue                      # one residue cannot serve two

            mutated = list(binder_sequence)
            for pick in picks:
                mutated[pick["binder_index"]] = "H"

            coupling = sum(p["handle_coupling"] for p in picks)
            variants.append({
                "design": path.stem,
                "n_substitutions": size,
                "substitutions": " ".join(
                    f"{p['current_aa']}{p['binder_position']}H" for p in picks),
                "engages": " ".join(
                    f"{p['handle_aa']}{p['handle_uniprot']}" for p in picks),
                "mean_distance_a": round(
                    float(np.mean([p["distance_a"] for p in picks])), 2),
                "min_orientation": round(
                    float(min(p["orientation"] for p in picks)), 3),
                "estimated_coupling_kcal_mol": round(coupling, 3),
                "estimated_fold_change": round(
                    float(np.exp(coupling / (0.0019872 * 298.15))), 1),
                "sequence": "".join(mutated),
            })

    variants.sort(key=lambda v: (-v["n_substitutions"],
                                 -v["estimated_coupling_kcal_mol"]))
    return {
        "design": path.stem,
        "binder_chain": binder_id,
        "target_chain": target_id,
        "binder_sequence": binder_sequence,
        "handles_present": sorted(present),
        "candidates": candidates,
        "variants": variants,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--designs", required=True,
                        help="directory of designed complex PDBs")
    parser.add_argument("--max-substitutions", type=int,
                        default=MAX_SUBSTITUTIONS,
                        help=f"histidines per variant (default {MAX_SUBSTITUTIONS})")
    parser.add_argument("--top", type=int, default=3,
                        help="variants to keep per design (default 3)")
    args = parser.parse_args()

    designs_dir = Path(args.designs)
    if not designs_dir.is_dir():
        raise SystemExit(f"not a directory: {designs_dir}")

    pdbs = sorted(designs_dir.glob("*.pdb"))
    if not pdbs:
        raise SystemExit(f"no PDB files in {designs_dir}")

    handles = load_handles()
    print(f"{len(handles)} usable target handles from step 03b:")
    for position, handle in sorted(handles.items()):
        print(f"  {handle['aa']}{position}  {handle['role']}, "
              f"pKa {handle['pka_free']}, {handle['coupling']} kcal/mol")

    print(f"\n{len(pdbs)} designs in {designs_dir}\n")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "sequences").mkdir(exist_ok=True)

    results, all_variants = [], []
    for path in pdbs:
        result = evaluate_design(path, handles, args.max_substitutions)
        if result is None:
            continue
        results.append(result)
        kept = result["variants"][:args.top]
        all_variants.extend(kept)

        handles_here = ", ".join(f"{h}" for h in result["handles_present"])
        print(f"  {result['design']}: binder chain {result['binder_chain']}, "
              f"{len(result['binder_sequence'])} aa, "
              f"handles nearby [{handles_here}], "
              f"{len(result['variants'])} variant(s)")
        for variant in kept:
            print(f"      {variant['substitutions']:<18} engages "
                  f"{variant['engages']:<18} "
                  f"{variant['estimated_coupling_kcal_mol']:>5.2f} kcal/mol  "
                  f"{variant['estimated_fold_change']:>7.1f}x")

    if not all_variants:
        print("\nNo variants proposed. Either no handle lies within reach of "
              "a substitutable binder position, or the designs do not engage "
              "the intended patch.")
        return 1

    csv_path = OUT_DIR / "05_variants.csv"
    columns = ["design", "n_substitutions", "substitutions", "engages",
               "mean_distance_a", "min_orientation",
               "estimated_coupling_kcal_mol", "estimated_fold_change",
               "sequence"]
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_variants)

    fasta_path = OUT_DIR / "sequences" / "ph_variants.fasta"
    with fasta_path.open("w") as fh:
        for variant in all_variants:
            name = f"{variant['design']}_{variant['substitutions'].replace(' ', '_')}"
            fh.write(f">{name} engages={variant['engages'].replace(' ', ',')} "
                     f"est_coupling={variant['estimated_coupling_kcal_mol']}\n")
            fh.write(variant["sequence"] + "\n")

    report_path = OUT_DIR / "05_report.md"
    with report_path.open("w") as fh:
        w = fh.write
        w("# Step 05 - Histidine engineering\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/05_engineer_ph_switch.py`.\n\n")
        w(f"{len(pdbs)} designs examined, {len(all_variants)} variants "
          f"proposed.\n\n")

        w("## Why substitutions come in sets\n\n")
        w("Wyman linkage caps the pH dependence of binding at about 1.23 "
          "kcal/mol per coupled protonation across the 0.9 pH units between "
          "7.4 and 6.5, which is roughly eight-fold in affinity. The "
          "challenge asks for binding at 6.5 and none detectable at 7.4, so "
          "a single histidine cannot deliver it. Variants engaging two or "
          "three handles are therefore ranked above single substitutions "
          "regardless of individual geometry.\n\n")

        w("## Target handles\n\n")
        w("| Handle | Role | pKa free | kcal/mol |\n")
        w("|--------|------|---------:|---------:|\n")
        for position, handle in sorted(handles.items()):
            w(f"| {handle['aa']}{position} | {handle['role']} "
              f"| {handle['pka_free']} | {handle['coupling']} |\n")

        w("\n## Proposed variants\n\n")
        w("| Design | Substitutions | Engages | Mean distance Å "
          "| kcal/mol | Fold |\n")
        w("|--------|---------------|---------|----------------"
          "|---------:|-----:|\n")
        for variant in all_variants:
            w(f"| {variant['design']} | {variant['substitutions']} "
              f"| {variant['engages']} | {variant['mean_distance_a']} "
              f"| {variant['estimated_coupling_kcal_mol']} "
              f"| {variant['estimated_fold_change']}× |\n")

        w("\n## Selection rules\n\n")
        w(f"- A binder position qualifies when its CB lies "
          f"{BRIDGE_RANGE[0]}–{BRIDGE_RANGE[1]} Å from the handle's "
          f"titratable atom, which is the range a histidine rotamer can "
          f"bridge, and its side chain points towards the handle.\n")
        w(f"- Cysteine, proline and glycine are not substituted: "
          f"{', '.join(f'{k} ({v})' for k, v in PROTECTED.items())}.\n")
        w(f"- Positions with burial above {MAX_BURIAL} are excluded. A charge "
          f"introduced into the binder core destabilises the fold instead of "
          f"forming a switch.\n")
        w("- One binder residue cannot serve two handles, so combinations "
          "reusing a position are discarded.\n")

        w("\n## Limitations\n\n")
        w("The coupling figures are carried over from step 03b's estimates "
          "for the target handles and assume the engineered histidine "
          "achieves the same pKa shift on binding. That shift is the "
          "quantity that actually decides the outcome and it cannot be "
          "computed before the variant complex exists. Nothing here "
          "establishes that a proposed substitution preserves binding at "
          "all: introducing a histidine changes the interface and may simply "
          "break it.\n\n")
        w("**Every variant must be refolded and rescored before it is "
          "believed.** This ranking is a prioritisation for that validation, "
          "not a result.\n")

    print(f"\n  -> {csv_path.relative_to(ROOT)}")
    print(f"  -> {fasta_path.relative_to(ROOT)}")
    print(f"  -> {report_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

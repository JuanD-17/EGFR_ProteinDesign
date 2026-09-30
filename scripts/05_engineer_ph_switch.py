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
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import freesasa
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "03_Design" / "ph_variants"
PKA_TABLE = ROOT / "02_Analysis" / "03b_pka_estimates.csv"
TARGET_PDB = ROOT / "03_Design" / "target" / "EGFR_domainIII.pdb"
HUMAN_FASTA = ROOT / "01_Target" / "human" / "P00533.fasta"

UNIPROT_MINUS_PDB = 24

freesasa.setVerbosity(freesasa.silent)

TITRATABLE_ATOMS = {
    "HIS": ["ND1", "NE2"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

# How far a histidine can actually span. Measured on the 23 histidines in
# 6ARU, CB-NE2 is 3.67 A; adding a 4.0 A charge-assisted contact puts the
# ceiling at 7.7 A from CB to the handle's titratable atom.
#
# An earlier version used 10.0 A, which proposed substitutions no rotamer
# can realise: of 33 proposals, 14 were later vetoed on reach alone and not
# one on steric clash. The bound is now derived from the structure rather
# than estimated.
BRIDGE_RANGE = (3.5, 7.7)

# Positions whose side chain points away from the handle cannot reach it
# however favourable the distance.
MIN_ORIENTATION = 0.0

# Residues that should not be mutated to histidine, and why.
PROTECTED = {
    "C": "may be disulfide-bonded",
    "P": "backbone conformation depends on proline",
    "G": "may occupy a position requiring positive phi",
}

# Relative SASA in the unbound binder below which a position counts as core.
# A charge introduced there destabilises the fold instead of forming a switch.
# Same 0.20 convention used for the target in step 01.
MIN_EXPOSURE = 0.20

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


def map_target_to_uniprot(target_chain) -> dict[int, object]:
    """Map UniProt positions to residues of a design's target chain.

    BindCraft renumbers the target in its output: the trimmed domain III
    arrives as 311-514 and comes back as 1-204. Assuming any fixed offset
    would break the moment the trim or the tool changes, so the
    correspondence is derived by aligning the chain against the human
    sequence, exactly as step 01 does for the crystal structure.
    """
    human = "".join(
        line.strip() for line in HUMAN_FASTA.read_text().splitlines()
        if not line.startswith(">")
    )

    letters, residues = [], []
    for residue in target_chain:
        if residue.id[0] != " ":
            continue
        try:
            letters.append(protein_letters_3to1[residue.get_resname()])
        except KeyError:
            continue
        residues.append(residue)
    sequence = "".join(letters)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    # Free end gaps: a 204-residue fragment is being placed inside a
    # 1210-residue sequence, so overhangs must not be penalised. Biopython
    # renamed these attributes; support both spellings.
    for name, value in (("end_insertion_score", 0.0),
                        ("end_deletion_score", 0.0)):
        try:
            setattr(aligner, name, value)
        except AttributeError:
            pass

    alignment = aligner.align(sequence, human)[0]
    mapping: dict[int, object] = {}
    for (s_start, s_end), (h_start, h_end) in zip(*alignment.aligned):
        for offset in range(s_end - s_start):
            mapping[int(h_start + offset) + 1] = residues[s_start + offset]
    return mapping


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


def binder_exposure(model, binder_id: str) -> dict[int, float]:
    """Relative SASA of each binder residue, computed on the binder alone.

    An earlier version approximated burial by counting heavy atoms within
    10 A and dividing by a constant. That does not discriminate: in a
    60-residue binder even a fully exposed position has upwards of a hundred
    neighbours, so the metric returned ~0.9 for everything and rejected 37 of
    52 candidate positions as "core".

    SASA is measured on the unbound binder deliberately. A position that is
    exposed in the free binder and contacts the target in the complex is
    exactly the position we want; measuring in the complex would score it as
    buried and discard it.
    """
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "binder.pdb"
        io = PDBIO()
        io.set_structure(model)
        io.save(str(path), select=SingleChain(binder_id))

        structure = freesasa.Structure(str(path))
        result = freesasa.calc(structure)

        exposure: dict[int, float] = {}
        for chain_id, residues in result.residueAreas().items():
            if chain_id != binder_id:
                continue
            for number, area in residues.items():
                try:
                    exposure[int(number)] = (area.relativeTotal
                                             if area.hasRelativeAreas else 1.0)
                except ValueError:
                    continue
    return exposure


class SingleChain(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return (residue.id[0] == " "
                and residue.get_resname() in protein_letters_3to1)

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


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

    exposure = binder_exposure(model, binder_id)

    # Locate the handles that are actually present in this complex, keyed by
    # UniProt position through the alignment rather than by residue number.
    uniprot_to_residue = map_target_to_uniprot(target_chain)

    present = {}
    for position, handle in handles.items():
        residue = uniprot_to_residue.get(position)
        if residue is None:
            continue
        if residue.get_resname() != handle["residue"]:
            # The alignment put a different residue here; trust the sequence
            # over the expectation and skip rather than mis-assign.
            continue
        atoms = [residue[n] for n in TITRATABLE_ATOMS[handle["residue"]]
                 if n in residue]
        if atoms:
            present[position] = {
                **handle,
                "coord": np.mean([a.coord for a in atoms], axis=0),
            }

    if not present:
        print(f"  {path.stem}: no target handle found in the complex "
              f"({len(uniprot_to_residue)} target residues mapped)")
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
            rel_sasa = exposure.get(residue.id[1], 1.0)

            blocked = PROTECTED.get(current)
            if rel_sasa < MIN_EXPOSURE:
                blocked = f"core of the binder (rel SASA {rel_sasa:.2f})"
            # A histidine already facing a handle needs no substitution: it
            # is an existing switch, and discarding it would throw away the
            # cheapest protonation event available.
            already = (current == "H" and blocked is None)

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
                "rel_sasa": round(rel_sasa, 2),
                "already_histidine": already,
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
                    (f"{p['current_aa']}{p['binder_position']}"
                     if p.get("already_histidine")
                     else f"{p['current_aa']}{p['binder_position']}H")
                    for p in picks),
                "n_new_histidines": sum(
                    0 if p.get("already_histidine") else 1 for p in picks),
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
        w(f"- Positions below {MIN_EXPOSURE} relative SASA in the unbound "
          f"binder are excluded as core: a charge introduced there "
          f"destabilises the fold instead of forming a switch. SASA is "
          f"measured on the binder alone, so a position that is exposed "
          f"when free and contacts the target when bound still qualifies.\n")
        w("- A position that is already histidine and faces a handle is kept "
          "as an existing switch rather than discarded, since it supplies a "
          "protonation event at no cost in substitutions.\n")
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

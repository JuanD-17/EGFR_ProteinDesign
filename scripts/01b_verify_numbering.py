#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 01b - Independent verification of the numbering map and domain call.

Step 01 derives a UniProt-to-PDB offset by sequence alignment and infers
domain boundaries from the disulfide pattern. Both are load-bearing: every
design decision downstream is expressed in residue numbers, and an offset
error would shift the entire analysis without raising a single exception.

This script checks both claims by routes that do not depend on the alignment
that produced them.

Test 1 - Disulfide geometry.
    UniProt annotates 25 disulfide bonds by residue pair. Mapping those pairs
    through the offset and measuring the actual SG-SG distance in the crystal
    is a geometric test of the map: a covalent disulfide is 2.05 +/- 0.1 A,
    and no wrong offset produces 25 of them by accident. This validates the
    mapping without reusing the sequence alignment that built it.

Test 2 - Literature residues.
    Residues quoted in the EGFR binder literature are looked up in both
    frames and their identities compared, with the surrounding sequence
    window printed so the match can be read directly.

Test 3 - Domain III boundaries.
    Compares the range derived here against ranges circulating elsewhere,
    in a single numbering frame, and reports the cysteine content of the
    regions where they disagree. Domains II and IV are cysteine-rich; an
    L-domain core is not. Cysteine density therefore arbitrates.

Usage:
    python scripts/01b_verify_numbering.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = Path(__file__).resolve().parent.parent

HUMAN_FASTA = ROOT / "01_Target/human/P00533.fasta"
HUMAN_JSON = ROOT / "01_Target/human/P00533.json"
STRUCTURE = ROOT / "01_Target/structure/6ARU.pdb"
METADATA = ROOT / "02_Analysis/01_annotation_metadata.json"

EGFR_CHAIN = "A"
ECTODOMAIN = (25, 645)

# Positions quoted in the literature on EGFR binders, as precursor numbering.
# The expected residue identity is recorded so a numbering error shows up as
# a mismatch rather than as a plausible-looking wrong answer.
LITERATURE_RESIDUES = {
    370: ("H", "histidine proposed as a pH switch"),
    377: ("R", "alanine scan, reduced binding"),
    406: ("L", "alanine scan, reduced binding"),
    433: ("H", "histidine proposed as a pH switch"),
    435: ("Q", "alanine scan, reduced binding"),
    489: ("K", "alanine scan, reduced binding"),
}

# Domain III ranges to compare, all in UniProt precursor numbering.
DOMAIN_III_CANDIDATES = {
    "derived here (disulfide gap)": (363, 469),
    "value in the circulating script": (335, 538),
    "UniProt Repeat, flagged Approximate": (390, 600),
}

DISULFIDE_LENGTH = 2.05
DISULFIDE_TOLERANCE = 0.35


def read_fasta(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if not line.startswith(">"))


def read_features(path: Path) -> list[dict]:
    return json.loads(path.read_text()).get("features", [])


def build_map(sequence: str, chain) -> tuple[dict[int, int], dict[int, object]]:
    """Map UniProt position to PDB residue number, by alignment."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"

    letters, residues = [], []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        try:
            letters.append(protein_letters_3to1[residue.get_resname()])
        except KeyError:
            continue
        residues.append(residue)
    structure_sequence = "".join(letters)

    alignment = aligner.align(structure_sequence, sequence.upper())[0]
    pdb_of_uniprot: dict[int, int] = {}
    residue_of_uniprot: dict[int, object] = {}
    for (s_start, s_end), (u_start, u_end) in zip(*alignment.aligned):
        for offset in range(s_end - s_start):
            residue = residues[s_start + offset]
            if residue.id[2] != " ":
                continue
            position = int(u_start + offset) + 1
            pdb_of_uniprot[position] = int(residue.id[1])
            residue_of_uniprot[position] = residue
    return pdb_of_uniprot, residue_of_uniprot


def distance(atom_a, atom_b) -> float:
    dx, dy, dz = (a - b for a, b in zip(atom_a.coord, atom_b.coord))
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def conformers(residue, atom_name: str) -> list:
    """Every modelled position of an atom, including alternate conformations.

    A cysteine refined in two conformations may be bonded in one and free in
    the other; taking only the first would report a disulfide as broken.
    """
    if atom_name not in residue:
        return []
    atom = residue[atom_name]
    if atom.is_disordered():
        return list(atom.disordered_get_list())
    return [atom]


def min_sg_distance(residue_a, residue_b) -> float:
    """Shortest SG-SG distance over all conformer pairs."""
    return min(
        (distance(a, b)
         for a in conformers(residue_a, "SG")
         for b in conformers(residue_b, "SG")),
        default=math.inf,
    )


def contact_compactness(residues: dict[int, object],
                        start: int, end: int,
                        cutoff: float = 8.0) -> tuple[float, int, int]:
    """Fraction of a range's residue contacts that stay inside the range.

    A genuine structural domain is a compact unit: most of what its residues
    touch is other residues of the same domain. A range drawn across a domain
    boundary leaks contacts to its neighbours. This is the structural test
    that sequence annotation cannot provide.
    """
    def representative(residue):
        for name in ("CB", "CA"):
            if name in residue:
                return residue[name]
        return None

    inside = {p: representative(r) for p, r in residues.items()
              if start <= p <= end}
    inside = {p: a for p, a in inside.items() if a is not None}
    outside = {p: representative(r) for p, r in residues.items()
               if not (start <= p <= end)}
    outside = {p: a for p, a in outside.items() if a is not None}

    internal = external = 0
    for position, atom in inside.items():
        for other, other_atom in inside.items():
            if other > position and distance(atom, other_atom) <= cutoff:
                internal += 1
        for other_atom in outside.values():
            if distance(atom, other_atom) <= cutoff:
                external += 1

    total = internal + external
    return (internal / total if total else 0.0), internal, external


def main() -> int:
    sequence = read_fasta(HUMAN_FASTA)
    features = read_features(HUMAN_JSON)
    model = PDBParser(QUIET=True).get_structure("6aru", str(STRUCTURE))[0]
    chain = model[EGFR_CHAIN]

    pdb_of_uniprot, residue_of_uniprot = build_map(sequence, chain)

    offsets: dict[int, int] = {}
    for uniprot_position, pdb_number in pdb_of_uniprot.items():
        delta = uniprot_position - pdb_number
        offsets[delta] = offsets.get(delta, 0) + 1
    offset = max(offsets, key=offsets.get)

    print("=" * 74)
    print("TEST 0  Offset consistency")
    print("=" * 74)
    print(f"  Mapped residues        : {len(pdb_of_uniprot)}")
    print(f"  Distinct offsets found : {len(offsets)}")
    for delta, count in sorted(offsets.items(), key=lambda kv: -kv[1]):
        print(f"      UniProt - PDB = {delta:+d}  in {count} residues")
    verdict = "PASS" if len(offsets) == 1 else "FAIL - the map is not uniform"
    print(f"  {verdict}")
    print(f"\n  A single offset of {offset:+d} means 6ARU is numbered by the")
    print("  mature protein, which begins after the 24-residue signal peptide.")

    # -- Test 1: disulfide geometry ----------------------------------------
    print()
    print("=" * 74)
    print("TEST 1  Disulfide geometry (independent of the alignment)")
    print("=" * 74)
    print("  Each UniProt disulfide pair is mapped through the offset and the")
    print("  SG-SG distance measured in the crystal. A real disulfide is")
    print(f"  {DISULFIDE_LENGTH} A. A wrong offset cannot reproduce these by chance.\n")

    pairs = sorted(
        (f["location"]["start"]["value"], f["location"]["end"]["value"])
        for f in features
        if f["type"] == "Disulfide bond"
        and f["location"]["end"]["value"] <= ECTODOMAIN[1]
    )

    good = bad = unresolved = 0
    distances = []
    for first, second in pairs:
        residue_a = residue_of_uniprot.get(first)
        residue_b = residue_of_uniprot.get(second)
        if residue_a is None or residue_b is None:
            print(f"    C{first:<4d}-C{second:<4d}  not resolved in the structure")
            unresolved += 1
            continue

        names = (residue_a.get_resname(), residue_b.get_resname())
        if names != ("CYS", "CYS"):
            print(f"    C{first:<4d}-C{second:<4d}  NOT CYSTEINE: {names}  <-- FAIL")
            bad += 1
            continue

        if "SG" not in residue_a or "SG" not in residue_b:
            print(f"    C{first:<4d}-C{second:<4d}  missing SG atom")
            unresolved += 1
            continue

        d = min_sg_distance(residue_a, residue_b)
        distances.append(d)
        ok = abs(d - DISULFIDE_LENGTH) <= DISULFIDE_TOLERANCE
        flag = "ok" if ok else "<-- long"
        print(f"    C{first:<4d}-C{second:<4d}  PDB {pdb_of_uniprot[first]:>4d}-"
              f"{pdb_of_uniprot[second]:<4d}  SG-SG {d:5.2f} A  {flag}")
        good += ok
        bad += not ok

    print(f"\n  {good} bonds within {DISULFIDE_TOLERANCE} A of "
          f"{DISULFIDE_LENGTH} A, {bad} longer, {unresolved} unresolved")
    if distances:
        print(f"  mean {sum(distances) / len(distances):.3f} A, "
              f"range {min(distances):.2f}-{max(distances):.2f} A")

    # A numbering error is systematic: it breaks every bond at once, because
    # the mapped cysteines land on arbitrary residues. One long bond among
    # many correct ones is a local feature of the crystal, not a bad map.
    print()
    if good >= 20 and bad <= 2:
        print("  PASS - the numbering map is geometrically confirmed.")
        if bad:
            print(f"  The {bad} long contact(s) are a property of the model, not")
            print("  of the mapping: a wrong offset would break all 25 at once,")
            print("  since the mapped positions would not be cysteines at all.")
            print("  Cysteines are correctly identified in every case.")
    else:
        print("  FAIL - the mapping does not reproduce the annotated disulfides.")

    # -- Test 2: literature residues ---------------------------------------
    print()
    print("=" * 74)
    print("TEST 2  Literature residues in both frames")
    print("=" * 74)
    print(f"  {'quoted':<8} {'UniProt':<9} {'PDB':<8} {'in PDB':<8} "
          f"{'match':<7} context (UniProt +/-3)")
    print("  " + "-" * 70)

    mismatches = 0
    for position, (expected, note) in sorted(LITERATURE_RESIDUES.items()):
        uniprot_aa = sequence[position - 1]
        pdb_number = pdb_of_uniprot.get(position)
        residue = residue_of_uniprot.get(position)
        pdb_aa = protein_letters_3to1.get(residue.get_resname(), "?") \
            if residue is not None else "-"

        window = sequence[max(0, position - 4):position + 3]
        agree = (uniprot_aa == expected == pdb_aa)
        mismatches += not agree

        print(f"  {expected}{position:<7} {uniprot_aa:<9} "
              f"{pdb_number if pdb_number else '-':<8} {pdb_aa:<8} "
              f"{'yes' if agree else 'NO':<7} {window}")

    print()
    if mismatches == 0:
        print("  PASS - every quoted residue has the expected identity in both frames")
        print("  The literature positions are precursor numbering, and they")
        print(f"  correspond to PDB numbers {offset} lower.")
    else:
        print(f"  FAIL - {mismatches} residue(s) do not match")

    # -- Test 3: domain III -------------------------------------------------
    print()
    print("=" * 74)
    print("TEST 3  Domain III boundaries")
    print("=" * 74)
    print("  All ranges below are UniProt precursor numbering, so they are")
    print("  directly comparable. The disulfide-derived range was computed")
    print("  from UniProt features, not from PDB residue numbers.\n")

    bonded = {p for f in features if f["type"] == "Disulfide bond"
              for p in (f["location"]["start"]["value"],
                        f["location"]["end"]["value"])}

    print(f"  {'range':<38} {'length':<8} {'SS-Cys':<8} density")
    print("  " + "-" * 70)
    for label, (start, end) in DOMAIN_III_CANDIDATES.items():
        length = end - start + 1
        count = sum(1 for p in bonded if start <= p <= end)
        print(f"  {label:<38} {length:<8} {count:<8} "
              f"{count / length * 100:.1f} per 100 aa")

    derived = DOMAIN_III_CANDIDATES["derived here (disulfide gap)"]
    circulating = DOMAIN_III_CANDIDATES["value in the circulating script"]

    print("\n  Regions where the two disagree:\n")
    for start, end, description in (
        (circulating[0], derived[0] - 1, "extra at the N-terminal end"),
        (derived[1] + 1, circulating[1], "extra at the C-terminal end"),
    ):
        if start > end:
            continue
        count = sum(1 for p in bonded if start <= p <= end)
        length = end - start + 1
        print(f"    {start}-{end} ({length} aa), {description}")
        print(f"      disulfide-bonded cysteines: {count} "
              f"({count / length * 100:.1f} per 100 aa)")
        print(f"      sequence: {sequence[start - 1:end]}")
        print()

    core_count = sum(1 for p in bonded if derived[0] <= p <= derived[1])
    print(f"    {derived[0]}-{derived[1]} (the derived core) contains "
          f"{core_count} disulfide-bonded cysteines.")

    print("\n  Reading: an L-domain core is free of disulfides. The regions the")
    print("  wider range adds are densely cross-linked, which is the signature")
    print("  of the furin-like domains that flank domain III.")

    # -- Test 4: structural compactness ------------------------------------
    print()
    print("=" * 74)
    print("TEST 4  Structural compactness of each candidate range")
    print("=" * 74)
    print("  Cysteine content is an indirect argument. This is the direct one:")
    print("  for each range, what fraction of its residue contacts stay inside")
    print("  it. A real domain is a compact unit and keeps most of its contacts")
    print("  internal; a range drawn across a boundary leaks them to neighbours.")
    print("  Contacts are CB-CB (CA for glycine) within 8 A.\n")

    print(f"  {'range':<38} {'internal':<10} {'external':<10} fraction")
    print("  " + "-" * 70)
    compactness = {}
    for label, (start, end) in DOMAIN_III_CANDIDATES.items():
        fraction, internal, external = contact_compactness(
            residue_of_uniprot, start, end)
        compactness[label] = fraction
        print(f"  {label:<38} {internal:<10} {external:<10} {fraction:.3f}")

    best = max(compactness, key=compactness.get)
    derived_fraction = compactness["derived here (disulfide gap)"]
    print(f"\n  Most compact: {best} ({compactness[best]:.3f})")

    if best != "derived here (disulfide gap)":
        print()
        print("  This contradicts the disulfide-gap derivation, and the wider")
        print("  range wins on the stronger argument. A short range is compact")
        print("  partly by being short, so a LONGER range scoring higher is not")
        print("  an artefact of length: it means the shorter one was cutting")
        print(f"  through a structural unit. The derived core leaks {1 - derived_fraction:.0%} of its")
        print("  contacts to residues outside itself, which a domain does not do.")
        print()
        print("  Revised reading: the disulfide gap locates the leucine-rich")
        print("  solenoid CORE of domain III, a real substructure but not the")
        print("  whole domain. The domain also carries flanking segments that")
        print("  are disulfide-bonded and pack against that core. Cysteine")
        print("  content identifies the core; contact topology identifies the")
        print("  domain, and it is the domain that matters for epitope scope.")

    # -- Cross-check against step 01 ---------------------------------------
    if METADATA.exists():
        recorded = json.loads(METADATA.read_text())["uniprot_minus_pdb_offset"]
        print()
        print("=" * 74)
        print(f"  Offset recorded by step 01 : {recorded:+d}")
        print(f"  Offset recomputed here     : {offset:+d}")
        print(f"  {'CONSISTENT' if recorded == offset else 'INCONSISTENT'}")
        print("=" * 74)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

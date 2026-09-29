#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 01 - Per-residue annotation of the EGFR ectodomain.

Builds the table that every downstream decision rests on: for each residue of
the human EGFR extracellular region, whether it is exposed, whether it is
conserved in mouse, whether it is buried under a glycan, and where it sits
structurally.

Four layers are computed and crossed:

  1. Conservation      human vs mouse, from a global pairwise alignment
  2. Accessibility     solvent-accessible surface area from the crystal
                       structure, absolute and relative to the residue maximum
  3. Glycan occlusion  N-glycosylation sites and the surface they shadow
  4. Architecture      domain assignment derived from the disulfide pattern

The output is a candidate epitope surface: residues that are exposed,
conserved across species, and not covered in sugar. Which domain that surface
falls in is a result of the analysis, not an input to it.

Two deliberate methodological choices are worth stating.

Domain boundaries are derived, not assumed. UniProt does not annotate the
ectodomain subdomains I-IV; its only Domain feature is the intracellular
kinase, and the two Repeat features covering the L-domains are flagged
"Approximate". Rather than hard-coding a range from the literature, domains
are inferred from the disulfide bond distribution: domains II and IV are
cysteine-rich furin-like modules, while domains I and III are leucine-rich
L-domains almost free of disulfides. The cysteine-free stretches are therefore
a structural signature that can be read off the data.

The cetuximab interface is computed but treated as an exclusion, not a target.
Cetuximab is human-specific and does not bind murine EGFR, so its epitope is
by definition poorly conserved. It is measured here to be avoided, and to
serve as a positive control that the geometry code is working.

Inputs:  01_Target/{human,mouse}/*.{fasta,json}, 01_Target/structure/6ARU.pdb
Outputs: 02_Analysis/01_residue_annotation.csv
         02_Analysis/01_annotation_summary.md
         02_Analysis/01_annotation_metadata.json

All positions are in UniProt precursor numbering unless stated otherwise.

Usage:
    python scripts/01_annotate_target.py
    python scripts/01_annotate_target.py --glycan-radius 15
"""

from __future__ import annotations

import argparse
import json
import math
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import freesasa
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

# freesasa prints a warning per unrecognised atom; the structure carries
# glycans and waters we deliberately exclude, so the stream is noise.
freesasa.setVerbosity(freesasa.silent)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "02_Analysis"

HUMAN_FASTA = ROOT / "01_Target/human/P00533.fasta"
HUMAN_JSON = ROOT / "01_Target/human/P00533.json"
MOUSE_FASTA = ROOT / "01_Target/mouse/Q01279.fasta"
MOUSE_JSON = ROOT / "01_Target/mouse/Q01279.json"
STRUCTURE = ROOT / "01_Target/structure/6ARU.pdb"

# Extracellular region, from the UniProt topological domain feature.
ECTODOMAIN = (25, 645)

EGFR_CHAIN = "A"
FAB_CHAINS = ("B", "C")

# A residue counts as surface-exposed above this fraction of its maximum
# accessible area. The 20-25% band is the usual convention; below it a side
# chain is buried and cannot be contacted by a binder.
EXPOSURE_THRESHOLD = 0.20

# Radius around a glycosylation site considered occluded. A complex N-glycan
# is a branched tree extending well beyond the asparagine it hangs from, and
# it is almost never resolved in a crystal structure, so the surface looks
# deceptively clear. 12 A is conservative; --glycan-radius overrides it.
GLYCAN_RADIUS = 12.0

# Buried surface area at which a residue counts as part of the Fab interface.
EPITOPE_DSASA = 1.0

# A gap this long between consecutive disulfide-bonded cysteines marks a
# cysteine-free L-domain core.
DISULFIDE_GAP = 80

GLYCAN_RESNAMES = {"NAG", "MAN", "BMA", "FUC", "GAL", "SIA", "GLC", "XYS"}


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def read_features(path: Path) -> list[dict]:
    return json.loads(path.read_text()).get("features", [])


def glycosylation_sites(features: list[dict], span: tuple[int, int]) -> dict[int, str]:
    """N-glycosylation positions within a span, mapped to their description."""
    sites = {}
    for feature in features:
        if feature["type"] != "Glycosylation":
            continue
        position = feature["location"]["start"]["value"]
        if span[0] <= position <= span[1]:
            sites[position] = feature.get("description", "")
    return sites


def disulfide_pairs(features: list[dict], span: tuple[int, int]) -> list[tuple[int, int]]:
    pairs = []
    for feature in features:
        if feature["type"] != "Disulfide bond":
            continue
        start = feature["location"]["start"]["value"]
        end = feature["location"]["end"]["value"]
        if span[0] <= start and end <= span[1]:
            pairs.append((start, end))
    return sorted(pairs)


# ---------------------------------------------------------------------------
# Layer 4: architecture from the disulfide pattern
# ---------------------------------------------------------------------------

def derive_domains(pairs: list[tuple[int, int]],
                   span: tuple[int, int],
                   min_gap: int) -> list[dict]:
    """Partition the ectodomain into cysteine-rich and cysteine-free regions.

    Domains II and IV of EGFR are furin-like modules packed with disulfides;
    domains I and III are L-domains that carry almost none. Long gaps in the
    disulfide positions therefore locate the L-domain cores without needing an
    external annotation.
    """
    positions = sorted({p for pair in pairs for p in pair})
    if not positions:
        return []

    # Find the gaps that separate disulfide clusters.
    boundaries: list[tuple[int, int]] = []
    for earlier, later in zip(positions, positions[1:]):
        if later - earlier >= min_gap:
            boundaries.append((earlier + 1, later - 1))

    regions: list[dict] = []
    cursor = span[0]
    for gap_start, gap_end in boundaries:
        if gap_start > cursor:
            regions.append({"start": cursor, "end": gap_start - 1,
                            "character": "cysteine-rich"})
        regions.append({"start": gap_start, "end": gap_end,
                        "character": "cysteine-free"})
        cursor = gap_end + 1
    if cursor <= span[1]:
        regions.append({"start": cursor, "end": span[1],
                        "character": "cysteine-rich"})

    # Label them. The ectodomain alternates L-domain / cysteine-rich starting
    # from domain I, so cysteine-free cores are domains I and III in order.
    l_domain_labels = ["I (L-domain core)", "III (L-domain core)"]
    rich_labels = ["I/II (N-terminal, cysteine-rich)",
                   "II (furin-like)", "IV (furin-like)"]
    l_index = rich_index = 0
    for region in regions:
        if region["character"] == "cysteine-free" and l_index < len(l_domain_labels):
            region["domain"] = l_domain_labels[l_index]
            l_index += 1
        elif region["character"] == "cysteine-rich" and rich_index < len(rich_labels):
            region["domain"] = rich_labels[rich_index]
            rich_index += 1
        else:
            region["domain"] = "unassigned"
        region["length"] = region["end"] - region["start"] + 1
    return regions


def domain_of(position: int, regions: list[dict]) -> str:
    for region in regions:
        if region["start"] <= position <= region["end"]:
            return region["domain"]
    return "outside"


# ---------------------------------------------------------------------------
# Layer 1: conservation
# ---------------------------------------------------------------------------

def build_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    return aligner


def aligned_pairs(aligner: Align.PairwiseAligner,
                  query: str, subject: str) -> dict[int, int]:
    """Map 0-based index in `query` to 0-based index in `subject`."""
    alignment = aligner.align(query.upper(), subject.upper())[0]
    mapping: dict[int, int] = {}
    for (q_start, q_end), (s_start, s_end) in zip(*alignment.aligned):
        for offset in range(q_end - q_start):
            mapping[q_start + offset] = s_start + offset
    return mapping


def classify_conservation(human_aa: str, mouse_aa: str | None,
                          matrix) -> str:
    if mouse_aa is None:
        return "unaligned"
    if human_aa == mouse_aa:
        return "identical"
    try:
        score = matrix[human_aa, mouse_aa]
    except (KeyError, IndexError):
        return "different"
    return "similar" if score >= 0 else "different"


# ---------------------------------------------------------------------------
# Structure handling
# ---------------------------------------------------------------------------

def load_structure(path: Path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("6aru", str(path))[0]


def chain_sequence(chain) -> tuple[str, list]:
    """Amino acid sequence of a chain plus the residue objects, in order."""
    letters = []
    residues = []
    for residue in chain:
        if residue.id[0] != " ":          # skip HETATM (water, glycans)
            continue
        try:
            letters.append(protein_letters_3to1[residue.get_resname()])
        except KeyError:
            continue
        residues.append(residue)
    return "".join(letters), residues


class ProteinChainSelect(Select):
    """Keep standard amino acid residues of the given chains, heavy atoms only.

    Writing PDB records by hand is a trap: the format is column-positional and
    a single-character offset silently corrupts residue names. PDBIO already
    gets this right.
    """

    def __init__(self, chains: tuple[str, ...]):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        return (residue.id[0] == " "
                and residue.get_resname() in protein_letters_3to1)

    def accept_atom(self, atom):
        if atom.element == "H":
            return False
        # Keep one conformer of alternate locations.
        return atom.get_altloc() in (" ", "A")


def write_subset_pdb(structure, chains: tuple[str, ...], destination: Path) -> None:
    """Write selected protein chains to a PDB file via Biopython."""
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(destination), select=ProteinChainSelect(chains))


def residue_sasa(pdb_path: Path) -> dict[tuple[str, int], tuple[float, float | None]]:
    """Per-residue SASA from a PDB file: (absolute, relative)."""
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    areas: dict[tuple[str, int], tuple[float, float | None]] = {}
    for chain_id, residues in result.residueAreas().items():
        for number, area in residues.items():
            try:
                key = (chain_id, int(number))
            except ValueError:
                continue
            relative = area.relativeTotal if area.hasRelativeAreas else None
            areas[key] = (area.total, relative)
    return areas


def glycan_anchor_atoms(model, chain_id: str,
                        glyco_positions: set[int],
                        pdb_of_uniprot: dict[int, int]) -> list[tuple]:
    """Coordinates that mark where a glycan tree hangs off the surface.

    Uses the modelled sugar residues when the crystal resolved them, and falls
    back to the ND2 atom of the glycosylated asparagine when it did not.
    """
    anchors = []
    chain = model[chain_id]

    for residue in chain:
        if residue.id[0].startswith("H_") and \
                residue.get_resname() in GLYCAN_RESNAMES:
            for atom in residue:
                anchors.append((atom.coord, f"{residue.get_resname()}{residue.id[1]}"))

    modelled = len(anchors)

    for uniprot_position in sorted(glyco_positions):
        pdb_number = pdb_of_uniprot.get(uniprot_position)
        if pdb_number is None:
            continue
        try:
            residue = chain[(" ", pdb_number, " ")]
        except KeyError:
            continue
        atom = residue["ND2"] if "ND2" in residue else \
            (residue["CB"] if "CB" in residue else None)
        if atom is not None:
            anchors.append((atom.coord, f"N{uniprot_position}"))

    return anchors, modelled


def minimum_distance(residue, anchors) -> float:
    best = math.inf
    for atom in residue:
        if atom.element == "H":
            continue
        ax, ay, az = atom.coord
        for coord, _ in anchors:
            dx, dy, dz = ax - coord[0], ay - coord[1], az - coord[2]
            distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            if distance < best:
                best = distance
    return best


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--glycan-radius", type=float, default=GLYCAN_RADIUS,
                        help=f"occlusion radius around glycans in A "
                             f"(default {GLYCAN_RADIUS})")
    parser.add_argument("--exposure", type=float, default=EXPOSURE_THRESHOLD,
                        help=f"relative SASA above which a residue counts as "
                             f"exposed (default {EXPOSURE_THRESHOLD})")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # -- Sequences and annotations -----------------------------------------
    print("Reading sequences and UniProt annotations")
    human_sequence = read_fasta(HUMAN_FASTA)
    mouse_sequence = read_fasta(MOUSE_FASTA)
    human_features = read_features(HUMAN_JSON)

    glycosites = glycosylation_sites(human_features, ECTODOMAIN)
    disulfides = disulfide_pairs(human_features, ECTODOMAIN)
    cysteines_in_bonds = {p for pair in disulfides for p in pair}
    print(f"  human {len(human_sequence)} aa, mouse {len(mouse_sequence)} aa")
    print(f"  {len(glycosites)} glycosylation sites, "
          f"{len(disulfides)} disulfide bonds in the ectodomain")

    # -- Layer 4: architecture ---------------------------------------------
    print("Deriving domain architecture from the disulfide pattern")
    regions = derive_domains(disulfides, ECTODOMAIN, DISULFIDE_GAP)
    for region in regions:
        print(f"  {region['start']:>4}-{region['end']:<4} "
              f"{region['length']:>4} aa  {region['character']:<14} "
              f"{region['domain']}")

    # -- Layer 1: conservation ---------------------------------------------
    print("Aligning human and mouse sequences")
    aligner = build_aligner()
    matrix = aligner.substitution_matrix
    human_to_mouse = aligned_pairs(aligner, human_sequence, mouse_sequence)

    # -- Structure ----------------------------------------------------------
    print(f"Reading structure {STRUCTURE.name}")
    model = load_structure(STRUCTURE)
    egfr_sequence, egfr_residues = chain_sequence(model[EGFR_CHAIN])
    print(f"  chain {EGFR_CHAIN}: {len(egfr_residues)} resolved residues")

    structure_to_human = aligned_pairs(aligner, egfr_sequence, human_sequence)
    identity = sum(
        1 for i, j in structure_to_human.items()
        if egfr_sequence[i] == human_sequence[j]
    )
    print(f"  {identity}/{len(egfr_sequence)} residues match P00533 "
          f"({100 * identity / len(egfr_sequence):.1f}%)")

    uniprot_of_pdb: dict[int, int] = {}
    pdb_of_uniprot: dict[int, int] = {}
    residue_by_uniprot: dict[int, object] = {}
    for structure_index, human_index in structure_to_human.items():
        residue = egfr_residues[structure_index]
        if residue.id[2] != " ":          # skip insertion codes
            continue
        # Biopython returns numpy integers from the alignment; cast so the
        # values stay JSON-serialisable downstream.
        uniprot_position = int(human_index) + 1
        pdb_number = int(residue.id[1])
        uniprot_of_pdb[pdb_number] = uniprot_position
        pdb_of_uniprot[uniprot_position] = pdb_number
        residue_by_uniprot[uniprot_position] = residue

    offsets = defaultdict(int)
    for pdb_number, uniprot_position in uniprot_of_pdb.items():
        offsets[uniprot_position - pdb_number] += 1
    dominant_offset, count = max(offsets.items(), key=lambda kv: kv[1])
    dominant_offset, count = int(dominant_offset), int(count)
    print(f"  numbering offset UniProt - PDB = {dominant_offset:+d} "
          f"({count}/{len(uniprot_of_pdb)} residues)")

    # -- Layer 2: accessibility --------------------------------------------
    print("Computing solvent accessibility")
    with tempfile.TemporaryDirectory() as tmp:
        free_pdb = Path(tmp) / "egfr_free.pdb"
        complex_pdb = Path(tmp) / "complex.pdb"
        write_subset_pdb(model, (EGFR_CHAIN,), free_pdb)
        write_subset_pdb(model, (EGFR_CHAIN,) + FAB_CHAINS, complex_pdb)
        sasa_free = residue_sasa(free_pdb)
        sasa_complex = residue_sasa(complex_pdb)
    print(f"  {len(sasa_free)} residues with SASA in the free receptor")

    # -- Layer 3: glycan occlusion -----------------------------------------
    print("Mapping glycan occlusion")
    anchors, modelled = glycan_anchor_atoms(
        model, EGFR_CHAIN, set(glycosites), pdb_of_uniprot)
    print(f"  {modelled} modelled sugar atoms, "
          f"{len(anchors) - modelled} asparagine anchors")

    # -- Assemble the table -------------------------------------------------
    print("Assembling per-residue table")
    rows = []
    for uniprot_position in range(ECTODOMAIN[0], ECTODOMAIN[1] + 1):
        human_aa = human_sequence[uniprot_position - 1]
        mouse_index = human_to_mouse.get(uniprot_position - 1)
        mouse_aa = mouse_sequence[mouse_index] if mouse_index is not None else None

        residue = residue_by_uniprot.get(uniprot_position)
        pdb_number = pdb_of_uniprot.get(uniprot_position)

        free = sasa_free.get((EGFR_CHAIN, pdb_number)) if pdb_number else None
        bound = sasa_complex.get((EGFR_CHAIN, pdb_number)) if pdb_number else None

        sasa_value = free[0] if free else None
        relative = free[1] if free else None
        delta = (free[0] - bound[0]) if (free and bound) else None

        glycan_distance = minimum_distance(residue, anchors) \
            if (residue is not None and anchors) else None

        exposed = relative is not None and relative >= args.exposure
        shadowed = glycan_distance is not None and glycan_distance <= args.glycan_radius
        conservation = classify_conservation(human_aa, mouse_aa, matrix)

        rows.append({
            "uniprot_pos": uniprot_position,
            "aa_human": human_aa,
            "aa_mouse": mouse_aa or "-",
            "mouse_pos": (mouse_index + 1) if mouse_index is not None else "",
            "conservation": conservation,
            "domain": domain_of(uniprot_position, regions),
            "pdb_resnum": pdb_number if pdb_number else "",
            "resolved": "yes" if residue is not None else "no",
            "sasa_free": round(sasa_value, 1) if sasa_value is not None else "",
            "rel_sasa": round(relative, 3) if relative is not None else "",
            "exposed": "yes" if exposed else "no",
            "dsasa_cetuximab": round(delta, 1) if delta is not None else "",
            "cetuximab_epitope": "yes" if (delta or 0) >= EPITOPE_DSASA else "no",
            "is_glycosite": "yes" if uniprot_position in glycosites else "no",
            "dist_to_glycan": round(glycan_distance, 1)
                              if glycan_distance is not None and
                              glycan_distance < math.inf else "",
            "glycan_shadow": "yes" if shadowed else "no",
            "ss_bonded_cys": "yes" if uniprot_position in cysteines_in_bonds else "no",
            "designable": "yes" if (
                exposed
                and not shadowed
                and conservation == "identical"
                and uniprot_position not in cysteines_in_bonds
            ) else "no",
        })

    # -- Write outputs ------------------------------------------------------
    import csv

    csv_path = OUT_DIR / "01_residue_annotation.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    resolved = [r for r in rows if r["resolved"] == "yes"]
    exposed_rows = [r for r in resolved if r["exposed"] == "yes"]
    designable = [r for r in rows if r["designable"] == "yes"]
    epitope = [r for r in rows if r["cetuximab_epitope"] == "yes"]

    conservation_counts = defaultdict(int)
    for row in rows:
        conservation_counts[row["conservation"]] += 1

    epitope_conservation = defaultdict(int)
    for row in epitope:
        epitope_conservation[row["conservation"]] += 1

    per_domain = defaultdict(lambda: {"total": 0, "designable": 0})
    for row in rows:
        per_domain[row["domain"]]["total"] += 1
        if row["designable"] == "yes":
            per_domain[row["domain"]]["designable"] += 1

    summary_path = OUT_DIR / "01_annotation_summary.md"
    with summary_path.open("w") as handle:
        w = handle.write
        w("# Step 01 - EGFR ectodomain annotation\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/01_annotate_target.py`.\n\n")
        w(f"Ectodomain residues {ECTODOMAIN[0]}-{ECTODOMAIN[1]} "
          f"({ECTODOMAIN[1] - ECTODOMAIN[0] + 1} aa), UniProt precursor "
          f"numbering. Structure 6ARU chain {EGFR_CHAIN}, "
          f"numbering offset UniProt - PDB = {dominant_offset:+d}.\n\n")

        w("## Domain architecture, derived from the disulfide pattern\n\n")
        w("UniProt does not annotate the ectodomain subdomains. These "
          "boundaries come from the distribution of disulfide bonds: "
          "cysteine-free stretches mark the L-domain cores.\n\n")
        w("| Range | Length | Character | Assignment |\n")
        w("|-------|--------|-----------|------------|\n")
        for region in regions:
            w(f"| {region['start']}-{region['end']} | {region['length']} aa "
              f"| {region['character']} | {region['domain']} |\n")

        w("\n## Coverage\n\n")
        w(f"- Resolved in the crystal structure: {len(resolved)}/{len(rows)}\n")
        w(f"- Exposed (relative SASA >= {args.exposure}): {len(exposed_rows)}\n")
        w(f"- Cetuximab interface: {len(epitope)} residues\n")
        w(f"- **Designable** (exposed, conserved, unglycosylated, not "
          f"disulfide-bonded): {len(designable)}\n")

        w("\n## Human/mouse conservation across the ectodomain\n\n")
        for label in ("identical", "similar", "different", "unaligned"):
            count = conservation_counts.get(label, 0)
            w(f"- {label}: {count} ({100 * count / len(rows):.1f}%)\n")

        w("\n## Conservation of the cetuximab epitope\n\n")
        w("This is the control on the decision to avoid it.\n\n")
        for label in ("identical", "similar", "different", "unaligned"):
            count = epitope_conservation.get(label, 0)
            if epitope:
                w(f"- {label}: {count} ({100 * count / len(epitope):.1f}%)\n")

        w("\n## Designable surface per domain\n\n")
        w("| Domain | Residues | Designable | Fraction |\n")
        w("|--------|----------|------------|----------|\n")
        for domain, counts in sorted(per_domain.items(),
                                     key=lambda kv: -kv[1]["designable"]):
            fraction = counts["designable"] / counts["total"] if counts["total"] else 0
            w(f"| {domain} | {counts['total']} | {counts['designable']} "
              f"| {100 * fraction:.1f}% |\n")

        w("\n## Glycosylation sites in the ectodomain\n\n")
        w("| Position | Domain | Description |\n")
        w("|----------|--------|-------------|\n")
        for position in sorted(glycosites):
            w(f"| N{position} | {domain_of(position, regions)} "
              f"| {glycosites[position]} |\n")

        w(f"\nSurface within {args.glycan_radius} A of a glycan is treated as "
          f"occluded. Glycan trees are rarely resolved crystallographically, "
          f"so this surface looks accessible in the structure but is not.\n")

        w("\n## Histidines on the designable surface\n\n")
        w("Relevant to objective 1. A target histidine near the interface can "
          "pair with a binder carboxylate, adding a second protonation-sensitive "
          "contact on top of the histidines engineered into the binder itself.\n\n")
        w("| Position | Domain | Relative SASA | Conservation | Glycan shadow |\n")
        w("|----------|--------|---------------|--------------|---------------|\n")
        for row in rows:
            if row["aa_human"] == "H" and row["exposed"] == "yes":
                w(f"| H{row['uniprot_pos']} | {row['domain']} "
                  f"| {row['rel_sasa']} | {row['conservation']} "
                  f"| {row['glycan_shadow']} |\n")

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "scripts/01_annotate_target.py",
        "numbering": "UniProt precursor (P00533)",
        "ectodomain": list(ECTODOMAIN),
        "structure": STRUCTURE.name,
        "egfr_chain": EGFR_CHAIN,
        "fab_chains": list(FAB_CHAINS),
        "uniprot_minus_pdb_offset": dominant_offset,
        "parameters": {
            "exposure_threshold": args.exposure,
            "glycan_radius_angstrom": args.glycan_radius,
            "epitope_dsasa_angstrom2": EPITOPE_DSASA,
            "disulfide_gap": DISULFIDE_GAP,
        },
        "domains": regions,
        "counts": {
            "ectodomain_residues": len(rows),
            "resolved": len(resolved),
            "exposed": len(exposed_rows),
            "cetuximab_epitope": len(epitope),
            "designable": len(designable),
        },
        "conservation": dict(conservation_counts),
        "glycosylation_sites": sorted(glycosites),
    }
    metadata_path = OUT_DIR / "01_annotation_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    print(f"\n  {len(rows)} residues annotated")
    print(f"  {len(designable)} designable "
          f"(exposed, conserved, unglycosylated, not disulfide-bonded)")
    print(f"\n  -> {csv_path.relative_to(ROOT)}")
    print(f"  -> {summary_path.relative_to(ROOT)}")
    print(f"  -> {metadata_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

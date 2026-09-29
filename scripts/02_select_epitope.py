#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 02 - Epitope patch selection.

Step 01 produced a pool of residues that a binder is structurally permitted
to touch: exposed, identical between human and mouse, clear of glycan, not
disulfide-bonded. That pool is not a set of hotspots. A binder does not
contact 45 scattered residues; it lands on a contiguous surface of roughly
600-900 A^2, which is on the order of 15-25 residues.

This script turns the pool into ranked candidate patches and emits the
hotspot specification that a binder design pipeline consumes.

Method
------
Every pool residue is taken in turn as a patch centre, and the patch is the
set of pool residues whose CB lies within a radius of it. Patches are then
scored, ranked, and greedily deduplicated so that the selection spans
distinct surfaces rather than 40 views of the same one.

Scoring follows the competition's own ranking order, in which pH selectivity
outranks cross-reactivity, which outranks affinity. Four terms, with the
weights stated in WEIGHTS so they can be argued with:

  ph_potential      Histidines in or adjacent to the patch, weighted by
                    exposure. A target histidine facing a binder carboxylate
                    is a protonation-sensitive contact that the binder does
                    not have to supply by itself.
  conservation      Fraction of the patch neighbourhood identical in mouse.
                    The patch core is identical by construction, but a binder
                    footprint spills past it, so the shell is what decides
                    cross-reactivity.
  functional        Overlap with the cetuximab interface. Treated as positive
                    evidence: the competition recommends a functional epitope
                    and names cetuximab's, and a surface with a therapeutic
                    antibody already bound to it is demonstrably druggable.
                    Non-conserved residues within it are already excluded by
                    the pool definition.
  geometry          How close the buriable area is to the target window, and
                    how far the patch sits from the nearest glycan.

Numbering
---------
Internally everything is UniProt precursor numbering. The BindCraft hotspot
file is emitted in PDB numbering, because BindCraft indexes into the PDB it
is given and 6ARU is numbered by the mature protein, 24 lower. Getting this
backwards would silently target a surface 24 residues away.

Inputs:  02_Analysis/01_residue_annotation.csv
         01_Target/structure/6ARU.pdb
Outputs: 02_Analysis/02_patches.csv
         02_Analysis/02_patch_summary.md
         02_Analysis/02_hotspots_bindcraft.json
         02_Analysis/02_inspect_patches.pml

Usage:
    python scripts/02_select_epitope.py
    python scripts/02_select_epitope.py --scope 157 362      # domain II
    python scripts/02_select_epitope.py --scope 25 645       # whole ectodomain
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "02_Analysis"
ANNOTATION = OUT_DIR / "01_residue_annotation.csv"
STRUCTURE = ROOT / "01_Target/structure/6ARU.pdb"

EGFR_CHAIN = "A"

# Domain III as scoped after the step 01b verification.
DEFAULT_SCOPE = (335, 538)

# Radius around a patch centre, CB to CB. 11 A gathers a surface of roughly
# the size an engineered binder buries.
PATCH_RADIUS = 11.0

# Buriable area window. Below roughly 600 A^2 an interface rarely reaches
# useful affinity; much above 900 A^2 is more surface than a small binder can
# cover.
AREA_WINDOW = (600.0, 900.0)

# A patch needs enough anchor points to be designable against.
MIN_PATCH_RESIDUES = 8

# Shell around the patch used to judge cross-reactivity, since a binder
# footprint extends past the pool residues it is aimed at.
SHELL_RADIUS = 13.0

# Two patches sharing more than this fraction of residues are the same
# surface, and only the better-scoring one is kept.
MAX_OVERLAP = 0.5

WEIGHTS = {
    "ph_potential": 0.40,
    "conservation": 0.30,
    "functional": 0.20,
    "geometry": 0.10,
}


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def load_annotation(path: Path) -> list[dict]:
    """Read step 01's table, coercing the numeric columns."""
    rows = []
    with path.open() as handle:
        for row in csv.DictReader(handle):
            row["uniprot_pos"] = int(row["uniprot_pos"])
            row["pdb_resnum"] = int(row["pdb_resnum"]) if row["pdb_resnum"] else None
            for key in ("sasa_free", "rel_sasa", "dsasa_cetuximab",
                        "dist_to_glycan"):
                row[key] = float(row[key]) if row[key] else None
            rows.append(row)
    return rows


def centroid_atoms(structure_chain, rows: list[dict]) -> dict[int, object]:
    """One representative atom per residue, keyed by UniProt position."""
    atoms = {}
    for row in rows:
        if row["pdb_resnum"] is None:
            continue
        try:
            residue = structure_chain[(" ", row["pdb_resnum"], " ")]
        except KeyError:
            continue
        atom = residue["CB"] if "CB" in residue else (
            residue["CA"] if "CA" in residue else None)
        if atom is not None:
            atoms[row["uniprot_pos"]] = atom
    return atoms


def distance(atom_a, atom_b) -> float:
    dx, dy, dz = (a - b for a, b in zip(atom_a.coord, atom_b.coord))
    return math.sqrt(dx * dx + dy * dy + dz * dz)


# ---------------------------------------------------------------------------
# Patch construction and scoring
# ---------------------------------------------------------------------------

def build_patch(centre: int, pool: list[int], atoms: dict[int, object],
                radius: float) -> list[int]:
    return sorted(
        position for position in pool
        if position in atoms
        and distance(atoms[centre], atoms[position]) <= radius
    )


def shell_of(patch: list[int], candidates: list[int],
             atoms: dict[int, object], radius: float) -> list[int]:
    """Residues near the patch but not in it: the rest of the footprint."""
    members = set(patch)
    shell = []
    for position in candidates:
        if position in members or position not in atoms:
            continue
        if any(distance(atoms[position], atoms[m]) <= radius for m in patch):
            shell.append(position)
    return sorted(shell)


def score_patch(patch: list[int], shell: list[int],
                by_position: dict[int, dict],
                atoms: dict[int, object]) -> dict:
    """Compute the four scoring terms and the composite."""
    members = [by_position[p] for p in patch]
    neighbourhood = members + [by_position[p] for p in shell]

    area = sum(r["sasa_free"] or 0.0 for r in members)

    # -- pH potential ------------------------------------------------------
    # Histidines count in proportion to how exposed they are; one in the
    # patch itself counts double one merely adjacent to it.
    histidines_in = [r for r in members if r["aa_human"] == "H"]
    histidines_near = [by_position[p] for p in shell
                       if by_position[p]["aa_human"] == "H"]
    ph_raw = (sum(2.0 * (r["rel_sasa"] or 0.0) for r in histidines_in)
              + sum(1.0 * (r["rel_sasa"] or 0.0) for r in histidines_near))
    ph_potential = min(1.0, ph_raw / 2.0)

    # -- Conservation of the whole footprint -------------------------------
    identical = sum(1 for r in neighbourhood if r["conservation"] == "identical")
    conservation = identical / len(neighbourhood) if neighbourhood else 0.0

    # -- Functional relevance ----------------------------------------------
    epitope_overlap = sum(1 for r in members if r["cetuximab_epitope"] == "yes")
    functional = min(1.0, epitope_overlap / 8.0)

    # -- Geometry ----------------------------------------------------------
    low, high = AREA_WINDOW
    if low <= area <= high:
        area_fit = 1.0
    elif area < low:
        area_fit = max(0.0, area / low)
    else:
        area_fit = max(0.0, 1.0 - (area - high) / high)

    glycan_distances = [r["dist_to_glycan"] for r in members
                        if r["dist_to_glycan"] is not None]
    nearest_glycan = min(glycan_distances) if glycan_distances else 0.0
    glycan_clearance = min(1.0, nearest_glycan / 20.0)
    geometry = 0.5 * area_fit + 0.5 * glycan_clearance

    terms = {
        "ph_potential": ph_potential,
        "conservation": conservation,
        "functional": functional,
        "geometry": geometry,
    }
    composite = sum(WEIGHTS[name] * value for name, value in terms.items())

    return {
        "residues": patch,
        "n_residues": len(patch),
        "area_a2": round(area, 1),
        "histidines": [r["uniprot_pos"] for r in histidines_in],
        "histidines_adjacent": [r["uniprot_pos"] for r in histidines_near],
        "cetuximab_overlap": epitope_overlap,
        "shell_size": len(shell),
        "nearest_glycan_a": round(nearest_glycan, 1),
        "conservation_identical_pct": round(100 * conservation, 1),
        **{f"score_{k}": round(v, 3) for k, v in terms.items()},
        "score": round(composite, 4),
    }


def deduplicate(patches: list[dict], max_overlap: float) -> list[dict]:
    """Keep the best patch from each group describing the same surface."""
    selected: list[dict] = []
    for patch in sorted(patches, key=lambda p: -p["score"]):
        members = set(patch["residues"])
        redundant = False
        for kept in selected:
            other = set(kept["residues"])
            union = members | other
            if union and len(members & other) / len(union) > max_overlap:
                redundant = True
                break
        if not redundant:
            selected.append(patch)
    return selected


# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--scope", nargs=2, type=int, default=DEFAULT_SCOPE,
                        metavar=("START", "END"),
                        help=f"UniProt range to search "
                             f"(default {DEFAULT_SCOPE[0]} {DEFAULT_SCOPE[1]})")
    parser.add_argument("--radius", type=float, default=PATCH_RADIUS,
                        help=f"patch radius in A (default {PATCH_RADIUS})")
    parser.add_argument("--top", type=int, default=5,
                        help="patches to report in detail (default 5)")
    args = parser.parse_args()

    scope = tuple(args.scope)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_annotation(ANNOTATION)
    by_position = {r["uniprot_pos"]: r for r in rows}

    model = PDBParser(QUIET=True).get_structure("6aru", str(STRUCTURE))[0]
    atoms = centroid_atoms(model[EGFR_CHAIN], rows)

    in_scope = [r for r in rows if scope[0] <= r["uniprot_pos"] <= scope[1]]
    pool = [r["uniprot_pos"] for r in in_scope if r["designable"] == "yes"]
    exposed = [r["uniprot_pos"] for r in in_scope if r["exposed"] == "yes"]

    print(f"Scope {scope[0]}-{scope[1]} ({len(in_scope)} residues)")
    print(f"  exposed            : {len(exposed)}")
    print(f"  designable pool    : {len(pool)}")
    print(f"  with coordinates   : {sum(1 for p in pool if p in atoms)}")

    # Report the literature histidines explicitly: a position the literature
    # treats as mechanistic is only usable if it is actually on the surface.
    print("\nLiterature histidines in scope:")
    for position in (370, 433):
        row = by_position.get(position)
        if row is None or not (scope[0] <= position <= scope[1]):
            print(f"  H{position}: outside the current scope")
            continue
        print(f"  H{position}: rel SASA {row['rel_sasa']}, "
              f"{row['conservation']}, exposed={row['exposed']}, "
              f"glycan_shadow={row['glycan_shadow']}, "
              f"in pool={'yes' if position in pool else 'no'}")

    usable = [p for p in pool if p in atoms]
    if len(usable) < MIN_PATCH_RESIDUES:
        print(f"\nOnly {len(usable)} usable pool residues; "
              f"need at least {MIN_PATCH_RESIDUES}.")
        return 1

    print(f"\nBuilding patches, radius {args.radius} A")
    candidates = []
    for centre in usable:
        patch = build_patch(centre, usable, atoms, args.radius)
        if len(patch) < MIN_PATCH_RESIDUES:
            continue
        shell = shell_of(patch, exposed, atoms, SHELL_RADIUS)
        scored = score_patch(patch, shell, by_position, atoms)
        scored["centre"] = centre
        candidates.append(scored)

    print(f"  {len(candidates)} patches of at least "
          f"{MIN_PATCH_RESIDUES} residues")

    selected = deduplicate(candidates, MAX_OVERLAP)
    print(f"  {len(selected)} distinct surfaces after deduplication")

    # -- CSV of every distinct patch ---------------------------------------
    csv_path = OUT_DIR / "02_patches.csv"
    columns = ["rank", "centre", "n_residues", "area_a2", "score",
               "score_ph_potential", "score_conservation", "score_functional",
               "score_geometry", "histidines", "histidines_adjacent",
               "cetuximab_overlap", "conservation_identical_pct",
               "nearest_glycan_a", "shell_size", "residues"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns,
                                extrasaction="ignore")
        writer.writeheader()
        for index, patch in enumerate(selected, start=1):
            record = dict(patch)
            record["rank"] = index
            record["histidines"] = " ".join(map(str, patch["histidines"]))
            record["histidines_adjacent"] = " ".join(
                map(str, patch["histidines_adjacent"]))
            record["residues"] = " ".join(map(str, patch["residues"]))
            writer.writerow(record)

    # -- BindCraft hotspot specification -----------------------------------
    hotspots = []
    for index, patch in enumerate(selected[:args.top], start=1):
        pdb_numbers = [by_position[p]["pdb_resnum"] for p in patch["residues"]
                       if by_position[p]["pdb_resnum"] is not None]
        hotspots.append({
            "rank": index,
            "centre_uniprot": patch["centre"],
            "score": patch["score"],
            "n_residues": patch["n_residues"],
            "area_a2": patch["area_a2"],
            "histidines_uniprot": patch["histidines"],
            "residues_uniprot": patch["residues"],
            "residues_pdb": sorted(pdb_numbers),
            "bindcraft_target_hotspot_residues":
                ",".join(f"{EGFR_CHAIN}{n}" for n in sorted(pdb_numbers)),
        })

    hotspot_path = OUT_DIR / "02_hotspots_bindcraft.json"
    hotspot_path.write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "scripts/02_select_epitope.py",
        "scope_uniprot": list(scope),
        "structure": STRUCTURE.name,
        "chain": EGFR_CHAIN,
        "numbering_note": (
            "residues_uniprot is precursor numbering; residues_pdb and the "
            "BindCraft string are PDB numbering, 24 lower, because BindCraft "
            "indexes into the structure file it is given"),
        "parameters": {
            "patch_radius_a": args.radius,
            "shell_radius_a": SHELL_RADIUS,
            "area_window_a2": list(AREA_WINDOW),
            "min_patch_residues": MIN_PATCH_RESIDUES,
            "max_overlap": MAX_OVERLAP,
            "weights": WEIGHTS,
        },
        "patches": hotspots,
    }, indent=2) + "\n")

    # -- Readable report ----------------------------------------------------
    summary_path = OUT_DIR / "02_patch_summary.md"
    with summary_path.open("w") as handle:
        w = handle.write
        w("# Step 02 - Epitope patch selection\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/02_select_epitope.py`.\n\n")
        w(f"Scope: UniProt {scope[0]}-{scope[1]}. "
          f"{len(exposed)} exposed residues, {len(pool)} in the designable "
          f"pool, {len(selected)} distinct candidate surfaces.\n\n")

        w("A pool residue is one a binder is *permitted* to touch. A patch is "
          "a contiguous surface a binder could actually land on. The pool is "
          "not a hotspot list.\n\n")

        w("## Ranked patches\n\n")
        w("| Rank | Centre | Residues | Area Å² | His | Cetuximab overlap "
          "| Conserved % | Score |\n")
        w("|------|--------|----------|---------|-----|-------------------"
          "|-------------|-------|\n")
        for index, patch in enumerate(selected, start=1):
            his = ", ".join(f"H{h}" for h in patch["histidines"]) or "—"
            w(f"| {index} | {patch['centre']} | {patch['n_residues']} "
              f"| {patch['area_a2']:.0f} | {his} | {patch['cetuximab_overlap']} "
              f"| {patch['conservation_identical_pct']} | {patch['score']} |\n")

        w(f"\n## Top {args.top} in detail\n\n")
        for index, patch in enumerate(selected[:args.top], start=1):
            w(f"### Patch {index} — centred on residue {patch['centre']}\n\n")
            w(f"- **Score {patch['score']}** "
              f"(pH {patch['score_ph_potential']}, "
              f"conservation {patch['score_conservation']}, "
              f"functional {patch['score_functional']}, "
              f"geometry {patch['score_geometry']})\n")
            w(f"- {patch['n_residues']} pool residues, "
              f"{patch['area_a2']:.0f} Å² buriable\n")
            w(f"- Histidines in patch: "
              f"{', '.join(f'H{h}' for h in patch['histidines']) or 'none'}\n")
            w(f"- Histidines adjacent: "
              f"{', '.join(f'H{h}' for h in patch['histidines_adjacent']) or 'none'}\n")
            w(f"- Overlaps the cetuximab interface at "
              f"{patch['cetuximab_overlap']} residues\n")
            w(f"- Nearest glycan {patch['nearest_glycan_a']} Å\n")
            w(f"- Footprint {patch['conservation_identical_pct']}% identical "
              f"in mouse, counting the surrounding shell\n\n")
            w("Residues, UniProt precursor numbering:\n\n")
            w("```\n" + " ".join(
                f"{by_position[p]['aa_human']}{p}"
                for p in patch["residues"]) + "\n```\n\n")
            pdb_numbers = sorted(
                by_position[p]["pdb_resnum"] for p in patch["residues"]
                if by_position[p]["pdb_resnum"] is not None)
            w("BindCraft `target_hotspot_residues`, PDB numbering:\n\n")
            w("```\n" + ",".join(f"{EGFR_CHAIN}{n}" for n in pdb_numbers) + "\n```\n\n")

        w("## Scoring\n\n")
        w("Weights follow the competition's ranking order, in which pH "
          "selectivity outranks cross-reactivity, which outranks affinity.\n\n")
        w("| Term | Weight | Meaning |\n|------|--------|---------|\n")
        w(f"| pH potential | {WEIGHTS['ph_potential']} | exposure-weighted "
          "histidines in and beside the patch |\n")
        w(f"| Conservation | {WEIGHTS['conservation']} | fraction of the "
          "footprint identical in mouse, shell included |\n")
        w(f"| Functional | {WEIGHTS['functional']} | overlap with the "
          "cetuximab interface, counted as positive evidence of "
          "druggability |\n")
        w(f"| Geometry | {WEIGHTS['geometry']} | buriable area in the "
          f"{AREA_WINDOW[0]:.0f}-{AREA_WINDOW[1]:.0f} Å² window, and "
          "clearance from glycans |\n")
        w("\nThe weights are a judgement, not a measurement. They are written "
          "in one place in the source so they can be changed and the ranking "
          "re-derived.\n")

    # -- PyMOL session for visual checking ---------------------------------
    pml_path = OUT_DIR / "02_inspect_patches.pml"
    colours = ["orange", "marine", "forest", "purple", "yellow"]
    with pml_path.open("w") as handle:
        w = handle.write
        w(f"# Visual check of the top patches. Numbering here is PDB.\n")
        w(f"load {Path('..') / STRUCTURE.relative_to(ROOT)}, egfr\n")
        w("hide everything\nshow cartoon\ncolor grey80\n")
        w(f"color palecyan, chain {EGFR_CHAIN}\n")
        w(f"select fab, not chain {EGFR_CHAIN}\ncolor salmon, fab\n")
        for index, patch in enumerate(selected[:args.top], start=1):
            pdb_numbers = sorted(
                by_position[p]["pdb_resnum"] for p in patch["residues"]
                if by_position[p]["pdb_resnum"] is not None)
            selector = "+".join(map(str, pdb_numbers))
            w(f"\nselect patch{index}, chain {EGFR_CHAIN} and resi {selector}\n")
            w(f"show surface, patch{index}\n")
            w(f"color {colours[(index - 1) % len(colours)]}, patch{index}\n")
        w("\nset transparency, 0.2\nbg_color white\norient\n")
        w(f"\n# UniProt position = PDB position + 24\n")

    print(f"\n  -> {csv_path.relative_to(ROOT)}")
    print(f"  -> {summary_path.relative_to(ROOT)}")
    print(f"  -> {hotspot_path.relative_to(ROOT)}")
    print(f"  -> {pml_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

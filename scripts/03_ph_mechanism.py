#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 03 - pH-switch mechanism analysis.

Objective 1 of the challenge asks for binding at pH 6.5 and no detectable
binding at pH 7.4. No binder design pipeline optimises for this, so the
mechanism has to be specified before any sequence is generated. This script
does the mechanistic groundwork: it asks what protonation changes between
the two pH values could plausibly translate into a binding free energy
difference, and whether the geometry of each candidate surface admits them.

Nothing here designs a binder. It defines what a binder would have to do.

The quantitative constraint
---------------------------
Wyman linkage relates the pH dependence of binding to the number of protons
taken up on binding:

    d(dG_bind)/d(pH) = 2.303 * R * T * dn(H+)

Over the 0.9 pH units between 7.4 and 6.5 the maximum achievable difference
is therefore about 1.23 kcal/mol per fully coupled ionizable group. A single
histidine, even ideally placed, buys roughly one order of magnitude in
affinity. "Binding versus no detectable binding" needs considerably more,
which means two to three coupled protonation events rather than one.

This is the analysis's central output: it rules out the intuitive
one-histidine design before any compute is spent on it.

Two mechanistic directions
--------------------------
Both must gain affinity on protonation, since the bound state is the acidic
one.

  Receptor-side   A target histidine protonates at 6.5 and pairs with a
                  carboxylate on the binder. Limited to the histidines the
                  target happens to present, and their pKa is not ours to
                  tune.

  Binder-side     A histidine on the binder protonates and pairs with a
                  carboxylate on the target. Fully under our control in
                  number and placement, so this is where Dn can be raised,
                  but it requires acidic residues on the target surface and
                  is opposed by nearby basic residues, which repel the
                  protonated histidine.

The script therefore maps the ionizable composition of each candidate
surface, tests whether there is physical room for a partner residue at each
titratable group, and reports which direction each surface can support.

Evidence is kept in separate tiers throughout, so that the methods report
can state which claims rest on experiment, which on this structure, and
which are hypotheses awaiting computational or experimental test.

Inputs:  02_Analysis/01_residue_annotation.csv, 02_patches.csv
         01_Target/structure/6ARU.pdb
Outputs: 02_Analysis/03_ph_handles.csv
         02_Analysis/03_mechanism_report.md
         02_Analysis/03_design_families.json

Usage:
    python scripts/03_ph_mechanism.py
    python scripts/03_ph_mechanism.py --patches 1 2 5
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "02_Analysis"
ANNOTATION = OUT_DIR / "01_residue_annotation.csv"
PATCHES = OUT_DIR / "02_patches.csv"
STRUCTURE = ROOT / "01_Target/structure/6ARU.pdb"

EGFR_CHAIN = "A"

# Thermodynamics.
R_KCAL = 0.0019872          # kcal/(mol*K)
TEMPERATURE = 298.15        # K
PH_HIGH, PH_LOW = 7.4, 6.5

# Titratable atoms per residue type, with the charge the group carries when
# protonated relative to its basic form.
TITRATABLE = {
    "HIS": (["ND1", "NE2"], "acquires +1 on protonation", "base"),
    "ASP": (["OD1", "OD2"], "carries -1 when deprotonated", "acid"),
    "GLU": (["OE1", "OE2"], "carries -1 when deprotonated", "acid"),
    "LYS": (["NZ"], "carries +1 across this pH range", "fixed_positive"),
    "ARG": (["NH1", "NH2", "NE"], "carries +1 across this pH range", "fixed_positive"),
}

# Residues within this distance of the patch are part of its electrostatic
# environment even if they are not in the pool.
ENVIRONMENT_RADIUS = 12.0

# Probe geometry for asking whether a partner residue could reach a group.
PROBE_DISTANCES = (3.0, 3.5, 4.0)
CLASH_RADIUS = 2.8
PROBE_DIRECTIONS = 60


def load_annotation() -> dict[int, dict]:
    rows = {}
    with ANNOTATION.open() as handle:
        for row in csv.DictReader(handle):
            position = int(row["uniprot_pos"])
            row["uniprot_pos"] = position
            row["pdb_resnum"] = int(row["pdb_resnum"]) if row["pdb_resnum"] else None
            row["rel_sasa"] = float(row["rel_sasa"]) if row["rel_sasa"] else None
            rows[position] = row
    return rows


def load_patches(wanted: list[int] | None) -> list[dict]:
    patches = []
    with PATCHES.open() as handle:
        for row in csv.DictReader(handle):
            rank = int(row["rank"])
            if wanted and rank not in wanted:
                continue
            patches.append({
                "rank": rank,
                "centre": int(row["centre"]),
                "score": float(row["score"]),
                "residues": [int(p) for p in row["residues"].split()],
            })
    return patches


def fibonacci_directions(count: int) -> np.ndarray:
    """Roughly uniform unit vectors on a sphere."""
    indices = np.arange(count) + 0.5
    phi = np.arccos(1 - 2 * indices / count)
    theta = np.pi * (1 + 5 ** 0.5) * indices
    return np.column_stack([
        np.cos(theta) * np.sin(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(phi),
    ])


DIRECTIONS = fibonacci_directions(PROBE_DIRECTIONS)


def partner_room(atom_coord: np.ndarray, protein_coords: np.ndarray,
                 outward: np.ndarray) -> int:
    """How many clash-free positions a partner atom could occupy.

    Probes are placed at hydrogen-bond and salt-bridge distances in every
    direction with an outward component, and rejected if any protein heavy
    atom comes within a clash radius. A titratable group with no clash-free
    probe position is buried against the protein and cannot pair with
    anything on a binder, whatever its accessibility score says.
    """
    outward_directions = DIRECTIONS[DIRECTIONS @ outward > 0.2]
    if len(outward_directions) == 0:
        return 0

    free = 0
    for radius in PROBE_DISTANCES:
        probes = atom_coord + outward_directions * radius
        deltas = probes[:, None, :] - protein_coords[None, :, :]
        nearest = np.sqrt((deltas ** 2).sum(axis=2)).min(axis=1)
        free += int((nearest >= CLASH_RADIUS).sum())
    return free


def linkage_table() -> list[tuple[int, float, float]]:
    """Maximum affinity change achievable for n coupled protonations."""
    rt = R_KCAL * TEMPERATURE
    slope = 2.303 * rt * (PH_HIGH - PH_LOW)
    rows = []
    for n in range(1, 5):
        ddg = slope * n
        fold = math.exp(ddg / rt)
        rows.append((n, ddg, fold))
    return rows


def required_protons(fold_change: float) -> float:
    rt = R_KCAL * TEMPERATURE
    return math.log(fold_change) / (2.303 * (PH_HIGH - PH_LOW))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--patches", nargs="*", type=int, default=None,
                        help="patch ranks to analyse (default: all)")
    args = parser.parse_args()

    annotation = load_annotation()
    patches = load_patches(args.patches)

    model = PDBParser(QUIET=True).get_structure("6aru", str(STRUCTURE))[0]
    chain = model[EGFR_CHAIN]

    protein_atoms = [a for r in chain if r.id[0] == " "
                     for a in r if a.element != "H"]
    protein_coords = np.array([a.coord for a in protein_atoms])
    protein_centre = protein_coords.mean(axis=0)

    pdb_to_uniprot = {row["pdb_resnum"]: position
                      for position, row in annotation.items()
                      if row["pdb_resnum"] is not None}

    print("Thermodynamic ceiling on pH selectivity")
    print(f"  between pH {PH_HIGH} and {PH_LOW}, "
          f"{PH_HIGH - PH_LOW:.1f} pH units\n")
    print(f"  {'protons':<10}{'max ddG kcal/mol':<20}{'max fold change'}")
    for n, ddg, fold in linkage_table():
        print(f"  {n:<10}{ddg:<20.2f}{fold:.0f}x")
    print(f"\n  A 100-fold switch needs "
          f"{required_protons(100):.1f} coupled protonation events.")
    print(f"  A 1000-fold switch needs {required_protons(1000):.1f}.")
    print("  One histidine is not enough.\n")

    records = []
    family_input = []

    for patch in patches:
        centroid = np.array([
            chain[(" ", annotation[p]["pdb_resnum"], " ")]["CA"].coord
            for p in patch["residues"]
            if annotation[p]["pdb_resnum"] is not None
        ]).mean(axis=0)

        # Every ionizable residue whose titratable atom lies near the patch.
        found = []
        for residue in chain:
            if residue.id[0] != " ":
                continue
            resname = residue.get_resname()
            if resname not in TITRATABLE:
                continue
            atom_names, description, kind = TITRATABLE[resname]
            atoms = [residue[n] for n in atom_names if n in residue]
            if not atoms:
                continue
            coords = np.array([a.coord for a in atoms])
            separation = float(np.linalg.norm(coords - centroid, axis=1).min())
            if separation > ENVIRONMENT_RADIUS:
                continue

            position = pdb_to_uniprot.get(residue.id[1])
            if position is None:
                continue
            row = annotation[position]

            best_room = 0
            for atom in atoms:
                coord = np.array(atom.coord)
                outward = coord - protein_centre
                outward /= np.linalg.norm(outward)
                best_room = max(best_room,
                                partner_room(coord, protein_coords, outward))

            found.append({
                "patch": patch["rank"],
                "uniprot_pos": position,
                "pdb_resnum": residue.id[1],
                "residue": resname,
                "aa": row["aa_human"],
                "kind": kind,
                "behaviour": description,
                "in_patch": "yes" if position in patch["residues"] else "no",
                "dist_to_centroid_a": round(separation, 1),
                "rel_sasa": row["rel_sasa"],
                "conservation": row["conservation"],
                "glycan_shadow": row["glycan_shadow"],
                "partner_positions": best_room,
                "reachable": "yes" if best_room >= 10 else "no",
            })

        records.extend(found)

        # Charge balance of the environment, which decides which mechanistic
        # direction the surface can support.
        positives = [f for f in found if f["kind"] == "fixed_positive"
                     and f["reachable"] == "yes"]
        acids = [f for f in found if f["kind"] == "acid"
                 and f["reachable"] == "yes"]
        histidines = [f for f in found if f["kind"] == "base"
                      and f["reachable"] == "yes"]
        net = len(positives) - len(acids)

        family_input.append({
            "patch": patch["rank"],
            "centre": patch["centre"],
            "score": patch["score"],
            "histidines": sorted(h["uniprot_pos"] for h in histidines),
            "acids": sorted(a["uniprot_pos"] for a in acids),
            "basics": sorted(p["uniprot_pos"] for p in positives),
            "net_charge": net,
        })

        print(f"Patch {patch['rank']} (centre {patch['centre']})")
        print(f"  reachable histidines : "
              f"{[h['uniprot_pos'] for h in histidines] or 'none'}")
        print(f"  reachable acids      : "
              f"{[a['uniprot_pos'] for a in acids] or 'none'}")
        print(f"  reachable basics     : "
              f"{[p['uniprot_pos'] for p in positives] or 'none'}")
        print(f"  net charge in the environment : {net:+d}")

        receptor_side = len(histidines)
        binder_side = len(acids)
        print(f"  receptor-side handles (target His) : {receptor_side}")
        print(f"  binder-side handles (target acids) : {binder_side}")
        if net > 1:
            print("  NOTE: basic environment opposes placing a protonated "
                  "histidine on the binder here")
        print(f"  ceiling on Dn if every handle couples : "
              f"{receptor_side + binder_side}")
        print()

    # -- Outputs ------------------------------------------------------------
    csv_path = OUT_DIR / "03_ph_handles.csv"
    columns = ["patch", "uniprot_pos", "pdb_resnum", "residue", "aa", "kind",
               "in_patch", "dist_to_centroid_a", "rel_sasa", "conservation",
               "glycan_shadow", "partner_positions", "reachable", "behaviour"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    families_path = OUT_DIR / "03_design_families.json"
    families_path.write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "scripts/03_ph_mechanism.py",
        "thermodynamics": {
            "ph_high": PH_HIGH,
            "ph_low": PH_LOW,
            "kcal_per_coupled_proton": round(
                2.303 * R_KCAL * TEMPERATURE * (PH_HIGH - PH_LOW), 3),
            "protons_for_100x": round(required_protons(100), 2),
            "protons_for_1000x": round(required_protons(1000), 2),
        },
        "patches": family_input,
    }, indent=2) + "\n")

    report_path = OUT_DIR / "03_mechanism_report.md"
    with report_path.open("w") as handle:
        w = handle.write
        w("# Step 03 - pH-switch mechanism\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/03_ph_mechanism.py`.\n\n")
        w("Nothing here designs a binder. It establishes what a binder would "
          "have to do, and whether each candidate surface can support it.\n\n")

        w("## How much selectivity is physically available\n\n")
        w("Wyman linkage gives the pH dependence of binding as\n\n")
        w("```\nd(dG_bind)/d(pH) = 2.303 * R * T * dn(H+)\n```\n\n")
        w(f"Over the {PH_HIGH - PH_LOW:.1f} pH units between {PH_HIGH} and "
          f"{PH_LOW}:\n\n")
        w("| Coupled protonations | Max ΔΔG kcal/mol | Max affinity change |\n")
        w("|---------------------:|-----------------:|--------------------:|\n")
        for n, ddg, fold in linkage_table():
            w(f"| {n} | {ddg:.2f} | {fold:.0f}× |\n")
        w(f"\nA 100-fold switch requires {required_protons(100):.1f} coupled "
          f"protonation events; 1000-fold requires "
          f"{required_protons(1000):.1f}.\n\n")
        w("**A single histidine cannot deliver the requested behaviour.** "
          "Ideally placed and ideally titrating, one gives about one order of "
          "magnitude. The challenge asks for binding at 6.5 and none "
          "detectable at 7.4. Designs must therefore couple two to three "
          "protonation events, and this rules out the intuitive "
          "one-histidine design before any compute is spent on it.\n\n")
        w("These are ceilings. Real coupling is lower: pKa values shift at "
          "interfaces, protonation is not all-or-nothing across this narrow "
          "window, and a neutral histidine can still hydrogen bond to a "
          "carboxylate, which erodes the difference between the two states.\n\n")

        w("## Two mechanistic directions\n\n")
        w("The bound state is the acidic one, so both directions must *gain* "
          "affinity on protonation.\n\n")
        w("**Receptor-side.** A target histidine protonates at 6.5 and pairs "
          "with a carboxylate placed on the binder. Limited to the "
          "histidines the target presents, and their pKa is not ours to "
          "tune.\n\n")
        w("**Binder-side.** A histidine on the binder protonates and pairs "
          "with a carboxylate on the target. Under our control in number and "
          "placement, so this is the direction that can raise Δn, but it "
          "requires acidic residues on the target surface and is opposed by "
          "nearby basic residues, which repel the protonated histidine.\n\n")

        w("## Ionizable environment of each surface\n\n")
        w("A group counts as reachable when a partner atom could occupy a "
          "clash-free position at hydrogen-bond or salt-bridge distance from "
          "it. Relative accessibility alone is not sufficient: a group can be "
          "exposed and still have no room for a partner.\n\n")
        w("| Patch | Target His | Target acids | Target basics | Net charge "
          "| Ceiling on Δn |\n")
        w("|-------|-----------|--------------|---------------|-----------"
          "|---------------|\n")
        for entry in family_input:
            his = ", ".join(f"H{p}" for p in entry["histidines"]) or "—"
            acids = ", ".join(str(p) for p in entry["acids"]) or "—"
            basics = ", ".join(str(p) for p in entry["basics"]) or "—"
            ceiling = len(entry["histidines"]) + len(entry["acids"])
            w(f"| {entry['patch']} | {his} | {acids} | {basics} "
              f"| {entry['net_charge']:+d} | {ceiling} |\n")

        w("\nNet charge is the count of reachable fixed-positive groups minus "
          "reachable acids. A strongly positive environment argues against "
          "the binder-side mechanism there, since a protonated histidine "
          "would be electrostatically opposed.\n")

        w("\n## Evidence tiers\n\n")
        w("Kept separate so the methods report can state what each claim "
          "rests on.\n\n")
        w("| Tier | Content |\n|------|---------|\n")
        w("| Experimental evidence | H433, Q435 and K489 appear in published "
          "alanine scanning as positions that reduce binding. H370 and H433 "
          "are described in the literature in connection with pH dependence. |\n")
        w("| Structural observation | Accessibility, conservation, glycan "
          "occlusion, ionizable environment and partner room, all measured "
          "from 6ARU and UniProt in steps 01 to 03. H370 and L406 are not "
          "accessible in this structure. |\n")
        w("| Mechanistic hypothesis | That a given arrangement of "
          "protonatable groups at this interface will couple to binding "
          "strongly enough to produce the requested switch. Untested. |\n")
        w("| Computational criterion | Metrics used to rank designs at the "
          "generation stage. Not yet defined. |\n")

        w("\n## Design families\n\n")
        w("Families are organised by *how many* protonation events they "
          "attempt to couple, not by which histidine they centre on. The "
          "thermodynamics above make the count the governing variable.\n\n")

        for entry in family_input:
            ceiling = len(entry["histidines"]) + len(entry["acids"])
            w(f"- **Patch {entry['patch']}** (centre {entry['centre']}): "
              f"ceiling Δn = {ceiling}, net charge "
              f"{entry['net_charge']:+d}. ")
            if entry["histidines"] and entry["acids"]:
                w("Supports both directions; a combined design is the only "
                  "route to Δn ≥ 2 here.\n")
            elif entry["histidines"]:
                w("Receptor-side only.\n")
            elif entry["acids"]:
                w("Binder-side only.\n")
            else:
                w("No usable handle; not a pH-switch candidate.\n")

        w("\nH370 is deliberately excluded as a design target. It is buried "
          "and glycan-shadowed in this structure, so it cannot be a direct "
          "contact. That does not exclude a mechanistic role in binding or "
          "pH response, but such a role would be allosteric or "
          "conformational and is not something a first-generation design can "
          "exploit rationally.\n")

    print(f"  -> {csv_path.relative_to(ROOT)}")
    print(f"  -> {report_path.relative_to(ROOT)}")
    print(f"  -> {families_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

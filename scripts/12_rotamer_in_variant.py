#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 12 - Can the histidine reach the handle in the re-predicted variant?

Step 11 found that only 2 of 16 engineered histidines lie within 4 A of a
target carboxylate in the re-predicted complexes, that the rest are buried
without a counter-charge, and that their pKa falls accordingly. That leaves
one question unresolved, and the interpretation of the whole project turns
on it:

    Did AlphaFold2 place the side chain badly, or can the position not
    reach the handle at all?

If some rotamer reaches, the placement rule was sound and the predictor
packed the interface wrongly. If none reaches, the rule's distance
criterion was measuring something that does not survive the substitution,
and the failure belongs to the method.

How this differs from step 07
-----------------------------
Step 07 ran the same probe on the *parent* backbone, before the
substitution, asking whether a histidine introduced there could reach. 18
of 23 passed. This runs it on the *variant complex as predicted*, with the
histidine already present, its own side chain removed and everything else -
including the target - left in place.

Same probe, same constants, different structure. The comparison is the
point.

Usage
-----
    python scripts/12_rotamer_in_variant.py
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
PREDICTIONS = ROOT / "03_Design" / "validation" / "predictions"
VARIANTS = ROOT / "03_Design" / "ph_variants" / "05_variants.csv"
BOUND_STATE = ROOT / "03_Design" / "validation" / "11_bound_state_pka.csv"
OUT_DIR = ROOT / "03_Design" / "validation"

# Step 07 owns the probe. Loading it by path keeps one implementation, so
# this step cannot answer a slightly different question by accident.
_spec = importlib.util.spec_from_file_location(
    "step07", ROOT / "scripts" / "07_validate_ph_variants.py")
step07 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step07)

TOKEN = re.compile(r"^([A-Z])(\d+)(H?)$")

# Predicted complexes renumber the target 1-204; the trim kept UniProt
# 335-538, so target residue n is UniProt n + 334.
TARGET_OFFSET = 334

# The atoms that carry the charge on each kind of handle.
HANDLE_ATOMS = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "HIS": ("ND1", "NE2"),
}


def histidine_positions(substitutions: str) -> list[int]:
    positions = []
    for token in substitutions.split():
        match = TOKEN.match(token)
        if not match:
            continue
        origin, number, engineered = match.groups()
        if engineered == "H" or origin == "H":
            positions.append(int(number))
    return sorted(set(positions))


def atoms_by_residue(model) -> dict:
    """{(chain, residue id): [(atom name, coord), ...]} over heavy atoms."""
    table = {}
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            table[(chain.id, residue.id)] = [
                (atom.get_name(), np.array(atom.coord))
                for atom in residue if atom.element != "H"
            ]
    return table


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", type=Path, default=PREDICTIONS)
    parser.add_argument("--label", default="",
                        help="suffix for the output files")
    parser.add_argument("--bound-state", type=Path, default=BOUND_STATE,
                        help="step 11 output to set beside the probe")
    args = parser.parse_args()

    if not VARIANTS.is_file():
        parser.error(f"run step 05 first: {VARIANTS} is missing")

    # The handles each variant was designed to engage, in UniProt numbering.
    intended: dict[tuple[str, str], list[str]] = {}
    for row in csv.DictReader(VARIANTS.open()):
        intended[(row["design"], row["substitutions"])] = \
            row["engages"].split()

    # Step 11's measured contact distances, to set beside the probe.
    measured: dict[tuple[str, str, int], dict] = {}
    if args.bound_state.is_file():
        for row in csv.DictReader(args.bound_state.open()):
            key = (row["design"], row["variant"], int(row["position"]))
            measured[key] = row

    pdb_parser = PDBParser(QUIET=True)
    rows = []

    complexes = sorted(p for p in args.predictions.glob("pred_*.pdb")
                       if not p.stem.endswith("__parent"))
    print(f"{len(complexes)} variant complexes, probe from step 07 "
          f"({step07.PROBE_DIRECTIONS} orientations, radii "
          f"{step07.HIS_TIP_DISTANCES[0]} and "
          f"{step07.HIS_TIP_DISTANCES[1]} A)\n")

    for path in complexes:
        design, _, variant = path.stem[len("pred_"):].partition("__")
        substitutions = variant.replace("_", " ")
        positions = histidine_positions(substitutions)
        handles = intended.get((design, substitutions), [])
        if not positions or not handles:
            continue

        model = pdb_parser.get_structure(path.stem, path)[0]
        lengths = {c.id: sum(1 for r in c if r.id[0] == " ") for c in model}
        binder_id = min(lengths, key=lengths.get)
        target_id = max(lengths, key=lengths.get)
        table = atoms_by_residue(model)

        print(f"  {design.split('_')[2]:<5} {substitutions}")

        for position in positions:
            try:
                residue = model[binder_id][position]
            except KeyError:
                continue
            if residue.get_resname() != "HIS":
                continue

            # Everything except this histidine's own side chain. The target
            # stays, which is the whole point: the probe has to fit in the
            # interface as predicted.
            environment = step07.environment_without(residue, table)

            best = None
            for handle in handles:
                number = int(re.sub(r"[^0-9]", "", handle)) - TARGET_OFFSET
                try:
                    handle_residue = model[target_id][number]
                except KeyError:
                    continue
                names = HANDLE_ATOMS.get(handle_residue.get_resname(), ())
                for name in names:
                    if name not in handle_residue:
                        continue
                    coord = np.array(handle_residue[name].coord)
                    reaches, closest, count = step07.rotamer_reaches(
                        residue, coord, environment)
                    if best is None or closest < best["closest"]:
                        best = {"handle": handle, "reaches": reaches,
                                "closest": closest, "orientations": count}
                    if reaches:
                        break
                if best and best["reaches"]:
                    break

            if best is None:
                print(f"      His{position:<4} no handle resolved")
                continue

            key = (design, substitutions, position)
            as_predicted = measured.get(key, {}).get("distance_a", "")

            rows.append({
                "design": design,
                "variant": substitutions,
                "position": position,
                "intended_handle": best["handle"],
                "rotamer_reaches_in_variant": "yes" if best["reaches"] else "no",
                "closest_possible_a": (round(best["closest"], 2)
                                       if math.isfinite(best["closest"])
                                       else ""),
                "n_contact_orientations": best["orientations"],
                "as_predicted_distance_a": as_predicted,
            })

            verdict = "REACHES" if best["reaches"] else "blocked"
            shown = ("inf" if not math.isfinite(best["closest"])
                     else f"{best['closest']:.2f}")
            print(f"      His{position:<4} -> {best['handle']:<6} "
                  f"best possible {shown:>5} A   "
                  f"AF2 built {as_predicted or '?':>5} A   {verdict}")

    if not rows:
        print("\nNothing probed.")
        return 1

    out_csv = OUT_DIR / f"12_rotamer_in_variant{args.label}.csv"
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    reaches = [r for r in rows if r["rotamer_reaches_in_variant"] == "yes"]

    report = OUT_DIR / f"12_rotamer_in_variant_report{args.label}.md"
    with report.open("w") as fh:
        w = fh.write
        stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        w("# Step 12 - Rotamer feasibility inside the re-predicted variant\n\n")
        w(f"Generated {stamp} by `scripts/12_rotamer_in_variant.py`.\n\n")
        w("Step 11 showed the engineered histidines are buried without a "
          "counter-charge and that their pKa falls. This asks whether the "
          "contact was reachable at all in the predicted variant, or "
          "whether AlphaFold2 simply packed the side chain away from it.\n\n")
        w("The probe is step 07's, unchanged: "
          f"{step07.PROBE_DIRECTIONS} orientations on spheres of "
          f"{step07.HIS_TIP_DISTANCES[0]} and "
          f"{step07.HIS_TIP_DISTANCES[1]} A around CB, clashes discarded, "
          f"a contact counted between {step07.CONTACT_RANGE[0]} and "
          f"{step07.CONTACT_RANGE[1]} A. Only the structure differs: the "
          "variant complex as predicted, with the histidine's own side "
          "chain removed and the target left in place.\n\n")

        w("## Result\n\n")
        w(f"**{len(reaches)} of {len(rows)}** histidines can reach their "
          "intended handle by some rotamer inside the predicted complex.\n\n")

        w("| Design | Variant | His | Handle | Best possible | AF2 built | "
          "Reaches |\n")
        w("|---|---|---:|---|---:|---:|---|\n")
        for row in sorted(rows, key=lambda r: (r["closest_possible_a"] == "",
                                               r["closest_possible_a"])):
            w(f"| {row['design'].split('_')[2]} | {row['variant']} "
              f"| {row['position']} | {row['intended_handle']} "
              f"| {row['closest_possible_a'] or '-'} "
              f"| {row['as_predicted_distance_a'] or '-'} "
              f"| {row['rotamer_reaches_in_variant']} |\n")

        w("\n## Reading this\n\n")
        w("*Best possible* is the closest approach any unclashed rotamer "
          "could make. *AF2 built* is the distance in the predicted "
          "structure. Where the two differ widely the predictor had a "
          "better option and did not take it; where both are large, no "
          "rotamer was available.\n\n")

        w("## Comparison with step 07\n\n")
        w("Step 07 ran this probe on the parent backbone and 18 of 23 "
          "substitutions passed. Step 12 runs it on the variant complex "
          f"and {len(reaches)} of {len(rows)} pass.\n\n")
        w("The gap between those two numbers is the finding. A position "
          "that can accept a reaching histidine before the substitution is "
          "not the same as a position that still can once the sequence has "
          "changed and the complex has been predicted again.\n")

    print(f"\n  {len(reaches)} of {len(rows)} reach their intended handle")
    print(f"  -> {out_csv.relative_to(ROOT)}")
    print(f"  -> {report.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

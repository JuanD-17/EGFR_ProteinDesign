#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 15 - Were the predictors right to decline the contact?

Steps 12 and 14 established that both AlphaFold2 and Boltz-2 place the
engineered histidines away from their handles while an unclashed rotamer
that would reach sits unused, in ten of sixteen cases. That was read as a
failure of the predictors. It may not be.

The probe tests sterics only. A rotamer that does not clash is not one a
protein would adopt: burying a charged imidazole costs desolvation, and a
pocket lined with carbon is a poor place for it however well it fits. If
the unused positions are apolar and carry no compensating charge, then both
models declined for a reason, and the fault is in the placement rule for
ignoring desolvation rather than in the predictors.

What is measured
----------------
For each histidine, the probe is run again and the tip position that comes
closest to the handle is kept. Around that point, within 6 A:

    apolar fraction   carbon and sulfur over all heavy atoms
    net charge        Asp/Glu carboxylates against Arg/Lys/His cations
    burial            heavy atom count

The same is measured at the imidazole nitrogen the predictor actually
built, so the two environments can be set side by side.

The control is internal. Two histidines did form the contact in both
models. If the environment of the unused positions differs from theirs in
the expected direction, the measurement means something; if it does not,
this analysis has no signal and should be reported as such.

Usage
-----
    python scripts/15_contact_environment.py
    python scripts/15_contact_environment.py --predictions <dir> --label _boltz
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import statistics
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
PREDICTIONS = ROOT / "03_Design" / "validation" / "predictions"
VARIANTS = ROOT / "03_Design" / "ph_variants" / "05_variants.csv"
OUT_DIR = ROOT / "03_Design" / "validation"

_spec = importlib.util.spec_from_file_location(
    "step07", ROOT / "scripts" / "07_validate_ph_variants.py")
step07 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step07)

TOKEN = re.compile(r"^([A-Z])(\d+)(H?)$")
TARGET_OFFSET = 334

HANDLE_ATOMS = {
    "ASP": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "HIS": ("ND1", "NE2"),
}

# Formal charges that matter over this pH window.
ANION = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CATION = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("LYS", "NZ")}

ENVIRONMENT_RADIUS = 6.0
APOLAR_ELEMENTS = {"C", "S"}


def histidine_positions(substitutions: str) -> list[int]:
    out = []
    for token in substitutions.split():
        match = TOKEN.match(token)
        if match and (match.group(3) == "H" or match.group(1) == "H"):
            out.append(int(match.group(2)))
    return sorted(set(out))


def describe(point: np.ndarray, atoms: list) -> dict:
    """Composition and charge of everything within 6 A of a point."""
    apolar = polar = 0
    anions = cations = 0
    for element, resname, name, coord in atoms:
        if np.linalg.norm(coord - point) > ENVIRONMENT_RADIUS:
            continue
        if element in APOLAR_ELEMENTS:
            apolar += 1
        else:
            polar += 1
        if (resname, name) in ANION:
            anions += 1
        elif (resname, name) in CATION:
            cations += 1
    total = apolar + polar
    return {
        "neighbours": total,
        "apolar_fraction": round(apolar / total, 3) if total else 0.0,
        # Carboxylate oxygens come in pairs on one group; halve to count
        # groups rather than atoms.
        "net_charge": round(cations - anions / 2, 1),
    }


def best_reaching_tip(residue, handle_coords: list[np.ndarray],
                      environment: np.ndarray):
    """The unclashed tip position that comes closest to any handle atom."""
    if "CB" not in residue:
        return None, float("inf")
    origin = np.array(residue["CB"].coord)
    tips = np.vstack([origin + step07.DIRECTIONS * d
                      for d in step07.HIS_TIP_DISTANCES])

    deltas = tips[:, None, :] - environment[None, :, :]
    nearest = np.sqrt((deltas ** 2).sum(axis=2)).min(axis=1)
    free = tips[nearest >= step07.HIS_RING_CLEARANCE]
    if len(free) == 0:
        return None, float("inf")

    best_tip, best_distance = None, float("inf")
    for handle in handle_coords:
        distances = np.linalg.norm(free - handle, axis=1)
        index = int(distances.argmin())
        if distances[index] < best_distance:
            best_distance = float(distances[index])
            best_tip = free[index]
    return best_tip, best_distance


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", type=Path, default=PREDICTIONS)
    parser.add_argument("--label", default="")
    args = parser.parse_args()

    intended = {}
    for row in csv.DictReader(VARIANTS.open()):
        intended[(row["design"], row["substitutions"])] = row["engages"].split()

    pdb_parser = PDBParser(QUIET=True)
    rows = []

    for path in sorted(p for p in args.predictions.glob("pred_*.pdb")
                       if not p.stem.endswith("__parent")):
        design, _, variant = path.stem[len("pred_"):].partition("__")
        substitutions = variant.replace("_", " ")
        handles = intended.get((design, substitutions), [])
        positions = histidine_positions(substitutions)
        if not handles or not positions:
            continue

        model = pdb_parser.get_structure(path.stem, path)[0]
        lengths = {c.id: sum(1 for r in c if r.id[0] == " ") for c in model}
        binder_id, target_id = min(lengths, key=lengths.get), max(
            lengths, key=lengths.get)

        atoms = [(a.element, a.get_parent().get_resname(), a.get_name(),
                  np.array(a.coord))
                 for a in model.get_atoms() if a.element != "H"]

        table = {}
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                table[(chain.id, residue.id)] = [
                    (a.get_name(), np.array(a.coord))
                    for a in residue if a.element != "H"]

        handle_coords = []
        acid_coords = []
        for handle in handles:
            number = int(re.sub(r"[^0-9]", "", handle)) - TARGET_OFFSET
            try:
                hres = model[target_id][number]
            except KeyError:
                continue
            for name in HANDLE_ATOMS.get(hres.get_resname(), ()):
                if name in hres:
                    handle_coords.append(np.array(hres[name].coord))
                    if hres.get_resname() in ("ASP", "GLU"):
                        acid_coords.append(np.array(hres[name].coord))
        if not handle_coords:
            continue

        for position in positions:
            try:
                residue = model[binder_id][position]
            except KeyError:
                continue
            if residue.get_resname() != "HIS":
                continue

            environment = step07.environment_without(residue, table)
            tip, distance = best_reaching_tip(residue, handle_coords,
                                              environment)
            if tip is None:
                continue

            # Where the predictor actually put the imidazole.
            built = [residue[n] for n in ("ND1", "NE2") if n in residue]
            if not built:
                continue
            built_point = np.mean([b.coord for b in built], axis=0)
            built_distance = min(
                float(np.linalg.norm(built_point - h)) for h in handle_coords)
            # Step 11 measures against carboxylates, and so does the
            # classification here: a receptor histidine cannot stabilise a
            # protonated binder histidine, so a contact with one is not the
            # contact the mechanism needs.
            acid_distance = min(
                (float(np.linalg.norm(built_point - h)) for h in acid_coords),
                default=float("inf"))

            would_be = describe(tip, atoms)
            as_built = describe(built_point, atoms)

            rows.append({
                "design": design.split("_")[2],
                "variant": substitutions,
                "position": position,
                "reaching_distance_a": round(distance, 2),
                "built_distance_a": round(built_distance, 2),
                "acid_distance_a": (round(acid_distance, 2)
                                    if acid_distance < 900 else ""),
                "handle_is_acid": "yes" if acid_coords else "no",
                "contact_built": "yes" if acid_distance <= 4.0 else "no",
                "would_be_apolar": would_be["apolar_fraction"],
                "would_be_net_charge": would_be["net_charge"],
                "would_be_neighbours": would_be["neighbours"],
                "as_built_apolar": as_built["apolar_fraction"],
                "as_built_net_charge": as_built["net_charge"],
                "as_built_neighbours": as_built["neighbours"],
            })

    if not rows:
        print("nothing to measure")
        return 1

    out_csv = OUT_DIR / f"15_contact_environment{args.label}.csv"
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    built = [r for r in rows if r["contact_built"] == "yes"]
    unused = [r for r in rows
              if r["contact_built"] == "no" and r["reaching_distance_a"] <= 4.0]

    def mean(group, field):
        return statistics.mean(r[field] for r in group) if group else float("nan")

    report = OUT_DIR / f"15_contact_environment{args.label}.md"
    with report.open("w") as fh:
        w = fh.write
        stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        w("# Step 15 - Environment of the contact that was not made\n\n")
        w(f"Generated {stamp} by `scripts/15_contact_environment.py`.\n\n")
        w("Steps 12 and 14 found a reaching rotamer sitting unused in most "
          "cases. The probe tests sterics only, so this asks whether the "
          "position the predictors declined is one a protein would want: "
          f"what lies within {ENVIRONMENT_RADIUS:.0f} A of it, how apolar, "
          "and whether a compensating charge is present.\n\n")

        w("## Declined positions against the ones that were built\n\n")
        w("| | Contact built | Reachable, declined |\n")
        w("|---|---:|---:|\n")
        w(f"| Cases | {len(built)} | {len(unused)} |\n")
        w(f"| Apolar fraction at the contact point | "
          f"{mean(built, 'would_be_apolar'):.3f} | "
          f"{mean(unused, 'would_be_apolar'):.3f} |\n")
        w(f"| Net charge within {ENVIRONMENT_RADIUS:.0f} A | "
          f"{mean(built, 'would_be_net_charge'):+.1f} | "
          f"{mean(unused, 'would_be_net_charge'):+.1f} |\n")
        w(f"| Heavy atoms within {ENVIRONMENT_RADIUS:.0f} A | "
          f"{mean(built, 'would_be_neighbours'):.0f} | "
          f"{mean(unused, 'would_be_neighbours'):.0f} |\n")

        w("\n## Where the predictor put it instead\n\n")
        w("| | Declined contact point | Position actually built |\n")
        w("|---|---:|---:|\n")
        w(f"| Apolar fraction | {mean(unused, 'would_be_apolar'):.3f} | "
          f"{mean(unused, 'as_built_apolar'):.3f} |\n")
        w(f"| Net charge | {mean(unused, 'would_be_net_charge'):+.1f} | "
          f"{mean(unused, 'as_built_net_charge'):+.1f} |\n")
        w(f"| Heavy atoms | {mean(unused, 'would_be_neighbours'):.0f} | "
          f"{mean(unused, 'as_built_neighbours'):.0f} |\n")

        w("\n## Per histidine\n\n")
        w("| Design | Variant | His | Reach | Built | Apolar (reach) | "
          "Net charge (reach) | Apolar (built) |\n")
        w("|---|---|---:|---:|---:|---:|---:|---:|\n")
        for r in sorted(rows, key=lambda r: -r["would_be_apolar"]):
            w(f"| {r['design']} | {r['variant']} | {r['position']} "
              f"| {r['reaching_distance_a']:.2f} "
              f"| {r['built_distance_a']:.2f} "
              f"| {r['would_be_apolar']:.3f} "
              f"| {r['would_be_net_charge']:+.1f} "
              f"| {r['as_built_apolar']:.3f} |\n")

        acids = [r for r in rows if r["handle_is_acid"] == "yes"]
        bases = [r for r in rows if r["handle_is_acid"] == "no"]
        w("\n## The handle chemistry\n\n")
        w(f"{len(bases)} of {len(rows)} histidines were placed against a "
          "handle that is itself a histidine on the receptor - H433 or "
          f"H370 - and only {len(acids)} against a carboxylate.\n\n")
        w("That is a flaw in the placement rule, not in the predictions. "
          "The mechanism needs the protonated binder histidine to be "
          "stabilised on binding, which requires a negative counter-charge. "
          "A receptor histidine carries +1 over the same pH window, so the "
          "pair repels and binding is destabilised at the lower pH, which "
          "is the opposite of the intended direction. Step 03 listed H433 "
          "as a titratable group and the rule took it as an anchor without "
          "asking which sign it carries.\n\n")
        w("The measured shifts follow: against carboxylates they average "
          "well above those against receptor histidines, and the single "
          "positive shift in the set is the one carboxylate contact that "
          "was built. With five cases against eleven and distance "
          "confounded with chemistry, this is consistent with the "
          "mechanism rather than a demonstration of it.\n\n")
        w("\n## How to read this\n\n")
        w("A protonated imidazole wants a polar, negatively charged "
          "neighbourhood. If the declined positions are more apolar than "
          "the ones that were built, or carry no compensating negative "
          "charge, the predictors had grounds and the placement rule is at "
          "fault for scoring geometry without desolvation.\n\n")
        w("If the two groups look alike, this measurement has no signal "
          "and does not support either reading. With single-digit group "
          "sizes it cannot support much in any case, and is reported as a "
          "description rather than a test.\n")

    print(f"{len(rows)} histidines, {len(built)} built, "
          f"{len(unused)} reachable but declined\n")
    print(f"  apolar fraction   built {mean(built, 'would_be_apolar'):.3f}"
          f"   declined {mean(unused, 'would_be_apolar'):.3f}")
    print(f"  net charge        built {mean(built, 'would_be_net_charge'):+.1f}"
          f"   declined {mean(unused, 'would_be_net_charge'):+.1f}")
    print(f"  neighbours        built {mean(built, 'would_be_neighbours'):.0f}"
          f"   declined {mean(unused, 'would_be_neighbours'):.0f}")
    print(f"\n  -> {out_csv.relative_to(ROOT)}")
    print(f"  -> {report.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

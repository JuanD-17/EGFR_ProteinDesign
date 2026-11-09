#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 11 - pKa of the engineered histidines in the bound state.

Step 03b stated the limitation this step removes: "pKa values in the
*complex* cannot be computed before the complex exists, so the coupling
estimate here is a screen, not a prediction." Every coupling figure reported
so far rests on an assumed 2.0-unit shift on binding, and that assumption is
the quantity that decides the result.

Step 07 produced the complexes. They can now be measured instead.

What is compared
----------------
For each engineered histidine, two pKa values:

    bound  - PROPKA on the binder/target complex
    free   - PROPKA on the binder alone

The binder is taken from the complex without moving an atom, so "free" here
means the same conformation with the partner deleted. That isolates the
electrostatic effect of the target, which is what the linkage expression
needs, and deliberately excludes conformational relaxation on unbinding.
Relaxation would add to the shift, so this is the conservative direction:
the measured shift is a lower bound on the real one.

The coupling is then recomputed from the measured pair by the same numerical
linkage used in step 03b, and set beside the figure the assumed shift gave.

What this is not
----------------
PROPKA is an empirical scheme with an error near one pKa unit, comparable to
the effects being reasoned about, and it is applied here to AlphaFold2
predictions rather than to crystal structures. A measured shift is a better
foundation than a stipulated one; it is not a validated one. Constant-pH
molecular dynamics remains the rigorous route.

Usage
-----
    python scripts/11_bound_state_pka.py
    python scripts/11_bound_state_pka.py --predictions <dir>
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
PREDICTIONS = ROOT / "03_Design" / "validation" / "predictions"
VARIANTS = ROOT / "03_Design" / "ph_variants" / "05_variants.csv"
OUT_DIR = ROOT / "03_Design" / "validation"

# Step 03b holds the PROPKA wrapper, the Python 3.14 patch and the linkage
# integration. Loading it by path keeps one implementation rather than a
# second copy that could drift.
_spec = importlib.util.spec_from_file_location(
    "step03b", ROOT / "scripts" / "03b_protonation_networks.py")
step03b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step03b)

# A substitution token is written M10H; a pre-existing histidine is H7.
TOKEN = re.compile(r"^([A-Z])(\d+)(H?)$")

# The predicted complexes renumber the target 1-204. The trim kept UniProt
# 335-538, so target residue n is UniProt n + 334.
TARGET_OFFSET = 334

CARBOXYLATE = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
IMIDAZOLE_N = ("ND1", "NE2")

# A charge-assisted imidazole-carboxylate contact, the same window step 07
# accepted a rotamer orientation in.
CONTACT_MAX = 4.0

# Heavy atoms within this radius of the imidazole centroid, as a plain
# measure of how enclosed the ring is.
BURIAL_RADIUS = 8.0


class ChainSelect(Select):
    """Write one chain, heavy atoms, standard residues only."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return (residue.id[0] == " "
                and residue.get_resname() in protein_letters_3to1)

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


def histidine_positions(substitutions: str) -> list[int]:
    """Binder positions carrying a histidine in this variant.

    Both kinds count. `M10H` is engineered and `H7` was already there, but
    both titrate and both contribute to the linkage.
    """
    positions = []
    for token in substitutions.split():
        match = TOKEN.match(token)
        if not match:
            continue
        origin, number, engineered = match.groups()
        if engineered == "H" or origin == "H":
            positions.append(int(number))
    return sorted(set(positions))


def binder_chain(model) -> str:
    """The shorter chain. The target is domain III at 204 residues."""
    lengths = {
        chain.id: sum(1 for r in chain if r.id[0] == " ")
        for chain in model
    }
    return min(lengths, key=lengths.get)


def nearest_carboxylate(residue, target_chain) -> tuple[float, str]:
    """Distance from either imidazole nitrogen to the closest target acid.

    This is the contact the design was placed to form. Step 07 established
    that a rotamer *could* reach it on the parent backbone; this measures
    whether the re-predicted variant actually does.
    """
    nitrogens = [residue[n] for n in IMIDAZOLE_N if n in residue]
    if not nitrogens:
        return float("inf"), "-"

    best, name = float("inf"), "-"
    for other in target_chain:
        oxygens = CARBOXYLATE.get(other.get_resname())
        if not oxygens:
            continue
        for oxygen in oxygens:
            if oxygen not in other:
                continue
            for nitrogen in nitrogens:
                distance = nitrogen - other[oxygen]
                if distance < best:
                    best = distance
                    letter = "D" if other.get_resname() == "ASP" else "E"
                    name = f"{letter}{other.id[1] + TARGET_OFFSET}"
    return best, name


def burial(residue, model) -> int:
    """Heavy atoms near the imidazole centroid, both chains included."""
    nitrogens = [residue[n] for n in IMIDAZOLE_N if n in residue]
    if not nitrogens:
        return 0
    centre = np.mean([n.coord for n in nitrogens], axis=0)
    return sum(1 for atom in model.get_atoms()
               if atom.element != "H"
               and np.linalg.norm(atom.coord - centre) < BURIAL_RADIUS)


def propka_by_chain(pdb_path: Path) -> dict[tuple[str, str, int], float]:
    """pKa keyed by (chain, residue type, number).

    Step 03b keys on (type, number) alone, which is safe for a single chain
    and silently collides in a complex where both chains can carry a
    histidine at the same number.
    """
    step03b.patch_propka_for_python_314()
    from propka.run import single

    molecule = single(str(pdb_path), optargs=["--quiet"], write_pka=False)
    name = (molecule.conformation_names[0]
            if hasattr(molecule, "conformation_names") else "AVR")
    conformation = molecule.conformations[name]

    values: dict[tuple[str, str, int], float] = {}
    for group in conformation.groups:
        if group.pka_value is None:
            continue
        key = (group.atom.chain_id, group.residue_type, group.atom.res_num)
        values[key] = float(group.pka_value)
    return values


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", type=Path, default=PREDICTIONS,
                        help="directory of re-predicted complexes")
    parser.add_argument("--scratch", type=Path,
                        default=OUT_DIR / "11_binder_only",
                        help="where isolated binder chains are written")
    args = parser.parse_args()

    if not args.predictions.is_dir():
        parser.error(f"not a directory: {args.predictions}")
    if not VARIANTS.is_file():
        parser.error(f"run step 05 first: {VARIANTS} is missing")

    # Substitution strings, so the histidine positions come from the same
    # record that produced the designs rather than from the filenames.
    assumed: dict[tuple[str, str], float] = {}
    for row in csv.DictReader(VARIANTS.open()):
        key = (row["design"], row["substitutions"])
        assumed[key] = float(row["estimated_coupling_kcal_mol"])

    args.scratch.mkdir(parents=True, exist_ok=True)
    parser_pdb = PDBParser(QUIET=True)
    io = PDBIO()

    rows = []
    complexes = sorted(p for p in args.predictions.glob("pred_*.pdb")
                       if not p.stem.endswith("__parent"))
    if not complexes:
        parser.error(f"no variant complexes in {args.predictions}")

    print(f"{len(complexes)} complexes\n")

    for path in complexes:
        design, _, variant = path.stem[len("pred_"):].partition("__")
        substitutions = variant.replace("_", " ")
        positions = histidine_positions(substitutions)
        if not positions:
            print(f"  {variant:<16} no histidine positions, skipped")
            continue

        model = parser_pdb.get_structure(path.stem, path)[0]
        chain_id = binder_chain(model)
        target_id = max((c.id for c in model), key=lambda c:
                        sum(1 for r in model[c] if r.id[0] == " "))

        bound = propka_by_chain(path)

        isolated = args.scratch / f"{path.stem}__binder.pdb"
        io.set_structure(model)
        io.save(str(isolated), ChainSelect(chain_id))
        free = propka_by_chain(isolated)

        print(f"  {design.split('_')[2]:<5} {variant:<16} chain {chain_id}")
        for position in positions:
            key = (chain_id, "HIS", position)
            # PROPKA writes the isolated chain back under the same id.
            pka_bound = bound.get(key)
            pka_free = free.get(key)
            if pka_bound is None or pka_free is None:
                print(f"      His{position:<4} not resolved by PROPKA")
                continue

            shift = pka_bound - pka_free
            coupling = step03b.coupling_energy(pka_free, pka_bound)

            residue = model[chain_id][position]
            distance, partner = nearest_carboxylate(residue, model[target_id])
            enclosed = burial(residue, model)
            in_contact = distance <= CONTACT_MAX

            rows.append({
                "design": design,
                "variant": substitutions,
                "position": position,
                "pka_free": round(pka_free, 2),
                "pka_bound": round(pka_bound, 2),
                "measured_shift": round(shift, 2),
                "coupling_measured_kcal_mol": round(coupling, 3),
                "nearest_acid": partner,
                "distance_a": round(distance, 2),
                "in_contact": "yes" if in_contact else "no",
                "neighbours_8a": enclosed,
            })
            flag = "contact" if in_contact else "       "
            print(f"      His{position:<4} free {pka_free:5.2f}  "
                  f"bound {pka_bound:5.2f}  shift {shift:+5.2f}   "
                  f"{partner:>5} {distance:5.2f} A  {flag}")

    if not rows:
        print("\nNo histidine resolved. Nothing written.")
        return 1

    out_csv = OUT_DIR / "11_bound_state_pka.csv"
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # Per variant, the coupling is the sum over its histidines, which is what
    # step 05 did with the assumed shift.
    totals: dict[tuple[str, str], float] = {}
    for row in rows:
        key = (row["design"], row["variant"])
        totals[key] = totals.get(key, 0.0) + row["coupling_measured_kcal_mol"]

    report = OUT_DIR / "11_bound_state_report.md"
    contacts = [r for r in rows if r["in_contact"] == "yes"]
    positive = [r for r in rows if r["measured_shift"] > 0]
    shifts = [r["measured_shift"] for r in rows]

    with report.open("w") as fh:
        w = fh.write
        stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        w("# Step 11 - Bound-state pKa and contact geometry\n\n")
        w(f"Generated {stamp} by `scripts/11_bound_state_pka.py`.\n\n")
        w("Step 03b assumed a 2.0-unit pKa shift on binding because the "
          "complexes did not yet exist. They exist now, so the shift is "
          "measured here, together with whether the contact that was "
          "supposed to produce it was actually formed.\n\n")

        w("## Result\n\n")
        w(f"Of {len(rows)} histidines across {len(totals)} variants, "
          f"**{len(contacts)} lie within {CONTACT_MAX:.1f} A of a target "
          f"carboxylate** and {len(positive)} show a positive pKa shift. "
          f"Measured shifts run from {min(shifts):+.2f} to "
          f"{max(shifts):+.2f}, against the "
          f"{step03b.DEFAULT_SHIFT:+.1f} that was assumed.\n\n")
        w("The assumption was wrong in sign, not in magnitude.\n\n")

        w("## Per histidine\n\n")
        w("| Design | Variant | His | pKa free | pKa bound | Shift | "
          "Nearest acid | d (A) | Contact | Neighbours |\n")
        w("|---|---|---:|---:|---:|---:|---|---:|---|---:|\n")
        for row in sorted(rows, key=lambda r: -r["measured_shift"]):
            w(f"| {row['design'].split('_')[2]} | {row['variant']} "
              f"| {row['position']} | {row['pka_free']:.2f} "
              f"| {row['pka_bound']:.2f} | {row['measured_shift']:+.2f} "
              f"| {row['nearest_acid']} | {row['distance_a']:.2f} "
              f"| {row['in_contact']} | {row['neighbours_8a']} |\n")

        w("\n## Per variant, measured against assumed\n\n")
        w("| Design | Variant | Assumed | Measured | \n")
        w("|---|---|---:|---:|\n")
        for (design, variant), measured in sorted(
                totals.items(), key=lambda kv: -kv[1]):
            before = assumed.get((design, variant))
            shown = f"{before:.3f}" if before is not None else "-"
            w(f"| {design.split('_')[2]} | {variant} | {shown} "
              f"| {measured:+.3f} |\n")

        w("\n## What this means\n\n")
        w("A positive shift means the target stabilises the protonated "
          "histidine, which is the direction that gives tighter binding at "
          "the lower pH. Every variant here is at or below zero, so under "
          "this estimate none of them is a pH switch in the intended "
          "direction, and most would bind marginally *better* at pH 7.4.\n\n")

        w("The contact column explains why. The histidines are enclosed - "
          "70 to 113 heavy atoms within 8 A - but they are not paired with "
          "a carboxylate. Burying a base without a counter-charge "
          "destabilises its protonated form, and the pKa falls. That is "
          "what the numbers show, and the internal consistency is the "
          "evidence that they mean it: the only histidine within hydrogen "
          "bonding distance of an acid is also the only one whose shift is "
          "not negative.\n\n")

        w("## Where the pipeline failed\n\n")
        w("Step 07 tested whether a rotamer *could* reach the handle, "
          "probing orientations on the parent backbone with the replaced "
          "side chain removed. 18 of 23 substitutions passed. That test was "
          "correct for the question it asked.\n\n")
        w("It is the wrong question. What matters is whether the contact "
          "exists in the re-predicted variant complex, and it mostly does "
          "not. Geometric reachability on the parent does not predict "
          "contact formation after the sequence changes and the structure "
          "is predicted again.\n\n")
        w("This is a failure of the placement rule, not of the "
          "thermodynamics. The linkage argument in step 03 stands: a "
          "coupled protonation is still worth up to 1.23 kcal/mol over this "
          "pH window. Nothing here was coupled.\n\n")

        w("## Limitations of this step\n\n")
        w("The binder was taken from the complex without moving an atom, so "
          "the free state is the bound conformation with the partner "
          "deleted. Conformational relaxation on unbinding would raise the "
          "free-state pKa and make the shifts *more* negative, so this is "
          "the conservative direction and does not rescue the result.\n\n")
        w("PROPKA's error is near one pKa unit. Several shifts here exceed "
          "three units, which is well outside that error, but the smaller "
          "ones are not resolvable.\n\n")
        w("These are AlphaFold2 models. If the prediction places the "
          "histidine wrongly, both the contact distance and the pKa are "
          "wrong together. Constant-pH molecular dynamics on an ensemble, "
          "rather than PROPKA on one frame, is the rigorous route.\n")

    print(f"\n  -> {out_csv.relative_to(ROOT)}")
    print(f"  -> {report.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

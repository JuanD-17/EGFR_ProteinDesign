#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 03b - pKa estimation and protonation network geometry.

Step 03 established a thermodynamic ceiling of about 1.23 kcal/mol per
coupled protonation, and concluded that two to three coupled events are
needed. That ceiling carries a condition which this step makes explicit and
then tests.

What actually produces pH selectivity
-------------------------------------
The linkage coefficient is the number of protons taken up on binding:

    dn(pH) = f_bound(pH) - f_free(pH)

If a group's pKa does not change on binding, f_bound equals f_free, dn is
zero, and there is no pH dependence at all, however many histidines the
interface contains. Selectivity comes from the pKa *shift* on binding, not
from the pKa itself, and the design requirement is therefore to place each
histidine so that binding stabilises its protonated form. Burying it against
a carboxylate is the mechanism that does this.

Two consequences follow, and both are testable now.

First, a group whose free pKa already lies far outside the 6.5-7.4 window
in the direction of the shift is useless: a histidine already protonated at
pH 7.4 cannot switch. Free-state pKa values are estimated here with PROPKA.

Second, the groups have to be simultaneously engageable. Two titratable
groups 30 A apart cannot both be contacted by one small binder, and two
groups 3 A apart will interact with each other rather than with the binder.
This step enumerates which subsets satisfy distance and orientation
constraints together, which is where a theoretical dn ceiling either becomes
a design mechanism or does not.

Scope and limitation
--------------------
pKa values in the *complex* cannot be computed before the complex exists, so
the coupling estimate here is a screen, not a prediction. PROPKA is an
empirical method with errors around one pKa unit, which is comparable to the
effects being reasoned about. Rigorous coupling free energies would need
constant-pH molecular dynamics, which is out of scope. What this step
delivers is a ranked set of geometrically viable protonation networks whose
free-state pKa values do not already rule them out.

Inputs:  01_Target/structure/6ARU.pdb, 02_Analysis/03_ph_handles.csv
Outputs: 02_Analysis/03b_pka_estimates.csv
         02_Analysis/03b_networks.csv
         02_Analysis/03b_network_report.md

Usage:
    python scripts/03b_protonation_networks.py
    python scripts/03b_protonation_networks.py --assumed-shift 1.5
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import tempfile
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "02_Analysis"
STRUCTURE = ROOT / "01_Target/structure/6ARU.pdb"
HANDLES = OUT_DIR / "03_ph_handles.csv"

EGFR_CHAIN = "A"
UNIPROT_MINUS_PDB = 24

R_KCAL = 0.0019872
TEMPERATURE = 298.15
PH_HIGH, PH_LOW = 7.4, 6.5

# A binder face of the size we are designing can span this much. Groups
# further apart than the upper bound cannot both be engaged; groups closer
# than the lower bound will titrate against each other rather than against
# the binder.
NETWORK_SPAN = (5.0, 22.0)

# Outward normals must point similarly enough for one face to reach all.
MAX_NORMAL_ANGLE = 80.0

TITRATABLE_ATOMS = {
    "HIS": ["ND1", "NE2"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}


def patch_propka_for_python_314() -> None:
    """Make PROPKA 3.5 run on Python 3.14.

    PROPKA dispatches its parameter-file parser on the annotation of each
    field, read as `self.__annotations__`. Python 3.14 evaluates annotations
    lazily (PEP 649) and no longer resolves `__annotations__` through an
    instance, so that lookup raises AttributeError and PROPKA cannot read its
    own parameter file. The annotations are still present on the class.

    This replaces the method with the identical logic reading the class
    annotations instead. It changes no behaviour on older interpreters and
    can be deleted once PROPKA supports 3.14.
    """
    import propka.parameters as parameters

    annotations = dict(parameters.Parameters.__annotations__)

    def parse_line(self, line):
        comment = line.find("#")
        if comment != -1:
            line = line[:comment]
        words = line.split()
        if not words:
            return
        annotation = annotations.get(words[0])
        if annotation is parameters._T_NUMBER_DICTIONARY:
            self.parse_to_number_dictionary(words)
        elif annotation is parameters._T_STRING_LIST:
            self.parse_to_string_list(words)
        elif annotation is parameters._T_STRING:
            self.parse_string(words)
        elif annotation is parameters._T_LIST_DICTIONARY:
            self.parse_to_list_dictionary(words)
        elif annotation in (parameters._T_MATRIX,
                            parameters._T_PAIR_WISE_MATRIX):
            self.parse_to_matrix(words)
        elif annotation is parameters._T_STRING_DICTIONARY:
            self.parse_to_string_dictionary(words)
        else:
            self.parse_parameter(words)

    parameters.Parameters.parse_line = parse_line

# Shift assumed to be achievable by burying a group against a complementary
# charge. Deliberately modest; --assumed-shift overrides it.
DEFAULT_SHIFT = 2.0

# Model pKa of a solvent-exposed histidine, used for the histidine an
# engineered binder would carry. Its shift on binding is what produces the
# binder-side switch.
BINDER_HIS_PKA = 6.5


class ChainOnly(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return (residue.id[0] == " "
                and residue.get_resname() in protein_letters_3to1)

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


# ---------------------------------------------------------------------------
# pKa
# ---------------------------------------------------------------------------

SUMMARY_LINE = re.compile(
    r"^\s*([A-Z]{2,3})\s+(\d+)\s+(\w)\s+([-\d.]+)\s+([-\d.]+)")


def run_propka(pdb_path: Path) -> dict[tuple[str, int], dict]:
    """Estimate pKa values for every ionizable group in a structure."""
    patch_propka_for_python_314()
    from propka.run import single

    molecule = single(str(pdb_path), optargs=["--quiet"], write_pka=False)
    conformation = molecule.conformations[
        molecule.conformation_names[0]
        if hasattr(molecule, "conformation_names") else "AVR"
    ]

    values: dict[tuple[str, int], dict] = {}
    for group in conformation.groups:
        if group.pka_value is None:
            continue
        key = (group.residue_type, group.atom.res_num)
        values[key] = {
            "pka": float(group.pka_value),
            "model_pka": float(group.model_pka),
            "chain": group.atom.chain_id,
        }
    return values


def fraction_protonated(pka: float, ph: float) -> float:
    return 1.0 / (1.0 + 10 ** (ph - pka))


def coupling_energy(pka_free: float, pka_bound: float,
                    steps: int = 200) -> float:
    """kcal/mol of pH selectivity from one group, by numerical linkage.

    Integrates dn over the pH window. Positive means tighter binding at the
    lower pH, which is the direction the challenge asks for.
    """
    grid = np.linspace(PH_LOW, PH_HIGH, steps)
    delta_n = np.array([
        fraction_protonated(pka_bound, ph) - fraction_protonated(pka_free, ph)
        for ph in grid
    ])
    integral = float(np.trapezoid(delta_n, grid))
    return 2.303 * R_KCAL * TEMPERATURE * integral


def fold_change(energy: float) -> float:
    return math.exp(energy / (R_KCAL * TEMPERATURE))


# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--assumed-shift", type=float, default=DEFAULT_SHIFT,
                        help=f"pKa shift assumed achievable on binding "
                             f"(default {DEFAULT_SHIFT})")
    args = parser.parse_args()

    # -- Which groups step 03 found reachable ------------------------------
    handles: dict[int, dict] = {}
    with HANDLES.open() as handle:
        for row in csv.DictReader(handle):
            if row["reachable"] != "yes":
                continue
            if row["residue"] not in TITRATABLE_ATOMS:
                continue
            position = int(row["uniprot_pos"])
            entry = handles.setdefault(position, {
                "uniprot_pos": position,
                "pdb_resnum": int(row["pdb_resnum"]),
                "residue": row["residue"],
                "aa": row["aa"],
                "conservation": row["conservation"],
                "glycan_shadow": row["glycan_shadow"],
                "patches": set(),
            })
            entry["patches"].add(int(row["patch"]))

    print(f"{len(handles)} reachable titratable groups from step 03\n")

    # -- pKa ----------------------------------------------------------------
    model = PDBParser(QUIET=True).get_structure("6aru", str(STRUCTURE))[0]
    chain = model[EGFR_CHAIN]

    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        free_pdb = workdir / "egfr_free.pdb"
        io = PDBIO()
        io.set_structure(model)
        io.save(str(free_pdb), select=ChainOnly(EGFR_CHAIN))

        print("Running PROPKA on the free receptor")
        pka_values = run_propka(free_pdb)
    print(f"  {len(pka_values)} ionizable groups predicted\n")

    print("pKa of the reachable groups in the free receptor")
    print(f"  {'group':<8}{'role':<36}{'pKa free':<10}"
          f"{'kcal/mol':<10}{'fold':<8}verdict")
    print("  " + "-" * 96)

    pka_rows = []
    for position in sorted(handles):
        entry = handles[position]
        key = (entry["residue"], entry["pdb_resnum"])
        record = pka_values.get(key)
        if record is None:
            continue

        pka_free = record["pka"]

        # The two mechanisms titrate different molecules, so they are scored
        # differently. Conflating them inverts the reading of an acid: a
        # carboxylate with a very low pKa is useless as a titrating group but
        # is the best possible anchor, because it is reliably ionised at both
        # pH values and therefore upshifts an engineered histidine's pKa
        # without itself changing state.
        if entry["residue"] == "HIS":
            role = "receptor-side switch"
            titrating_pka = pka_free
            pka_bound = pka_free + args.assumed_shift
            energy = coupling_energy(titrating_pka, pka_bound)

            if pka_free > PH_HIGH + 1.5:
                verdict = "protonated at both pH, cannot switch"
            elif pka_free < PH_LOW - 2.5:
                verdict = "too acidic for the assumed shift"
            elif energy >= 0.35:
                verdict = "usable"
            else:
                verdict = "weak"
        else:
            role = "anchor for an engineered histidine"
            # The binder's histidine titrates, starting from its model pKa
            # in solvent and shifted up by burial against this anchor.
            titrating_pka = BINDER_HIS_PKA
            pka_bound = BINDER_HIS_PKA + args.assumed_shift
            energy = coupling_energy(titrating_pka, pka_bound)

            # Anchor quality is how reliably the carboxylate stays ionised
            # across the window. Lower pKa is better.
            if pka_free <= PH_LOW - 2.0:
                verdict = "excellent anchor, fully ionised"
            elif pka_free <= PH_LOW - 0.5:
                verdict = "good anchor"
            else:
                verdict = "weak anchor, partly neutral in range"

        fold = fold_change(energy)

        pka_rows.append({
            "uniprot_pos": position,
            "pdb_resnum": entry["pdb_resnum"],
            "residue": entry["residue"],
            "aa": entry["aa"],
            "role": role,
            "pka_free": round(pka_free, 2),
            "model_pka": round(record["model_pka"], 2),
            "titrating_pka": round(titrating_pka, 2),
            "assumed_pka_bound": round(pka_bound, 2),
            "coupling_kcal_mol": round(energy, 3),
            "fold_change": round(fold, 2),
            "verdict": verdict,
            "conservation": entry["conservation"],
            "glycan_shadow": entry["glycan_shadow"],
            "patches": " ".join(map(str, sorted(entry["patches"]))),
        })

        print(f"  {entry['aa']}{position:<7}{role:<36}{pka_free:<10.2f}"
              f"{energy:<10.3f}{fold:<8.2f}{verdict}")

    usable = [r for r in pka_rows
              if r["verdict"] in ("usable", "excellent anchor, fully ionised",
                                  "good anchor")]
    switches = [r for r in usable if r["residue"] == "HIS"]
    anchors = [r for r in usable if r["residue"] != "HIS"]
    print(f"\n  {len(switches)} receptor-side switch(es), "
          f"{len(anchors)} usable anchor(s), "
          f"under a {args.assumed_shift} unit shift")

    # -- Network geometry ---------------------------------------------------
    print("\nEnumerating protonation networks")

    protein_coords = np.array([
        a.coord for r in chain if r.id[0] == " " for a in r if a.element != "H"
    ])
    protein_centre = protein_coords.mean(axis=0)

    geometry: dict[int, dict] = {}
    for row in usable:
        try:
            residue = chain[(" ", row["pdb_resnum"], " ")]
        except KeyError:
            continue
        names = TITRATABLE_ATOMS[row["residue"]]
        atoms = [residue[n] for n in names if n in residue]
        if not atoms:
            continue
        centre = np.mean([a.coord for a in atoms], axis=0)
        outward = centre - protein_centre
        outward /= np.linalg.norm(outward)
        geometry[row["uniprot_pos"]] = {
            "coord": centre, "normal": outward, "row": row,
        }

    networks = []
    positions = sorted(geometry)
    for size in (2, 3):
        for combo in combinations(positions, size):
            coords = [geometry[p]["coord"] for p in combo]
            normals = [geometry[p]["normal"] for p in combo]

            distances = [float(np.linalg.norm(a - b))
                         for a, b in combinations(coords, 2)]
            angles = [
                math.degrees(math.acos(float(np.clip(np.dot(a, b), -1, 1))))
                for a, b in combinations(normals, 2)
            ]

            if min(distances) < NETWORK_SPAN[0]:
                continue
            if max(distances) > NETWORK_SPAN[1]:
                continue
            if max(angles) > MAX_NORMAL_ANGLE:
                continue

            rows = [geometry[p]["row"] for p in combo]
            total = sum(r["coupling_kcal_mol"] for r in rows)
            shared = set.intersection(
                *[set(r["patches"].split()) for r in rows])

            receptor_side = sum(1 for r in rows if r["residue"] == "HIS")
            binder_side = size - receptor_side

            networks.append({
                "groups": " ".join(f"{r['aa']}{r['uniprot_pos']}" for r in rows),
                "size": size,
                "receptor_side_his": receptor_side,
                "engineered_his_needed": binder_side,
                "max_span_a": round(max(distances), 1),
                "min_separation_a": round(min(distances), 1),
                "max_normal_angle_deg": round(max(angles), 1),
                "total_coupling_kcal_mol": round(total, 3),
                "predicted_fold_change": round(fold_change(total), 1),
                "shared_patches": " ".join(sorted(shared)) or "none",
            })

    networks.sort(key=lambda n: -n["total_coupling_kcal_mol"])
    print(f"  {len(networks)} geometrically viable networks")
    for network in networks[:8]:
        print(f"    {network['groups']:<22} span "
              f"{network['max_span_a']:>5.1f} A  "
              f"{network['total_coupling_kcal_mol']:>6.2f} kcal/mol  "
              f"{network['predicted_fold_change']:>7.1f}x  "
              f"patches {network['shared_patches']}")

    # -- Outputs ------------------------------------------------------------
    pka_path = OUT_DIR / "03b_pka_estimates.csv"
    with pka_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pka_rows[0]))
        writer.writeheader()
        writer.writerows(pka_rows)

    network_path = OUT_DIR / "03b_networks.csv"
    if networks:
        with network_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(networks[0]))
            writer.writeheader()
            writer.writerows(networks)

    report_path = OUT_DIR / "03b_network_report.md"
    with report_path.open("w") as handle:
        w = handle.write
        w("# Step 03b - pKa and protonation networks\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/03b_protonation_networks.py`.\n\n")

        w("## Selectivity comes from the pKa shift, not the pKa\n\n")
        w("The linkage coefficient is the number of protons taken up on "
          "binding, `Δn(pH) = f_bound(pH) − f_free(pH)`. A group whose pKa "
          "does not change on binding contributes exactly nothing, however "
          "favourable its interaction. The design requirement is therefore "
          "to place each histidine so that **binding stabilises its "
          "protonated form** — burying it against a carboxylate is the "
          "mechanism that achieves this.\n\n")
        w(f"Estimates below assume a shift of {args.assumed_shift} pKa units "
          f"on binding, which is modest for a buried salt bridge. Coupling "
          f"is integrated numerically over the {PH_LOW}–{PH_HIGH} window "
          f"rather than taken from the two-state ceiling.\n\n")

        w("## Free-state pKa of the reachable groups\n\n")
        w("| Group | Role | pKa free | Titrating pKa | Bound | kcal/mol "
          "| Fold | Verdict |\n")
        w("|-------|------|---------:|--------------:|------:|---------:"
          "|-----:|---------|\n")
        for row in pka_rows:
            w(f"| {row['aa']}{row['uniprot_pos']} | {row['role']} "
              f"| {row['pka_free']} | {row['titrating_pka']} "
              f"| {row['assumed_pka_bound']} "
              f"| {row['coupling_kcal_mol']} | {row['fold_change']}× "
              f"| {row['verdict']} |\n")

        w("\n## Geometrically viable networks\n\n")
        w(f"A network is viable when its groups lie between "
          f"{NETWORK_SPAN[0]} and {NETWORK_SPAN[1]} Å apart — close enough "
          f"for one binder face to reach all of them, far enough that they "
          f"titrate against the binder rather than each other — and their "
          f"outward normals fall within {MAX_NORMAL_ANGLE}°.\n\n")

        if networks:
            w("| Groups | Receptor His | Engineered His | Span Å | Angle ° "
              "| kcal/mol | Fold | Patches |\n")
            w("|--------|-------------:|---------------:|-------:|--------:"
              "|---------:|-----:|---------|\n")
            for network in networks:
                w(f"| {network['groups']} | {network['receptor_side_his']} "
                  f"| {network['engineered_his_needed']} "
                  f"| {network['max_span_a']} "
                  f"| {network['max_normal_angle_deg']} "
                  f"| {network['total_coupling_kcal_mol']} "
                  f"| {network['predicted_fold_change']}× "
                  f"| {network['shared_patches']} |\n")
        else:
            w("None. No combination satisfies the distance and orientation "
              "constraints simultaneously.\n")

        w("\n## Limitations\n\n")
        w("pKa values in the complex cannot be computed before the complex "
          "exists, so the assumed shift is a stand-in for the quantity that "
          "actually decides the outcome. PROPKA is empirical with errors "
          "near one pKa unit, comparable to the effects reasoned about here. "
          "Rigorous coupling free energies would require constant-pH "
          "molecular dynamics, which is out of scope.\n\n")
        w("Every anchor is scored with the same engineered-histidine model, "
          "so all anchors return an identical coupling energy and the "
          "ranking among anchor-only networks is driven by how many groups "
          "they contain rather than by their quality. In this model anchor "
          "pKa governs *reliability* — whether the carboxylate is dependably "
          "ionised across the window — not the magnitude of the shift it "
          "induces. A more careful treatment would let anchor pKa and burial "
          "depth modulate the achievable shift, and would separate the "
          "networks that here score identically.\n\n")
        w("These numbers are a screen that removes groups whose free-state "
          "pKa already rules them out, and a geometric filter on which "
          "combinations a single binder could engage. They are not "
          "predictions of experimental pH selectivity.\n")

    print(f"\n  -> {pka_path.relative_to(ROOT)}")
    print(f"  -> {network_path.relative_to(ROOT)}")
    print(f"  -> {report_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

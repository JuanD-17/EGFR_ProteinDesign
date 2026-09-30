#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 07 - Structural validation of the pH variants.

Steps 05 and 06 established that the designed histidine placements engage
target handles more often than chance does. Neither says whether the
variant still folds, or still binds. Introducing three charged side chains
into a packed interface can do the opposite of what is intended:

    +His  ->  new contact  ->  repacking  ->  clash  ->  weaker binder

A variant with a good protonation network and a broken interface is worth
nothing, so the network score is never merged with the structural scores
here. Three layers are reported separately:

    A  Folding plausibility     does the binder still fold as designed?
    B  Interface quality        does it still bind, by AlphaFold2?
    C  pH-network geometry      is the designed network still formed?

Keeping them apart is what allows the eventual claim to be "these variants
retain a structurally plausible interface *and* the proposed protonation
network", rather than the much weaker "the score improved".

Two stages
----------
`--prepare` runs on CPU and does two things. It checks, for every proposed
substitution, whether a histidine side chain can physically reach the
handle from that backbone position without clashing into the binder or the
target: a geometric veto that AlphaFold2 will not give and that can kill a
variant before any GPU time is spent on it. It then writes the surviving
variants as complex FASTA records ready for re-prediction.

`--score` runs after the GPU re-prediction and reads the predicted
complexes back, reporting the three layers against the parent design.

Inputs:  03_Design/ph_variants/05_variants.csv
         03_Design/accepted/*.pdb
Outputs: 03_Design/validation/07_rotamer_check.csv
         03_Design/validation/07_variants_for_af2.fasta
         03_Design/validation/07_validation_report.md

Usage:
    python scripts/07_validate_ph_variants.py --prepare
    python scripts/07_validate_ph_variants.py --score --predictions <dir>
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "03_Design" / "validation"
DESIGNS = ROOT / "03_Design" / "accepted"
VARIANTS = ROOT / "03_Design" / "ph_variants" / "05_variants.csv"

_spec = importlib.util.spec_from_file_location(
    "step05", ROOT / "scripts" / "05_engineer_ph_switch.py")
step05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step05)

# Histidine geometry, measured from the 23 histidines resolved in 6ARU
# rather than recalled: CB-ND1 is 2.53 +/- 0.01 A and CB-NE2 3.67 +/- 0.01 A.
# Both nitrogens titrate, so both are probed; using a single averaged
# distance discards the shorter reach entirely.
HIS_TIP_DISTANCES = (2.53, 3.67)
HIS_RING_CLEARANCE = 2.9

# Contact between an imidazole nitrogen and the handle. The upper bound is
# 4.0 rather than a strict hydrogen-bond 3.2, because the interaction being
# designed is between a protonated histidine and a carboxylate, and a
# charge-assisted contact remains productive at that separation.
CONTACT_RANGE = (2.6, 4.0)

PROBE_DIRECTIONS = 200


def fibonacci_sphere(count: int) -> np.ndarray:
    indices = np.arange(count) + 0.5
    phi = np.arccos(1 - 2 * indices / count)
    theta = np.pi * (1 + 5 ** 0.5) * indices
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)])


DIRECTIONS = fibonacci_sphere(PROBE_DIRECTIONS)


BACKBONE = {"N", "CA", "C", "O", "CB", "OXT"}


def environment_without(residue, atoms_by_residue) -> np.ndarray:
    """All heavy atoms except the side chain being replaced.

    A substitution removes the existing side chain, so testing a histidine
    against it asks whether two side chains can occupy one position. They
    cannot, which is why an earlier version vetoed every substitution.
    Backbone atoms stay: they are not going anywhere.
    """
    key = (residue.get_parent().id, residue.id)
    coords = []
    for (chain_id, res_id), residue_atoms in atoms_by_residue.items():
        if (chain_id, res_id) == key:
            coords.extend(c for name, c in residue_atoms if name in BACKBONE)
        else:
            coords.extend(c for _, c in residue_atoms)
    return np.array(coords)


def rotamer_reaches(residue, handle_coord: np.ndarray,
                    environment: np.ndarray) -> tuple[bool, float, int]:
    """Can a histidine at this backbone position reach the handle?

    The imidazole tip is placed at every orientation on a sphere of
    HIS_TIP_DISTANCE around CB, positions clashing with the surroundings are
    discarded, and what remains is checked for one that lands in
    hydrogen-bonding range of the handle.

    This is a rotamer feasibility test, not a rotamer library search: it
    asks whether any placement of the tip is geometrically possible, which
    is the question that decides whether the substitution is worth
    predicting at all.
    """
    if "CB" not in residue:
        return False, math.inf, 0
    origin = np.array(residue["CB"].coord)

    tips = np.vstack([origin + DIRECTIONS * d for d in HIS_TIP_DISTANCES])

    # Reject tips colliding with the surroundings. `environment` must
    # already exclude the side chain being replaced.
    deltas = tips[:, None, :] - environment[None, :, :]
    nearest = np.sqrt((deltas ** 2).sum(axis=2)).min(axis=1)
    free = tips[nearest >= HIS_RING_CLEARANCE]

    if len(free) == 0:
        return False, math.inf, 0

    distances = np.linalg.norm(free - handle_coord, axis=1)
    in_range = int(((distances >= CONTACT_RANGE[0])
                    & (distances <= CONTACT_RANGE[1])).sum())
    return in_range > 0, float(distances.min()), in_range


def prepare() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    handles_raw = step05.load_handles()
    with VARIANTS.open() as fh:
        variants = list(csv.DictReader(fh))

    by_design: dict[str, list[dict]] = {}
    for row in variants:
        by_design.setdefault(row["design"], []).append(row)

    checks, records = [], []
    print(f"{len(variants)} variants across {len(by_design)} designs\n")

    for path in sorted(DESIGNS.glob("*.pdb")):
        name = path.stem
        if name not in by_design:
            continue

        model = step05.PDBParser(QUIET=True).get_structure(name, str(path))[0]
        target_id, binder_id = step05.identify_chains(model, 204)
        target_sequence, _ = step05.chain_sequence(model[target_id])
        binder_sequence, binder_residues = step05.chain_sequence(
            model[binder_id])

        atoms_by_residue = {
            (chain.id, r.id): [(a.get_name(), a.coord) for a in r
                               if a.element != "H"]
            for chain in model for r in chain if r.id[0] == " "
        }

        mapping = step05.map_target_to_uniprot(model[target_id])
        handles = {}
        for uniprot, handle in handles_raw.items():
            residue = mapping.get(uniprot)
            if residue is None or residue.get_resname() != handle["residue"]:
                continue
            atoms = [residue[n]
                     for n in step05.TITRATABLE_ATOMS[handle["residue"]]
                     if n in residue]
            if atoms:
                handles[uniprot] = np.mean([a.coord for a in atoms], axis=0)

        print(f"  {name}")
        for variant in by_design[name]:
            subs = variant["substitutions"].split()
            targets = [int(t[1:]) for t in variant["engages"].split()]

            feasible, details = True, []
            for substitution, uniprot in zip(subs, targets):
                position = int("".join(c for c in substitution[1:]
                                       if c.isdigit()))
                residue = binder_residues[position - 1]
                reaches, closest, count = rotamer_reaches(
                    residue, handles[uniprot],
                    environment_without(residue, atoms_by_residue))
                feasible &= reaches
                details.append((substitution, uniprot, reaches, closest, count))

                checks.append({
                    "design": name,
                    "variant": variant["substitutions"],
                    "substitution": substitution,
                    "engages": uniprot,
                    "rotamer_reaches": "yes" if reaches else "no",
                    "closest_approach_a": round(closest, 2)
                                          if closest != math.inf else "",
                    "n_contact_orientations": count,
                })

            flag = "OK " if feasible else "VETO"
            print(f"      [{flag}] {variant['substitutions']:<20} "
                  + "  ".join(
                      f"{s}->{u}:{'y' if r else 'n'}"
                      for s, u, r, _, _ in details))

            if feasible:
                mutated = list(binder_sequence)
                for substitution in subs:
                    position = int("".join(c for c in substitution[1:]
                                           if c.isdigit()))
                    mutated[position - 1] = "H"
                records.append({
                    "name": f"{name}__{variant['substitutions'].replace(' ', '_')}",
                    "design": name,
                    "substitutions": variant["substitutions"],
                    "engages": variant["engages"],
                    "coupling": variant["estimated_coupling_kcal_mol"],
                    "target": target_sequence,
                    "binder": "".join(mutated),
                })

        # The unmodified parent is predicted too: without it there is no
        # baseline to say whether the substitutions cost anything.
        records.append({
            "name": f"{name}__parent",
            "design": name,
            "substitutions": "none (parent)",
            "engages": "",
            "coupling": "0",
            "target": target_sequence,
            "binder": binder_sequence,
        })
        print()

    check_path = OUT_DIR / "07_rotamer_check.csv"
    with check_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(checks[0]))
        writer.writeheader()
        writer.writerows(checks)

    fasta_path = OUT_DIR / "07_variants_for_af2.fasta"
    with fasta_path.open("w") as fh:
        for record in records:
            fh.write(f">{record['name']} coupling={record['coupling']} "
                     f"engages={record['engages'].replace(' ', ',')}\n")
            fh.write(f"{record['target']}:{record['binder']}\n")

    index_path = OUT_DIR / "07_prediction_index.csv"
    with index_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["name", "design", "substitutions", "engages",
                            "coupling"], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    vetoed = sum(1 for c in checks if c["rotamer_reaches"] == "no")
    parents = sum(1 for r in records if r["substitutions"] == "none (parent)")

    report = OUT_DIR / "07_validation_report.md"
    with report.open("w") as fh:
        w = fh.write
        w("# Step 07 - Structural validation, preparation\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/07_validate_ph_variants.py --prepare`.\n\n")

        w("## Why the scores stay separate\n\n")
        w("Introducing charged side chains into a packed interface can make "
          "it worse. A variant with an elegant protonation network and a "
          "broken interface is worthless, so the network score is never "
          "combined with the structural scores. Three layers are reported "
          "independently: folding plausibility, interface quality, and "
          "whether the designed network survives.\n\n")

        w("## Rotamer feasibility\n\n")
        w(f"For each substitution, an imidazole nitrogen is placed at "
          f"{PROBE_DIRECTIONS} orientations on spheres of "
          f"{HIS_TIP_DISTANCES[0]} and {HIS_TIP_DISTANCES[1]} Å around CB, "
          f"the CB-ND1 and CB-NE2 distances measured from the 23 histidines "
          f"resolved in 6ARU. Orientations clashing with the "
          f"binder or the target within {HIS_RING_CLEARANCE} Å are "
          f"discarded, and the substitution passes only if a surviving "
          f"orientation lands {CONTACT_RANGE[0]}–{CONTACT_RANGE[1]} Å from "
          f"the handle.\n\n")
        w("This is a feasibility test rather than a rotamer library search. "
          "It answers whether the contact is geometrically possible at all, "
          "which is a veto AlphaFold2 does not provide and which is worth "
          "applying before spending GPU time.\n\n")
        w(f"{len(checks) - vetoed} of {len(checks)} substitutions can reach "
          f"their handle; {vetoed} cannot.\n\n")

        w("## Prepared for prediction\n\n")
        w(f"{len(records) - parents} variants and {parents} unmodified "
          f"parents, as complex records in "
          f"`07_variants_for_af2.fasta`.\n\n")
        w("The parents matter as much as the variants. Without predicting "
          "the unmodified design under identical settings there is no "
          "baseline, and any change in interface confidence could not be "
          "attributed to the substitutions rather than to prediction "
          "noise.\n\n")

        w("## Next\n\n")
        w("Re-predict every record on GPU, then:\n\n")
        w("```bash\n"
          "python scripts/07_validate_ph_variants.py --score "
          "--predictions <dir>\n"
          "```\n")

    print(f"  {len(checks) - vetoed}/{len(checks)} substitutions feasible")
    print(f"  {len(records)} records written "
          f"({parents} parents as baseline)")
    print(f"\n  -> {check_path.relative_to(ROOT)}")
    print(f"  -> {fasta_path.relative_to(ROOT)}")
    print(f"  -> {index_path.relative_to(ROOT)}")
    print(f"  -> {report.relative_to(ROOT)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--prepare", action="store_true",
                        help="CPU stage: rotamer veto and FASTA for AlphaFold2")
    parser.add_argument("--score", action="store_true",
                        help="read predicted complexes back and report layers")
    parser.add_argument("--predictions", type=str,
                        help="directory of predicted complexes, for --score")
    args = parser.parse_args()

    if args.prepare:
        return prepare()
    if args.score:
        raise SystemExit(
            "--score is not implemented yet: it needs the predicted "
            "complexes to exist so their format can be read rather than "
            "guessed. Run --prepare, predict on GPU, then come back.")
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

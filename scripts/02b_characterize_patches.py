#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 02b - Geometric characterisation of the candidate patches.

Step 02 groups pool residues by proximity to a seed. Proximity in a distance
matrix is not the same thing as a binding face: a set of residues can be
mutually close and still be scattered around a curved surface, or split
across a cleft, and a score computed from composition would not notice.

Nor does the residue-overlap deduplication in step 02 distinguish two
genuinely different faces from one face found twice from different seeds.
Two patches can share few residues and still occupy the same surface.

This script measures what the composition scores cannot:

  Compactness    diameter, radius of gyration, and how far the patch
                 deviates from a plane. A binder lands on something roughly
                 flat or gently curved, not on a 30 A sprawl.

  Continuity     whether the patch is one connected surface. Residues are
                 linked when any heavy atoms come within 5 A, and the number
                 of connected components is counted. More than one means the
                 patch is not a single face.

  Orientation    the outward normal of each patch, from a best-fit plane
                 oriented away from the protein centre. Two patches on the
                 same face have nearby centroids and roughly parallel
                 normals; two patches on different faces do not.

  Fab relation   graded rather than binary: direct contact, close
                 proximity, same region, or distant.

The output groups the patches into distinct faces, which is the thing needed
before committing GPU time to a design campaign.

Inputs:  02_Analysis/02_patches.csv, 02_Analysis/01_residue_annotation.csv
         01_Target/structure/6ARU.pdb
Outputs: 02_Analysis/02b_patch_geometry.csv
         02_Analysis/02b_geometry_report.md

Usage:
    python scripts/02b_characterize_patches.py
    python scripts/02b_characterize_patches.py --top 8
"""

from __future__ import annotations

import argparse
import csv
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "02_Analysis"
PATCHES = OUT_DIR / "02_patches.csv"
ANNOTATION = OUT_DIR / "01_residue_annotation.csv"
STRUCTURE = ROOT / "01_Target/structure/6ARU.pdb"

EGFR_CHAIN = "A"
FAB_CHAINS = ("B", "C")

# Heavy atoms within this distance make two residues neighbours on the
# surface. 5 A is about a van der Waals contact plus a water.
CONTACT_CUTOFF = 5.0

# Two patches count as the same face when their centroids are within this
# distance and their outward normals are within this angle.
SAME_FACE_DISTANCE = 12.0
SAME_FACE_ANGLE = 45.0

# Graded distance bands to the Fab, in Angstrom.
FAB_CONTACT = 4.0
FAB_PROXIMAL = 8.0
FAB_REGIONAL = 15.0


def load_patches(path: Path, limit: int) -> list[dict]:
    patches = []
    with path.open() as handle:
        for row in csv.DictReader(handle):
            patches.append({
                "rank": int(row["rank"]),
                "centre": int(row["centre"]),
                "score": float(row["score"]),
                "residues": [int(p) for p in row["residues"].split()],
            })
            if len(patches) >= limit:
                break
    return patches


def load_pdb_numbers(path: Path) -> dict[int, int]:
    mapping = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            if row["pdb_resnum"]:
                mapping[int(row["uniprot_pos"])] = int(row["pdb_resnum"])
    return mapping


def heavy_atoms(residue) -> list:
    return [a for a in residue if a.element != "H"]


def best_fit_normal(points: np.ndarray) -> np.ndarray:
    """Unit normal of the least-squares plane through a point cloud."""
    centred = points - points.mean(axis=0)
    _, _, vh = np.linalg.svd(centred, full_matrices=False)
    return vh[2] / np.linalg.norm(vh[2])


def planarity(points: np.ndarray, normal: np.ndarray) -> float:
    """RMS deviation from the best-fit plane, in Angstrom."""
    centred = points - points.mean(axis=0)
    return float(np.sqrt(np.mean((centred @ normal) ** 2)))


def connected_components(residues: list, cutoff: float) -> list[list[int]]:
    """Group residues linked by heavy-atom contacts."""
    atom_sets = [np.array([a.coord for a in heavy_atoms(r)]) for r in residues]
    count = len(residues)
    adjacency: dict[int, set[int]] = {i: set() for i in range(count)}

    for i in range(count):
        for j in range(i + 1, count):
            deltas = atom_sets[i][:, None, :] - atom_sets[j][None, :, :]
            if np.sqrt((deltas ** 2).sum(axis=2)).min() <= cutoff:
                adjacency[i].add(j)
                adjacency[j].add(i)

    seen: set[int] = set()
    groups = []
    for start in range(count):
        if start in seen:
            continue
        stack, group = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            group.append(node)
            for neighbour in adjacency[node] - seen:
                seen.add(neighbour)
                stack.append(neighbour)
        groups.append(sorted(group))
    return sorted(groups, key=len, reverse=True)


def fab_relation(residues: list, fab_coords: np.ndarray) -> tuple[float, str]:
    coords = np.array([a.coord for r in residues for a in heavy_atoms(r)])
    deltas = coords[:, None, :] - fab_coords[None, :, :]
    minimum = float(np.sqrt((deltas ** 2).sum(axis=2)).min())

    if minimum <= FAB_CONTACT:
        band = "direct contact"
    elif minimum <= FAB_PROXIMAL:
        band = "proximal"
    elif minimum <= FAB_REGIONAL:
        band = "same region"
    else:
        band = "distant"
    return minimum, band


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--top", type=int, default=6,
                        help="patches to characterise (default 6)")
    args = parser.parse_args()

    patches = load_patches(PATCHES, args.top)
    pdb_of_uniprot = load_pdb_numbers(ANNOTATION)

    model = PDBParser(QUIET=True).get_structure("6aru", str(STRUCTURE))[0]
    chain = model[EGFR_CHAIN]

    fab_coords = np.array([
        a.coord for chain_id in FAB_CHAINS for residue in model[chain_id]
        if residue.id[0] == " " for a in heavy_atoms(residue)
    ])

    protein_centre = np.array([
        a.coord for residue in chain if residue.id[0] == " "
        for a in heavy_atoms(residue)
    ]).mean(axis=0)

    print(f"Characterising the top {len(patches)} patches\n")

    results = []
    for patch in patches:
        residues, cb_points = [], []
        for position in patch["residues"]:
            pdb_number = pdb_of_uniprot.get(position)
            if pdb_number is None:
                continue
            try:
                residue = chain[(" ", pdb_number, " ")]
            except KeyError:
                continue
            residues.append(residue)
            atom = residue["CB"] if "CB" in residue else residue["CA"]
            cb_points.append(atom.coord)

        points = np.array(cb_points)
        centroid = points.mean(axis=0)

        # Compactness.
        deltas = points[:, None, :] - points[None, :, :]
        pairwise = np.sqrt((deltas ** 2).sum(axis=2))
        diameter = float(pairwise.max())
        gyration = float(np.sqrt(((points - centroid) ** 2).sum(axis=1).mean()))

        normal = best_fit_normal(points)
        # Orient outward, away from the protein centre.
        if np.dot(normal, centroid - protein_centre) < 0:
            normal = -normal
        flatness = planarity(points, normal)

        # Continuity. Knowing *which* residues detach matters: a patch whose
        # core is continuous apart from one stray residue is repairable by
        # dropping it, while one that splits down the middle is not.
        components = connected_components(residues, CONTACT_CUTOFF)
        largest = len(components[0]) if components else 0
        present = [p for p in patch["residues"] if pdb_of_uniprot.get(p)
                   and (" ", pdb_of_uniprot[p], " ") in chain]
        core = sorted(present[i] for i in components[0]) if components else []
        orphans = sorted(set(present) - set(core))

        # Fab.
        fab_distance, fab_band = fab_relation(residues, fab_coords)

        record = {
            "rank": patch["rank"],
            "centre": patch["centre"],
            "score": patch["score"],
            "n_residues": len(residues),
            "diameter_a": round(diameter, 1),
            "radius_gyration_a": round(gyration, 1),
            "planarity_rms_a": round(flatness, 2),
            "components": len(components),
            "largest_component": largest,
            "fragmented": "yes" if len(components) > 1 else "no",
            "fab_min_dist_a": round(fab_distance, 1),
            "fab_relation": fab_band,
            "core_residues": " ".join(map(str, core)),
            "detached_residues": " ".join(map(str, orphans)),
            "_centroid": centroid,
            "_normal": normal,
            "_residues": patch["residues"],
        }
        results.append(record)

        print(f"  Patch {patch['rank']:>2} (centre {patch['centre']}): "
              f"{len(residues)} residues, diameter {diameter:.1f} A, "
              f"Rg {gyration:.1f} A, planarity {flatness:.2f} A")
        print(f"       {len(components)} component(s), largest {largest}; "
              f"Fab {fab_distance:.1f} A, {fab_band}")
        if orphans:
            sizes = ", ".join(str(len(c)) for c in components)
            print(f"       fragmented into {sizes}; detached: "
                  f"{', '.join(map(str, orphans))}")

    # -- Face grouping ------------------------------------------------------
    print("\nComparing patches pairwise\n")
    comparisons = []
    for i, a in enumerate(results):
        for b in results[i + 1:]:
            separation = float(np.linalg.norm(a["_centroid"] - b["_centroid"]))
            cosine = float(np.clip(np.dot(a["_normal"], b["_normal"]), -1, 1))
            angle = math.degrees(math.acos(cosine))
            shared = set(a["_residues"]) & set(b["_residues"])
            union = set(a["_residues"]) | set(b["_residues"])
            jaccard = len(shared) / len(union) if union else 0.0
            same = (separation <= SAME_FACE_DISTANCE and angle <= SAME_FACE_ANGLE)
            # A verdict that flips under a small change of threshold is not a
            # verdict. Flag those instead of reporting them as settled.
            borderline = (separation <= SAME_FACE_DISTANCE
                          and abs(angle - SAME_FACE_ANGLE) <= 10.0)
            comparisons.append({
                "patch_a": a["rank"], "patch_b": b["rank"],
                "centroid_sep_a": round(separation, 1),
                "normal_angle_deg": round(angle, 1),
                "shared_residues": len(shared),
                "jaccard": round(jaccard, 3),
                "same_face": "yes" if same else "no",
                "borderline": "yes" if borderline else "no",
            })
            verdict = "SAME FACE" if same else "different faces"
            if borderline:
                verdict += "  <-- BORDERLINE, within 10 deg of the cutoff"
            print(f"  {a['rank']} vs {b['rank']}: centroids {separation:5.1f} A, "
                  f"normals {angle:5.1f} deg, shared {len(shared):>2}, {verdict}")

    # Group into faces by transitive closure of the same-face relation.
    ranks = [r["rank"] for r in results]
    parent = {rank: rank for rank in ranks}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for comparison in comparisons:
        if comparison["same_face"] == "yes":
            a, b = find(comparison["patch_a"]), find(comparison["patch_b"])
            if a != b:
                parent[b] = a

    faces: dict[int, list[int]] = {}
    for rank in ranks:
        faces.setdefault(find(rank), []).append(rank)

    print(f"\n  {len(faces)} distinct face(s) among {len(results)} patches")
    for index, (_, members) in enumerate(sorted(faces.items()), start=1):
        best = min(members)
        print(f"    Face {index}: patches {members}  (best ranked: {best})")

    # -- Outputs ------------------------------------------------------------
    csv_path = OUT_DIR / "02b_patch_geometry.csv"
    columns = [k for k in results[0] if not k.startswith("_")]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    report_path = OUT_DIR / "02b_geometry_report.md"
    with report_path.open("w") as handle:
        w = handle.write
        w("# Step 02b - Patch geometry\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/02b_characterize_patches.py`.\n\n")
        w("Step 02 groups residues by proximity and scores them by "
          "composition. Neither operation can tell whether a patch is one "
          "usable binding face, and neither can tell whether two patches are "
          "the same face found from different seeds. This step measures "
          "both.\n\n")

        w("## Geometry\n\n")
        w("| Patch | Centre | Residues | Diameter Å | Rg Å | Planarity Å "
          "| Components | Fab | Relation |\n")
        w("|-------|--------|----------|------------|------|-------------"
          "|------------|-----|----------|\n")
        for record in results:
            w(f"| {record['rank']} | {record['centre']} "
              f"| {record['n_residues']} | {record['diameter_a']} "
              f"| {record['radius_gyration_a']} | {record['planarity_rms_a']} "
              f"| {record['components']} | {record['fab_min_dist_a']} "
              f"| {record['fab_relation']} |\n")

        w("\nPlanarity is the RMS deviation of the CB atoms from their "
          "best-fit plane: small values mean a flat face, large values a "
          "patch wrapped around curvature. Components counts connected "
          "groups under a 5 Å heavy-atom contact rule; more than one means "
          "the residues do not form a single continuous surface.\n")

        w("\n## Are these distinct surfaces?\n\n")
        w("| A | B | Centroid separation Å | Normal angle ° | Shared residues "
          "| Jaccard | Same face |\n")
        w("|---|---|----------------------|----------------|-----------------"
          "|---------|----------|\n")
        for comparison in comparisons:
            note = " (borderline)" if comparison["borderline"] == "yes" else ""
            w(f"| {comparison['patch_a']} | {comparison['patch_b']} "
              f"| {comparison['centroid_sep_a']} "
              f"| {comparison['normal_angle_deg']} "
              f"| {comparison['shared_residues']} | {comparison['jaccard']} "
              f"| {comparison['same_face']}{note} |\n")

        flagged = [c for c in comparisons if c["borderline"] == "yes"]
        if flagged:
            w(f"\n{len(flagged)} pair(s) fall within 10° of the "
              f"{SAME_FACE_ANGLE}° cutoff and their classification would flip "
              f"under a modest change of threshold. They are reported as "
              f"unresolved rather than decided:\n\n")
            for comparison in flagged:
                w(f"- Patches {comparison['patch_a']} and "
                  f"{comparison['patch_b']}: centroids "
                  f"{comparison['centroid_sep_a']} Å apart, normals "
                  f"{comparison['normal_angle_deg']}° apart. Most likely "
                  f"adjacent regions of one curved surface rather than two "
                  f"independent faces.\n")

        w(f"\nTwo patches are called the same face when their centroids are "
          f"within {SAME_FACE_DISTANCE} Å and their outward normals within "
          f"{SAME_FACE_ANGLE}°. Residue overlap alone cannot decide this: two "
          f"patches can share few residues and still sit on one surface, "
          f"which is the failure mode the step 02 deduplication has.\n")

        w(f"\n**{len(faces)} distinct face(s) among the top "
          f"{len(results)} patches.**\n\n")
        for index, (_, members) in enumerate(sorted(faces.items()), start=1):
            w(f"- Face {index}: patches {', '.join(map(str, sorted(members)))}"
              f" — best ranked is patch {min(members)}\n")

        w("\nA design campaign needs one target per distinct face. Patches "
          "on the same face are alternative framings of one surface and do "
          "not diversify the portfolio.\n")

    print(f"\n  -> {csv_path.relative_to(ROOT)}")
    print(f"  -> {report_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

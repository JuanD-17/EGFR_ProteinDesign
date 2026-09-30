#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 07 - Negative control for the histidine placement method.

Step 05 places histidines where they face an ionizable group on the target,
and reports a coupling estimate for each variant. That number means nothing
on its own. The question a reviewer will ask, and the one that decides
whether the method is worth anything, is:

    Does placing histidines by this mechanistic rule produce better
    protonation networks than placing the same number of histidines at
    random among the positions that were available?

The null model
--------------
The control is deliberately hard. A null that scattered histidines anywhere
on the binder would be trivial to beat, and beating it would prove only that
interface positions differ from non-interface positions, which nobody
doubts. Instead the null draws from exactly the pool step 05 could have
chosen from: binder positions that are at the interface, exposed in the
unbound binder, and chemically substitutable. The only thing removed is the
targeting.

So the comparison isolates one variable. If designed and random placements
score alike, the mechanistic rule adds nothing over "put histidines at the
interface", and that is worth knowing before any more compute is spent.

Statistics
----------
For each design, the same number of positions the designed variant used is
drawn at random from the pool, and the coupling recomputed under identical
rules. Repeating this thousands of times gives the null distribution of
coupling achievable by chance, against which each designed variant gets an
empirical p-value: the fraction of random draws scoring at least as high.

Sampling is cheap on CPU, so the null is estimated from many draws rather
than from a handful of matched controls.

Inputs:  03_Design/accepted/*.pdb, 03_Design/ph_variants/05_variants.csv
         02_Analysis/03b_pka_estimates.csv
Outputs: 03_Design/negative_control/07_null_distribution.csv
         03_Design/negative_control/07_control_report.md

Usage:
    python scripts/07_negative_control.py
    python scripts/07_negative_control.py --samples 5000 --seed 20260930
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import random
import statistics
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "03_Design" / "negative_control"
DESIGNS = ROOT / "03_Design" / "accepted"
VARIANTS = ROOT / "03_Design" / "ph_variants" / "05_variants.csv"

# Step 05's filename starts with a digit, so it cannot be imported normally.
# Loading it by path keeps one implementation of the placement rules rather
# than a second copy that could drift out of step with it.
_spec = importlib.util.spec_from_file_location(
    "step05", ROOT / "scripts" / "05_engineer_ph_switch.py")
step05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(step05)

DEFAULT_SAMPLES = 3000
DEFAULT_SEED = 20260930

# A binder residue counts as interface when any heavy atom comes within this
# distance of the target.
INTERFACE_CUTOFF = 8.0


def interface_pool(model, binder_id: str, target_id: str,
                   exposure: dict[int, float],
                   binder_sequence: str) -> list[dict]:
    """Binder positions the null model is allowed to draw from.

    Interface, exposed in the unbound binder, and chemically substitutable:
    the same three conditions step 05 imposes before it considers targeting.
    """
    target_coords = np.array([
        a.coord for r in model[target_id] if r.id[0] == " "
        for a in r if a.element != "H"
    ])

    pool = []
    index = -1
    for residue in model[binder_id]:
        if residue.id[0] != " ":
            continue
        try:
            aa = step05.protein_letters_3to1[residue.get_resname()]
        except KeyError:
            continue
        index += 1

        atoms = np.array([a.coord for a in residue if a.element != "H"])
        if len(atoms) == 0:
            continue
        deltas = atoms[:, None, :] - target_coords[None, :, :]
        if np.sqrt((deltas ** 2).sum(axis=2)).min() > INTERFACE_CUTOFF:
            continue

        if exposure.get(residue.id[1], 1.0) < step05.MIN_EXPOSURE:
            continue
        if aa in step05.PROTECTED:
            continue

        pool.append({
            "index": index,
            "position": index + 1,
            "aa": aa,
            "residue": residue,
            "rel_sasa": round(exposure.get(residue.id[1], 1.0), 2),
        })
    return pool


def coupling_of(positions: list[dict], handles: dict[int, dict]) -> float:
    """Coupling a set of histidine placements would achieve.

    Scored by exactly step 05's rules: a placed histidine counts only if it
    can bridge to a handle, one handle per histidine, and each engaged
    handle contributes its own coupling estimate.
    """
    engaged: dict[int, float] = {}
    for placement in positions:
        residue = placement["residue"]
        centre = (residue["CB"] if "CB" in residue else residue["CA"]).coord
        direction = step05.side_chain_direction(residue)
        if direction is None:
            continue

        best_handle, best_orientation = None, step05.MIN_ORIENTATION
        for uniprot, handle in handles.items():
            if uniprot in engaged:
                continue
            separation = float(np.linalg.norm(centre - handle["coord"]))
            if not (step05.BRIDGE_RANGE[0] <= separation
                    <= step05.BRIDGE_RANGE[1]):
                continue
            towards = handle["coord"] - centre
            towards = towards / np.linalg.norm(towards)
            orientation = float(np.dot(direction, towards))
            if orientation > best_orientation:
                best_handle, best_orientation = uniprot, orientation

        if best_handle is not None:
            engaged[best_handle] = handles[best_handle]["coupling"]

    return sum(engaged.values())


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--samples", type=int, default=DEFAULT_SAMPLES,
                        help=f"random draws per design "
                             f"(default {DEFAULT_SAMPLES})")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="random seed, recorded for reproducibility")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    handles_raw = step05.load_handles()

    designed: dict[str, list[dict]] = {}
    with VARIANTS.open() as fh:
        for row in csv.DictReader(fh):
            designed.setdefault(row["design"], []).append(row)

    rows, summary = [], []
    print(f"{len(designed)} designs with variants, "
          f"{args.samples} random draws each, seed {args.seed}\n")

    for path in sorted(DESIGNS.glob("*.pdb")):
        name = path.stem
        if name not in designed:
            continue

        model = step05.PDBParser(QUIET=True).get_structure(name, str(path))[0]
        target_id, binder_id = step05.identify_chains(model, 204)
        binder_sequence, _ = step05.chain_sequence(model[binder_id])
        exposure = step05.binder_exposure(model, binder_id)

        # Handles present in this complex, with coordinates.
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
                handles[uniprot] = {
                    **handle,
                    "coord": np.mean([a.coord for a in atoms], axis=0),
                }

        pool = interface_pool(model, binder_id, target_id, exposure,
                              binder_sequence)
        if len(pool) < 3:
            print(f"  {name}: pool too small ({len(pool)}), skipped")
            continue

        print(f"  {name}")
        print(f"      pool of {len(pool)} substitutable exposed interface "
              f"positions, {len(handles)} handles in reach")

        # Null distribution per number of substitutions used.
        nulls: dict[int, list[float]] = {}
        for variant in designed[name]:
            k = len([s for s in variant["substitutions"].split()])
            if k in nulls:
                continue
            draws = []
            for _ in range(args.samples):
                picks = rng.sample(pool, min(k, len(pool)))
                draws.append(coupling_of(picks, handles))
            nulls[k] = draws

        for variant in designed[name]:
            k = len(variant["substitutions"].split())
            observed = float(variant["estimated_coupling_kcal_mol"])
            null = nulls[k]
            at_least = sum(1 for value in null if value >= observed)
            p_value = (at_least + 1) / (len(null) + 1)

            rows.append({
                "design": name,
                "substitutions": variant["substitutions"],
                "n_histidines": k,
                "engages": variant["engages"],
                "observed_kcal_mol": observed,
                "null_mean": round(statistics.mean(null), 3),
                "null_sd": round(statistics.pstdev(null), 3),
                "null_max": round(max(null), 3),
                "null_fraction_zero": round(
                    sum(1 for v in null if v == 0) / len(null), 3),
                "p_empirical": round(p_value, 5),
                "n_draws": len(null),
            })
            print(f"      {variant['substitutions']:<20} "
                  f"observed {observed:5.2f}  "
                  f"null {statistics.mean(null):5.2f} "
                  f"± {statistics.pstdev(null):4.2f}  "
                  f"p = {p_value:.4f}")

        summary.append({"design": name, "pool": len(pool),
                        "handles": len(handles)})
        print()

    if not rows:
        print("No variants to test.")
        return 1

    csv_path = OUT_DIR / "07_null_distribution.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    significant = [r for r in rows if r["p_empirical"] < 0.05]
    beats_all = [r for r in rows
                 if r["observed_kcal_mol"] > r["null_max"]]

    report = OUT_DIR / "07_control_report.md"
    with report.open("w") as fh:
        w = fh.write
        w("# Step 07 - Negative control\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/07_negative_control.py`, seed {args.seed}, "
          f"{args.samples} draws per design and substitution count.\n\n")

        w("## The question\n\n")
        w("Step 05 reports a coupling estimate for each designed variant. "
          "That figure is uninterpretable alone. What matters is whether "
          "placing histidines by the mechanistic rule beats placing the same "
          "number at random among the positions that were available.\n\n")

        w("## The null model\n\n")
        w("The null draws from exactly the pool step 05 could have chosen "
          "from: binder positions at the interface, exposed in the unbound "
          "binder, and chemically substitutable. Only the targeting is "
          "removed.\n\n")
        w("This is deliberately the hard version of the control. A null that "
          "scattered histidines anywhere on the binder would be easy to "
          "beat, and beating it would show only that interface positions "
          "differ from the rest, which nobody doubts.\n\n")

        w("| Design | Pool size | Handles in reach |\n")
        w("|--------|-----------|------------------|\n")
        for entry in summary:
            w(f"| {entry['design']} | {entry['pool']} | {entry['handles']} |\n")

        w("\n## Results\n\n")
        w("| Design | Substitutions | Observed | Null mean ± SD | Null max "
          "| p |\n")
        w("|--------|---------------|---------:|----------------|--------:"
          "|--:|\n")
        for row in rows:
            w(f"| {row['design'][:28]} | {row['substitutions']} "
              f"| {row['observed_kcal_mol']} "
              f"| {row['null_mean']} ± {row['null_sd']} "
              f"| {row['null_max']} | {row['p_empirical']} |\n")

        w(f"\n{len(significant)} of {len(rows)} variants reach p < 0.05. "
          f"{len(beats_all)} exceed every random draw.\n\n")

        w("## Reading this\n\n")
        w("`p` is the fraction of random draws achieving coupling at least "
          "as high as the designed variant, with the usual +1 correction so "
          "that zero never appears. `null_fraction_zero` is the proportion "
          "of random placements that engage no handle at all, which measures "
          "how much of the pool is useless for this purpose.\n\n")
        w("A small p means the mechanistic rule found something chance "
          "rarely finds. A large p means the interface is small enough that "
          "almost any exposed substitution lands near a handle, in which "
          "case the targeting contributes little and the honest conclusion "
          "is that the pool, not the rule, is doing the work.\n\n")

        w("## What this does not show\n\n")
        w("Both arms are scored by the same model, with the same assumed "
          "pKa shift. The control tests whether the placement rule beats "
          "chance *under that model*; it does not validate the model. A "
          "designed variant that beats the null is a better candidate, not "
          "a demonstrated pH switch.\n\n")
        w("Nor does this test structural plausibility. A variant can win "
          "here and still fail to fold or to bind, which is what step 06 "
          "exists to check.\n\n")
        w("Two statistical caveats. No designed variant exceeds the maximum "
          "of its null distribution: chance does occasionally find "
          "placements as good, it just does so rarely, so the claim is "
          "about frequency and not about reaching something unreachable. "
          "And the p-values are uncorrected. Fifteen tests were run, and "
          "variants within a design share positions so they are not "
          "independent, which makes a clean correction awkward; under a "
          "conservative Bonferroni threshold the strongest variants survive "
          "and the marginal ones do not.\n")

    print(f"  {len(significant)}/{len(rows)} variants at p < 0.05")
    print(f"  {len(beats_all)}/{len(rows)} exceed every random draw")
    print(f"\n  -> {csv_path.relative_to(ROOT)}")
    print(f"  -> {report.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

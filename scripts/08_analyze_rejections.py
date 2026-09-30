#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 08 - Where the design campaign loses its candidates.

A campaign that accepts 5 designs out of 56 trajectories has thrown away 51,
and the useful question is not how many survived but *where* the others
died. The answer decides whether more compute is worth spending, and on
what.

Three places a candidate can be lost, and they mean different things:

  Trajectory     the hallucination never converged. Wasted compute, but
                 cheap, because these fail early.
  MPNN           a backbone was produced but no redesigned sequence
                 survived re-prediction. The expensive failure.
  Filters        a complex was produced and scored, and one or more
                 thresholds rejected it. The informative failure, because
                 it says which property is limiting.

That last category is the one that could justify changing something. If
rejections cluster on a peripheral criterion, such as surface
hydrophobicity, relaxing it costs little. If they cluster on AlphaFold2's
confidence in the interface, relaxing it means accepting designs the model
does not believe in, and the honest response is to accept the rate and run
more trajectories instead.

The script reads a run directory as BindCraft leaves it, so it also works
for comparing one run against another: a rejection profile that shifts
between runs means something changed, and one that holds means the rate is
a property of the target rather than of luck.

Inputs:  a BindCraft run directory (or several)
Outputs: 03_Design/rejection_analysis/08_rejection_by_filter.csv
         03_Design/rejection_analysis/08_run_summary.csv
         03_Design/rejection_analysis/08_rejection_report.md

Usage:
    python scripts/08_analyze_rejections.py --runs 03_Design/run1/gpu0 \\
                                                   03_Design/run1/gpu1
    python scripts/08_analyze_rejections.py --runs 03_Design/run*/gpu* \\
                                            --label-by parent
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "03_Design" / "rejection_analysis"

# Filters that reject on AlphaFold2's own confidence in the prediction.
# Separating these matters: they are not tunable in the way a composition
# threshold is, because relaxing them means accepting designs the model
# does not believe in.
CONFIDENCE_FILTERS = {
    "i_pAE", "i_pTM", "pLDDT", "pAE", "pTM", "i_pLDDT", "ss_pLDDT",
    "Binder_pLDDT", "Binder_pTM", "Binder_pAE",
    "Trajectory_logits_pLDDT", "Trajectory_softmax_pLDDT",
    "Trajectory_one-hot_pLDDT", "Trajectory_final_pLDDT",
}

TRAJECTORY_LINE = re.compile(r"Starting trajectory: (\S+)")
SUCCESS_LINE = re.compile(r"Trajectory successful")
NO_MPNN_LINE = re.compile(r"No accepted MPNN designs found")
DURATION_LINE = re.compile(
    r"took: (\d+) hours?, (\d+) minutes?, (\d+) seconds?")


def read_failures(run: Path) -> Counter:
    """Rejection counts per filter, from BindCraft's failure table."""
    path = run / "failure_csv.csv"
    counts: Counter = Counter()
    if not path.exists():
        return counts
    with path.open() as fh:
        for row in csv.DictReader(fh):
            for key, value in row.items():
                try:
                    number = int(float(value))
                except (TypeError, ValueError):
                    continue
                if number > 0:
                    counts[key] += number
    return counts


def read_log(path: Path) -> dict:
    """Trajectory counts and durations from a run log."""
    if not path.exists():
        return {}
    text = path.read_text(errors="ignore")
    durations = [int(h) * 3600 + int(m) * 60 + int(s)
                 for h, m, s in DURATION_LINE.findall(text)]
    return {
        "started": len(TRAJECTORY_LINE.findall(text)),
        "successful": len(SUCCESS_LINE.findall(text)),
        "no_mpnn_accepted": len(NO_MPNN_LINE.findall(text)),
        "durations_s": durations,
    }


def count_pdbs(run: Path, folder: str) -> int:
    directory = run / folder
    if not directory.exists():
        return 0
    return len(list(directory.glob("*.pdb")))


def find_log(run: Path) -> Path | None:
    """The log belonging to a run, searched by the run's directory name."""
    stem = run.name
    parent = run.parent.name
    # Logs from different runs share a device name, so the run directory's
    # parent disambiguates them: run1/gpu0 and run2/gpu0 both end in gpu0.
    for candidate in [
        ROOT / "03_Design" / "logs" / f"{parent}_bc_{stem}.log",
        ROOT / "03_Design" / "logs" / f"bc_{stem}.log",
        run.parent / f"bc_{stem}.log",
        run / "log.txt",
    ]:
        if candidate.exists():
            return candidate
    matches = list((ROOT / "03_Design" / "logs").glob(f"*{stem}*.log"))
    return matches[0] if matches else None


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--runs", nargs="+", required=True,
                        help="BindCraft run directories")
    parser.add_argument("--label-by", choices=["name", "parent"],
                        default="name",
                        help="group runs by their own name or by their "
                             "parent directory (default: name)")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    runs = [Path(r) for r in args.runs]
    missing = [r for r in runs if not r.is_dir()]
    if missing:
        raise SystemExit("not directories: "
                         + ", ".join(str(m) for m in missing))

    per_run, all_filters = [], Counter()
    grouped: dict[str, Counter] = defaultdict(Counter)

    print(f"{len(runs)} run directories\n")

    for run in runs:
        label = run.name if args.label_by == "name" else run.parent.name
        failures = read_failures(run)
        log = find_log(run)
        stats = read_log(log) if log else {}

        accepted = count_pdbs(run, "Accepted")
        rejected = count_pdbs(run, "Rejected")
        trajectories = count_pdbs(run, "Trajectory")
        mpnn = count_pdbs(run, "MPNN")

        started = stats.get("started", 0)
        successful = stats.get("successful", 0)
        durations = stats.get("durations_s", [])

        record = {
            "run": str(run.relative_to(ROOT)) if run.is_absolute()
                   else str(run),
            "group": label,
            "trajectories_started": started,
            "trajectories_successful": successful,
            "trajectory_success_pct": round(
                100 * successful / started, 1) if started else "",
            "no_mpnn_accepted": stats.get("no_mpnn_accepted", 0),
            "accepted_pdb": accepted,
            "rejected_pdb": rejected,
            "trajectory_pdb": trajectories,
            "mpnn_pdb": mpnn,
            "acceptance_per_trajectory_pct": round(
                100 * accepted / started, 2) if started else "",
            "median_trajectory_s": round(
                sorted(durations)[len(durations) // 2]) if durations else "",
            "total_rejections_logged": sum(failures.values()),
            "log": str(log.name) if log else "not found",
        }
        per_run.append(record)
        all_filters.update(failures)
        grouped[label].update(failures)

        print(f"  {run}")
        print(f"      started {started}, successful {successful}, "
              f"accepted {accepted}"
              + (f", median {record['median_trajectory_s']}s"
                 if durations else ""))
        if not log:
            print("      no log found; trajectory counts unavailable")

    if not all_filters:
        print("\nNo failure tables found. Nothing to analyse.")
        return 1

    total = sum(all_filters.values())
    confidence = sum(v for k, v in all_filters.items()
                     if k in CONFIDENCE_FILTERS)
    other = total - confidence

    print(f"\n{total} rejections logged across all runs")
    print(f"  {confidence} ({100 * confidence / total:.0f}%) on AlphaFold2 "
          f"confidence")
    print(f"  {other} ({100 * other / total:.0f}%) on other criteria\n")
    print("  top filters:")
    for name, count in all_filters.most_common(8):
        kind = "confidence" if name in CONFIDENCE_FILTERS else "other"
        print(f"      {name:<34} {count:>5}  ({kind})")

    # -- Outputs ------------------------------------------------------------
    by_filter = OUT_DIR / "08_rejection_by_filter.csv"
    with by_filter.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["filter", "category", "rejections",
                         "percent_of_total"])
        for name, count in all_filters.most_common():
            writer.writerow([
                name,
                "confidence" if name in CONFIDENCE_FILTERS else "other",
                count, round(100 * count / total, 2),
            ])

    summary = OUT_DIR / "08_run_summary.csv"
    with summary.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(per_run[0]))
        writer.writeheader()
        writer.writerows(per_run)

    report = OUT_DIR / "08_rejection_report.md"
    with report.open("w") as fh:
        w = fh.write
        w("# Step 08 - Where the campaign loses candidates\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/08_analyze_rejections.py` over "
          f"{len(runs)} run directories.\n\n")

        w("## Per run\n\n")
        w("| Run | Started | Successful | Accepted | Success % "
          "| Accept % | Median s |\n")
        w("|-----|--------:|-----------:|---------:|----------:"
          "|---------:|---------:|\n")
        for record in per_run:
            w(f"| {record['run']} | {record['trajectories_started']} "
              f"| {record['trajectories_successful']} "
              f"| {record['accepted_pdb']} "
              f"| {record['trajectory_success_pct']} "
              f"| {record['acceptance_per_trajectory_pct']} "
              f"| {record['median_trajectory_s']} |\n")

        started = sum(r["trajectories_started"] for r in per_run)
        accepted = sum(r["accepted_pdb"] for r in per_run)
        if started:
            w(f"\nOverall, {accepted} of {started} trajectories yielded an "
              f"accepted design, {100 * accepted / started:.1f} percent.\n")

        w("\n## What rejects them\n\n")
        w(f"{confidence} of {total} rejections ({100 * confidence / total:.0f}"
          f"%) are on AlphaFold2's confidence in the re-predicted complex "
          f"rather than on composition or geometry.\n\n")
        w("| Filter | Category | Rejections | % |\n")
        w("|--------|----------|-----------:|--:|\n")
        for name, count in all_filters.most_common(15):
            w(f"| {name} "
              f"| {'confidence' if name in CONFIDENCE_FILTERS else 'other'} "
              f"| {count} | {100 * count / total:.1f} |\n")

        w("\n## Reading this\n\n")
        w("The split between confidence filters and the rest is what decides "
          "whether the acceptance rate can be improved cheaply.\n\n")
        w("Rejections concentrated on peripheral criteria, such as surface "
          "hydrophobicity or unsatisfied hydrogen bonds, would suggest a "
          "threshold worth revisiting: those measure properties of a design "
          "that is otherwise believed to bind.\n\n")
        w("Rejections concentrated on i_pAE, i_pTM and pLDDT do not. Those "
          "are the model's confidence that the complex forms at all, and "
          "relaxing them means accepting designs AlphaFold2 does not believe "
          "in, which moves the failure from the filter to the assay. The "
          "honest response to a low rate of this kind is to accept it and "
          "run more trajectories.\n")

        if len(grouped) > 1:
            w("\n## Between runs\n\n")
            w("A rejection profile that holds across runs indicates the "
              "acceptance rate is a property of the target and protocol "
              "rather than of chance in any one campaign.\n\n")
            names = sorted(grouped)
            w("| Filter | " + " | ".join(names) + " |\n")
            w("|--------|" + "|".join(["---:"] * len(names)) + "|\n")
            for name, _ in all_filters.most_common(10):
                cells = " | ".join(str(grouped[n][name]) for n in names)
                w(f"| {name} | {cells} |\n")

    metadata = OUT_DIR / "08_rejection_analysis.json"
    metadata.write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        "runs": [r["run"] for r in per_run],
        "totals": {
            "trajectories_started": sum(
                r["trajectories_started"] for r in per_run),
            "accepted": sum(r["accepted_pdb"] for r in per_run),
            "rejections_logged": total,
            "confidence_rejections": confidence,
            "other_rejections": other,
        },
        "by_filter": dict(all_filters),
        "per_run": per_run,
    }, indent=2) + "\n")

    print(f"\n  -> {by_filter.relative_to(ROOT)}")
    print(f"  -> {summary.relative_to(ROOT)}")
    print(f"  -> {report.relative_to(ROOT)}")
    print(f"  -> {metadata.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

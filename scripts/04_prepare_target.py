#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 04 - Prepare the trimmed target and BindCraft configuration.

BindCraft's own guidance is to trim the target PDB to the smallest structure
that still presents the epitope, because AlphaFold2 backpropagation scales
badly with target size and GPU memory is the binding constraint on the free
tiers this project is limited to. The 6ARU ectodomain is 609 resolved
residues; domain III alone is 204, which is the difference between a job
that fits on a T4 and one that does not.

What this produces
------------------
  A trimmed PDB carrying only domain III of chain A, with the Fab and all
  heteroatoms removed. Residue numbering is preserved, so every hotspot
  number derived in steps 02 and 03 remains valid.

  A BindCraft target JSON for each selected patch, with hotspots written in
  the numbering BindCraft will actually see.

  A record of what was cut, so the trim can be audited rather than trusted.

Numbering
---------
The trimmed file keeps 6ARU's own numbering, which is the mature protein and
runs 24 below UniProt. Domain III at UniProt 335-538 is therefore PDB
311-514. Hotspots are emitted in PDB numbering; the JSON records both.

Inputs:  01_Target/structure/6ARU.pdb
         02_Analysis/02_hotspots_bindcraft.json
Outputs: 03_Design/target/EGFR_domainIII.pdb
         03_Design/target/bindcraft_patch<N>.json
         03_Design/target/trim_report.md

Usage:
    python scripts/04_prepare_target.py
    python scripts/04_prepare_target.py --range 335 538 --patches 1 2 5
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "03_Design" / "target"
STRUCTURE = ROOT / "01_Target/structure/6ARU.pdb"
HOTSPOTS = ROOT / "02_Analysis/02_hotspots_bindcraft.json"
GEOMETRY = ROOT / "02_Analysis/02b_patch_geometry.csv"

EGFR_CHAIN = "A"
UNIPROT_MINUS_PDB = 24
DEFAULT_RANGE = (335, 538)          # domain III, UniProt precursor numbering

# Binder lengths to sample. The competition allows 10-250 residues; 55-100
# is the minibinder band, which is where AlphaFold2-based pipelines have
# their best documented success rate and where a 600-900 A^2 interface fits.
BINDER_LENGTHS = [55, 65, 75, 85, 100]
DESIGNS_PER_PATCH = 100


class DomainSelect(Select):
    """Keep one chain's standard residues within a PDB-numbered range."""

    def __init__(self, chain_id: str, first: int, last: int):
        self.chain_id = chain_id
        self.first = first
        self.last = last

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return (residue.id[0] == " "
                and residue.get_resname() in protein_letters_3to1
                and self.first <= residue.id[1] <= self.last)

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--range", nargs=2, type=int, default=DEFAULT_RANGE,
                        metavar=("START", "END"),
                        help="UniProt range to keep "
                             f"(default {DEFAULT_RANGE[0]} {DEFAULT_RANGE[1]})")
    parser.add_argument("--patches", nargs="*", type=int, default=[1, 2, 5],
                        help="patch ranks to write configurations for")
    parser.add_argument("--keep-detached", action="store_true",
                        help="keep hotspots that step 02b found detached "
                             "from the patch surface")
    args = parser.parse_args()

    uniprot_first, uniprot_last = args.range
    pdb_first = uniprot_first - UNIPROT_MINUS_PDB
    pdb_last = uniprot_last - UNIPROT_MINUS_PDB

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    model = PDBParser(QUIET=True).get_structure("6aru", str(STRUCTURE))[0]

    before = sum(1 for r in model[EGFR_CHAIN]
                 if r.id[0] == " " and r.get_resname() in protein_letters_3to1)
    other_chains = [c.id for c in model if c.id != EGFR_CHAIN]

    trimmed_path = OUT_DIR / "EGFR_domainIII.pdb"
    io = PDBIO()
    io.set_structure(model)
    io.save(str(trimmed_path),
            select=DomainSelect(EGFR_CHAIN, pdb_first, pdb_last))

    # Read back what was actually written, rather than assuming.
    check = PDBParser(QUIET=True).get_structure("trim", str(trimmed_path))[0]
    kept = [r for r in check[EGFR_CHAIN] if r.id[0] == " "]
    numbers = [r.id[1] for r in kept]
    atoms = sum(len(list(r)) for r in kept)

    print(f"Trimmed target")
    print(f"  UniProt range      : {uniprot_first}-{uniprot_last}")
    print(f"  PDB range          : {pdb_first}-{pdb_last}")
    print(f"  residues kept      : {len(kept)} of {before}")
    print(f"  actual PDB span    : {min(numbers)}-{max(numbers)}")
    print(f"  heavy atoms        : {atoms}")
    print(f"  chains removed     : {', '.join(other_chains)}")
    print(f"  file               : {trimmed_path.relative_to(ROOT)} "
          f"({trimmed_path.stat().st_size / 1024:.0f} KB)")

    # Gaps matter: an unresolved loop inside the trim leaves a break that
    # AlphaFold2 will see, and it is better to know now than to wonder later.
    gaps = [(a, b) for a, b in zip(numbers, numbers[1:]) if b - a > 1]
    if gaps:
        print(f"  chain breaks       : {len(gaps)}")
        for a, b in gaps:
            print(f"      {a} -> {b}  ({b - a - 1} residues missing, "
                  f"UniProt {a + UNIPROT_MINUS_PDB + 1}-"
                  f"{b + UNIPROT_MINUS_PDB - 1})")
    else:
        print("  chain breaks       : none")

    # -- BindCraft configurations ------------------------------------------
    hotspot_data = json.loads(HOTSPOTS.read_text())
    patches = {p["rank"]: p for p in hotspot_data["patches"]}

    # Detached residues identified by step 02b, translated to PDB numbering.
    geometry_cores: dict[int, dict] = {}
    if GEOMETRY.exists():
        import csv as _csv
        with GEOMETRY.open() as handle:
            for row in _csv.DictReader(handle):
                detached = [int(v) - UNIPROT_MINUS_PDB
                            for v in row["detached_residues"].split()]
                geometry_cores[int(row["rank"])] = {"detached": detached}

    written = []
    print("\nBindCraft target configurations")
    for rank in args.patches:
        patch = patches.get(rank)
        if patch is None:
            print(f"  patch {rank}: not in the hotspot file, skipped")
            continue

        # Step 02b found that some patches carry residues which are close in
        # a distance matrix but detached from the patch surface. Handing one
        # to BindCraft as a hotspot points the design at a residue on the
        # wrong face, so the connected core is used where 02b computed one.
        candidates = patch["residues_pdb"]
        detached = geometry_cores.get(rank, {}).get("detached", [])
        if detached and not args.keep_detached:
            candidates = [n for n in candidates if n not in detached]
            names = ", ".join(f"{n} (UniProt {n + UNIPROT_MINUS_PDB})"
                              for n in sorted(detached))
            print(f"  patch {rank}: dropping {len(detached)} detached "
                  f"residue(s) per step 02b: {names}")

        inside = [n for n in candidates if pdb_first <= n <= pdb_last]
        outside = [n for n in candidates if n not in inside]
        if outside:
            print(f"  patch {rank}: {len(outside)} hotspot(s) fall outside "
                  f"the trim and were dropped: {outside}")
        if not inside:
            print(f"  patch {rank}: no hotspots survive the trim, skipped")
            continue

        config = {
            "design_path": f"./designs_patch{rank}/",
            "binder_name": f"EGFRd3_p{rank}",
            "starting_pdb": f"./EGFR_domainIII.pdb",
            "chains": EGFR_CHAIN,
            "target_hotspot_residues": ",".join(map(str, sorted(inside))),
            "lengths": [min(BINDER_LENGTHS), max(BINDER_LENGTHS)],
            "number_of_final_designs": DESIGNS_PER_PATCH,
        }
        path = OUT_DIR / f"bindcraft_patch{rank}.json"
        path.write_text(json.dumps(config, indent=4) + "\n")
        written.append((rank, config, inside, patch))

        print(f"  patch {rank}: {len(inside)} hotspots -> "
              f"{path.relative_to(ROOT)}")
        print(f"      {config['target_hotspot_residues']}")

    # -- Report -------------------------------------------------------------
    report = OUT_DIR / "trim_report.md"
    with report.open("w") as handle:
        w = handle.write
        w("# Step 04 - Trimmed target and BindCraft configuration\n\n")
        w(f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
          f"by `scripts/04_prepare_target.py`.\n\n")
        w("## What was cut\n\n")
        w(f"- Source: `{STRUCTURE.name}`, chain {EGFR_CHAIN}, "
          f"{before} resolved residues\n")
        w(f"- Kept: domain III, UniProt {uniprot_first}–{uniprot_last} "
          f"= PDB {pdb_first}–{pdb_last}\n")
        w(f"- Result: {len(kept)} residues, {atoms} heavy atoms\n")
        w(f"- Removed: chains {', '.join(other_chains)} (the cetuximab Fab), "
          f"all heteroatoms including glycans and water, hydrogens, and "
          f"alternate conformations beyond the first\n")
        if gaps:
            w(f"- Chain breaks inside the trim: {len(gaps)}\n")
            for a, b in gaps:
                w(f"  - PDB {a}→{b}, {b - a - 1} residues unresolved\n")
        else:
            w("- No chain breaks inside the trim\n")

        w("\nTrimming is BindCraft's own recommendation: AlphaFold2 "
          "backpropagation scales badly with target size, and GPU memory is "
          "the binding constraint on the free tiers available here. Going "
          f"from {before} to {len(kept)} residues is the difference between "
          "a job that fits on a T4 and one that does not.\n")

        w("\n## Numbering\n\n")
        w("The trimmed file preserves 6ARU's numbering, which is the mature "
          "protein and runs 24 below UniProt. Every hotspot number derived "
          "in steps 02 and 03 therefore remains valid without translation. "
          "Hotspots below are in PDB numbering, which is what BindCraft "
          "reads.\n")

        w("\n## Configurations\n\n")
        for rank, config, inside, patch in written:
            w(f"### Patch {rank}\n\n")
            w(f"- File: `bindcraft_patch{rank}.json`\n")
            w(f"- Hotspots, PDB numbering: "
              f"`{config['target_hotspot_residues']}`\n")
            w(f"- Same residues, UniProt: "
              f"{', '.join(str(n + UNIPROT_MINUS_PDB) for n in sorted(inside))}\n")
            w(f"- Binder lengths: {config['lengths'][0]}–{config['lengths'][1]} aa\n")
            w(f"- Target accepted designs: {config['number_of_final_designs']}\n\n")

        w("## Running\n\n")
        w("Upload `EGFR_domainIII.pdb` and the JSON files to the GPU "
          "environment, then for each patch:\n\n")
        w("```bash\n")
        w("python -u ./bindcraft.py \\\n")
        w("    --settings './bindcraft_patch1.json' \\\n")
        w("    --filters './settings_filters/default_filters.json' \\\n")
        w("    --advanced './settings_advanced/default_4stage_multimer.json'\n")
        w("```\n\n")
        w("Run a handful of trajectories first and time them before "
          "committing the batch, since the per-trajectory cost sets the "
          "whole budget and it is cheaper to learn it early.\n")

    print(f"\n  -> {report.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 13 - Write Boltz-2 inputs for the variant complexes.

Steps 11 and 12 together leave one question open. AlphaFold2 declines
contacts that geometry allows: 12 of 16 engineered histidines can reach
their handle by an unclashed rotamer, and the prediction places 11 of them
5 to 11 A away. Every electrostatic number downstream rests on those
coordinates.

An independent predictor is the way to test that. Boltz-2 is a different
architecture trained differently, so if it builds the contacts the
placement was sound and AlphaFold2 is the problem; if it does not, the
behaviour is general and the finding is stronger. Either answer is worth
having, which is the mark of a test worth running.

This step only writes the inputs. The prediction needs a GPU and runs
elsewhere; `--notebook` prints the cells for it.

MSA
---
The target gets a real alignment from the MSA server. The binder gets none:
these sequences were hallucinated and have no homologues, so an alignment
for them would be noise. That matches how they were generated.

Usage
-----
    python scripts/13_prepare_boltz.py
    python scripts/13_prepare_boltz.py --notebook
"""

from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FASTA = ROOT / "03_Design" / "validation" / "07_variants_for_af2.fasta"
OUT_DIR = ROOT / "03_Design" / "validation" / "boltz_input"

NOTEBOOK = '''
# ---------------------------------------------------------------- cell 1
# Kaggle: Settings -> Accelerator = GPU T4 x2, Internet = On
!pip install -q boltz
!boltz --help | head -5


# ---------------------------------------------------------------- cell 2
# Upload boltz_input/ as a Kaggle dataset first, then point at it.
import os, glob
IN = "/kaggle/input/egfr-boltz-input"          # adjust to the dataset slug
print(len(glob.glob(f"{IN}/*.yaml")), "inputs")


# ---------------------------------------------------------------- cell 3
# ~19 complexes of 272 residues. Expect roughly 2-4 minutes each on a T4.
!boltz predict {IN} \\
    --out_dir /kaggle/working/boltz_out \\
    --use_msa_server \\
    --output_format pdb \\
    --diffusion_samples 1


# ---------------------------------------------------------------- cell 4
# Keep only the structures and the confidence files. Do not zip the whole
# working directory: the model weights are in there.
!cd /kaggle/working/boltz_out && \\
  find . -name "*.pdb" -o -name "confidence*.json" | \\
  zip -q /kaggle/working/boltz_results.zip -@
!ls -lh /kaggle/working/boltz_results.zip
'''


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", type=Path, default=FASTA)
    parser.add_argument("--out", type=Path, default=OUT_DIR)
    parser.add_argument("--notebook", action="store_true",
                        help="print the Kaggle cells and exit")
    args = parser.parse_args()

    if args.notebook:
        print(NOTEBOOK)
        return 0

    if not args.fasta.is_file():
        parser.error(f"run step 07 --prepare first: {args.fasta} is missing")

    args.out.mkdir(parents=True, exist_ok=True)
    for stale in args.out.glob("*.yaml"):
        stale.unlink()

    written = 0
    name = None
    for line in args.fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            continue
        if not line.strip() or name is None:
            continue

        target, _, binder = line.strip().partition(":")
        if not binder:
            print(f"  {name}: no chain separator, skipped")
            continue

        # Chain ids match the AlphaFold2 predictions: A is the 204-residue
        # target, B the binder. Keeping them identical means steps 11 and 12
        # run on Boltz output without a translation layer.
        (args.out / f"{name}.yaml").write_text(
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {target}\n"
            "  - protein:\n"
            "      id: B\n"
            f"      sequence: {binder}\n"
            "      msa: empty\n"
        )
        written += 1
        name = None

    print(f"{written} Boltz-2 inputs -> {args.out.relative_to(ROOT)}")
    print("\nNext:")
    print(f"  1. zip -r boltz_input.zip {args.out.relative_to(ROOT)}")
    print("  2. upload as a Kaggle dataset")
    print("  3. python scripts/13_prepare_boltz.py --notebook")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

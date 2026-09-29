#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate data-provenance text for a manuscript from MANIFEST.json.

The point is to never write provenance by hand. A Materials and Methods
paragraph typed from memory drifts silently the moment the data is
re-downloaded and a database release changes; nobody notices until peer
review. Generating the text from the manifest means the manuscript and the
files on disk cannot disagree.

Three output formats:

    --format paragraph   Materials and Methods prose (default)
    --format table       Markdown provenance table for supplementary material
    --format bibtex      BibTeX entries for the database records

Usage:
    python scripts/citation.py
    python scripts/citation.py --format table
    python scripts/citation.py --format bibtex > docs/data_sources.bib

Dependencies: Python 3.9+ standard library only.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "01_Target" / "MANIFEST.json"

MONTHS = ["January", "February", "March", "April", "May", "June", "July",
          "August", "September", "October", "November", "December"]


def load_manifest() -> dict:
    if not MANIFEST.exists():
        raise SystemExit(
            f"ERROR: {MANIFEST.relative_to(ROOT)} not found.\n"
            "Run: python scripts/00_download_data.py"
        )
    return json.loads(MANIFEST.read_text())


def format_date(iso_timestamp: str) -> str:
    """2026-09-29T15:05:08+00:00 -> '29 September 2026'."""
    moment = datetime.fromisoformat(iso_timestamp)
    return f"{moment.day} {MONTHS[moment.month - 1]} {moment.year}"


def organism_from_header(header: str) -> str:
    """Pull the species out of a UniProt FASTA header (OS= field)."""
    match = re.search(r"OS=(.+?)\s+OX=", header)
    return match.group(1) if match else "unknown organism"


def group_by_source(manifest: dict) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {}
    for record in manifest["files"]:
        groups.setdefault(record["source"], []).append(record)
    return groups


# ---------------------------------------------------------------------------
# Output formats
# ---------------------------------------------------------------------------

def render_paragraph(manifest: dict) -> str:
    """Materials and Methods prose."""
    files = manifest["files"]
    access_date = format_date(manifest["generated_utc"])

    sentences: list[str] = []

    # Sequences, grouped so that one sentence covers all UniProt records that
    # share a release.
    fasta_records = [r for r in files if r["file"].endswith(".fasta")]
    if fasta_records:
        releases = {r.get("uniprot_release") for r in fasta_records
                    if r.get("uniprot_release")}
        descriptors = []
        for record in fasta_records:
            organism = organism_from_header(record.get("fasta_header", ""))
            # The matching .json entry carries the entry version.
            sibling = next(
                (r for r in files
                 if r["accession"] == record["accession"]
                 and r["file"].endswith(".json")),
                None,
            )
            version = sibling.get("entry_version") if sibling else None
            version_text = f", entry version {version}" if version else ""
            descriptors.append(
                f"{organism} EGFR (UniProt {record['accession']}{version_text})"
            )

        joined = " and ".join(descriptors) if len(descriptors) == 2 \
            else ", ".join(descriptors[:-1]) + f", and {descriptors[-1]}" \
            if len(descriptors) > 2 else descriptors[0]
        release_text = f" release {sorted(releases)[0]}" if len(releases) == 1 else ""
        sentences.append(
            f"Reference sequences for {joined} were retrieved from "
            f"UniProtKB{release_text} on {access_date}."
        )

    # Structures.
    structures = {r["accession"] for r in files if r["source"] == "RCSB PDB"}
    for accession in sorted(structures):
        record = next(r for r in files
                      if r["accession"] == accession and r["file"].endswith(".pdb"))
        chains = record.get("atoms_per_chain", {})
        chain_text = ""
        if chains:
            chain_text = (f" The asymmetric unit contains chains "
                          f"{', '.join(sorted(chains))}.")
        sentences.append(
            f"The crystal structure of the human EGFR extracellular region in "
            f"complex with a cetuximab Fab mutant (PDB {accession}) was "
            f"retrieved from the RCSB Protein Data Bank on {access_date}."
            f"{chain_text}"
        )

    sentences.append(
        "SHA-256 checksums, download timestamps and source URLs for all input "
        "files are recorded in the accompanying repository "
        "(01_Target/MANIFEST.json) and can be re-verified with "
        "`python scripts/00_download_data.py --verify`."
    )

    return " ".join(sentences)


def render_table(manifest: dict) -> str:
    """Markdown provenance table for supplementary material."""
    lines = [
        "### Data provenance",
        "",
        f"Generated from `01_Target/MANIFEST.json` "
        f"({manifest['generated_utc']}).",
        "",
        "| File | Source | Accession | Version | Retrieved | SHA-256 |",
        "|------|--------|-----------|---------|-----------|---------|",
    ]

    for record in manifest["files"]:
        version = (record.get("entry_version")
                   or record.get("uniprot_release")
                   or "—")
        lines.append(
            f"| `{Path(record['file']).name}` "
            f"| {record['source']} "
            f"| {record['accession']} "
            f"| {version} "
            f"| {format_date(record['downloaded_utc'])} "
            f"| `{record['sha256'][:16]}…` |"
        )

    lines += [
        "",
        "Checksums are truncated to 16 characters for readability; full values "
        "are in the manifest.",
    ]
    return "\n".join(lines)


def render_bibtex(manifest: dict) -> str:
    """BibTeX entries for the database records."""
    entries: list[str] = []
    seen: set[str] = set()

    for record in manifest["files"]:
        accession = record["accession"]
        if accession in seen:
            continue
        seen.add(accession)

        year = datetime.fromisoformat(record["downloaded_utc"]).year
        access_date = format_date(record["downloaded_utc"])

        if record["source"] == "UniProtKB":
            release = record.get("uniprot_release", "")
            note = f"UniProtKB release {release}" if release else "UniProtKB"
            entries.append(
                f"@misc{{uniprot_{accession.lower()},\n"
                f"  title        = {{{{UniProtKB entry {accession}}}}},\n"
                f"  author       = {{{{The UniProt Consortium}}}},\n"
                f"  year         = {{{year}}},\n"
                f"  howpublished = {{\\url{{https://www.uniprot.org/uniprotkb/{accession}}}}},\n"
                f"  note         = {{{note}. Accessed {access_date}}}\n"
                f"}}"
            )
        else:
            entries.append(
                f"@misc{{pdb_{accession.lower()},\n"
                f"  title        = {{{{RCSB PDB entry {accession}}}}},\n"
                f"  year         = {{{year}}},\n"
                f"  howpublished = {{\\url{{https://www.rcsb.org/structure/{accession}}}}},\n"
                f"  note         = {{Accessed {access_date}}}\n"
                f"}}"
            )

    return "\n\n".join(entries)


# ---------------------------------------------------------------------------

RENDERERS = {
    "paragraph": render_paragraph,
    "table": render_table,
    "bibtex": render_bibtex,
}


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--format", choices=sorted(RENDERERS), default="paragraph",
                        help="output format (default: paragraph)")
    args = parser.parse_args()

    print(RENDERERS[args.format](load_manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

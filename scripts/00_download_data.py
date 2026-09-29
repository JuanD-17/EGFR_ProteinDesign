#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 00 - Reproducible download of EGFR target data.

Downloads the reference sequences and the crystal structure that serve as the
starting point for binder design, and records the exact provenance of every
file in a manifest (MANIFEST.json): URL, download timestamp, size, SHA-256
checksum and, where the server exposes it, the database release version.

This matters for reproducibility. UniProt and the RCSB PDB revise their
entries over time, so "we downloaded P00533" does not identify a file. A
checksum does.

Usage:
    python scripts/00_download_data.py            # download whatever is missing
    python scripts/00_download_data.py --verify   # only check against manifest
    python scripts/00_download_data.py --force    # re-download everything

Dependencies: Python 3.9+ standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Declarative source definitions.
# Adding a new dataset means adding one entry here and nothing else.
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "01_Target" / "MANIFEST.json"

SOURCES = [
    {
        "id": "egfr_human_fasta",
        "url": "https://rest.uniprot.org/uniprotkb/P00533.fasta",
        "path": "01_Target/human/P00533.fasta",
        "description": (
            "Human EGFR, canonical precursor sequence (isoform P00533-1, "
            "1210 aa, includes the 1-24 signal peptide). This is the "
            "challenge target."
        ),
        "source": "UniProtKB",
        "accession": "P00533",
    },
    {
        "id": "egfr_mouse_fasta",
        "url": "https://rest.uniprot.org/uniprotkb/Q01279.fasta",
        "path": "01_Target/mouse/Q01279.fasta",
        "description": (
            "Mouse (Mus musculus) EGFR precursor. Required for challenge "
            "objective 2, cross-species reactivity."
        ),
        "source": "UniProtKB",
        "accession": "Q01279",
    },
    {
        "id": "egfr_human_json",
        "url": "https://rest.uniprot.org/uniprotkb/P00533.json",
        "path": "01_Target/human/P00533.json",
        "description": (
            "Full UniProt entry for human EGFR in JSON. Carries the domain "
            "annotations (I-IV), glycosylation sites, disulfide bonds and "
            "sequence variants. Domain III boundaries are read from here "
            "rather than assumed."
        ),
        "source": "UniProtKB",
        "accession": "P00533",
    },
    {
        "id": "egfr_mouse_json",
        "url": "https://rest.uniprot.org/uniprotkb/Q01279.json",
        "path": "01_Target/mouse/Q01279.json",
        "description": "Full UniProt entry for mouse EGFR in JSON.",
        "source": "UniProtKB",
        "accession": "Q01279",
    },
    {
        "id": "structure_6aru_pdb",
        "url": "https://files.rcsb.org/download/6ARU.pdb",
        "path": "01_Target/structure/6ARU.pdb",
        "description": (
            "Human EGFR extracellular region in complex with a cetuximab Fab "
            "mutant. Chain A = EGFR, chain B = Fab light chain, chain C = Fab "
            "heavy chain. This is the reference structure named on the "
            "challenge page."
        ),
        "source": "RCSB PDB",
        "accession": "6ARU",
    },
    {
        "id": "structure_6aru_cif",
        "url": "https://files.rcsb.org/download/6ARU.cif",
        "path": "01_Target/structure/6ARU.cif",
        "description": (
            "The same 6ARU structure in mmCIF format. Kept because the legacy "
            "PDB format can truncate residue numbering and chain identifiers; "
            "mmCIF is the authoritative version."
        ),
        "source": "RCSB PDB",
        "accession": "6ARU",
    },
]

USER_AGENT = "EGFR-ProteinDesign/1.0 (Anthropic x Adaptyv 2026 competition)"
MAX_ATTEMPTS = 3

# HTTP headers worth preserving: they pin the database release.
HEADERS_OF_INTEREST = {
    "x-uniprot-release",
    "x-uniprot-release-date",
    "last-modified",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sha256(path: Path) -> str:
    """SHA-256 checksum of a file, read in 1 MiB blocks."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url: str) -> tuple[bytes, dict[str, str]]:
    """Download a URL, returning its body and the headers worth recording."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None

    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                body = response.read()
                headers = {
                    key: value
                    for key, value in response.headers.items()
                    if key.lower() in HEADERS_OF_INTEREST
                }
                return body, headers
        except (urllib.error.URLError, TimeoutError) as error:
            last_error = error
            if attempt < MAX_ATTEMPTS:
                print(f"      attempt {attempt} failed ({error}); retrying...")

    raise RuntimeError(f"could not download {url}: {last_error}")


def validate(path: Path) -> None:
    """Sanity-check a downloaded file.

    A failing server can return an HTML error page with status 200, which
    would otherwise be written to disk and silently corrupt the analysis.
    """
    if path.stat().st_size == 0:
        raise ValueError(f"{path.name} is empty")

    with path.open("rb") as fh:
        head = fh.read(512).lstrip()

    if head.startswith(b"<"):
        raise ValueError(f"{path.name} looks like HTML, not a data file")

    suffix = path.suffix
    if suffix == ".fasta" and not head.startswith(b">"):
        raise ValueError(f"{path.name} does not start with '>' as FASTA requires")
    if suffix == ".pdb" and not (head.startswith(b"HEADER") or b"ATOM" in head):
        raise ValueError(f"{path.name} does not look like a PDB file")
    if suffix == ".cif" and not head.startswith(b"data_"):
        raise ValueError(f"{path.name} does not look like an mmCIF file")
    if suffix == ".json":
        try:
            json.loads(path.read_text())
        except json.JSONDecodeError as error:
            raise ValueError(f"{path.name} is not valid JSON: {error}") from None


def summarize(path: Path) -> dict:
    """Extract a few human-readable facts so the manifest can be audited at a
    glance, without opening the data files themselves."""
    summary: dict = {}

    if path.suffix == ".fasta":
        lines = path.read_text().splitlines()
        sequence = "".join(line.strip() for line in lines[1:]
                           if not line.startswith(">"))
        summary["fasta_header"] = lines[0] if lines else ""
        summary["length_aa"] = len(sequence)

    elif path.suffix == ".pdb":
        atoms_per_chain: dict[str, int] = {}
        title_parts: list[str] = []
        for line in path.read_text(errors="ignore").splitlines():
            if line.startswith("TITLE"):
                title_parts.append(line[10:].strip())
            elif line.startswith("ATOM"):
                chain = line[21]
                atoms_per_chain[chain] = atoms_per_chain.get(chain, 0) + 1
        summary["title"] = " ".join(title_parts)
        summary["atoms_per_chain"] = atoms_per_chain

    elif path.suffix == ".json":
        entry = json.loads(path.read_text())
        audit = entry.get("entryAudit", {})
        summary["entry_version"] = audit.get("entryVersion")
        summary["last_sequence_update"] = audit.get("lastSequenceUpdateDate")

    return summary


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def download_all(force: bool) -> int:
    records = []
    print(f"Project root: {ROOT}\n")

    for source in SOURCES:
        destination = ROOT / source["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        headers: dict[str, str] = {}

        if destination.exists() and not force:
            print(f"  [cached] {source['path']}")
        else:
            print(f"  [http]   {source['path']}")
            body, headers = fetch(source["url"])
            destination.write_bytes(body)

        validate(destination)

        record = {
            "id": source["id"],
            "file": source["path"],
            "url": source["url"],
            "source": source["source"],
            "accession": source["accession"],
            "description": source["description"],
            "bytes": destination.stat().st_size,
            "sha256": sha256(destination),
            "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        if headers.get("X-UniProt-Release"):
            record["uniprot_release"] = headers["X-UniProt-Release"]
        if headers.get("Last-Modified"):
            record["server_last_modified"] = headers["Last-Modified"]
        record.update(summarize(destination))
        records.append(record)

        print(f"           {record['bytes']:>9,} bytes  "
              f"sha256:{record['sha256'][:16]}...")

    manifest = {
        "project": "Conditional EGFR binder design",
        "competition": "Anthropic x Adaptyv 2026, Challenge 01 (EGFR)",
        "generated_by": "scripts/00_download_data.py",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "files": records,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"\n  {len(records)} files downloaded and validated.")
    print(f"  Manifest -> {MANIFEST.relative_to(ROOT)}")
    print("\n  To check integrity later:")
    print("      python scripts/00_download_data.py --verify")
    return 0


def verify() -> int:
    """Recompute checksums and compare them against the manifest."""
    if not MANIFEST.exists():
        print(f"ERROR: {MANIFEST} not found. Run the script without --verify first.")
        return 1

    manifest = json.loads(MANIFEST.read_text())
    print(f"Manifest generated {manifest['generated_utc']}\n")

    failures = 0
    for record in manifest["files"]:
        path = ROOT / record["file"]
        if not path.exists():
            print(f"  MISSING    {record['file']}")
            failures += 1
            continue

        observed = sha256(path)
        if observed != record["sha256"]:
            print(f"  ALTERED    {record['file']}")
            print(f"             expected {record['sha256']}")
            print(f"             observed {observed}")
            failures += 1
        else:
            print(f"  OK         {record['file']}")

    if failures:
        print(f"\n  {failures} problem(s). Run with --force to restore.")
        return 1

    print(f"\n  All {len(manifest['files'])} files match the manifest.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--force", action="store_true",
                        help="re-download even if the file already exists")
    parser.add_argument("--verify", action="store_true",
                        help="do not download; check files against MANIFEST.json")
    args = parser.parse_args()

    return verify() if args.verify else download_all(args.force)


if __name__ == "__main__":
    raise SystemExit(main())

"""Download the documented source snapshots, verifying known checksums."""

import hashlib
import json
import urllib.request

from pipelines.build_database import ROOT, SOURCES


# Download missing source snapshots and verify each file against its recorded checksum.
def fetch():
    manifest = json.loads((ROOT / "data/provenance.json").read_text())
    for key, (_, url, _, name) in SOURCES.items():
        dest = ROOT / "data/raw" / name
        expected = manifest[key]["sha256"]
        if dest.exists():
            with dest.open("rb") as f:
                existing = hashlib.file_digest(f, "sha256").hexdigest()
            if existing == expected:
                print("Verified:", name)
                continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        part = dest.with_suffix(".part")
        with (
            urllib.request.urlopen(url, timeout=90) as response,
            part.open("wb") as out,
        ):
            while block := response.read(1024 * 1024):
                out.write(block)
        with part.open("rb") as f:
            actual = hashlib.file_digest(f, "sha256").hexdigest()
        if actual != expected:
            raise SystemExit(
                f"Snapshot changed for {name}. Inspect the source before importing; downloaded file retained as {part.name}."
            )
        part.replace(dest)
        print("Downloaded and verified:", name)


if __name__ == "__main__":
    fetch()

"""Restore verified deployment data after an ephemeral host restarts."""

import gzip
import hashlib
import json
import os
import sqlite3
import tempfile
from contextlib import closing
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKED = ROOT / "deployment/data"


# Compute a bounded-memory digest for compressed and restored snapshot verification.
def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


# Restore only expected SQLite snapshots and refuse to overwrite the research databases.
def restore(destination):
    destination = Path(destination).resolve()
    if destination == (ROOT / "data").resolve():
        raise ValueError(
            "Use a separate hosting directory, not the research data directory."
        )
    destination.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((PACKED / "manifest.json").read_text())
    for record in manifest["files"]:
        name = record["filename"]
        if Path(name).name != name or not name.endswith(".sqlite"):
            raise ValueError("Invalid snapshot filename.")
        metadata = ROOT / "data" / record["metadata_file"]
        if name != "rice.sqlite" and not metadata.exists():
            continue
        if digest(metadata) != record["metadata_sha256"]:
            raise ValueError(
                "Study metadata changed; regenerate the deployment snapshot."
            )
        packed = PACKED / (name + ".gz")
        if digest(packed) != record["compressed_sha256"]:
            raise ValueError("Compressed deployment snapshot checksum failed.")
        target = destination / name
        if target.exists() and digest(target) == record["sha256"]:
            continue
        with tempfile.NamedTemporaryFile(dir=destination, delete=False) as output:
            temporary = Path(output.name)
            try:
                with gzip.open(packed, "rb") as source:
                    total = 0
                    while block := source.read(1024 * 1024):
                        total += len(block)
                        if total > record["bytes"]:
                            raise ValueError(
                                "Deployment snapshot exceeds recorded size."
                            )
                        output.write(block)
                output.close()
                if total != record["bytes"] or digest(temporary) != record["sha256"]:
                    raise ValueError("Restored deployment snapshot checksum failed.")
                with closing(
                    sqlite3.connect(f"file:{temporary}?mode=ro", uri=True)
                ) as db:
                    if db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                        raise ValueError("Restored database integrity failed.")
                temporary.chmod(0o444)
                temporary.replace(target)
            finally:
                temporary.unlink(missing_ok=True)
        print(f"Restored {name}", flush=True)


if __name__ == "__main__":
    configured = os.environ.get("RICELENS_DATA_DIR")
    if not configured:
        raise SystemExit("Set RICELENS_DATA_DIR to a separate hosting data directory.")
    restore(configured)

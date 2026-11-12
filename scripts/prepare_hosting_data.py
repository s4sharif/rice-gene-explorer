"""Package compact read-only databases for the free reviewer deployment."""

import gzip
import hashlib
import json
import shutil
import sqlite3
import tempfile
from contextlib import closing
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESTINATION = ROOT / "deployment/data"
FILES = {
    "rice.sqlite": "provenance.json",
    "expression.sqlite": "expression_manifest.json",
    "expression_wang.sqlite": "expression_wang_manifest.json",
    "expression_li.sqlite": "expression_li_manifest.json",
    "expression_li_author.sqlite": "expression_li_author_manifest.json",
}


# Hash data incrementally so snapshot creation and restoration have bounded memory use.
def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


# Preserve all serving tables while omitting only offline per-sample expression summaries.
def prepare():
    DESTINATION.mkdir(parents=True, exist_ok=True)
    records = []
    with tempfile.TemporaryDirectory(prefix="ricelens-serving-") as directory:
        for name, metadata_name in FILES.items():
            source = ROOT / "data/databases" / name
            target = Path(directory) / name
            with closing(
                sqlite3.connect(f"file:{source}?mode=ro", uri=True)
            ) as original:
                with closing(sqlite3.connect(target)) as serving:
                    original.backup(serving)
                    excluded = [] if name == "rice.sqlite" else ["sample_summary"]
                    for table in excluded:
                        serving.execute(f"DROP TABLE {table}")
                    serving.commit()
                    serving.execute("VACUUM")
                    tables = [
                        r[0]
                        for r in serving.execute(
                            "SELECT name FROM sqlite_master WHERE type='table'"
                        )
                    ]
                    counts = {}
                    for table in tables:
                        counts[table] = serving.execute(
                            f'SELECT COUNT(*) FROM "{table}"'
                        ).fetchone()[0]
                        if (
                            counts[table]
                            != original.execute(
                                f'SELECT COUNT(*) FROM "{table}"'
                            ).fetchone()[0]
                        ):
                            raise ValueError("Serving snapshot lost rows.")
                    if serving.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                        raise ValueError("Serving snapshot integrity failed.")
            packed = DESTINATION / (name + ".gz")
            with target.open("rb") as incoming, packed.open("wb") as raw:
                with gzip.GzipFile(
                    filename="", fileobj=raw, mode="wb", mtime=0, compresslevel=6
                ) as outgoing:
                    shutil.copyfileobj(incoming, outgoing, 1024 * 1024)
            records.append(
                {
                    "filename": name,
                    "bytes": target.stat().st_size,
                    "sha256": digest(target),
                    "compressed_bytes": packed.stat().st_size,
                    "compressed_sha256": digest(packed),
                    "source_sha256": digest(source),
                    "metadata_file": metadata_name,
                    "metadata_sha256": digest(ROOT / "data" / metadata_name),
                    "retained_tables": counts,
                    "excluded_tables": excluded,
                }
            )
            print(
                f"Packaged {name}: {packed.stat().st_size:,} compressed bytes",
                flush=True,
            )
    (DESTINATION / "manifest.json").write_text(
        json.dumps(
            {
                "purpose": "Read-only reviewer snapshot; raw matrices and offline sample summaries excluded",
                "files": records,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    prepare()

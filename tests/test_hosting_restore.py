"""Verify snapshot integrity failures cannot silently replace deployed research data."""

import gzip
import json
import sqlite3
import tempfile
import unittest
from contextlib import closing
from pathlib import Path
from unittest.mock import patch

from scripts import restore_hosting_data as restore


class RestoreTests(unittest.TestCase):
    # Create a tiny valid snapshot and its independently checked provenance fixture.
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.packed = self.root / "packed"
        self.packed.mkdir()
        (self.root / "data").mkdir()
        metadata = self.root / "data/provenance.json"
        metadata.write_text('{"fixture":true}')
        database = self.root / "original.sqlite"
        with closing(sqlite3.connect(database)) as db:
            db.execute("CREATE TABLE evidence(value INTEGER)")
            db.execute("INSERT INTO evidence VALUES(7)")
            db.commit()
        compressed = self.packed / "rice.sqlite.gz"
        compressed.write_bytes(gzip.compress(database.read_bytes(), mtime=0))
        self.manifest = {
            "files": [
                {
                    "filename": "rice.sqlite",
                    "metadata_file": "provenance.json",
                    "metadata_sha256": restore.digest(metadata),
                    "bytes": database.stat().st_size,
                    "sha256": restore.digest(database),
                    "compressed_sha256": restore.digest(compressed),
                }
            ]
        }
        self.write_manifest()
        self.root_patch = patch.object(restore, "ROOT", self.root)
        self.packed_patch = patch.object(restore, "PACKED", self.packed)
        self.root_patch.start()
        self.packed_patch.start()
        self.addCleanup(self.root_patch.stop)
        self.addCleanup(self.packed_patch.stop)

    # Write a fixture manifest without relying on the snapshot generator under test.
    def write_manifest(self):
        (self.packed / "manifest.json").write_text(json.dumps(self.manifest))

    # Restore a correct database read-only and keep repeated startup restoration safe.
    def test_valid_restore_and_restart(self):
        destination = self.root / "served"
        restore.restore(destination)
        restore.restore(destination)
        target = destination / "rice.sqlite"
        with closing(sqlite3.connect(f"file:{target}?mode=ro", uri=True)) as db:
            self.assertEqual(db.execute("SELECT value FROM evidence").fetchone()[0], 7)
        self.assertEqual(target.stat().st_mode & 0o222, 0)

    # Reject corrupt compressed data before touching a previously restored database.
    def test_corruption_does_not_replace_data(self):
        destination = self.root / "served"
        restore.restore(destination)
        target = destination / "rice.sqlite"
        before = restore.digest(target)
        (self.packed / "rice.sqlite.gz").write_bytes(b"invalid snapshot")
        with self.assertRaises(ValueError):
            restore.restore(destination)
        self.assertEqual(restore.digest(target), before)

    # Refuse stale provenance, traversal filenames and the original research data path.
    def test_invalid_snapshot_targets(self):
        with self.assertRaises(ValueError):
            restore.restore(self.root / "data")
        self.manifest["files"][0]["filename"] = "../outside.sqlite"
        self.write_manifest()
        with self.assertRaises(ValueError):
            restore.restore(self.root / "served")
        self.manifest["files"][0]["filename"] = "rice.sqlite"
        self.manifest["files"][0]["metadata_sha256"] = "wrong"
        self.write_manifest()
        with self.assertRaises(ValueError):
            restore.restore(self.root / "served")


if __name__ == "__main__":
    unittest.main()

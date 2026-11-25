"""Check sparse source extraction, context separation and expression calculations."""

import gzip
import json
import math
import struct
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from backend import expression
from pipelines.prepare_wang_atlas import (
    ORGANS,
    group_cells,
    open_counts,
    summarize_blocks,
)
from scripts.stream_wang_atlas import AtlasReader


class AtlasImportTests(unittest.TestCase):
    # Extract original counts rather than the differently scaled data slot from an independent R fixture.
    def test_streamed_rna_counts_and_metadata(self):
        fixture = Path(__file__).parent / "fixtures/wang_atlas_counts.rds"
        with tempfile.TemporaryDirectory() as directory:
            with gzip.open(fixture, "rb") as source:
                report = AtlasReader(source, directory).export()
            path = Path(directory)
            self.assertEqual(report["shape"], [3, 2])
            self.assertEqual(report["assay_names"], ["RNA", "ATAC"])
            self.assertEqual(
                struct.unpack(">3d", (path / "x.bin").read_bytes()), (2, 3, 5)
            )
            self.assertEqual(
                struct.unpack(">3i", (path / "p.bin").read_bytes()), (0, 1, 3)
            )
            self.assertIn("epidermis", (path / "cell_metadata.tsv").read_text())
            spec, arrays = open_counts(path)
            detected, total, sizes = summarize_blocks(
                spec, arrays, [0, 1], 2, block_cells=1
            )
            self.assertEqual(sizes.tolist(), [2, 8])
            self.assertEqual(detected.tolist(), [[0, 1], [1, 0], [0, 1]])
            self.assertAlmostEqual(total[1, 0], math.log1p(10000))
            self.assertAlmostEqual(total[0, 1], math.log1p(3750))
            self.assertAlmostEqual(total[2, 1], math.log1p(6250))
            self.assertEqual(total[0, 0], 0)

    # Reject an ATAC assay even when its sparse matrix has the same dimensions as RNA.
    def test_atac_is_not_substituted_for_rna(self):
        fixture = Path(__file__).parent / "fixtures/wang_atlas_counts.rds"
        with tempfile.TemporaryDirectory() as directory:
            with (
                gzip.open(fixture, "rb") as source,
                self.assertRaisesRegex(ValueError, "not verified RNA"),
            ):
                AtlasReader(source, directory, assay_index=1).export()

    # Reject a truncated sparse vector before any expression values are computed.
    def test_truncated_counts_fail(self):
        fixture = Path(__file__).parent / "fixtures/wang_atlas_counts.rds"
        with tempfile.TemporaryDirectory() as directory:
            with gzip.open(fixture, "rb") as source:
                AtlasReader(source, directory).export()
            (Path(directory) / "x.bin").write_bytes(b"bad")
            with self.assertRaisesRegex(ValueError, "Truncated"):
                open_counts(Path(directory))

    # Keep the same source cell-type label distinct in every organ and retain unknown labels.
    def test_organ_contexts_and_missing_annotations(self):
        rows = [
            {
                "barcode": f"{tissue}-{index}",
                "sample": f"{tissue}_{index}",
                "tissue": tissue,
                "cluster_names": "unknown",
            }
            for tissue in ORGANS
            for index in [1, 2]
        ]
        barcodes = [row["barcode"] for row in reversed(rows)]
        aligned, counts = group_cells(rows, barcodes, expected_cells=16)
        self.assertEqual(len(counts), 16)
        self.assertEqual(len({key[1] for key in aligned}), 8)
        self.assertTrue(all(key[3] == "unknown" for key in aligned))
        rows[0]["cluster_names"] = ""
        with self.assertRaisesRegex(ValueError, "Missing"):
            group_cells(rows, barcodes, expected_cells=16)

    # Require matching validation provenance before enabling the new dataset.
    def test_stale_validation_disables_atlas(self):
        dataset = "wang_2025_eight_organ_atlas"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "data").mkdir()
            manifest = root / "data/expression_wang_atlas_manifest.json"
            validation = root / "data/expression_wang_atlas_validation.json"
            manifest.write_text(json.dumps({"source_sha256": "source-one"}))
            with patch.object(
                expression, "__file__", str(root / "backend/expression.py")
            ):
                self.assertFalse(expression.dataset_enabled(dataset))
                validation.write_text(
                    json.dumps(
                        {
                            "all_deposited_cells_and_features": "pass",
                            "source_sha256": "different-source",
                        }
                    )
                )
                self.assertFalse(expression.dataset_enabled(dataset))
                validation.write_text(
                    json.dumps(
                        {
                            "all_deposited_cells_and_features": "pass",
                            "source_sha256": "source-one",
                        }
                    )
                )
                self.assertTrue(expression.dataset_enabled(dataset))

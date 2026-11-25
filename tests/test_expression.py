"""Check expression semantics independently of the large source dataset."""

import gzip
import json
import math
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from backend import expression
from scripts.extract_rds_metadata import SlotReader, attributes, column_values


class ExpressionTests(unittest.TestCase):
    # Compare a bounded-memory metadata extraction with an independently serialized R fixture.
    def test_metadata_slot_reader(self):
        source = Path(__file__).parent / "fixtures/metadata_slot.rds"
        with gzip.open(source, "rb") as handle:
            metadata = SlotReader(handle).extract()
        attrs = attributes(metadata["attrs"])
        self.assertEqual(
            column_values(attrs["names"]), ["sample", "CellAnnotation", "nCount_RNA"]
        )
        self.assertEqual(column_values(attrs["row.names"]), ["cell-A", "cell-B"])
        self.assertEqual(
            column_values(metadata["value"][1]), ["Mesophyll", "Epidermis"]
        )
        self.assertEqual(column_values(metadata["value"][2]), [3, 7])

    # Keep an organ distinct from the cell's tissue-system family and preserve original labels.
    def test_organ_and_tissue_system_are_distinct(self):
        annotation = expression.context_annotation(
            {"organ": "Developing leaf"}, "Xylem"
        )
        self.assertEqual(
            annotation, {"organ": "Developing leaf", "source_cell_type": "Xylem"}
        )
        self.assertEqual(
            expression.context_annotation({}, "Xylem")["organ"], "Not specified"
        )
        atlas = {
            "organ": "Eight organs",
            "context_annotations": {
                "Root | Cortex": {"organ": "Crown root", "source_cell_type": "Cortex"}
            },
        }
        self.assertEqual(
            expression.context_annotation(atlas, "Root | Cortex")["organ"], "Crown root"
        )

    # Confirm that unavailable datasets return an explicit state rather than fabricated zeros.
    def test_missing_dataset(self):
        with patch.object(expression, "EXPRESSION_DB", Path("/does-not-exist.sqlite")):
            self.assertFalse(
                expression.expression_evidence(["Os01g0100100"])["available"]
            )

    # Check detection percentages, tied ranks and missing features against a known count example.
    def test_summary_semantics(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "expression.sqlite"
            with sqlite3.connect(path) as db:
                db.executescript(
                    """
                    CREATE TABLE metadata(key,value);
                    CREATE TABLE feature(feature_id,gene_id,mapping_status);
                    CREATE TABLE cell_group(cell_type,tissue_system,n_cells);
                    CREATE TABLE summary(gene_id,cell_type,n_cells,n_expressing,mean_log,n_samples);
                """
                )
                db.execute(
                    "INSERT INTO metadata VALUES('dataset',?)",
                    (json.dumps({"citation": "test fixture"}),),
                )
                db.execute(
                    "INSERT INTO feature VALUES('Os01g0100100','Os01g0100100','exact RAP locus')"
                )
                for name, detected, mean in [
                    ("A", 1, math.log1p(10000) / 2),
                    ("B", 0, 0),
                    ("C", 0, 0),
                ]:
                    db.execute("INSERT INTO cell_group VALUES(?,?,2)", (name, "Leaf"))
                    db.execute(
                        "INSERT INTO summary VALUES(?,?,2,?,?,1)",
                        ("Os01g0100100", name, detected, mean),
                    )
            db.close()
            with patch.object(expression, "EXPRESSION_DB", path):
                data = expression.expression_evidence(["Os01g0100100", "Os01g0100200"])
            rows = data["genes"][0]["groups"]
            self.assertEqual(rows[0]["percent_expressing"], 50)
            self.assertAlmostEqual(rows[0]["mean_log"], math.log1p(10000) / 2)
            self.assertEqual([r["rank"] for r in rows], [1, 2, 2])
            self.assertEqual(rows[1]["percent_expressing"], 0)
            self.assertEqual(data["genes"][1]["status"], "not_available")
            self.assertEqual(data["genes"][1]["groups"], [])

    # Keep measurements and metadata separate when the same gene appears in two studies.
    def test_studies_are_not_pooled(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / "leaf.sqlite", Path(directory) / "root.sqlite"]
            ids = ["robertson_2026_control", "wang_2025_normal_root"]
            for index, path in enumerate(paths, 1):
                with sqlite3.connect(path) as db:
                    db.executescript(
                        "CREATE TABLE metadata(key,value); CREATE TABLE feature(feature_id,gene_id,mapping_status); CREATE TABLE cell_group(cell_type,tissue_system,n_cells); CREATE TABLE summary(gene_id,cell_type,n_cells,n_expressing,mean_log,n_samples);"
                    )
                    metadata = {
                        "id": ids[index - 1],
                        "citation": ids[index - 1],
                        "organ": "test",
                        "condition": "control",
                        "n_cells": index * 10,
                        "n_samples": 1,
                    }
                    db.execute(
                        "INSERT INTO metadata VALUES('dataset',?)",
                        (json.dumps(metadata),),
                    )
                    db.execute("INSERT INTO feature VALUES('gene','gene','exact')")
                    db.execute(
                        "INSERT INTO cell_group VALUES('celltype','test',?)",
                        (index * 10,),
                    )
                    db.execute(
                        "INSERT INTO summary VALUES('gene','celltype',?,1,?,1)",
                        (index * 10, index),
                    )
                db.close()
            with (
                patch.object(expression, "EXPRESSION_DB", paths[0]),
                patch.object(expression, "WANG_DB", paths[1]),
                patch.object(
                    expression, "WANG_ATLAS_DB", Path(directory) / "atlas_absent.sqlite"
                ),
                patch.object(expression, "dataset_enabled", return_value=True),
                patch.object(
                    expression, "LI_DB", Path(directory) / "not_imported.sqlite"
                ),
                patch.object(
                    expression, "LI_AUTHOR_DB", Path(directory) / "author_absent.sqlite"
                ),
            ):
                leaf = expression.expression_evidence(["gene"], ids[0])
                root = expression.expression_evidence(["gene"], ids[1])
                self.assertEqual(leaf["genes"][0]["groups"][0]["n_cells"], 10)
                self.assertEqual(root["genes"][0]["groups"][0]["n_cells"], 20)
                self.assertEqual(root["genes"][0]["groups"][0]["mean_log"], 2)
                self.assertEqual(len(expression.expression_datasets()), 2)

    # Hide a reverted study even when its ignored downloaded database still exists.
    def test_reverted_study_is_disabled(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "downloaded.sqlite"
            path.touch()
            with (
                patch.object(expression, "WANG_DB", path),
                patch.object(expression, "dataset_enabled", return_value=False),
            ):
                self.assertFalse(
                    expression.expression_evidence(["gene"], "wang_2025_normal_root")[
                        "available"
                    ]
                )
                self.assertEqual(expression.expression_datasets(), [])

    # Reject unknown study names rather than interpreting them as database paths.
    def test_unknown_study_rejected(self):
        with self.assertRaises(ValueError):
            expression.expression_evidence(["Os01g0100100"], "../private")

    # Bound comparison requests to prevent accidental full-atlas web responses.
    def test_comparison_limit(self):
        for genes in [[], ["x" * 41], list(map(str, range(13)))]:
            with self.assertRaises(ValueError):
                expression.expression_evidence(genes)

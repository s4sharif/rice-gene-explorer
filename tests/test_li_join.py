"""Check sample matching and non-overlapping exclusions for curated pistil annotations."""

import contextlib
import csv
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import prepare_li


class PistilJoinTests(unittest.TestCase):
    # Keep exact diploid matches while accounting for haploid, unknown, residue and unmatched nuclei.
    def test_curated_annotation_partition(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "data/raw/li"
            raw.mkdir(parents=True)
            config = {
                "prefixes": {"GSM1": "one", "GSM2": "two"},
                "experiments": {"SRX1": "GSM1", "SRX2": "GSM2"},
                "excluded_labels": ["Unknow", "Stamen residue"],
            }
            (root / "data/expression_li_selection.json").write_text(json.dumps(config))
            for sample, experiment, prefix in [
                ("GSM1", "SRX1", "one"),
                ("GSM2", "SRX2", "two"),
            ]:
                (raw / (sample + ".soft")).write_text(
                    "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term="
                    + experiment
                )
                with gzip.open(raw / (prefix + "_barcodes.tsv.gz"), "wt") as handle:
                    handle.write("A-1\nB-1\nC-1\n")
            records = [
                ("SRX1", "A-1", "Style"),
                ("SRX2", "A-1", "Nucellus"),
                ("SRX1", "B-1", "Unknow"),
                ("SRX1", "C-1", "Stamen residue"),
                ("SRX1", "missing-1", "Style"),
                ("SRX16335166", "A-1", "Style"),
            ]
            with (raw / "cell_metadata.tsv").open("w") as handle:
                fields = [
                    "Dataset",
                    "Condition",
                    "Tissue",
                    "Genotype",
                    "Orig.ident",
                    "barcode",
                    "Celltype",
                ]
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for experiment, barcode, label in records:
                    writer.writerow(
                        {
                            "Dataset": "SRP386976",
                            "Condition": "Normal",
                            "Tissue": "Pistil",
                            "Genotype": "japonica cv. Nipponbare",
                            "Orig.ident": experiment,
                            "barcode": experiment + "@@_" + barcode,
                            "Celltype": label,
                        }
                    )
            with (
                patch.object(prepare_li, "ROOT", root),
                patch.object(prepare_li, "RAW", raw),
                contextlib.redirect_stdout(io.StringIO()),
            ):
                prepare_li.main()
            report = json.loads(
                (root / "data/expression_li_join_validation.json").read_text()
            )
            self.assertEqual(report["included_nuclei"], 2)
            self.assertEqual(report["included_groups"], 2)
            self.assertEqual(
                report["included_nuclei"] + sum(report["exclusions"].values()),
                len(records),
            )
            self.assertEqual(
                report["exclusions"]["curated_barcode_absent_from_GEO_matrix"], 1
            )

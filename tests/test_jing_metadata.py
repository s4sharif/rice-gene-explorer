"""Reject ambiguous gene transfers and keep metadata outside quantitative expression."""

import unittest

from prepare_jing_metadata import map_genes
from source_coverage import source_coverage


class JingMetadataTests(unittest.TestCase):
    # Reject author multi-target mappings and collisions even when local aliases resolve uniquely.
    def test_conservative_crosswalk(self):
        records = map_genes(
            {
                "ZH1": {"msu1"},
                "ZH2": {"msu2"},
                "ZH3": {"msu3", "ChrUn"},
                "ZH4": {"msu4"},
                "ZH5": {"msu5"},
            },
            {
                "msu1": {"rap1"},
                "msu2": {"rap1"},
                "msu3": {"rap3"},
                "msu4": {"rap4"},
                "msu5": {"rap5", "rap6"},
            },
        )
        self.assertEqual(
            [r["rap_id"] for r in records], [None, None, None, "rap4", None]
        )

    # Keep reviewed publications out of the expression data model and reject malformed gene keys.
    def test_coverage_has_no_expression_values(self):
        coverage = source_coverage("Os03g0407400")
        for source in coverage["sources"]:
            self.assertNotIn("genes", source)
            self.assertNotIn("mean_log", source)
        with self.assertRaises(ValueError):
            source_coverage("../../private")


if __name__ == "__main__":
    unittest.main()

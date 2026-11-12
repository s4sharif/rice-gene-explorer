"""Protect exact author barcode joins and preserve unresolved cluster identities."""

import unittest

from pipelines.prepare_li_author import join_annotations


class LiAuthorTests(unittest.TestCase):
    # Retain original author clusters while excluding explicitly identified stamen residue.
    def test_author_labels_and_residue(self):
        rows = [
            {"A": "rep1_A-1", "B": "ricerep1", "C": "12", "D": "C12", "E": "C12"},
            {"A": "rep2_B-1", "B": "ricerep2", "C": "11", "D": "sta", "E": "sta"},
        ]
        matched, excluded = join_annotations(
            rows, {"GSM6346349": {"A-1"}, "GSM6346350": {"B-1"}}
        )
        self.assertEqual(matched[0]["cell_type"], "Cluster 12 (C12)")
        self.assertEqual(matched[0]["source_cluster"], "12")
        self.assertEqual(matched[0]["sample"], "GSM6346349")
        self.assertEqual(excluded, {"stamen_residue": 1})

    # Fail on missing, repeated or mislabeled replicate keys rather than silently losing nuclei.
    def test_bad_join_is_rejected(self):
        row = {"A": "rep1_A-1", "B": "ricerep1", "C": "1", "D": "ep", "E": "ca"}
        known = {"GSM6346349": {"A-1"}, "GSM6346350": {"B-1"}}
        with self.assertRaises(ValueError):
            join_annotations([row, row], known)
        with self.assertRaises(ValueError):
            join_annotations([{**row, "A": "rep1_missing-1"}], known)
        with self.assertRaises(ValueError):
            join_annotations([{**row, "B": "ricerep2"}], known)


if __name__ == "__main__":
    unittest.main()

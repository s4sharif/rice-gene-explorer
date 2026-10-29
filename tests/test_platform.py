import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import app


class PlatformTests(unittest.TestCase):
    # Verify identifier mappings while retaining every candidate for an ambiguous symbol.
    def test_gene_alias_and_symbol_resolve_same_locus(self):
        result = app.search("GS3")["matches"]
        gene = "Os03g0407400"
        self.assertIn(gene, [r["id"] for r in result])
        # Shared symbols must remain ambiguous rather than selecting the first gene.
        self.assertGreater(len(result), 1)
        with app.connect() as db:
            aliases = [
                x[0]
                for x in db.execute(
                    "SELECT name FROM alias WHERE gene_id=? AND kind='MSU locus'",
                    ("Os01g0100100",),
                )
            ]
        self.assertTrue(aliases)
        for alias in aliases:
            self.assertIn(
                "Os01g0100100", [x["id"] for x in app.search(alias)["matches"]]
            )
        self.assertEqual(app.search(gene.lower())["matches"][0]["id"], gene)

    # Check inclusive gene boundaries and strand-aware upstream and downstream labels.
    def test_coordinate_boundaries_and_strand(self):
        g = {"start": 100, "end": 200, "strand": "+"}
        self.assertEqual(app.relation(g, 100), ("overlaps gene", 0))
        self.assertEqual(app.relation(g, 200), ("overlaps gene", 0))
        self.assertEqual(app.relation(g, 99), ("upstream", 1))
        self.assertEqual(app.relation(g, 201), ("downstream", 1))
        g["strand"] = "-"
        self.assertEqual(app.relation(g, 99), ("downstream", 1))
        self.assertEqual(app.relation(g, 201), ("upstream", 1))

    # Compare regional results and complete exports with an independent SQL count.
    def test_region_and_export_match_direct_sql(self):
        d = app.region("Os03g0407400", 25000)
        with app.connect() as db:
            n = db.execute(
                "SELECT count(*) FROM association WHERE assembly_id=? AND chromosome=? AND position BETWEEN ? AND ?",
                (d["gene"]["assembly_id"], 3, d["start"], d["end"]),
            ).fetchone()[0]
        self.assertGreater(n, 0)
        self.assertEqual(n, d["total_associations"])
        self.assertEqual(len(app.region("Os03g0407400", 25000, export=True)), n)
        self.assertTrue(
            all(d["start"] <= a["position"] <= d["end"] for a in d["associations"])
        )
        for g in d["neighbors"]:
            self.assertTrue(g["start"] <= d["end"] and g["end"] >= d["start"])

    # Verify that larger windows preserve evidence and trait filters restrict results.
    def test_window_expansion_and_trait_filter(self):
        small = app.region("Os03g0407400", 10000)
        large = app.region("Os03g0407400", 100000)
        self.assertGreaterEqual(
            large["total_associations"], small["total_associations"]
        )
        filtered = app.region("Os03g0407400", 100000, "grain")
        self.assertTrue(
            all("grain" in a["trait"].lower() for a in filtered["associations"])
        )
        self.assertEqual(
            app.region("Os03g0407400", 10000, "not-a-real-trait")["total_associations"],
            0,
        )

    # Confirm that unresolved references are quarantined and database relationships remain valid.
    def test_ambiguous_references_excluded(self):
        with app.connect() as db:
            self.assertEqual(
                db.execute(
                    "SELECT count(*) FROM association WHERE source_reference='1'"
                ).fetchone()[0],
                0,
            )
            self.assertGreater(
                db.execute(
                    "SELECT count(*) FROM quarantine WHERE reason LIKE 'Unresolved reference%'"
                ).fetchone()[0],
                0,
            )
            self.assertEqual(db.execute("PRAGMA integrity_check").fetchone()[0], "ok")
            self.assertEqual(db.execute("PRAGMA foreign_key_check").fetchall(), [])

    # Reject invalid windows and treat query text as literal data.
    def test_validation_and_literal_queries(self):
        for w in [-1, 1000001]:
            with self.assertRaises(ValueError):
                app.region("Os03g0407400", w)
        with self.assertRaises(ValueError):
            app.region("fake")
        self.assertEqual(app.search("' OR 1=1 --")["matches"], [])
        self.assertEqual(app.search("%")["matches"], [])


if __name__ == "__main__":
    unittest.main()

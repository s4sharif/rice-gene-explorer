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

    # Preserve shared SD1 aliases while prioritizing the exact display symbol and returning descriptions.
    def test_sd1_ambiguity_is_explained(self):
        with app.connect() as db:
            duplicate = db.execute(
                "SELECT MIN(name) FROM alias GROUP BY gene_id,lower(name) HAVING COUNT(*)>1 LIMIT 1"
            ).fetchone()[0]
        matches = app.search(duplicate)["matches"]
        self.assertEqual(len(matches), len({g["id"] for g in matches}))
        result = app.search("SD1")
        self.assertEqual(result["match_type"], "exact")
        self.assertEqual(result["matches"][0]["id"], "Os01g0883800")
        self.assertIn("Os03g0130100", [g["id"] for g in result["matches"]])
        for gene in result["matches"]:
            self.assertEqual(gene["matched_identifier"].lower(), "sd1")
            self.assertTrue(gene["description"])

    # Keep unverified adjustment status distinct from verified FDR values in plots and exports.
    def test_adjustment_status_and_coordinate_provenance(self):
        data = app.region("Os03g0407400", 25000)
        rows = data["plot_associations"]
        self.assertTrue(
            any(a["adjustment_status"] == "verified_fdr_adjusted" for a in rows)
        )
        self.assertTrue(any(a["adjustment_status"] == "unverified" for a in rows))
        for row in app.region("Os03g0407400", 25000, export=True):
            self.assertIn("IRGSP-1.0", row["coordinate_status"])
            if row["adjustment_status"] == "verified_fdr_adjusted":
                self.assertEqual(row["statistic_kind"], "FDR-adjusted p-value")
            else:
                self.assertEqual(row["adjustment_status"], "unverified")

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
                "SELECT count(*) FROM association_display WHERE assembly_id=? AND chromosome=? AND position BETWEEN ? AND ?",
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

    # Verify eight-row pages are distinct and keep every matching record available for export.
    def test_eight_record_pagination(self):
        first = app.region("Os03g0407400", 25000, page=1)
        second = app.region("Os03g0407400", 25000, page=2)
        self.assertEqual(first["page_size"], 8)
        self.assertEqual(len(first["associations"]), 8)
        self.assertEqual(len(second["associations"]), 8)
        first_ids = {row["id"] for row in first["associations"]}
        second_ids = {row["id"] for row in second["associations"]}
        self.assertTrue(first_ids.isdisjoint(second_ids))
        exported = app.region("Os03g0407400", 25000, export=True)
        self.assertEqual(len(exported), first["total_associations"])

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

    # Verify study, publication, p-value and strand-aware position filters agree with exports.
    def test_evidence_filters(self):
        baseline = app.region("Os03g0407400", 25000)
        publication = "33367255"
        study = baseline["study_options"][0]["id"]
        for filters in [
            {"publication": publication},
            {"study": study},
            {"maximum_p": 1e-10},
            {"position_filter": "upstream"},
            {"position_filter": "overlap"},
        ]:
            result = app.region("Os03g0407400", 25000, **filters)
            exported = app.region("Os03g0407400", 25000, export=True, **filters)
            self.assertEqual(result["total_associations"], len(exported))
            for row in exported:
                if "publication" in filters:
                    self.assertEqual(row["pmid"], publication)
                if "maximum_p" in filters:
                    self.assertLessEqual(row["p_value"], filters["maximum_p"])
                if "position_filter" in filters:
                    expected = (
                        "overlaps gene"
                        if filters["position_filter"] == "overlap"
                        else filters["position_filter"]
                    )
                    self.assertEqual(row["relation"], expected)

    # Verify repeat grouping preserves all source rows and exposes an explicit expansion option.
    def test_repeat_grouping(self):
        grouped = app.region("Os03g0407400", 25000)
        raw = app.region("Os03g0407400", 25000, include_repeats=True)
        self.assertGreater(raw["total_associations"], grouped["total_associations"])
        exported = app.region("Os03g0407400", 25000, export=True)
        self.assertEqual(
            sum(row["source_record_count"] for row in exported),
            raw["total_associations"],
        )

    # Verify out-of-range pages clamp safely and chromosome edges bound every search window.
    def test_page_and_chromosome_edges(self):
        last = app.region("Os03g0407400", 25000, page=9999)
        self.assertEqual(last["page"], last["page_count"])
        self.assertLessEqual(len(last["associations"]), 8)
        first_gene = app.region("Os01g0100100", 1000000)
        self.assertEqual(first_gene["start"], 1)
        with app.connect() as db:
            last_gene = db.execute(
                "SELECT id FROM gene WHERE chromosome=1 ORDER BY end DESC LIMIT 1"
            ).fetchone()[0]
        last_gene_region = app.region(last_gene, 1000000)
        self.assertEqual(last_gene_region["end"], 43270923)
        empty = app.region("Os03g0407400", trait="no-such-trait", page=9999)
        self.assertEqual(
            (empty["page"], empty["page_count"], empty["associations"]), (1, 1, [])
        )

    # Verify original-paper conflicts and adjusted statistics remain visible in record exports.
    def test_publication_review_flags(self):
        records = app.region("Os03g0407400", 25000, publication="33367255", export=True)
        conflicts = [row for row in records if row["review_status"] == "conflict"]
        self.assertTrue(conflicts)
        self.assertTrue(
            any("Trait conflict" in row["curation_note"] for row in conflicts)
        )
        self.assertTrue(
            any(row["statistic_kind"] == "FDR-adjusted p-value" for row in records)
        )

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

"""Verify shared genomic searches and the revised GWAS presentation contracts."""

import csv
import io
import json
import os
import unittest
from unittest.mock import patch

from test_hosting import PASSWORD, USERNAME, request

from backend import app


class GenomicLayoutTests(unittest.TestCase):
    # Preserve an exact interval search even when it contains several genes or no genes.
    def test_coordinate_search_and_empty_interval(self):
        search = app.search("Chr3:16,729,501–16,735,109")
        self.assertEqual(search["match_type"], "locus")
        self.assertEqual(search["selector"], "chr3:16729501-16735109")
        locus = app.region(search["selector"], 25000)
        gene = app.region("Os03g0407400", 25000)
        self.assertEqual((locus["start"], locus["end"]), (gene["start"], gene["end"]))
        self.assertEqual(locus["total_associations"], gene["total_associations"])
        self.assertEqual(locus["selection_type"], "locus")
        self.assertEqual(locus["gene"]["id"], "")
        self.assertEqual(app.locus_genes("chr1:1", 0)["genes"], [])
        for invalid in [
            "chr13:1",
            "chr1:0",
            "chr3:100-99",
            "chr1:999999999",
            "chr1:1-3000000",
            "chrX:1",
        ]:
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                app.search(invalid)

    # Sort every filtered association before pagination and keep the chart's strongest-record policy.
    def test_sorting_summary_and_exports(self):
        selector = "Os03g0407400"
        ascending = app.region(selector, export=True)
        descending = app.region(selector, export=True, sort="p_desc")
        self.assertEqual(
            [r["p_value"] for r in ascending], sorted(r["p_value"] for r in ascending)
        )
        self.assertEqual(
            [r["p_value"] for r in descending],
            sorted((r["p_value"] for r in descending), reverse=True),
        )
        self.assertEqual({r["id"] for r in ascending}, {r["id"] for r in descending})
        second = app.region(selector, page=2, sort="p_desc")
        self.assertEqual(
            [r["id"] for r in second["associations"]],
            [r["id"] for r in descending[8:16]],
        )
        self.assertEqual(second["summary"]["smallest_p"], ascending[0]["p_value"])
        self.assertEqual(second["summary"]["records"], len(ascending))
        for top in second["summary"]["top_traits"]:
            self.assertEqual(
                top["count"], sum(r["trait"] == top["trait"] for r in ascending)
            )
        empty = app.region(selector, trait="no-such-trait")
        self.assertIsNone(empty["summary"]["smallest_p"])
        self.assertEqual(empty["summary"]["top_traits"], [])
        with self.assertRaises(ValueError):
            app.region(selector, sort="p_value;DROP TABLE gene")

    # Include aliases and only overlapping source-observed sites without guessing reference alleles.
    def test_gene_table_and_qtl_provenance(self):
        data = app.locus_genes("Os03g0407400")
        self.assertEqual(len(data["genes"]), 3)
        for gene in data["genes"]:
            self.assertTrue(gene["aliases"])
            self.assertIsNone(gene["reference_allele"])
            self.assertIsNone(gene["alternate_allele"])
            for variant in gene["variants"]:
                self.assertLessEqual(gene["start"], variant["position"])
                self.assertLessEqual(variant["position"], gene["end"])
        with app.connect() as db:
            row = db.execute(
                "SELECT chromosome,position,json_extract(raw_json,'$.QTL') FROM association WHERE trim(coalesce(json_extract(raw_json,'$.QTL'),''))!='' LIMIT 1"
            ).fetchone()
        records = app.region(
            f"chr{row[0]}:{row[1]}", 0, export=True, include_repeats=True
        )
        self.assertTrue(any(r["qtl_label"] == row[2] for r in records))
        self.assertTrue(all(r["evidence_type"] == "GWAS" for r in records))
        self.assertTrue(all("experimental_condition" in r for r in records))

    # Derive homepage statistics from imported data and exclude overlapping annotation versions.
    def test_overview_counts(self):
        data = app.overview()
        self.assertEqual(data["study_contexts"], 828)
        self.assertEqual(data["publications"], 152)
        self.assertEqual(data["traits"], 357)
        self.assertEqual(
            data["atlas_cells"], sum(d["n_cells"] for d in data["atlas_datasets"])
        )
        self.assertNotIn(
            "li_2023_pistil_2n_scplantdb", [d["id"] for d in data["atlas_datasets"]]
        )


class GenomicRouteTests(unittest.TestCase):
    # Configure disposable credentials so API routes are tested behind the hosted access gate.
    def setUp(self):
        environment = patch.dict(
            os.environ, RICELENS_USERNAME=USERNAME, RICELENS_PASSWORD=PASSWORD
        )
        environment.start()
        self.addCleanup(environment.stop)

    # Serve new pages and exports through the same mandatory reviewer authentication.
    def test_new_routes_and_full_gene_csv(self):
        for route in [
            "/locus",
            "/homolog",
            "/download",
            "/tutorial",
            "/contact",
            "/api/overview",
            "/api/genes?gene=Os03g0407400",
            "/api/genes/export?gene=Os03g0407400",
        ]:
            with self.subTest(route=route):
                self.assertTrue(request(route)["status"].startswith("401"))
                self.assertTrue(
                    request(route, USERNAME + ":" + PASSWORD)["status"].startswith(
                        "200"
                    )
                )
        response = request(
            "/api/genes/export?gene=Os03g0407400", USERNAME + ":" + PASSWORD
        )
        rows = list(csv.DictReader(io.StringIO(response["body"].decode())))
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(r["reference_allele"] == "" for r in rows))
        response = request(
            "/api/region?gene=chr3:16729501-16735109&sort=p_desc",
            USERNAME + ":" + PASSWORD,
        )
        self.assertEqual(json.loads(response["body"])["sort"], "p_desc")

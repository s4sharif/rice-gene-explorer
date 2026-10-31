"""Local, read-only gene search and genomic-neighborhood API."""

import argparse
import csv
import io
import json
import math
import sqlite3
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from build_database import CHROMOSOME_LENGTHS
from expression import expression_datasets, expression_evidence

ROOT = Path(__file__).resolve().parent
DB = ROOT / "data/rice.sqlite"
PAGE_SIZE = 8
PUBLICATION_REVIEW = json.loads((ROOT / "data/publication_review.json").read_text())
PUBLICATION_REVIEWS = PUBLICATION_REVIEW["records"]


# Open a read-only database connection and close it after the query finishes.
@contextmanager
def connect():
    db = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    try:
        yield db
    finally:
        db.close()


# Resolve an identifier or symbol to matching genes without choosing among ambiguous matches.
def search(query):
    query = query.strip()
    if not query or len(query) > 120:
        raise ValueError("Enter a gene identifier or symbol (maximum 120 characters).")
    with connect() as db:
        exact = db.execute(
            "SELECT g.*,MIN(a.name) AS matched_identifier,MIN(a.kind) AS matched_identifier_kind FROM gene g JOIN alias a ON a.gene_id=g.id WHERE a.name=? COLLATE NOCASE GROUP BY g.id ORDER BY CASE WHEN g.symbol=? COLLATE NOCASE THEN 0 ELSE 1 END,g.id",
            (query, query),
        ).fetchall()
        if exact:
            return {
                "matches": [dict(r) for r in exact],
                "match_type": "exact",
                "truncated": False,
            }
        # Literal prefix search; SQL wildcards typed by users are escaped.
        prefix = (
            query.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "%"
        )
        rows = db.execute(
            "SELECT DISTINCT g.* FROM gene g JOIN alias a ON a.gene_id=g.id WHERE a.name LIKE ? ESCAPE '\\' COLLATE NOCASE ORDER BY g.id LIMIT 51",
            (prefix,),
        ).fetchall()
        return {
            "matches": [dict(r) for r in rows[:50]],
            "match_type": "prefix",
            "truncated": len(rows) > 50,
        }


# Classify a variant relative to the gene strand and measure its distance from the gene span.
def relation(gene, position):
    if gene["start"] <= position <= gene["end"]:
        return "overlaps gene", 0
    left = position < gene["start"]
    distance = gene["start"] - position if left else position - gene["end"]
    upstream = left if gene["strand"] == "+" else not left
    return ("upstream" if upstream else "downstream"), distance


# Retrieve a gene neighborhood and its matching published association evidence.
def region(
    gene_id,
    window=25000,
    trait="",
    page=1,
    export=False,
    study="",
    publication="",
    maximum_p=1.0,
    position_filter="",
    include_repeats=False,
):
    if not 0 <= window <= 1000000:
        raise ValueError("Window must be between 0 and 1,000,000 bp on each side.")
    if page < 1 or page > 100000:
        raise ValueError("Invalid results page.")
    if len(trait) > 200:
        raise ValueError("Trait filter too long.")
    if len(study) > 100 or (publication and not publication.isdigit()):
        raise ValueError("Select a valid study or publication.")
    if not math.isfinite(maximum_p) or not 0 < maximum_p <= 1:
        raise ValueError("Maximum p-value must be greater than zero and at most one.")
    if position_filter not in {"", "overlap", "upstream", "downstream"}:
        raise ValueError("Select a recognized position filter.")
    with connect() as db:
        row = db.execute("SELECT * FROM gene WHERE id=?", (gene_id,)).fetchone()
        if not row:
            raise ValueError("Select a recognized RAP gene ID from the search results.")
        gene = dict(row)
        start = max(1, gene["start"] - window)
        end = min(CHROMOSOME_LENGTHS[gene["chromosome"]], gene["end"] + window)
        genes = [
            dict(r)
            for r in db.execute(
                "SELECT * FROM gene WHERE assembly_id=? AND chromosome=? AND start<=? AND end>=? ORDER BY start,id",
                (gene["assembly_id"], gene["chromosome"], end, start),
            )
        ]
        args = [gene["assembly_id"], gene["chromosome"], start, end]
        where = "a.assembly_id=? AND a.chromosome=? AND a.position BETWEEN ? AND ?"
        base_args = args.copy()
        base_where = where
        if trait:
            where += " AND instr(lower(t.label),lower(?))>0"
            args.append(trait)
        if study:
            where += " AND s.id=?"
            args.append(study)
        if publication:
            where += " AND p.pmid=?"
            args.append(publication)
        where += " AND a.p_value<=?"
        args.append(maximum_p)
        if position_filter == "overlap":
            where += " AND a.position BETWEEN ? AND ?"
            args.extend([gene["start"], gene["end"]])
        elif position_filter:
            left = (position_filter == "upstream") == (gene["strand"] == "+")
            where += " AND a.position" + ("<?" if left else ">?")
            args.append(gene["start"] if left else gene["end"])
        table = "association" if include_repeats else "association_display"
        join = (
            " FROM "
            + table
            + " a JOIN trait t ON t.id=a.trait_id JOIN study s ON s.id=a.study_id LEFT JOIN publication p ON p.pmid=s.publication_id WHERE "
            + where
        )
        total = db.execute("SELECT count(*)" + join, args).fetchone()[0]
        source_total = db.execute(
            "SELECT count(*)"
            + join.replace(" FROM " + table + " a ", " FROM association a "),
            args,
        ).fetchone()[0]
        options_join = (
            " FROM association a JOIN study s ON s.id=a.study_id LEFT JOIN publication p ON p.pmid=s.publication_id WHERE "
            + base_where
        )
        study_options = [
            dict(r)
            for r in db.execute(
                "SELECT DISTINCT s.id,s.source_study_id,s.model"
                + options_join
                + " ORDER BY s.source_study_id,s.id",
                base_args,
            )
        ]
        publication_options = [
            dict(r)
            for r in db.execute(
                "SELECT DISTINCT p.pmid,p.title"
                + options_join
                + " AND p.pmid IS NOT NULL ORDER BY p.pmid",
                base_args,
            )
        ]
        page_count = max(1, math.ceil(total / PAGE_SIZE))
        page = min(page, page_count)
        repeat_count = "1" if include_repeats else "a.source_record_count"
        cols = (
            "SELECT "
            + repeat_count
            + " AS source_record_count,a.id,a.chromosome,a.position,a.p_value,a.p_value_text,a.effect_allele,a.effect_allele_frequency,a.other_allele,a.effect_size,a.r2_percent,a.reported_location,a.source_reference,a.source_row,t.label AS trait,t.ontology_id,s.source_study_id,s.model,s.sample_size,s.population,s.environment,s.tissue,p.pmid,p.title,p.url AS publication_url"
        )

        # Add gene-relative distances, plot values and source links to association records.
        def enrich(rows):
            result = []
            for row in rows:
                a = dict(row)
                a["relation"], a["distance_bp"] = relation(gene, a["position"])
                a["neg_log10_p"] = -math.log10(a["p_value"])
                a["assembly"] = gene["assembly_id"]
                a["source_url"] = "https://ngdc.cncb.ac.cn/gwas/download"
                a["source_record_url"] = "/api/association?id=" + a["id"]
                review = PUBLICATION_REVIEWS.get(a["id"], {})
                a["review_status"] = (
                    "conflict"
                    if review.get("conflicts")
                    else (
                        "coordinate/statistic checked"
                        if review
                        else "not independently reviewed"
                    )
                )
                a["statistic_kind"] = review.get(
                    "statistic_kind", "Adjustment status not independently checked"
                )
                a["adjustment_status"] = (
                    "verified_fdr_adjusted"
                    if a["statistic_kind"] == "FDR-adjusted p-value"
                    else "unverified"
                )
                a["coordinate_status"] = (
                    "Source-labelled IRGSP-1.0; original-study remapping not independently verified"
                )
                a["curation_note"] = " ".join(review.get("conflicts", []))
                a["original_table_trait"] = review.get("original_trait")
                a["original_table_url"] = review.get("original_source_url")
                a["original_major_allele"] = review.get("major_allele")
                a["original_minor_allele"] = review.get("minor_allele")
                a["original_minor_allele_frequency"] = review.get(
                    "minor_allele_frequency"
                )
                result.append(a)
            return result

        order = " ORDER BY a.p_value,a.position,a.id"
        if export:
            return enrich(db.execute(cols + join + order, args))
        associations = enrich(
            db.execute(
                cols + join + order + " LIMIT ? OFFSET ?",
                args + [PAGE_SIZE, (page - 1) * PAGE_SIZE],
            )
        )
        plot = enrich(db.execute(cols + join + order + " LIMIT 500", args))
        aliases = [
            dict(r)
            for r in db.execute(
                "SELECT name,kind FROM alias WHERE gene_id=? ORDER BY kind,name",
                (gene_id,),
            )
        ]
        return {
            "gene": gene,
            "aliases": aliases,
            "start": start,
            "end": end,
            "window": window,
            "neighbors": genes,
            "associations": associations,
            "plot_associations": plot,
            "total_associations": total,
            "plot_count": len(plot),
            "page": page,
            "page_size": PAGE_SIZE,
            "trait_filter": trait,
            "page_count": page_count,
            "source_record_total": source_total,
            "repeated_rows_hidden": source_total - total,
            "study_options": study_options,
            "publication_options": publication_options,
            "filters": {
                "study": study,
                "publication": publication,
                "maximum_p": maximum_p,
                "position": position_filter,
                "include_repeats": include_repeats,
            },
            "interpretation": "Regional association evidence does not establish gene causality. Overlap refers to the gene span, including introns. Upstream/downstream follow gene strand.",
        }


class Handler(BaseHTTPRequestHandler):
    # Route a browser request to the appropriate read-only API or static resource.
    def do_GET(self):
        try:
            parsed = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(parsed.query).items()}
            if parsed.path == "/api/expression/datasets":
                return self.send_json({"datasets": expression_datasets()})
            if parsed.path in {"/api/expression", "/api/expression/export"}:
                genes = q.get("genes", "").split(",")
                with connect() as db:
                    for gene in genes:
                        if not db.execute(
                            "SELECT 1 FROM gene WHERE id=?", (gene,)
                        ).fetchone():
                            raise ValueError("Select recognized RAP gene identifiers.")
                data = expression_evidence(
                    genes, q.get("dataset", "robertson_2026_control")
                )
                if parsed.path.endswith("/export"):
                    if not data["available"]:
                        raise ValueError(data["message"])
                    out = io.StringIO()
                    fields = [
                        "gene_id",
                        "status",
                        "source_feature_id",
                        "mapping_status",
                        "cell_type",
                        "tissue_system",
                        "rank",
                        "mean_log",
                        "n_expressing",
                        "n_cells",
                        "percent_expressing",
                        "n_samples",
                        "dataset",
                        "publication_url",
                        "metric",
                        "annotation_source",
                        "annotation_source_url",
                        "annotation_publication_url",
                        "annotation_note",
                    ]
                    writer = csv.DictWriter(out, fieldnames=fields)
                    writer.writeheader()
                    for gene in data["genes"]:
                        for group in gene["groups"] or [{}]:
                            writer.writerow(
                                {
                                    **group,
                                    "gene_id": gene["gene_id"],
                                    "status": gene["status"],
                                    "source_feature_id": (gene["feature"] or {}).get(
                                        "feature_id"
                                    ),
                                    "mapping_status": (gene["feature"] or {}).get(
                                        "mapping_status"
                                    ),
                                    "dataset": data["dataset"]["id"],
                                    "publication_url": data["dataset"][
                                        "publication_url"
                                    ],
                                    "metric": data["dataset"]["metric_description"],
                                    "annotation_source": data["dataset"].get(
                                        "annotation_source_name",
                                        "Original deposited annotations",
                                    ),
                                    "annotation_source_url": data["dataset"].get(
                                        "annotation_source_url",
                                        data["dataset"]["repository_url"],
                                    ),
                                    "annotation_publication_url": data["dataset"].get(
                                        "annotation_publication_url",
                                        data["dataset"]["publication_url"],
                                    ),
                                    "annotation_note": data["dataset"].get(
                                        "annotation_note", ""
                                    ),
                                }
                            )
                    return self.respond(
                        out.getvalue().encode(),
                        "text/csv; charset=utf-8",
                        attachment="ricelens_baseline_expression.csv",
                    )
                return self.send_json(data)
            if parsed.path == "/api/search":
                return self.send_json(search(q.get("q", "")))
            if parsed.path in {"/api/region", "/api/export"}:
                result = region(
                    q.get("gene", ""),
                    int(q.get("window", 25000)),
                    q.get("trait", ""),
                    int(q.get("page", 1)),
                    parsed.path == "/api/export",
                    study=q.get("study", ""),
                    publication=q.get("publication", ""),
                    maximum_p=float(q.get("maximum_p", 1)),
                    position_filter=q.get("position", ""),
                    include_repeats=q.get("include_repeats", "false") == "true",
                )
                if parsed.path == "/api/export":
                    out = io.StringIO()
                    fields = (
                        list(result[0])
                        if result
                        else [
                            "id",
                            "position",
                            "trait",
                            "p_value",
                            "relation",
                            "distance_bp",
                        ]
                    )
                    writer = csv.DictWriter(out, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(result)
                    return self.respond(
                        out.getvalue().encode(),
                        "text/csv; charset=utf-8",
                        attachment="rice_region_associations.csv",
                    )
                return self.send_json(result)
            if parsed.path == "/api/association":
                with connect() as db:
                    record = db.execute(
                        "SELECT * FROM association WHERE id=?", (q.get("id", ""),)
                    ).fetchone()
                    if not record:
                        return self.send_json(
                            {"error": "Association record not found."}, 404
                        )
                    record = dict(record)
                    record["primary_table_review"] = PUBLICATION_REVIEWS.get(
                        record["id"]
                    )
                    record["original_source_row"] = json.loads(record.pop("raw_json"))
                    return self.send_json(record)
            if parsed.path == "/api/status":
                manifest = json.loads((ROOT / "data/provenance.json").read_text())
                return self.send_json(manifest)
            files = {"/": "index.html", "/app.js": "app.js", "/style.css": "style.css"}
            if parsed.path not in files:
                return self.respond(b"Not found", "text/plain", 404)
            name = files[parsed.path]
            mime = (
                "text/html"
                if name.endswith("html")
                else "application/javascript" if name.endswith("js") else "text/css"
            )
            return self.respond(
                (ROOT / "static" / name).read_bytes(), mime + "; charset=utf-8"
            )
        except (ValueError, OverflowError) as e:
            self.send_json({"error": str(e)}, 400)
        except Exception:
            import traceback

            traceback.print_exc()
            self.send_json({"error": "The request could not be completed."}, 500)

    # Serialize a response as JSON without allowing nonfinite numeric values.
    def send_json(self, data, status=200):
        self.respond(
            json.dumps(data, allow_nan=False).encode(),
            "application/json; charset=utf-8",
            status,
        )

    # Send a response with explicit content headers and an optional download filename.
    def respond(self, body, mime, status=200, attachment=None):
        self.send_response(status)
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        if attachment:
            self.send_header(
                "Content-Disposition", f'attachment; filename="{attachment}"'
            )
        self.end_headers()
        self.wfile.write(body)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    if not DB.exists():
        raise SystemExit("Run python3 build_database.py first.")
    print(f"RiceLENS: http://127.0.0.1:{args.port}", flush=True)
    ThreadingHTTPServer(("127.0.0.1", args.port), Handler).serve_forever()

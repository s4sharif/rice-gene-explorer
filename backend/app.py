"""Local, read-only gene search and genomic-neighborhood API."""

import argparse
import csv
import io
import json
import math
import re
import sqlite3
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from backend.expression import expression_datasets, expression_evidence
from backend.runtime_data import DATA_DIR
from backend.source_coverage import source_coverage
from pipelines.build_database import CHROMOSOME_LENGTHS

ROOT = Path(__file__).resolve().parents[1]
DB = DATA_DIR / "rice.sqlite"
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


# Validate a point or inclusive genomic interval against the loaded reference assembly.
def parse_locus(query):
    normalized = re.sub(r"[\s,]", "", query.lower())
    match = re.fullmatch(r"(?:chr)?(\d+):(\d+)(?:[-–](\d+))?", normalized)
    if not match:
        if ":" in normalized:
            raise ValueError(
                "Use a genomic position such as chr3:16729501 or chr3:16729501-16735109."
            )
        return None
    chromosome, start = int(match[1]), int(match[2])
    end = int(match[3] or start)
    if (
        chromosome not in CHROMOSOME_LENGTHS
        or not 1 <= start <= end <= CHROMOSOME_LENGTHS[chromosome]
    ):
        raise ValueError(
            "Coordinates must be within rice chromosomes 1–12 on IRGSP-1.0, with start no greater than end."
        )
    if end - start > 2000000:
        raise ValueError(
            "Search intervals must span at most 2 Mb before adding the flanking window."
        )
    return chromosome, start, end


# Format a validated coordinate search so bookmarks preserve its exact interval.
def canonical_locus(locus):
    chromosome, start, end = locus
    return f"chr{chromosome}:{start}-{end}"


# Report live database counts with an explicit canonical annotation choice for overlapping nuclei.
def overview():
    with connect() as db:
        result = {
            "study_contexts": db.execute(
                "SELECT count(DISTINCT study_id) FROM association"
            ).fetchone()[0],
            "publications": db.execute(
                "SELECT count(DISTINCT s.publication_id) FROM association a JOIN study s ON s.id=a.study_id"
            ).fetchone()[0],
            "traits": db.execute("SELECT count(*) FROM trait").fetchone()[0],
            "genes": db.execute("SELECT count(*) FROM gene").fetchone()[0],
            "association_rows": db.execute(
                "SELECT count(*) FROM association"
            ).fetchone()[0],
            "qtl_named_rows": db.execute(
                "SELECT count(*) FROM association WHERE trim(coalesce(json_extract(raw_json,'$.QTL'),''))!=''"
            ).fetchone()[0],
        }
    datasets = expression_datasets()
    ids = {d["id"] for d in datasets}
    selected = [
        d
        for d in datasets
        if d["id"] != "li_2023_pistil_2n_scplantdb"
        or "li_2023_pistil_2n_author" not in ids
    ]
    result["atlas_cells"] = sum(d["n_cells"] for d in selected)
    result["atlas_datasets"] = selected
    result["cell_count_note"] = (
        "Cells and nuclei in the selected annotation versions. Li author annotations are counted once; the overlapping scPlantDB annotation version is excluded. The Cell Expression preview currently exposes the validated full Wang atlas; other imported versions remain separate."
    )
    result["qtl_note"] = (
        "Source-reported QTL names are retained on GWAS records; this export does not independently identify linkage-QTL studies."
    )
    return result


# Return all regional genes and overlapping source-observed SNP positions without inferring allele orientation.
def locus_genes(selector, window=25000):
    data = region(selector, window)
    with connect() as db:
        aliases = db.execute(
            "SELECT gene_id,name FROM alias WHERE gene_id IN (SELECT id FROM gene WHERE chromosome=? AND start<=? AND end>=?) ORDER BY name",
            (data["gene"]["chromosome"], data["end"], data["start"]),
        ).fetchall()
        variants = [
            dict(r)
            for r in db.execute(
                "SELECT DISTINCT position,effect_allele,other_allele FROM association WHERE assembly_id=? AND chromosome=? AND position BETWEEN ? AND ? ORDER BY position",
                (
                    data["gene"]["assembly_id"],
                    data["gene"]["chromosome"],
                    data["start"],
                    data["end"],
                ),
            )
        ]
    alias_map = {}
    for row in aliases:
        alias_map.setdefault(row["gene_id"], set()).add(row["name"])
    genes = []
    for gene in data["neighbors"]:
        genes.append(
            {
                **gene,
                "aliases": sorted(alias_map.get(gene["id"], set())),
                "variants": [
                    v for v in variants if gene["start"] <= v["position"] <= gene["end"]
                ],
                "reference_allele": None,
                "alternate_allele": None,
            }
        )
    return {
        "genes": genes,
        "start": data["start"],
        "end": data["end"],
        "chromosome": data["gene"]["chromosome"],
        "selector": selector,
    }


# Resolve an identifier or symbol to matching genes without choosing among ambiguous matches.
def search(query):
    query = query.strip()
    if not query or len(query) > 120:
        raise ValueError("Enter a gene identifier or symbol (maximum 120 characters).")
    locus = parse_locus(query)
    if locus:
        return {
            "matches": [],
            "match_type": "locus",
            "selector": canonical_locus(locus),
            "truncated": False,
        }
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
    variant_query="",
    suggestions=False,
    sort="p_asc",
):
    if sort not in {"p_asc", "p_desc"}:
        raise ValueError("Select ascending or descending p-value order.")
    if not 0 <= window <= 1000000:
        raise ValueError("Window must be between 0 and 1,000,000 bp on each side.")
    if page < 1 or page > 100000:
        raise ValueError("Invalid results page.")
    if len(variant_query) > 100:
        raise ValueError("Variant search must contain at most 100 characters.")
    variant_query = "".join(variant_query.lower().replace(",", "").split())
    if len(trait) > 200:
        raise ValueError("Trait filter too long.")
    if len(study) > 100 or (publication and not publication.isdigit()):
        raise ValueError("Select a valid study or publication.")
    if not math.isfinite(maximum_p) or not 0 < maximum_p <= 1:
        raise ValueError("Maximum p-value must be greater than zero and at most one.")
    if position_filter not in {"", "overlap", "upstream", "downstream"}:
        raise ValueError("Select a recognized position filter.")
    with connect() as db:
        locus = parse_locus(gene_id)
        if locus:
            chromosome, locus_start, locus_end = locus
            gene = {
                "id": "",
                "assembly_id": "Os-Nipponbare-Reference-IRGSP-1.0",
                "chromosome": chromosome,
                "start": locus_start,
                "end": locus_end,
                "strand": "+",
                "symbol": "Selected region",
                "description": "Genomic interval search",
            }
        else:
            row = db.execute("SELECT * FROM gene WHERE id=?", (gene_id,)).fetchone()
            if not row:
                raise ValueError(
                    "Select a recognized RAP gene ID or enter chr3:16729501-16735109."
                )
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
        if re.fullmatch(r"\d+:\d+", variant_query):
            where += " AND a.id=?"
            args.append(variant_query)
        elif variant_query:
            where += " AND (instr(lower('chr'||a.chromosome||':'||a.position),?)>0 OR instr(lower(a.id),?)>0)"
            args.extend([variant_query, variant_query])
        table = "association" if include_repeats else "association_display"
        join = (
            " FROM "
            + table
            + " a JOIN trait t ON t.id=a.trait_id JOIN study s ON s.id=a.study_id LEFT JOIN publication p ON p.pmid=s.publication_id WHERE "
            + where
        )
        total = db.execute("SELECT count(*)" + join, args).fetchone()[0]
        if suggestions:
            matches = [
                dict(row)
                for row in db.execute(
                    "SELECT a.id,a.chromosome,a.position,t.label AS trait"
                    + join
                    + " ORDER BY a.position,a.id LIMIT 8",
                    args,
                )
            ]
            return {"suggestions": matches, "total": total}

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
            + " AS source_record_count,a.id,a.chromosome,a.position,a.p_value,a.p_value_text,a.effect_allele,a.effect_allele_frequency,a.other_allele,a.effect_size,a.r2_percent,a.reported_location,a.source_reference,a.source_row,json_extract(a.raw_json,'$.QTL') AS qtl_label,json_extract(a.raw_json,'$.\"Sampling condition\"') AS experimental_condition,t.label AS trait,t.ontology_id,s.source_study_id,s.model,s.sample_size,s.population,s.environment,s.tissue,p.pmid,p.title,p.url AS publication_url"
        )

        # Add gene-relative distances, plot values and source links to association records.
        def enrich(rows):
            result = []
            for row in rows:
                a = dict(row)
                a["relation"], a["distance_bp"] = relation(gene, a["position"])
                if locus:
                    a["relation"] = (
                        "within searched interval"
                        if a["distance_bp"] == 0
                        else (
                            "lower coordinate"
                            if a["position"] < gene["start"]
                            else "higher coordinate"
                        )
                    )
                a["evidence_type"] = "GWAS"
                a["reference_allele"] = None
                a["alternate_allele"] = None
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

        top_traits = [
            dict(row)
            for row in db.execute(
                "SELECT t.label AS trait,count(*) AS count"
                + join
                + " GROUP BY t.label ORDER BY count DESC,t.label",
                args,
            )
        ]
        smallest_p = db.execute("SELECT MIN(a.p_value)" + join, args).fetchone()[0]
        trait_options = [
            r[0]
            for r in db.execute(
                "SELECT DISTINCT t.label FROM association a JOIN trait t ON t.id=a.trait_id WHERE "
                + base_where
                + " ORDER BY t.label",
                base_args,
            )
        ]
        order = (
            " ORDER BY a.p_value"
            + (" DESC" if sort == "p_desc" else " ASC")
            + ",a.position,a.id"
        )

        if export:
            return enrich(db.execute(cols + join + order, args))
        associations = enrich(
            db.execute(
                cols + join + order + " LIMIT ? OFFSET ?",
                args + [PAGE_SIZE, (page - 1) * PAGE_SIZE],
            )
        )
        plot = enrich(
            db.execute(
                cols + join + " ORDER BY a.p_value,a.position,a.id LIMIT 500", args
            )
        )
        aliases = [
            dict(r)
            for r in db.execute(
                "SELECT name,kind FROM alias WHERE gene_id=? ORDER BY kind,name",
                (gene_id,),
            )
        ]
        return {
            "gene": gene,
            "selector": gene_id,
            "selection_type": "locus" if locus else "gene",
            "sort": sort,
            "trait_options": trait_options,
            "summary": {
                "records": total,
                "smallest_p": smallest_p,
                "top_traits": (
                    [r for r in top_traits if r["count"] == top_traits[0]["count"]]
                    if top_traits
                    else []
                ),
            },
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
            if parsed.path == "/api/expression/sources":
                return self.send_json(source_coverage(q.get("gene", "")))
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
                        "organ",
                        "source_cell_type",
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
            if parsed.path == "/api/overview":
                return self.send_json(overview())
            if parsed.path in {"/api/genes", "/api/genes/export"}:
                data = locus_genes(q.get("gene", ""), int(q.get("window", 25000)))
                if parsed.path.endswith("/export"):
                    out = io.StringIO()
                    fields = [
                        "id",
                        "symbol",
                        "chromosome",
                        "start",
                        "end",
                        "strand",
                        "aliases",
                        "description",
                        "assembly_id",
                        "snp_positions",
                        "reference_allele",
                        "alternate_allele",
                    ]
                    writer = csv.DictWriter(
                        out, fieldnames=fields, extrasaction="ignore"
                    )
                    writer.writeheader()
                    for gene in data["genes"]:
                        writer.writerow(
                            {
                                **gene,
                                "aliases": "; ".join(gene["aliases"]),
                                "snp_positions": "; ".join(
                                    f"chr{gene['chromosome']}:{v['position']}"
                                    for v in gene["variants"]
                                ),
                            }
                        )
                    return self.respond(
                        out.getvalue().encode(),
                        "text/csv; charset=utf-8",
                        attachment="ricelens_region_genes.csv",
                    )
                return self.send_json(data)
            if parsed.path == "/api/search":
                return self.send_json(search(q.get("q", "")))
            if parsed.path in {
                "/api/region",
                "/api/export",
                "/api/variant-suggestions",
            }:
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
                    variant_query=q.get("variant_query", ""),
                    suggestions=parsed.path == "/api/variant-suggestions",
                    sort=q.get("sort", "p_asc"),
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
            files = {
                "/": "index.html",
                "/variation": "index.html",
                "/locus": "index.html",
                "/homolog": "index.html",
                "/download": "index.html",
                "/tutorial": "index.html",
                "/contact": "index.html",
                "/gwas": "index.html",
                "/expression": "index.html",
                "/app.js": "app.js",
                "/style.css": "style.css",
            }
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
        raise SystemExit("Run python3 -m pipelines.build_database first.")
    print(f"RiceLENS: http://127.0.0.1:{args.port}", flush=True)
    ThreadingHTTPServer(("127.0.0.1", args.port), Handler).serve_forever()

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

ROOT = Path(__file__).resolve().parent
DB = ROOT / "data/rice.sqlite"
PAGE_SIZE = 8


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
            "SELECT DISTINCT g.* FROM gene g JOIN alias a ON a.gene_id=g.id WHERE a.name=? COLLATE NOCASE ORDER BY g.id",
            (query,),
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
def region(gene_id, window=25000, trait="", page=1, export=False):
    if not 0 <= window <= 1000000:
        raise ValueError("Window must be between 0 and 1,000,000 bp on each side.")
    if page < 1 or page > 100000:
        raise ValueError("Invalid results page.")
    if len(trait) > 200:
        raise ValueError("Trait filter too long.")
    with connect() as db:
        row = db.execute("SELECT * FROM gene WHERE id=?", (gene_id,)).fetchone()
        if not row:
            raise ValueError("Select a recognized RAP gene ID from the search results.")
        gene = dict(row)
        start = max(1, gene["start"] - window)
        end = gene["end"] + window
        genes = [
            dict(r)
            for r in db.execute(
                "SELECT * FROM gene WHERE assembly_id=? AND chromosome=? AND start<=? AND end>=? ORDER BY start,id",
                (gene["assembly_id"], gene["chromosome"], end, start),
            )
        ]
        args = [gene["assembly_id"], gene["chromosome"], start, end]
        where = "a.assembly_id=? AND a.chromosome=? AND a.position BETWEEN ? AND ?"
        if trait:
            where += " AND instr(lower(t.label),lower(?))>0"
            args.append(trait)
        join = (
            " FROM association a JOIN trait t ON t.id=a.trait_id JOIN study s ON s.id=a.study_id LEFT JOIN publication p ON p.pmid=s.publication_id WHERE "
            + where
        )
        total = db.execute("SELECT count(*)" + join, args).fetchone()[0]
        cols = "SELECT a.id,a.chromosome,a.position,a.p_value,a.p_value_text,a.effect_allele,a.effect_allele_frequency,a.other_allele,a.effect_size,a.r2_percent,a.reported_location,a.source_reference,a.source_row,t.label AS trait,t.ontology_id,s.source_study_id,s.model,s.sample_size,s.population,s.environment,s.tissue,p.pmid,p.title,p.url AS publication_url"

        # Add gene-relative distances, plot values and source links to association records.
        def enrich(rows):
            result = []
            for row in rows:
                a = dict(row)
                a["relation"], a["distance_bp"] = relation(gene, a["position"])
                a["neg_log10_p"] = -math.log10(a["p_value"])
                a["assembly"] = gene["assembly_id"]
                a["source_url"] = "https://ngdc.cncb.ac.cn/gwas/download"
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
            "interpretation": "Regional association evidence does not establish gene causality. Overlap refers to the gene span, including introns. Upstream/downstream follow gene strand.",
        }


class Handler(BaseHTTPRequestHandler):
    # Route a browser request to the appropriate read-only API or static resource.
    def do_GET(self):
        try:
            parsed = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(parsed.query).items()}
            if parsed.path == "/api/search":
                return self.send_json(search(q.get("q", "")))
            if parsed.path in {"/api/region", "/api/export"}:
                result = region(
                    q.get("gene", ""),
                    int(q.get("window", 25000)),
                    q.get("trait", ""),
                    int(q.get("page", 1)),
                    parsed.path == "/api/export",
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
    print(f"Rice Gene Explorer: http://127.0.0.1:{args.port}", flush=True)
    ThreadingHTTPServer(("127.0.0.1", args.port), Handler).serve_forever()

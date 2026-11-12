"""Build a compact read-only expression database from validated baseline summaries."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import sqlite3
from collections import defaultdict
from contextlib import closing
from pathlib import Path

GENE_DB = Path(__file__).resolve().parents[1] / "data/databases/rice.sqlite"

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/expression"
OUTPUT = ROOT / "data/databases/expression.sqlite"


# Hash a source file incrementally without loading its complete contents into memory.
def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


# Load validated per-sample summaries and retain unmapped features explicitly.
def build_expression(study="robertson"):
    wang = study == "wang"
    li = study in {"li", "li_author"}
    suffix = "" if study == "robertson" else "_" + study
    raw = ROOT / "data/raw" / ("expression" if study == "robertson" else study)
    output = ROOT / f"data/databases/expression{suffix}.sqlite"
    output.parent.mkdir(parents=True, exist_ok=True)
    inputs_path = ROOT / f"data/expression{suffix}_inputs.json"
    inputs = json.loads(inputs_path.read_text())
    for source in inputs["sources"]:
        if file_hash(raw / source["filename"]) != source["sha256"]:
            raise ValueError(f"Source checksum mismatch: {source['filename']}")
    report = json.loads((raw / "extraction.json").read_text())
    expected_samples = {
        "robertson": ["Well.Watered.Rep1", "Well.Watered.Rep2", "Well.Watered.Rep3"],
        "wang": ["NIP_G_1", "NIP_G_2"],
        "li": ["GSM6346349", "GSM6346350"],
        "li_author": ["GSM6346349", "GSM6346350"],
    }[study]
    if report["eligible_samples"] != expected_samples:
        raise ValueError("Unexpected baseline sample selection.")
    features = list(csv.DictReader((raw / "features.tsv").open(), delimiter="\t"))
    with closing(sqlite3.connect(f"file:{GENE_DB}?mode=ro", uri=True)) as db:
        recognized = {r[0] for r in db.execute("SELECT id FROM gene")}
        aliases = defaultdict(set)
        for name, gene_id in db.execute(
            "SELECT name,gene_id FROM alias WHERE kind IN ('MSU locus','MSU transcript')"
        ):
            aliases[name.lower()].add(gene_id)
    mapped = {}
    candidates = defaultdict(list)
    for feature in features:
        original = feature["feature_id"]
        matches = (
            {original}
            if original in recognized
            else (aliases.get(original.lower(), set()) if (wang or li) else set())
        )
        if len(matches) == 1:
            candidates[next(iter(matches))].append(original)
    for canonical, originals in candidates.items():
        if len(originals) == 1:
            mapped[originals[0]] = canonical
    temporary = output.with_suffix(".building.sqlite")
    temporary.unlink(missing_ok=True)
    db = sqlite3.connect(temporary)
    try:
        db.executescript(
            """
            PRAGMA foreign_keys=ON;
            CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
            CREATE TABLE feature(feature_id TEXT PRIMARY KEY,gene_id TEXT UNIQUE,mapping_status TEXT NOT NULL);
            CREATE TABLE cell_group(cell_type TEXT PRIMARY KEY,tissue_system TEXT NOT NULL,n_cells INTEGER NOT NULL);
            CREATE TABLE sample_summary(feature_id TEXT REFERENCES feature(feature_id),sample TEXT,cell_type TEXT REFERENCES cell_group(cell_type),n_cells INTEGER CHECK(n_cells>0),n_expressing INTEGER CHECK(n_expressing>=0 AND n_expressing<=n_cells),sum_log REAL CHECK(sum_log>=0),PRIMARY KEY(feature_id,sample,cell_type));
            CREATE TABLE summary(gene_id TEXT,cell_type TEXT REFERENCES cell_group(cell_type),n_cells INTEGER,n_expressing INTEGER,mean_log REAL,n_samples INTEGER,PRIMARY KEY(gene_id,cell_type));
        """
        )
        for feature in features:
            original = feature["feature_id"]
            db.execute(
                "INSERT INTO feature VALUES(?,?,?)",
                (
                    original,
                    mapped.get(original),
                    (
                        "exact RAP locus"
                        if original in mapped
                        else "unmapped or ambiguous"
                    ),
                ),
            )
        for group in report["groups"]:
            db.execute(
                "INSERT INTO cell_group VALUES(?,?,?)",
                (group["cell_type"], group["tissue_system"], group["n_cells"]),
            )
        expected_cells = {
            (g["sample"], g["cell_type"]): g["n_cells"] for g in report["sample_groups"]
        }
        with gzip.open(raw / "summaries.tsv.gz", "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                n_cells, n_detected = int(row["n_cells"]), int(row["n_expressing"])
                value = float(row["sum_log"])
                if (
                    expected_cells.get((row["sample"], row["cell_type"])) != n_cells
                    or not math.isfinite(value)
                    or value < 0
                ):
                    raise ValueError(
                        "Expression summary does not match eligible cell counts."
                    )
                db.execute(
                    "INSERT INTO sample_summary VALUES(?,?,?,?,?,?)",
                    (
                        row["feature_id"],
                        row["sample"],
                        row["cell_type"],
                        n_cells,
                        n_detected,
                        value,
                    ),
                )
        actual = db.execute("SELECT count(*) FROM sample_summary").fetchone()[0]
        expected_rows = sum(
            report["features_per_sample"][sample] for sample, group in expected_cells
        )
        if actual != expected_rows:
            raise ValueError("Missing feature-by-sample-group rows.")
        db.execute(
            """INSERT INTO summary SELECT f.gene_id,s.cell_type,SUM(n_cells),SUM(n_expressing),SUM(sum_log)/SUM(n_cells),COUNT(*) FROM sample_summary s JOIN feature f USING(feature_id) WHERE f.gene_id IS NOT NULL GROUP BY f.gene_id,s.cell_type"""
        )
        dataset = {
            "id": "robertson_2026_control",
            "citation": "Robertson et al. (2026)",
            "publication_url": "https://doi.org/10.1111/nph.71378",
            "repository_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313726",
            "organ": "Developing third leaf (~6 cm)",
            "condition": "Well-watered controls only",
            "assay": "single-nucleus RNA-seq",
            "reference": "IRGSP1.0",
            "n_cells": report["n_cells"],
            "n_samples": len(expected_samples),
            "n_features": len(features),
            "n_mapped_features": len(mapped),
            "metric_description": "Mean log1p(10,000 × RNA UMI count / processed-matrix RNA library size), including zero-count nuclei. Detection = RNA count > 0. Replicates are pooled with nuclei-weighted means; ranks use that mean within this study.",
            "coverage_note": "All available processed control-matrix features are included, not only markers or variable genes. Features absent from a replicate are excluded from its denominator, not imputed zero. Source processing and annotation filtering limit coverage.",
            "source_sha256": file_hash(
                raw
                / {
                    "robertson": "GSE313726_OsLeafStressIntegrated.rds",
                    "wang": "GSE273875_LTP_NIP_label_240622.rds.gz",
                    "li": "SRP386976.rds.gz",
                    "li_author": "mmc3.xlsx",
                }[study]
            ),
            "extraction": report,
        }
        if wang:
            dataset.update(
                {
                    "id": "wang_2025_normal_root",
                    "citation": "Wang et al. (2025)",
                    "publication_url": "https://doi.org/10.1038/s41586-025-09251-0",
                    "repository_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273875",
                    "organ": "Root (normal Nipponbare subset)",
                    "condition": "Normal wild-type Nipponbare only; mutants excluded",
                    "assay": "single-cell RNA-seq",
                    "cell_unit": "cells",
                    "reference": "Source alignment MSU7; deposited features already RAP-labelled; exact ID matching, no coordinate transfer",
                    "metric_description": "Mean log1p(10,000 × RNA UMI count / retained RNA-count library size), including zero-count cells. Detection = RNA count > 0. Two normal replicates are pooled with cell-weighted means; ranks are within this study only.",
                    "coverage_note": "All features in the original RNA count assay are retained, not only markers or variable genes. Unmapped and ambiguous gene mappings are unavailable in the RAP search. This is the paper's normal-root experiment, not the full eight-organ atlas.",
                }
            )
        if li:
            dataset.update(
                {
                    "id": "li_2023_pistil_2n_scplantdb",
                    "citation": "Li et al. (2023) · scPlantDB-curated",
                    "publication_url": "https://doi.org/10.1016/j.devcel.2023.03.004",
                    "repository_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208432",
                    "organ": "Pistil · diploid (2N) nuclei",
                    "condition": "Wild-type Nipponbare, before flowering; 1N library excluded",
                    "assay": "single-nucleus RNA-seq",
                    "reference": "IRGSP1.0; exact gene IDs or unambiguous frozen MSU crosswalk",
                    "metric_description": "Mean log1p(10,000 × original GEO RNA count / per-nucleus library size), including zeros. Detection = count > 0. Two 2N replicates are pooled with nuclei-weighted means; ranks are within this curated dataset only.",
                    "coverage_note": "All features in the original GEO matrices are retained for nuclei with matched scPlantDB annotations. No marker-only or variable-gene restriction; absent features are unavailable rather than zero.",
                    "annotation_source_name": "scPlantDB · SRP386976",
                    "annotation_source_url": "https://biobigdata.nju.edu.cn/scplantdb/",
                    "annotation_publication_url": "https://doi.org/10.1093/nar/gkad706",
                    "annotation_note": "Labels and retained cells come from the scPlantDB-curated version, not a verified copy of the authors' final barcode annotation. Sampling is enriched for 2N nuclei; expression results do not estimate natural cell-type proportions.",
                }
            )
        if study == "li_author":
            dataset.update(
                {
                    "id": "li_2023_pistil_2n_author",
                    "citation": "Li et al. (2023) · author annotations",
                    "annotation_source_name": "Author Table S2 · Rice_Pistil_index",
                    "annotation_source_url": "https://doi.org/10.1016/j.devcel.2023.03.004",
                    "annotation_publication_url": "https://doi.org/10.1016/j.devcel.2023.03.004",
                    "annotation_note": "Exact original barcodes and author cell-type labels from Table S2; 96 stamen-residue nuclei excluded. C12 and C13 remain author cluster labels. This shares original libraries with the scPlantDB version and is not independent replication. Sampling is enriched for 2N nuclei.",
                    "metric_description": "Mean log1p(10,000 × original GEO RNA count / per-nucleus library size), including zeros. Detection = count > 0. Two 2N replicates are pooled with nuclei-weighted means; ranks are within this author-annotated dataset only.",
                    "coverage_note": "All 38,726 original count-matrix features are retained for author-annotated pistil nuclei. No marker or variable-gene restriction; source mappings remain explicit. This is the published author cell selection, rather than the scPlantDB re-curation.",
                }
            )
        db.execute("INSERT INTO metadata VALUES('dataset',?)", (json.dumps(dataset),))
        if (
            db.execute("PRAGMA integrity_check").fetchone()[0] != "ok"
            or db.execute("PRAGMA foreign_key_check").fetchall()
        ):
            raise ValueError("Expression database integrity check failed.")
        db.commit()
    finally:
        db.close()
    temporary.replace(output)
    (ROOT / f"data/expression{suffix}_manifest.json").write_text(
        json.dumps(dataset, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "cells": report["n_cells"],
                "features": len(features),
                "mapped": len(mapped),
                "sample_summary_rows": actual,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--study", choices=["robertson", "wang", "li", "li_author"], default="robertson"
    )
    build_expression(parser.parse_args().study)

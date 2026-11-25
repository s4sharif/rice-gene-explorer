"""Independently verify full-atlas calculations from disk-backed original counts."""

import csv
import json
import math
import sqlite3
from collections import Counter, defaultdict
from contextlib import closing
from pathlib import Path

from pipelines.build_expression import file_hash
from pipelines.prepare_wang_atlas import ROOT, group_cells, open_counts

TARGETS = {"Os03g0407400", "Os12g0528400", "Os01g0883800", "Os01g0100100"}


# Recompute selected genes cell by cell with scalar log1p rather than the summarizer's block reduction.
def independent_checks(spec, arrays, features, aligned, targets):
    import numpy as np

    wanted = [(features.index(feature), feature) for feature in targets]
    sums = defaultdict(float)
    detected = Counter()
    cells = Counter()
    for column, key in enumerate(aligned):
        sample, context = key[:2]
        cells[(sample, context)] += 1
        left, right = int(arrays["p"][column]), int(arrays["p"][column + 1])
        rows = arrays["i"][left:right]
        values = arrays["x"][left:right]
        library = float(values.sum(dtype=np.float64))
        if not math.isfinite(library) or library <= 0:
            raise ValueError("Invalid original count library.")
        for index, feature in wanted:
            offset = int(np.searchsorted(rows, index))
            count = (
                float(values[offset])
                if offset < len(rows) and rows[offset] == index
                else 0.0
            )
            group = (feature, sample, context)
            detected[group] += int(count > 0)
            sums[group] += math.log1p(10000 * count / library)
        if column % 25000 == 0:
            print(
                f"Independently checked {column:,}/{spec['shape'][1]:,} cells",
                flush=True,
            )
    return cells, detected, sums


# Verify provenance, all-feature coverage, group bounds and independently recalculated targets before enabling serving.
def validate(raw=None, database=None):
    raw = Path(raw or ROOT / "data/raw/wang_atlas")
    database = Path(database or ROOT / "data/databases/expression_wang_atlas.sqlite")
    report = json.loads((raw / "extraction.json").read_text())
    inputs = json.loads((ROOT / "data/expression_wang_atlas_inputs.json").read_text())
    for source in inputs["sources"]:
        if file_hash(raw / source["filename"]) != source["sha256"]:
            raise ValueError("Original full-atlas source checksum failed.")
    for name, digest in report["streamed_files"].items():
        if file_hash(raw / name) != digest:
            raise ValueError(f"Extracted count/annotation checksum failed: {name}")
    spec, arrays = open_counts(raw)
    with (raw / "features.tsv").open() as handle:
        features = [r["feature_id"] for r in csv.DictReader(handle, delimiter="\t")]
    with (raw / "barcodes.tsv").open() as handle:
        barcodes = [r["barcode"] for r in csv.DictReader(handle, delimiter="\t")]
    with (raw / "cell_metadata.tsv").open() as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    aligned, expected_groups = group_cells(metadata, barcodes)
    with closing(sqlite3.connect(f"file:{database}?mode=ro", uri=True)) as db:
        db.row_factory = sqlite3.Row
        if (
            db.execute("PRAGMA integrity_check").fetchone()[0] != "ok"
            or db.execute("PRAGMA foreign_key_check").fetchall()
        ):
            raise ValueError("Full-atlas database integrity failed.")
        if db.execute("SELECT count(*) FROM feature").fetchone()[0] != len(features):
            raise ValueError("An original RNA feature was omitted.")
        if db.execute(
            "SELECT count(*) FROM summary WHERE n_cells<=0 OR n_expressing<0 OR n_expressing>n_cells OR mean_log<0 OR mean_log>? OR n_samples NOT BETWEEN 1 AND 2",
            (math.log1p(10000) + 1e-10,),
        ).fetchone()[0]:
            raise ValueError("Invalid expression summary bounds.")
        disagreements = db.execute(
            "SELECT count(*) FROM summary a JOIN (SELECT f.gene_id,s.cell_type,SUM(s.n_cells) AS n,SUM(s.n_expressing) AS d,SUM(s.sum_log)/SUM(s.n_cells) AS m,COUNT(*) AS samples FROM sample_summary s JOIN feature f USING(feature_id) WHERE f.gene_id IS NOT NULL GROUP BY f.gene_id,s.cell_type) b USING(gene_id,cell_type) WHERE a.n_cells!=b.n OR a.n_expressing!=b.d OR abs(a.mean_log-b.m)>1e-10 OR a.n_samples!=b.samples"
        ).fetchone()[0]
        if disagreements:
            raise ValueError(
                "Study-level aggregation differs from original sample summaries."
            )
        mapped = db.execute(
            "SELECT count(*) FROM feature WHERE gene_id IS NOT NULL"
        ).fetchone()[0]
        contexts = db.execute("SELECT count(*) FROM cell_group").fetchone()[0]
        if (
            db.execute("SELECT count(*) FROM summary").fetchone()[0]
            != mapped * contexts
        ):
            raise ValueError(
                "Mapped zero-expression gene/context combinations were lost."
            )
        if db.execute("SELECT SUM(n_cells) FROM cell_group").fetchone()[0] != len(
            barcodes
        ):
            raise ValueError(
                "Cells are double-counted or missing from the atlas contexts."
            )
        target_rows = db.execute(
            "SELECT feature_id,gene_id FROM feature WHERE gene_id IN (?,?,?,?)",
            tuple(sorted(TARGETS)),
        ).fetchall()
        targets = {row["feature_id"] for row in target_rows}
        if len(targets) < 2:
            raise ValueError(
                "Too few independent validation genes are measured in this atlas."
            )
        cells, detected, sums = independent_checks(
            spec, arrays, features, aligned, targets
        )
        if cells != Counter(
            {(key[0], key[1]): value for key, value in expected_groups.items()}
        ):
            raise ValueError("Independent annotation grouping disagrees.")
        checks = 0
        for feature in targets:
            for (sample, context), n_cells in cells.items():
                row = db.execute(
                    "SELECT n_cells,n_expressing,sum_log FROM sample_summary WHERE feature_id=? AND sample=? AND cell_type=?",
                    (feature, sample, context),
                ).fetchone()
                key = (feature, sample, context)
                if (
                    row is None
                    or (row[0], row[1]) != (n_cells, detected[key])
                    or not math.isclose(row[2], sums[key], rel_tol=1e-10, abs_tol=1e-8)
                ):
                    raise ValueError(f"Independent count calculation differs: {key}")
                checks += 1
        dataset = json.loads(
            db.execute("SELECT value FROM metadata WHERE key='dataset'").fetchone()[0]
        )
    result = {
        "integrity": "pass",
        "source_checksums": "pass",
        "all_deposited_cells_and_features": "pass",
        "source_sha256": dataset["source_sha256"],
        "n_cells": len(barcodes),
        "n_features": len(features),
        "n_organs": report["n_organs"],
        "n_contexts": contexts,
        "n_mapped_features": mapped,
        "independent_gene_sample_context_checks": checks,
        "genes": [r["gene_id"] for r in target_rows],
        "cell_count_discrepancy": report["cell_count_discrepancy"],
        "unknown_labels": "Original unknown cell labels retained; not relabelled or omitted",
    }
    (ROOT / "data/expression_wang_atlas_validation.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2), flush=True)
    return result


if __name__ == "__main__":
    validate()

"""Validate expression summaries against independent calculations from source MEX counts."""

import argparse
import csv
import gzip
import json
import math
import re
import sqlite3
from collections import defaultdict
from contextlib import closing
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data/raw/expression"
TARGETS = {"Os03g0407400", "Os12g0528400", "Os07g0141400", "Os01g0100100"}


# Recalculate selected genes directly from each original control count matrix.
def main(study="robertson"):
    wang = study == "wang"
    li = study == "li"
    suffix = "" if study == "robertson" else "_" + study
    raw = ROOT / "data/raw" / ("expression" if study == "robertson" else study)
    if not wang:
        selection = json.loads((raw / "selection.json").read_text())
        with (
            raw / ("matched_cells.tsv" if li else "cell_metadata.tsv")
        ).open() as handle:
            metadata = list(csv.DictReader(handle, delimiter="\t"))
    checks = 0
    with closing(sqlite3.connect(ROOT / f"data/expression{suffix}.sqlite")) as db:
        db.row_factory = sqlite3.Row
        if db.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise ValueError("Expression database integrity failed.")
        invalid = db.execute(
            "SELECT count(*) FROM summary WHERE n_cells<=0 OR n_expressing<0 OR n_expressing>n_cells OR mean_log<0 OR n_samples NOT BETWEEN 1 AND 3"
        ).fetchone()[0]
        if invalid:
            raise ValueError("Invalid expression summary bounds.")
        target_rows = db.execute(
            "SELECT feature_id,gene_id FROM feature WHERE gene_id IN (?,?,?,?)",
            tuple(sorted(TARGETS)),
        ).fetchall()
        targets = {r[0] for r in target_rows}
        target_genes = {r[1] for r in target_rows}
        for index in range(1, 3 if (wang or li) else 4):
            if li:
                sample = ["GSM6346349", "GSM6346350"][index - 1]
                prefix = raw / selection["prefixes"][sample]
                with gzip.open(str(prefix) + "_features.tsv.gz", "rt") as handle:
                    features = [line.split("\t")[0] for line in handle]
                with gzip.open(str(prefix) + "_barcodes.tsv.gz", "rt") as handle:
                    barcodes = [line.strip() for line in handle]
                groups = {
                    row["count_barcode"]: row["cell_type"]
                    for row in metadata
                    if row["sample"] == sample
                }
            elif wang:
                sample = f"NIP_G_{index}"
                prefix = raw / sample
                with gzip.open(str(prefix) + ".features.tsv.gz", "rt") as handle:
                    features = [line.strip() for line in handle]
                with (raw / f"{sample}.cells.tsv").open() as handle:
                    cells = list(csv.DictReader(handle, delimiter="\t"))
                barcodes = [row["barcode"] for row in cells]
                groups = {row["barcode"]: row["cell_type"] for row in cells}
                if len(groups) != len(cells):
                    raise ValueError("Duplicate exported cell barcode.")
            else:
                sample = f"Well.Watered.Rep{index}"
                prefix = raw / f"GSM{9373071 + index}_{sample}"
                with gzip.open(str(prefix) + ".genes.tsv.gz", "rt") as handle:
                    features = [line.split("\t")[0] for line in handle]
                with gzip.open(str(prefix) + ".barcodes.tsv.gz", "rt") as handle:
                    barcodes = [line.strip() for line in handle]
                groups = {}
                for row in metadata:
                    if (
                        row[selection["sample_column"]]
                        != selection["sample_values"][sample]
                    ):
                        continue
                    barcode = (
                        row["barcode"]
                        if selection["barcode_column"] == "rownames"
                        else row[selection["barcode_column"]]
                    )
                    barcode = re.sub(selection["barcode_remove_prefix"], "", barcode)
                    if barcode in groups:
                        raise ValueError("Duplicate annotation barcode.")
                    groups[barcode] = row["CellAnnotation"]
            libraries = [0] * len(barcodes)
            values = defaultdict(dict)
            target_indices = {
                i + 1: feature
                for i, feature in enumerate(features)
                if feature in targets
            }
            with gzip.open(
                str(prefix) + ("_matrix.mtx.gz" if li else ".matrix.mtx.gz"), "rt"
            ) as handle:
                dimensions = None
                entries = 0
                for line in handle:
                    if line.startswith("%"):
                        continue
                    if dimensions is None:
                        dimensions = tuple(map(int, line.split()))
                        continue
                    parts = line.split()
                    gene_index, column = int(parts[0]), int(parts[1])
                    count = float(parts[2])
                    if not math.isfinite(count) or count < 0 or count != int(count):
                        raise ValueError(
                            "RNA counts must be finite nonnegative integers."
                        )
                    libraries[column - 1] += count
                    if gene_index in target_indices:
                        values[target_indices[gene_index]][column - 1] = count
                    entries += 1
            if dimensions != (len(features), len(barcodes), entries):
                raise ValueError("MEX dimensions or entry count failed.")
            group_cells = defaultdict(list)
            for column, barcode in enumerate(barcodes):
                if barcode in groups:
                    group_cells[groups[barcode]].append(column)
            if sum(map(len, group_cells.values())) != len(groups):
                raise ValueError("Annotation join lost cells.")
            for gene in target_indices.values():
                for group, columns in group_cells.items():
                    detected = sum(values[gene].get(col, 0) > 0 for col in columns)
                    total_log = sum(
                        math.log1p(10000 * values[gene].get(col, 0) / libraries[col])
                        for col in columns
                    )
                    row = db.execute(
                        "SELECT * FROM sample_summary WHERE feature_id=? AND sample=? AND cell_type=?",
                        (gene, sample, group),
                    ).fetchone()
                    if (
                        not row
                        or (row["n_cells"], row["n_expressing"])
                        != (len(columns), detected)
                        or not math.isclose(
                            row["sum_log"], total_log, rel_tol=1e-10, abs_tol=1e-9
                        )
                    ):
                        raise ValueError(
                            f"Independent expression calculation disagrees: {gene} {sample} {group}"
                        )
                    checks += 1
            print(f"Validated independent MEX calculations: {sample}", flush=True)
        report = {
            "integrity": "pass",
            "summary_bounds": "pass",
            "independent_gene_sample_group_checks": checks,
            "genes": sorted(target_genes),
            "requested_genes": sorted(TARGETS),
            "source_features_checked": sorted(targets),
            "normalization": "log1p CP10K recalculated independently of R summaries from "
            + (
                "uncorrected RNA counts exported from the source object"
                if wang
                else "original source MEX counts"
            ),
        }
    (ROOT / f"data/expression{suffix}_validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--study", choices=["robertson", "wang", "li"], default="robertson"
    )
    main(parser.parse_args().study)

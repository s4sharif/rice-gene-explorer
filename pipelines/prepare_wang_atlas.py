"""Summarize every deposited RNA feature by original sample, organ and cell label."""

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

from pipelines.build_expression import file_hash

ROOT = Path(__file__).resolve().parents[1]
ORGANS = {
    "Root": "Crown root",
    "ST": "Culm",
    "leaf": "Third leaf",
    "Flag": "Flag leaf",
    "SAM": "Shoot apical meristem",
    "Bud": "Tiller bud",
    "SP": "Young spike",
    "Seed": "Spikelet",
}


# Match each deposited sample to GEO's original normal Nipponbare RNA library record.
def baseline_samples(path):
    source = {}
    sample = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("^SAMPLE = "):
                sample = {"accession": line.strip().split(" = ", 1)[1]}
            elif sample is not None and line.startswith("!Sample_title = "):
                sample["title"] = line.strip().split(" = ", 1)[1]
                source[sample["title"]] = sample
            elif sample is not None and line.startswith(
                "!Sample_characteristics_ch1 = "
            ):
                field = line.strip().split(" = ", 1)[1]
                if ": " in field:
                    key, value = field.split(": ", 1)
                    sample[key.lower()] = value
    sample_codes = {
        "Flag": ("Flag_Leaf", ["Flag_1", "Flag_2"]),
        "leaf": ("Leaf", ["leaf_1", "leaf_2"]),
        "Root": ("Root", ["Root_1", "Root_2"]),
        "SAM": ("SAM", ["SAM_1", "SAM_2"]),
        "SP": ("SP", ["SP_1", "SP_2"]),
        "ST": ("Culm", ["ST_1", "ST_2"]),
        "Bud": ("Tiller_Bud", ["Bud_21", "Bud_22"]),
        "Seed": ("Spikelet", ["Seed_31", "Seed_32"]),
    }
    result = {}
    for tissue, (title, codes) in sample_codes.items():
        for replicate, code in enumerate(codes, 1):
            record = source.get(f"Rice-scRNA-{title}-Rep{replicate}")
            if (
                not record
                or record.get("genotype") != "Nipponbare"
                or record.get("treatment") != "normal"
            ):
                raise ValueError(
                    f"Normal wild-type source metadata could not be verified: {code}"
                )
            result[code] = {**record, "tissue_code": tissue}
    return result


# Open validated disk-backed sparse vectors without materializing the whole RNA matrix.
def open_counts(raw):
    import numpy as np

    specification = json.loads((raw / "matrix.json").read_text())
    arrays = {}
    if specification["assay"] != "RNA" or specification["slot"] != "counts":
        raise ValueError("Only verified original RNA counts may be normalized.")
    for name, spec in specification["arrays"].items():
        path = raw / spec["file"]
        if path.stat().st_size != spec["length"] * np.dtype(spec["dtype"]).itemsize:
            raise ValueError(f"Truncated sparse array: {name}")
        arrays[name] = np.memmap(
            path, mode="r", dtype=spec["dtype"], shape=(spec["length"],)
        )
    n_features, n_cells = specification["shape"]
    if (
        len(arrays["p"]) != n_cells + 1
        or arrays["p"][0] != 0
        or arrays["p"][-1] != len(arrays["x"])
        or len(arrays["i"]) != len(arrays["x"])
        or np.any(np.diff(arrays["p"].astype("int64")) < 0)
    ):
        raise ValueError("Invalid sparse column pointers.")
    return specification, arrays


# Require complete original annotations and keep identical cell labels separate between organs.
def group_cells(rows, barcodes, expected_cells=115395):
    if len(rows) != expected_cells or len(barcodes) != expected_cells:
        raise ValueError("The cell count must match the deposited eight-organ atlas.")
    lookup = {row["barcode"]: row for row in rows}
    if len(lookup) != len(rows) or set(lookup) != set(barcodes):
        raise ValueError("Metadata-to-count-cell join is not one-to-one.")
    if any(not {"sample", "tissue", "cluster_names"}.issubset(row) for row in rows):
        raise ValueError(
            "Original sample, tissue and cluster_names annotations are required."
        )
    observed = {row["tissue"] for row in rows}
    if observed != set(ORGANS):
        raise ValueError(
            f"Expected all eight original tissue codes, found {sorted(observed)}"
        )
    groups = Counter()
    aligned = []
    source_samples = defaultdict(set)
    for barcode in barcodes:
        row = lookup[barcode]
        sample, tissue, label = row["sample"], row["tissue"], row["cluster_names"]
        if not sample or not label or label.lower() in {"na", "nan"}:
            raise ValueError("Missing original cell annotation; do not invent labels.")
        context = ORGANS[tissue] + " | " + label
        key = (sample, context, ORGANS[tissue], label)
        groups[key] += 1
        aligned.append(key)
        source_samples[tissue].add(sample)
    if any(len(samples) != 2 for samples in source_samples.values()):
        raise ValueError(
            f"Expected two deposited samples per organ: {dict(source_samples)}"
        )
    return aligned, groups


# Aggregate sparse blocks while including every annotated zero cell in each denominator.
def summarize_blocks(specification, arrays, group_indices, n_groups, block_cells=512):
    import numpy as np

    n_features, n_cells = specification["shape"]
    detected = np.zeros((n_features, n_groups), dtype=np.int64)
    totals = np.zeros((n_features, n_groups), dtype=np.float64)
    libraries = np.zeros(n_cells, dtype=np.float64)
    for start in range(0, n_cells, block_cells):
        end = min(start + block_cells, n_cells)
        pointers = arrays["p"][start : end + 1].astype("int64")
        left, right = int(pointers[0]), int(pointers[-1])
        counts = arrays["x"][left:right].astype("float64")
        rows = arrays["i"][left:right].astype("int64")
        columns = np.repeat(np.arange(end - start), np.diff(pointers))
        if (
            np.any(~np.isfinite(counts))
            or np.any(counts < 0)
            or np.any(counts != np.floor(counts))
            or np.any(rows < 0)
            or np.any(rows >= n_features)
        ):
            raise ValueError(
                "Counts must be finite nonnegative integers with valid feature indices."
            )
        if np.any((np.diff(columns) == 0) & (np.diff(rows) <= 0)):
            raise ValueError("Duplicate or unsorted features in a count column.")
        sizes = np.bincount(columns, weights=counts, minlength=end - start)
        if np.any(sizes <= 0):
            raise ValueError("Annotated cells must have a positive RNA library size.")
        libraries[start:end] = sizes
        groups = np.asarray(group_indices[start:end])[columns]
        positive = counts > 0
        np.add.at(detected, (rows[positive], groups[positive]), 1)
        values = np.log1p(10000 * counts / sizes[columns])
        np.add.at(totals, (rows, groups), values)
        if start % (block_cells * 40) == 0:
            print(f"Summarized {end:,}/{n_cells:,} cells", flush=True)
    return detected, totals, libraries


# Write reproducible sample summaries, context annotations and source checksums for the full atlas.
def prepare(raw):
    import numpy as np

    raw = Path(raw)
    spec, arrays = open_counts(raw)
    sample_evidence = baseline_samples(raw / "GSE232863_family.soft.gz")
    with (raw / "cell_metadata.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    with (raw / "barcodes.tsv").open() as handle:
        barcodes = [r["barcode"] for r in csv.DictReader(handle, delimiter="\t")]
    with (raw / "features.tsv").open() as handle:
        features = [r["feature_id"] for r in csv.DictReader(handle, delimiter="\t")]
    if len(features) != spec["shape"][0]:
        raise ValueError("Feature names disagree with count dimensions.")
    aligned, counts = group_cells(rows, barcodes)
    if {row["sample"] for row in rows} != set(sample_evidence):
        raise ValueError(
            "Deposited sample names do not match the verified normal-sample metadata."
        )
    keys = sorted(counts)
    group_index = {key: index for index, key in enumerate(keys)}
    indices = [group_index[key] for key in aligned]
    detected, sums, libraries = summarize_blocks(spec, arrays, indices, len(keys))
    expected_libraries = (
        {row["barcode"]: float(row["nCount_RNA"]) for row in rows}
        if "nCount_RNA" in rows[0]
        else None
    )
    if expected_libraries is not None and not np.array_equal(
        libraries, [expected_libraries[b] for b in barcodes]
    ):
        raise ValueError(
            "Computed RNA library sizes disagree with source nCount_RNA metadata."
        )
    with (raw / "library_sizes.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["barcode", "library_size"])
        writer.writerows(zip(barcodes, libraries))
    with gzip.open(raw / "summaries.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["feature_id", "sample", "cell_type", "n_cells", "n_expressing", "sum_log"]
        )
        for index, key in enumerate(keys):
            sample, context, organ, label = key
            for feature_index, feature in enumerate(features):
                writer.writerow(
                    [
                        feature,
                        sample,
                        context,
                        counts[key],
                        int(detected[feature_index, index]),
                        format(sums[feature_index, index], ".17g"),
                    ]
                )
    contexts = {}
    for key in keys:
        sample, context, organ, label = key
        if context not in contexts:
            contexts[context] = {
                "cell_type": context,
                "tissue_system": organ,
                "organ": organ,
                "source_cell_type": label,
                "n_cells": 0,
            }
        contexts[context]["n_cells"] += counts[key]
    samples = sorted({k[0] for k in keys})
    report = {
        "eligible_samples": samples,
        "n_cells": len(barcodes),
        "n_organs": len(ORGANS),
        "organs": sorted(ORGANS.values()),
        "groups": list(contexts.values()),
        "sample_groups": [
            {
                "sample": k[0],
                "cell_type": k[1],
                "tissue_system": k[2],
                "n_cells": counts[k],
            }
            for k in keys
        ],
        "features_per_sample": {sample: len(features) for sample in samples},
        "metadata_columns": list(rows[0]),
        "source_cell_count": len(rows),
        "published_cell_count": 116564,
        "cell_count_discrepancy": "Paper reports 116,564 cells; deposited RNA count matrix and metadata contain 115,395. All deposited cells are retained; the reason for the 1,169-cell difference is not established.",
        "sample_evidence": sample_evidence,
        "source_library_sizes_agree": expected_libraries is not None,
        "normalization": "log1p CP10K from uncorrected RNA counts; zeros included",
        "selection": "All deposited atlas cells, retaining the original organ, sample and cluster_names labels; normal Nipponbare source samples",
        "feature_missingness": "All original RNA count features; no marker or highly-variable-gene restriction",
        "streamed_files": {
            name: file_hash(raw / name)
            for name in [
                "i.bin",
                "p.bin",
                "x.bin",
                "matrix.json",
                "cell_metadata.tsv",
                "barcodes.tsv",
                "features.tsv",
            ]
        },
    }
    (raw / "extraction.json").write_text(json.dumps(report, indent=2) + "\n")
    sources = [
        {
            "filename": "GSE232863_scRNA_omics.Rds",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232863/suppl/GSE232863_scRNA_omics.Rds",
            "bytes": (raw / "GSE232863_scRNA_omics.Rds").stat().st_size,
            "sha256": file_hash(raw / "GSE232863_scRNA_omics.Rds"),
        },
        {
            "filename": "GSE232863_family.soft.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232863/soft/GSE232863_family.soft.gz",
            "sha256": file_hash(raw / "GSE232863_family.soft.gz"),
        },
    ]
    (ROOT / "data/expression_wang_atlas_inputs.json").write_text(
        json.dumps(
            {
                "accession": "GSE232863",
                "scope": "Full eight-organ normal Nipponbare single-cell multiomics RNA atlas",
                "sources": sources,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            {
                "cells": len(barcodes),
                "features": len(features),
                "organs": len(ORGANS),
                "contexts": len(contexts),
                "sample_groups": len(keys),
            },
            indent=2,
        ),
        flush=True,
    )


# Summarize the full atlas locally or inside a scheduled Fir compute job.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=ROOT / "data/raw/wang_atlas")
    args = parser.parse_args()
    prepare(args.raw)


if __name__ == "__main__":
    main()

"""Join Li's published Table S2 annotations to the original diploid count matrices."""

import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from build_expression import file_hash
from scripts.read_supplement import worksheet_rows

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data/raw/li_author"


# Match the authors' explicit replicate prefixes and preserve their original label codes.
def join_annotations(rows, known):
    labels = {
        "ep": "Epidermis (ep)",
        "ow": "Ovary wall (ow)",
        "nu": "Nucellus (nu)",
        "ov": "Ovule (ov)",
        "sti": "Stigma (sti)",
        "ch": "Chalaza (ch)",
        "C12": "Cluster 12 (C12)",
        "C13": "Cluster 13 (C13)",
        "sta": "Stamen residue (sta)",
    }
    samples = {"rep1": "GSM6346349", "rep2": "GSM6346350"}
    matched, seen, excluded = [], set(), Counter()
    for row in rows:
        prefix, barcode = row["A"].split("_", 1)
        if prefix not in samples or row["B"] != "rice" + prefix:
            raise ValueError("Unexpected author replicate or sample label.")
        sample = samples[prefix]
        key = sample, barcode
        if key in seen or barcode not in known[sample]:
            raise ValueError("Duplicate or unmatched author sample/barcode key.")
        seen.add(key)
        if row["D"] not in labels or not 0 <= int(row["C"]) <= 13:
            raise ValueError("Unrecognized author cluster or cell-type annotation.")
        if row["D"] == "sta":
            excluded["stamen_residue"] += 1
            continue
        matched.append(
            {
                "sample": sample,
                "count_barcode": barcode,
                "source_barcode": row["A"],
                "cell_type": labels[row["D"]],
                "source_cell_type": row["D"],
                "source_cluster": row["C"],
                "source_temporal_label": row["E"],
            }
        )
    return matched, dict(excluded)


# Verify frozen sources and write the exact author-cell join for sparse summarization.
def main():
    inputs = json.loads((ROOT / "data/expression_li_author_inputs.json").read_text())
    for source in inputs["sources"]:
        if file_hash(RAW / source["filename"]) != source["sha256"]:
            raise ValueError("Source checksum mismatch: " + source["filename"])
    config = json.loads((ROOT / "data/expression_li_selection.json").read_text())
    known = {}
    for sample, prefix in config["prefixes"].items():
        with gzip.open(RAW / (prefix + "_barcodes.tsv.gz"), "rt") as handle:
            known[sample] = {line.strip() for line in handle}
    rows = worksheet_rows(RAW / "mmc3.xlsx", "Rice_Pistil_index")
    if rows[0] != {
        "A": "cellcode",
        "B": "samplename",
        "C": "cluster",
        "D": "celltype",
        "E": "youngmodeltype",
    }:
        raise ValueError("Unexpected Table S2 worksheet header.")
    matched, exclusions = join_annotations(rows[1:], known)
    if len(rows) - 1 != 8173 or len(matched) + sum(exclusions.values()) != 8173:
        raise ValueError("Author atlas cell partition disagrees with the paper.")
    with (RAW / "matched_cells.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(matched[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(matched)
    config["scope"] = "Author Table S2 diploid pistil nuclei; stamen residue excluded"
    (RAW / "selection.json").write_text(json.dumps(config, indent=2) + "\n")
    report = {
        "source_author_nuclei": 8173,
        "included_nuclei": len(matched),
        "included_groups": len({r["cell_type"] for r in matched}),
        "included_by_sample": dict(Counter(r["sample"] for r in matched)),
        "exclusions": exclusions,
        "unmatched_barcodes": 0,
        "annotation_source": "Li et al. (2023), Table S2, Rice_Pistil_index",
        "source_sha256": file_hash(RAW / "mmc3.xlsx"),
    }
    (ROOT / "data/expression_li_author_join_validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

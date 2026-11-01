"""Match curated pistil labels to original diploid GEO matrices with explicit exclusions."""

import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data/raw/li"


# Join the curated annotations to exact GEO sample/barcode keys and record every exclusion.
def main():
    config = json.loads((ROOT / "data/expression_li_selection.json").read_text())
    (RAW / "selection.json").write_text(json.dumps(config, indent=2) + "\n")
    known = {}
    for sample, prefix in config["prefixes"].items():
        text = (RAW / (sample + ".soft")).read_text()
        experiment = re.search(r"term=(SRX\d+)", text).group(1)
        if config["experiments"].get(experiment) != sample:
            raise ValueError(
                "Curated experiment does not match the original GEO sample."
            )
        with gzip.open(RAW / (prefix + "_barcodes.tsv.gz"), "rt") as handle:
            known[sample] = {line.strip() for line in handle}
    matched = []
    seen = set()
    excluded = Counter()
    unmatched = []
    with (RAW / "cell_metadata.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if (
            row["Dataset"] != "SRP386976"
            or row["Condition"] != "Normal"
            or row["Tissue"] != "Pistil"
            or row["Genotype"] != "japonica cv. Nipponbare"
        ):
            raise ValueError("Unexpected dataset or baseline sample metadata.")
        experiment = row["Orig.ident"]
        if experiment not in config["experiments"]:
            if experiment != "SRX16335166":
                raise ValueError("Unrecognized source experiment.")
            excluded["haploid_1N_library"] += 1
            continue
        sample = config["experiments"][experiment]
        prefix = experiment + "@@_"
        if not row["barcode"].startswith(prefix):
            raise ValueError("Unexpected curated barcode prefix.")
        barcode = row["barcode"][len(prefix) :]
        key = (sample, barcode)
        if key in seen:
            raise ValueError("Duplicate curated sample/barcode pair.")
        seen.add(key)
        if row["Celltype"] in config["excluded_labels"]:
            excluded["label: " + row["Celltype"]] += 1
            continue
        if barcode not in known[sample]:
            excluded["curated_barcode_absent_from_GEO_matrix"] += 1
            unmatched.append(
                {
                    "sample": sample,
                    "source_barcode": row["barcode"],
                    "label": row["Celltype"],
                }
            )
            continue
        matched.append(
            {
                "sample": sample,
                "count_barcode": barcode,
                "source_barcode": row["barcode"],
                "cell_type": row["Celltype"],
            }
        )
    if not matched or len(matched) + sum(excluded.values()) != len(rows):
        raise ValueError("Curated cell partition is incomplete.")
    with (RAW / "matched_cells.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "count_barcode", "source_barcode", "cell_type"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(matched)
    report = {
        "source_curated_nuclei": len(rows),
        "included_nuclei": len(matched),
        "included_groups": len({r["cell_type"] for r in matched}),
        "included_by_sample": dict(Counter(r["sample"] for r in matched)),
        "exclusions": dict(excluded),
        "unmatched_barcodes": unmatched,
        "annotation_source": "scPlantDB SRP386976; original labels retained",
    }
    (ROOT / "data/expression_li_join_validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

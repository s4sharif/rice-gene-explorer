"""Create descriptive local dataset links and a calculation provenance inventory."""

import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "datasets"
STUDIES = [
    ("robertson_2026_well_watered_leaf", "expression", "expression", "robertson"),
    ("wang_2025_normal_root", "wang", "expression_wang", "wang"),
    (
        "wang_2025_eight_organ_atlas",
        "wang_atlas",
        "expression_wang_atlas",
        "wang_atlas",
    ),
    ("li_2023_diploid_pistil_scplantdb_annotations", "li", "expression_li", "li"),
    (
        "li_2023_diploid_pistil_author_annotations",
        "li_author",
        "expression_li_author",
        "li_author",
    ),
]


# Create a relative descriptive link while preserving the original source filename.
def link_file(folder, name, target):
    link = folder / name
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise FileExistsError(f"Refusing to replace a regular file: {link}")
    link.symlink_to(os.path.relpath(target, folder))


# Inventory each study's inputs, summaries, databases and validation calculations.
def organize():
    CATALOG.mkdir(exist_ok=True)
    rows = []
    for label, raw_name, stem, study in STUDIES:
        folder = CATALOG / label
        folder.mkdir(exist_ok=True)
        raw = ROOT / "data/raw" / raw_name
        pipeline = f"python3 -m pipelines.fetch_expression --study {study}; python3 -m pipelines.build_expression --study {study}; python3 -m pipelines.validate_expression --study {study}"
        if study == "wang_atlas":
            pipeline = "python3 -m pipelines.fetch_expression --study wang_atlas; python3 -m scripts.stream_wang_atlas data/raw/wang_atlas/GSE232863_scRNA_omics.Rds data/raw/wang_atlas; python3 -m pipelines.prepare_wang_atlas; python3 -m pipelines.build_expression --study wang_atlas; python3 -m pipelines.validate_wang_atlas"
        files = [
            (p, "source_or_intermediate") for p in sorted(raw.glob("*")) if p.is_file()
        ]
        files += [
            (p, "metadata_or_validation")
            for suffix in [
                "inputs",
                "manifest",
                "selection",
                "validation",
                "join_validation",
                "source_review",
            ]
            if (p := ROOT / "data" / f"{stem}_{suffix}.json").exists()
        ]
        files += [(ROOT / "data/databases" / (stem + ".sqlite"), "expression_database")]
        for source, role in files:
            if not source.exists():
                continue
            alias = (
                label
                + "_"
                + (
                    source.name
                    if source.parent == raw
                    else source.name.replace(stem, "expression", 1)
                )
            )
            link_file(folder, alias, source)
            calculation = "Source selection, extraction or provenance; consult the study input manifest."
            if source.name == "summaries.tsv.gz":
                calculation = "Per-sample/cell-type n_cells, n_expressing and sum_log; input to pipelines/build_expression.py."
            elif source.name.endswith(".sqlite"):
                calculation = "Cell-weighted mean log1p(10000*count/library_size), percent expressing=100*n_expressing/n_cells and within-study ranking; served by expression.py."
            elif "validation" in source.name:
                calculation = "Independent checks of counts, detection and normalized expression; not expression input."
            rows.append(
                [
                    label,
                    role,
                    str((folder / alias).relative_to(ROOT)),
                    str(source.relative_to(ROOT)),
                    calculation,
                    pipeline,
                ]
            )
    for label, filename, calculation in [
        (
            "jing_2026_gene_mapping_metadata",
            "expression_jing_metadata.json",
            "prepare_jing_metadata.py: conservative ZH11/MSU/RAP mapping and supporting annotations only; no quantitative expression.",
        ),
        (
            "liu_2026_seed_data_review",
            "expression_liu_review.json",
            "Source availability review only; no quantitative expression or cell-type calculations.",
        ),
        (
            "rapdb_gene_coordinates",
            "rice.sqlite",
            "build_database.py: RAP gene/alias import plus GWAS Atlas association filtering, coordinate checks and deduplication; app.py: neighboring-gene windows and plots.",
        ),
    ]:
        folder = CATALOG / label
        folder.mkdir(exist_ok=True)
        source = (
            ROOT
            / ("data/databases" if filename.endswith(".sqlite") else "data")
            / filename
        )
        alias = label + source.suffix
        link_file(folder, alias, source)
        rows.append(
            [
                label,
                (
                    "supporting_metadata"
                    if source.suffix == ".json"
                    else "gene_and_gwas_database"
                ),
                str((folder / alias).relative_to(ROOT)),
                str(source.relative_to(ROOT)),
                calculation,
                "See calculation column",
            ]
        )
    for label, sources, usage in [
        (
            "rapdb_gene_coordinates",
            list((ROOT / "data/raw").glob("rap_*")),
            "Gene coordinates, annotations and MSU/RAP mapping; build_database.py.",
        ),
        (
            "gwas_atlas_published_associations",
            [
                ROOT / "data/raw/rice_associations.txt.gz",
                ROOT / "data/provenance.json",
                ROOT / "data/validation.json",
                ROOT / "data/publication_review.json",
            ],
            "Published association import and selected-source audits; build_database.py, validate_data.py and audit_publications.py.",
        ),
        (
            "liu_2026_seed_data_review",
            list((ROOT / "data/raw/liu").glob("*")),
            "Downloaded source under review; not used in RiceLENS expression calculations.",
        ),
    ]:
        folder = CATALOG / label
        folder.mkdir(exist_ok=True)
        for source in sources:
            if not source.is_file():
                continue
            alias = label + "_" + source.name
            link_file(folder, alias, source)
            rows.append(
                [
                    label,
                    "source_or_review",
                    str((folder / alias).relative_to(ROOT)),
                    str(source.relative_to(ROOT)),
                    usage,
                    "See calculation column",
                ]
            )
    with (CATALOG / "calculation_index.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "dataset",
                "role",
                "descriptive_path",
                "original_path",
                "calculation_or_usage",
                "pipeline",
            ]
        )
        writer.writerows(rows)
    (CATALOG / "inventory.json").write_text(
        json.dumps(
            [
                dict(
                    zip(
                        [
                            "dataset",
                            "role",
                            "descriptive_path",
                            "original_path",
                            "calculation_or_usage",
                            "pipeline",
                        ],
                        row,
                    )
                )
                for row in rows
            ],
            indent=2,
        )
        + "\n"
    )
    print(f"Indexed {len(rows)} files without copying dataset contents.")


if __name__ == "__main__":
    organize()

"""Read precomputed baseline expression without loading single-cell matrices."""

import json
import sqlite3
from contextlib import closing
from pathlib import Path

from runtime_data import DATA_DIR

EXPRESSION_DB = DATA_DIR / "expression.sqlite"
WANG_DB = DATA_DIR / "expression_wang.sqlite"
LI_AUTHOR_DB = DATA_DIR / "expression_li_author.sqlite"
LI_DB = DATA_DIR / "expression_li.sqlite"


# Resolve a study identifier to its own database without accepting arbitrary paths.
def dataset_path(dataset_id):
    if dataset_id == "robertson_2026_control":
        return EXPRESSION_DB
    if dataset_id == "wang_2025_normal_root":
        return WANG_DB
    if dataset_id == "li_2023_pistil_2n_author":
        return LI_AUTHOR_DB
    if dataset_id == "li_2023_pistil_2n_scplantdb":
        return LI_DB
    raise ValueError("Select a recognized expression dataset.")


# Enable a study only while its versioned manifest remains present in the checkout.
def dataset_enabled(dataset_id):
    filenames = {
        "robertson_2026_control": "expression_manifest.json",
        "wang_2025_normal_root": "expression_wang_manifest.json",
        "li_2023_pistil_2n_scplantdb": "expression_li_manifest.json",
        "li_2023_pistil_2n_author": "expression_li_author_manifest.json",
    }
    return (Path(__file__).resolve().parent / "data" / filenames[dataset_id]).exists()


# List only imported studies so separate normalizations are never combined silently.
def expression_datasets():
    datasets = []
    for dataset_id in [
        "robertson_2026_control",
        "wang_2025_normal_root",
        "li_2023_pistil_2n_author",
        "li_2023_pistil_2n_scplantdb",
    ]:
        path = dataset_path(dataset_id)
        if path.exists() and dataset_enabled(dataset_id):
            with closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as db:
                metadata = json.loads(
                    db.execute(
                        "SELECT value FROM metadata WHERE key='dataset'"
                    ).fetchone()[0]
                )
            datasets.append(
                {
                    key: metadata[key]
                    for key in [
                        "id",
                        "citation",
                        "organ",
                        "condition",
                        "n_cells",
                        "n_samples",
                    ]
                }
            )
    return datasets


# Return measured features and within-study cell-type summaries for selected RAP genes.
def expression_evidence(gene_ids, dataset_id="robertson_2026_control"):
    genes = list(dict.fromkeys(gene_ids))
    if not genes or len(genes) > 12 or any(len(g) > 40 for g in genes):
        raise ValueError("Select between one and twelve gene identifiers.")
    path = dataset_path(dataset_id)
    if not path.exists() or not dataset_enabled(dataset_id):
        return {
            "available": False,
            "message": "This expression study is not enabled or its data have not been imported.",
        }
    with closing(sqlite3.connect(f"file:{path}?mode=ro", uri=True)) as db:
        db.row_factory = sqlite3.Row
        dataset = json.loads(
            db.execute("SELECT value FROM metadata WHERE key='dataset'").fetchone()[0]
        )
        results = []
        for gene in genes:
            feature = db.execute(
                "SELECT * FROM feature WHERE gene_id=?", (gene,)
            ).fetchone()
            rows = [
                dict(r)
                for r in db.execute(
                    "SELECT s.*,g.tissue_system FROM summary s JOIN cell_group g USING(cell_type) WHERE gene_id=? ORDER BY mean_log DESC,cell_type",
                    (gene,),
                )
            ]
            previous = None
            rank = 0
            for index, row in enumerate(rows, 1):
                if row["mean_log"] != previous:
                    rank = index
                row["rank"] = rank
                previous = row["mean_log"]
                row["percent_expressing"] = 100 * row["n_expressing"] / row["n_cells"]
            results.append(
                {
                    "gene_id": gene,
                    "status": "measured" if feature else "not_available",
                    "feature": dict(feature) if feature else None,
                    "groups": rows,
                }
            )
        return {
            "available": True,
            "dataset": dataset,
            "genes": results,
            "cell_groups": [
                dict(r)
                for r in db.execute("SELECT * FROM cell_group ORDER BY cell_type")
            ],
        }

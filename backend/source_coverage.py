"""Expose reviewed source coverage separately from measured expression results."""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


# Report supporting evidence without advertising pending datasets as imported expression.
def source_coverage(gene_id):
    if not re.fullmatch(r"Os\d{2}g\d{7}", gene_id):
        raise ValueError("Select a recognized RAP locus identifier.")
    sources = []
    jing = ROOT / "data/expression_jing_metadata.json"
    if jing.exists():
        metadata = json.loads(jing.read_text())
        sources.append(
            {
                "citation": metadata["citation"],
                "publication_url": metadata["publication_url"],
                "status": "Supporting metadata integrated; expression import pending",
                "message": f"The published supplement provides sample and cell-type definitions and {metadata['n_unique_rap_mappings']:,} unique RAP cross-references. It does not contain an all-gene expression matrix. Marker statistics are not used as baseline measurements.",
                "gene_mappings": [
                    r for r in metadata["gene_crossrefs"] if r["rap_id"] == gene_id
                ],
                "mapping_note": "A missing match means no unique cross-reference in this limited table, not absent expression. ZH11 coordinates are not transferred to Nipponbare.",
                "links": [
                    {"label": "Original atlas", "url": metadata["atlas_url"]},
                    {"label": "Author code", "url": metadata["code_url"]},
                ],
            }
        )
    liu = ROOT / "data/expression_liu_review.json"
    if liu.exists():
        sources.append(json.loads(liu.read_text()))
    return {"gene_id": gene_id, "sources": sources}

"""Import author sample definitions and conservative gene cross-references from Jing Table S1."""

import argparse
import json
import re
import sqlite3
from collections import defaultdict
from contextlib import closing
from pathlib import Path

from build_expression import file_hash
from scripts.read_supplement import worksheet_rows

ROOT = Path(__file__).resolve().parent


# Keep only unique author ZH11-to-MSU-to-RAP paths and reject many-to-one collisions.
def map_genes(pairs, aliases):
    records, candidates = [], defaultdict(list)
    for source_id, msu_ids in sorted(pairs.items()):
        targets = set()
        for msu in msu_ids:
            targets.update(aliases.get(msu.lower(), set()))
        record = {
            "source_gene_id": source_id,
            "msu_ids": sorted(msu_ids),
            "candidate_rap_ids": sorted(targets),
            "rap_id": None,
            "mapping_status": "unmapped or ambiguous",
        }
        records.append(record)
        if len(msu_ids) == 1 and len(targets) == 1:
            candidates[next(iter(targets))].append(record)
    for gene, matches in candidates.items():
        if len(matches) == 1:
            matches[0]["rap_id"] = gene
            matches[0]["mapping_status"] = "unique author MSU cross-reference"
    return records


# Extract supporting metadata without treating marker statistics as baseline expression.
def prepare(path):
    rows = worksheet_rows(path, "Table S1f")
    if rows[3].get("K") != "MSU.MSU" or rows[3].get("AA") != "MSU.MSU":
        raise ValueError("Unexpected author cross-reference table headers.")
    pairs = defaultdict(set)
    for row in rows[4:]:
        for source_col, target_col in [("H", "K"), ("Z", "AA")]:
            source = row.get(source_col, "").strip()
            targets = row.get(target_col, "").strip()
            if source in {"", "NA"} or targets in {"", "NA"}:
                continue
            if not re.fullmatch(r"ZH\d{2}G\d{5}", source):
                raise ValueError("Unexpected source gene identifier.")
            msu_ids = targets.split("|")
            if any(
                x != "ChrUn" and not re.fullmatch(r"LOC_Os\d{2}g\d{5}", x)
                for x in msu_ids
            ):
                raise ValueError("Unexpected author MSU cross-reference.")
            pairs[source].update(msu_ids)
    aliases = defaultdict(set)
    with closing(
        sqlite3.connect(f"file:{ROOT / 'data/rice.sqlite'}?mode=ro", uri=True)
    ) as db:
        for name, gene in db.execute(
            "SELECT name,gene_id FROM alias WHERE kind='MSU locus'"
        ):
            aliases[name.lower()].add(gene)
    mappings = map_genes(pairs, aliases)
    types = worksheet_rows(path, "Table S1c")
    cell_types = [
        {
            "project": r.get("A", ""),
            "cell_type": r.get("B", ""),
            "abbreviation": r.get("C", ""),
            "assay": r.get("D", ""),
        }
        for r in types[3:]
        if r.get("B")
    ]
    sample_rows = worksheet_rows(path, "Table S1b")
    output = {
        "id": "jing_metadata",
        "citation": "Jing et al. (online 2026; Cell 2027)",
        "publication_url": "https://doi.org/10.1016/j.cell.2026.09.023",
        "atlas_url": "https://db.genomics.cn/stomics/rice/",
        "code_url": "https://github.com/RiceSpatioTemporalAtlas/RSTA",
        "source_file": path.name,
        "source_sha256": file_hash(path),
        "status": "supporting_metadata_only",
        "expression_imported": False,
        "mapping_scope": "Only cross-references present in Table S1f marker tables; not a complete genome crosswalk. Marker effect sizes and detection percentages are not imported as baseline measurements.",
        "n_source_genes_with_crossrefs": len(mappings),
        "n_unique_rap_mappings": sum(r["rap_id"] is not None for r in mappings),
        "cell_types": cell_types,
        "sample_table_rows": sample_rows,
        "gene_crossrefs": mappings,
    }
    (ROOT / "data/expression_jing_metadata.json").write_text(
        json.dumps(output, indent=2) + "\n"
    )
    print(json.dumps({k: v for k, v in output.items() if k.startswith("n_")}, indent=2))


# Accept a local publisher supplement while keeping the original workbook unchanged.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook", type=Path)
    prepare(parser.parse_args().workbook)


if __name__ == "__main__":
    main()

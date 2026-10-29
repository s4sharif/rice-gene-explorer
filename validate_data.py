"""Audit imported coordinates and statistics against original export rows."""

import csv
import gzip
import json
import sqlite3

from build_database import ROOT, sha


# Check source hashes, imported coordinates, statistics and database integrity.
def validate():
    manifest = json.loads((ROOT / "data/provenance.json").read_text())
    for key in ["rap_loci", "rap_annotation", "rap_msu", "gwas_atlas"]:
        assert sha(ROOT / "data/raw" / manifest[key]["file"]) == manifest[key]["sha256"]
    db = sqlite3.connect(ROOT / "data/rice.sqlite")
    imported = {
        r[0]: r[1:]
        for r in db.execute(
            "SELECT source_row,chromosome,position,p_value,p_value_text FROM association"
        )
    }
    held = {r[0] for r in db.execute("SELECT source_row FROM quarantine")}
    checked = 0
    unicode_minus = 0
    with gzip.open(ROOT / "data/raw/rice_associations.txt.gz", "rt") as f:
        for line, row in enumerate(csv.DictReader(f, delimiter="\t"), 2):
            assert (line in imported) != (line in held)
            if line not in imported:
                continue
            chrom, pos, p, original = imported[line]
            assert chrom == int(row["Chr"].removeprefix("chr")) and pos == int(
                row["Pos"]
            )
            assert original == row["P-value"] and p == float(original.replace("−", "-"))
            unicode_minus += "−" in original
            checked += 1
    assert db.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
    assert not db.execute("PRAGMA foreign_key_check").fetchall()
    example = db.execute(
        "SELECT count(*) FROM association WHERE chromosome=3 AND position BETWEEN 16704501 AND 16760109"
    ).fetchone()[0]
    duplicate_ids = db.execute(
        "SELECT count(*) FROM (SELECT json_extract(raw_json,'$.\"#GaP_id\"') FROM association GROUP BY 1 HAVING count(*)>1)"
    ).fetchone()[0]
    db.close()
    result = {
        "source_checksums": "pass",
        "row_partition": "pass",
        "imported_coordinates_and_p_values_checked": checked,
        "unicode_minus_values_normalized": unicode_minus,
        "quarantined_rows": len(held),
        "integrity_check": "ok",
        "foreign_key_check": "ok",
        "reused_source_ids_retained_as_separate_rows": duplicate_ids,
        "GS3_25kb_association_records": example,
    }
    (ROOT / "data/validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    validate()

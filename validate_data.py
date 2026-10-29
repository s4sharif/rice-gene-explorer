"""Audit imported coordinates and statistics against original export rows."""

import csv
import gzip
import json
import sqlite3

from build_database import CHROMOSOME_LENGTHS, ROOT, sha


# Check source hashes, imported coordinates, statistics and database integrity.
def validate():
    manifest = json.loads((ROOT / "data/provenance.json").read_text())
    for key in ["rap_loci", "rap_annotation", "rap_msu", "gwas_atlas"]:
        assert sha(ROOT / "data/raw" / manifest[key]["file"]) == manifest[key]["sha256"]
    review = json.loads((ROOT / "data/publication_review.json").read_text())
    assert review["association_source_sha256"] == manifest["gwas_atlas"]["sha256"]
    assert (
        sha(ROOT / "data/chromosomes.json")
        == manifest["chromosome_reference"]["sha256"]
    )
    db = sqlite3.connect(ROOT / "data/rice.sqlite")
    for record in review["records"]:
        assert (
            db.execute(
                "SELECT count(*) FROM association WHERE id=?", (record,)
            ).fetchone()[0]
            == 1
        )
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
    invalid_genes = sum(
        db.execute(
            "SELECT count(*) FROM gene WHERE chromosome=? AND (start<1 OR end>?)",
            (chromosome, length),
        ).fetchone()[0]
        for chromosome, length in CHROMOSOME_LENGTHS.items()
    )
    invalid_variants = sum(
        db.execute(
            "SELECT count(*) FROM association WHERE chromosome=? AND (position<1 OR position>?)",
            (chromosome, length),
        ).fetchone()[0]
        for chromosome, length in CHROMOSOME_LENGTHS.items()
    )
    assert invalid_genes == invalid_variants == 0
    unique_evidence = db.execute("SELECT count(*) FROM evidence_summary").fetchone()[0]
    repeat_rows = len(imported) - unique_evidence
    repeat_groups = db.execute(
        "SELECT count(*) FROM evidence_summary WHERE source_record_count>1"
    ).fetchone()[0]
    quarantine_reasons = dict(
        db.execute("SELECT reason,count(*) FROM quarantine GROUP BY reason")
    )
    quarantine_publications = dict(
        db.execute(
            "SELECT json_extract(raw_json,'$.PMID'),count(*) FROM quarantine GROUP BY 1"
        )
    )
    allele_values = dict(
        db.execute(
            "SELECT COALESCE(effect_allele,'not supplied'),count(*) FROM association GROUP BY 1"
        )
    )
    assert not db.execute(
        "SELECT count(*) FROM association WHERE evidence_key NOT IN (SELECT evidence_key FROM evidence_summary)"
    ).fetchone()[0]
    duplicate_ids = db.execute(
        "SELECT count(*) FROM (SELECT json_extract(raw_json,'$.\"#GaP_id\"') FROM association GROUP BY 1 HAVING count(*)>1)"
    ).fetchone()[0]
    availability = {
        field: db.execute(
            f"SELECT count(*) FROM association WHERE {field} IS NOT NULL"
        ).fetchone()[0]
        for field in ["effect_allele_frequency", "effect_size", "other_allele"]
    }
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
        "chromosome_bounds": "pass",
        "unique_evidence_records": unique_evidence,
        "exact_repeated_source_rows": repeat_rows,
        "repeated_evidence_groups": repeat_groups,
        "quarantine_reasons": quarantine_reasons,
        "quarantine_by_publication": quarantine_publications,
        "effect_allele_values": allele_values,
        "effect_frequency_nonmissing": availability["effect_allele_frequency"],
        "effect_size_nonmissing": availability["effect_size"],
        "other_allele_nonmissing": availability["other_allele"],
    }
    (ROOT / "data/validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    validate()

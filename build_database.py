"""Import versioned public rice data; standard-library Python only."""

import csv
import gzip
import hashlib
import io
import json
import math
import re
import sqlite3
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ASSEMBLY = "Os-Nipponbare-Reference-IRGSP-1.0"
CHROMOSOMES = json.loads((ROOT / "data/chromosomes.json").read_text())
CHROMOSOME_LENGTHS = {int(k): v for k, v in CHROMOSOMES["lengths"].items()}
RAP = "https://rapdb.dna.naro.go.jp/download/archive/"
SOURCES = {
    "rap_loci": (
        "RAP-DB supported gene loci",
        RAP + "irgsp1/IRGSP-1.0_representative_2026-02-05.tar.gz",
        "2026-02-05",
        "rap_representative.tar.gz",
    ),
    "rap_annotation": (
        "RAP-DB gene annotations",
        RAP + "irgsp1/IRGSP-1.0_representative_annotation_2026-02-05.tsv.gz",
        "2026-02-05",
        "rap_annotation.tsv.gz",
    ),
    "rap_msu": (
        "RAP-DB RAP/MSU mapping",
        RAP + "RAP-MSU_2026-02-05.txt.gz",
        "2026-02-05",
        "rap_msu.txt.gz",
    ),
    "gwas_atlas": (
        "CNCB GWAS Atlas rice association export",
        "https://download.cncb.ac.cn/gwas_atlas/gwas_association_result_for_rice.txt.gz",
        "Unversioned export; identified by SHA-256",
        "rice_associations.txt.gz",
    ),
}
# Only explicit, equivalent assembly names. Numeric "1" is intentionally excluded.
ACCEPTED_REFERENCE = {"IRGSP-1.0", "IRGSP 1.0", "IRGSP1.0", "Nipponbare (IRGSP 1.0)"}


# Calculate a file checksum in blocks without loading the entire file into memory.
def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(1024 * 1024):
            h.update(b)
    return h.hexdigest()


# Convert empty source placeholders into explicit missing values.
def clean(value):
    value = (value or "").strip()
    return None if value in {"", "-", "NA", "None"} else value


# Import the documented rice snapshots into a validated replacement database.
def build():
    path = ROOT / "data/rice.sqlite"
    temporary = path.with_suffix(".building.sqlite")
    temporary.unlink(missing_ok=True)
    db = sqlite3.connect(temporary)
    db.executescript((ROOT / "schema.sql").read_text())
    db.execute(
        "INSERT INTO assembly VALUES (?,?,?,?)",
        (ASSEMBLY, "Oryza sativa", "Nipponbare", "1-based inclusive"),
    )
    manifest = {}
    for key, (name, url, release, file) in SOURCES.items():
        digest = sha(ROOT / "data/raw" / file)
        db.execute(
            "INSERT INTO source VALUES (?,?,?,?,?)", (key, name, url, release, digest)
        )
        manifest[key] = {
            "name": name,
            "url": url,
            "release": release,
            "file": file,
            "sha256": digest,
        }
    with tarfile.open(ROOT / "data/raw/rap_representative.tar.gz") as archive:
        with io.TextIOWrapper(
            archive.extractfile("IRGSP-1.0_representative/locus.gff")
        ) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                chrom, _, kind, start, end, _, strand, _, attrs = line.rstrip().split(
                    "\t"
                )
                if kind != "gene":
                    continue
                a = dict(x.split("=", 1) for x in attrs.split(";") if "=" in x)
                gene = a["ID"]
                chromosome = int(chrom.removeprefix("chr"))
                if not 1 <= int(start) <= int(end) <= CHROMOSOME_LENGTHS[chromosome]:
                    raise ValueError(
                        f"Gene coordinates outside chromosome bounds: {gene}"
                    )
                db.execute(
                    "INSERT INTO gene VALUES (?,?,?,?,?,?,?,?,?)",
                    (
                        gene,
                        ASSEMBLY,
                        "rap_loci",
                        int(chrom.removeprefix("chr")),
                        int(start),
                        int(end),
                        strand,
                        None,
                        a.get("Note", ""),
                    ),
                )
                db.execute(
                    "INSERT INTO alias VALUES (?,?,?,?)",
                    (gene, gene, "RAP locus", "rap_loci"),
                )
    genes = {x[0] for x in db.execute("SELECT id FROM gene")}
    with gzip.open(ROOT / "data/raw/rap_annotation.tsv.gz", "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gene = row["Locus_ID"]
            if gene not in genes:
                continue
            db.execute(
                "INSERT OR IGNORE INTO alias VALUES (?,?,?,?)",
                (gene, row["Transcript_ID"], "RAP transcript", "rap_annotation"),
            )
            symbol = None
            for field in [
                "CGSNL Gene Symbol",
                "RAP-DB Gene Symbol Synonym(s)",
                "Oryzabase Gene Symbol Synonym(s)",
            ]:
                for item in re.split(r"[,;|]", row[field]):
                    item = item.strip()
                    if not clean(item):
                        continue
                    symbol = symbol or item
                    db.execute(
                        "INSERT OR IGNORE INTO alias VALUES (?,?,?,?)",
                        (gene, item, "gene symbol", "rap_annotation"),
                    )
            db.execute(
                "UPDATE gene SET symbol=COALESCE(symbol,?),description=? WHERE id=?",
                (symbol, row["Description"], gene),
            )
    with gzip.open(ROOT / "data/raw/rap_msu.txt.gz", "rt") as f:
        for line in f:
            if not line.strip():
                continue
            gene, mapping = line.rstrip("\n").split("\t")
            if gene not in genes:
                continue
            for alias in mapping.split(","):
                if not clean(alias):
                    continue
                for value, kind in [
                    (alias, "MSU transcript"),
                    (alias.split(".")[0], "MSU locus"),
                ]:
                    db.execute(
                        "INSERT OR IGNORE INTO alias VALUES (?,?,?,?)",
                        (gene, value, kind, "rap_msu"),
                    )
    total = 0
    with gzip.open(ROOT / "data/raw/rice_associations.txt.gz", "rt") as f:
        for line, row in enumerate(csv.DictReader(f, delimiter="\t"), 2):
            total += 1
            reason = None
            try:
                if row["Ref ver"] not in ACCEPTED_REFERENCE:
                    raise ValueError("Unresolved reference label: " + row["Ref ver"])
                chrom = int(row["Chr"].removeprefix("chr"))
                pos = int(row["Pos"])
                p = float(row["P-value"].replace("−", "-"))
                if (
                    chrom not in CHROMOSOME_LENGTHS
                    or not 1 <= pos <= CHROMOSOME_LENGTHS[chrom]
                ):
                    raise ValueError(
                        "Genomic coordinate outside reference chromosome bounds"
                    )
                if not math.isfinite(p) or not 0 < p <= 1:
                    raise ValueError("P-value outside (0,1]")
            except (ValueError, TypeError) as exc:
                reason = str(exc)
            raw = json.dumps(row, ensure_ascii=False, sort_keys=True)
            evidence_fields = {
                key: value for key, value in row.items() if key != "#GaP_id"
            }
            evidence_key = hashlib.sha256(
                json.dumps(evidence_fields, sort_keys=True).encode()
            ).hexdigest()
            if reason:
                db.execute(
                    "INSERT INTO quarantine VALUES (?,?,?,?)",
                    ("gwas_atlas", line, reason, raw),
                )
                continue
            pmid = clean(row["PMID"])
            if pmid:
                db.execute(
                    "INSERT OR IGNORE INTO publication VALUES (?,?,?,?,?,?)",
                    (
                        pmid,
                        row["Title"],
                        row["Journal"],
                        row["First author"],
                        row["Published time"],
                        f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    ),
                )
            study_values = (
                row["StudyId"],
                pmid,
                row["GWAS model"],
                row["Species/Population"],
                row["Sample size"],
                row["Environment"],
                row["Tissue"],
            )
            # Retain distinct study contexts even if the source reuses a StudyId.
            study_id = (
                row["StudyId"]
                + "_"
                + hashlib.sha256(json.dumps(study_values).encode()).hexdigest()[:16]
            )
            db.execute(
                "INSERT OR IGNORE INTO study VALUES (?,?,?,?,?,?,?,?)",
                (study_id, *study_values),
            )
            trait_id = hashlib.sha256(
                (row["Trait"] + "|" + row["Trait accession"]).encode()
            ).hexdigest()[:20]
            db.execute(
                "INSERT OR IGNORE INTO trait VALUES (?,?,?)",
                (trait_id, row["Trait"], clean(row["Trait accession"])),
            )
            values = (
                f"{row['#GaP_id']}:{line}",
                "gwas_atlas",
                line,
                study_id,
                trait_id,
                ASSEMBLY,
                chrom,
                pos,
                p,
                row["P-value"],
                clean(row["Effect allele"]),
                clean(row["Effect allele frequence"]),
                None,
                None,
                clean(row["R2(%)"]),
                row["Ref ver"],
                row["Reported location"],
                clean(row["Reported gene(S)"]),
                raw,
                evidence_key,
            )
            db.execute(
                "INSERT INTO association VALUES (" + ",".join("?" * len(values)) + ")",
                values,
            )
    db.execute(
        "INSERT INTO evidence_summary SELECT evidence_key,MIN(source_row),COUNT(*) "
        "FROM association GROUP BY evidence_key"
    )
    db.commit()
    assert db.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
    assert not db.execute("PRAGMA foreign_key_check").fetchall()
    counts = {
        t: db.execute(f"SELECT count(*) FROM {t}").fetchone()[0]
        for t in [
            "gene",
            "alias",
            "association",
            "publication",
            "study",
            "trait",
            "quarantine",
        ]
    }
    assert counts["association"] + counts["quarantine"] == total
    counts["unique_evidence"] = db.execute(
        "SELECT count(*) FROM evidence_summary"
    ).fetchone()[0]
    manifest["chromosome_reference"] = {
        **CHROMOSOMES,
        "sha256": sha(ROOT / "data/chromosomes.json"),
    }
    manifest["audit"] = {
        "counts": counts,
        "source_association_rows": total,
        "accepted_reference_labels": sorted(ACCEPTED_REFERENCE),
        "quarantine_reasons": dict(
            db.execute("SELECT reason,count(*) FROM quarantine GROUP BY reason")
        ),
        "p_value_parsing": "Unicode minus U+2212 normalized to ASCII minus for numeric parsing only; original text preserved",
        "integrity_check": "ok",
        "foreign_key_check": "ok",
    }
    db.close()
    temporary.replace(path)
    (ROOT / "data/provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["audit"], indent=2))


if __name__ == "__main__":
    build()

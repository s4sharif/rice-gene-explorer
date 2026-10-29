"""Compare selected original-paper tables with imported association records."""

import argparse
import hashlib
import json
import math
import re
import urllib.request
import xml.etree.ElementTree as ET

from app import connect
from build_database import ROOT

PAPERS = {"33367255": "PMC7748012", "29616055": "PMC5864933"}


# Parse a publication's scientific-notation p-value while retaining its reported precision.
def parse_publication_p(text):
    text = text.replace("−", "-").replace(" ", "")
    if "×10" in text:
        coefficient, exponent = text.split("×10")
        return float(coefficient) * 10 ** int(exponent)
    return float(text)


# Extract the first primary-results table into coordinate and statistic comparisons.
def extract_table(pmid, path):
    root = ET.parse(path)
    table = root.findall(".//table-wrap")[0]
    result = []
    for row in table.findall(".//tr")[1:]:
        values = [" ".join("".join(cell.itertext()).split()) for cell in row]
        if pmid == "33367255":
            chromosome, position = int(values[1]), int(values[2].replace(" ", ""))
            trait = values[0].lstrip("∗*#$").strip().lower()
            result.append(
                {
                    "chromosome": chromosome,
                    "position": position,
                    "trait": trait,
                    "p": parse_publication_p(values[7]),
                    "statistic_kind": "FDR-adjusted p-value",
                    "model": None,
                    "locator": "Table 1, " + values[0],
                    "major_allele": values[3],
                    "minor_allele": values[4],
                    "minor_allele_frequency": values[5],
                }
            )
        else:
            chromosome, position = map(int, re.findall(r"\d+", values[2]))
            for column, model in [
                (3, "Compressed Mixed Linear Model"),
                (4, "General Linear Model"),
                (5, "Compressed Mixed Linear Model"),
            ]:
                result.append(
                    {
                        "chromosome": chromosome,
                        "position": position,
                        "trait": "mesocotyl length",
                        "p": 10 ** -float(values[column]),
                        "statistic_kind": "p-value reconstructed from rounded -log10(p)",
                        "model": model,
                        "locator": f"Table 1, statistic column {column - 2}",
                    }
                )
    return result


# Download only bounded open-access XML documents needed for the primary-table audit.
def download_papers():
    folder = ROOT / "data/raw/validation"
    folder.mkdir(parents=True, exist_ok=True)
    for pmid, pmcid in PAPERS.items():
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
        with urllib.request.urlopen(url, timeout=30) as response:
            body = response.read(5_000_001)
        if len(body) > 5_000_000:
            raise ValueError("Unexpectedly large article XML response.")
        ET.fromstring(body)
        (folder / f"{pmid}.xml").write_bytes(body)


# Record matched statistics and source conflicts without silently correcting the imported evidence.
def audit():
    reviews = {}
    sources = {}
    with connect() as db:
        for pmid, pmcid in PAPERS.items():
            path = ROOT / f"data/raw/validation/{pmid}.xml"
            rows = extract_table(pmid, path)
            sources[pmid] = {
                "url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/",
                "xml_url": f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "table": "Table 1",
            }
            associations = db.execute(
                "SELECT a.id,a.chromosome,a.position,a.p_value,t.label,s.model FROM association a JOIN trait t ON t.id=a.trait_id JOIN study s ON s.id=a.study_id WHERE s.publication_id=?",
                (pmid,),
            )
            for association in associations:
                matches = [
                    row
                    for row in rows
                    if row["chromosome"] == association["chromosome"]
                    and row["position"] == association["position"]
                    and math.isclose(row["p"], association["p_value"], rel_tol=0.02)
                ]
                if not matches:
                    continue

                # Prefer a compatible table column when rounded statistics admit several matches.
                def match_priority(row):
                    return (
                        row["model"] is not None
                        and row["model"] != association["model"],
                        row["trait"] not in association["label"].lower(),
                    )

                matches.sort(key=match_priority)
                match = matches[0]
                conflicts = []
                original_trait = match["trait"]
                if original_trait not in association["label"].lower():
                    conflicts.append(
                        f"Trait conflict: original table reports {original_trait}; export reports {association['label']}."
                    )
                if match["model"] and association["model"] != match["model"]:
                    conflicts.append(
                        f"Model conflict: original table identifies {match['model']}; export reports {association['model']}."
                    )
                reviews[association["id"]] = {
                    "pmid": pmid,
                    "locator": match["locator"],
                    "statistic_kind": match["statistic_kind"],
                    "original_trait": original_trait,
                    "original_model": match["model"],
                    "major_allele": match.get("major_allele"),
                    "minor_allele": match.get("minor_allele"),
                    "minor_allele_frequency": match.get("minor_allele_frequency"),
                    "coordinate_and_statistic_match": True,
                    "comparison_tolerance": "2% relative tolerance for rounded table values",
                    "conflicts": conflicts,
                    "original_source_url": sources[pmid]["url"],
                }
    manifest = json.loads((ROOT / "data/provenance.json").read_text())
    report = {
        "association_source_sha256": manifest["gwas_atlas"]["sha256"],
        "scope": "Selected records matching Table 1 in two original publications; unmatched records are not validated by this audit.",
        "sources": sources,
        "matched_records": len(reviews),
        "records_with_conflicts": sum(
            bool(row["conflicts"]) for row in reviews.values()
        ),
        "records": reviews,
    }
    (ROOT / "data/publication_review.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(
        json.dumps(
            {key: value for key, value in report.items() if key != "records"}, indent=2
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    if args.download:
        download_papers()
    audit()

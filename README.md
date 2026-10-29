# Rice Gene Explorer

A working local prototype for searching rice genes, exploring neighboring genes and viewing published GWAS evidence. Single-cell expression is reserved for a later stage. Searches query preloaded data; they do not run GWAS or CATCH.

## Clone and run

```sh
git clone https://github.com/s4sharif/rice-gene-explorer.git
cd rice-gene-explorer
python3 fetch_data.py
python3 build_database.py
python3 validate_data.py
python3 -m unittest discover -s tests -v
python3 app.py --port 8766
```

Open [Rice Gene Explorer](http://127.0.0.1:8766) in a browser after starting the server. This address runs on your own computer. Keep the terminal open; Ctrl+C stops the server. After the first database build, subsequent launches need only `python3 app.py --port 8766`. If that port is occupied, use `python3 app.py --port 8767` and open [port 8767](http://127.0.0.1:8767) instead. On macOS, `Launch.command` starts the default port, [8765](http://127.0.0.1:8765).

Use Python 3.11 or newer. The application and data importer use only the Python standard library; no Python packages, Node.js, Docker or HPC service are required to run the prototype. Internet access is needed for the initial source downloads, approximately 21 MB compressed. The database occupies approximately 224 MB after building; allow additional space while rebuilding. Development tools are described in [CONTRIBUTING.md](CONTRIBUTING.md).

The repository contains source code, documentation, a preview and provenance manifests. Downloaded third-party source files, the derived SQLite database and local caches are excluded from Git. The data scripts reproduce the database from the recorded source URLs and checksums. A private repository must be accessible to your GitHub account before cloning.

![Rice Gene Explorer interface](docs/design_preview.png)

## What is implemented

- Search RAP gene IDs, RAP transcripts, mapped MSU identifiers, and annotated gene symbols; case-insensitive exact matching followed by literal prefix matching.
- Show all matching genes when an identifier is ambiguous. For example, `GS3` matches two loci in this annotation; the prototype does not choose one automatically.
- Retrieve reference coordinates, strand, description, aliases and annotation release.
- Choose 10, 25, 50, 100, 250 kb or 1 Mb on each side of the gene span.
- Plot gene positions and orientations with GWAS association strength; select a neighboring gene to investigate its region.
- Filter associations by trait text, study, publication, reported p-value and strand-aware position. Show source p-values, available effect-allele information and R², study details, publication links and original reported locations.
- Distinguish overlap, upstream and downstream using the selected gene's strand, with distance to the nearest gene boundary.
- Paginate evidence in groups of 8 and export all matching associations as CSV. Plots show the 500 strongest records when a region has more than 500.

## Data loaded

| Data                                                      |   Scope |
| --------------------------------------------------------- | ------: |
| RAP-DB supported gene loci, annotation release 2026-02-05 |  36,061 |
| Identifier/alias mappings                                 | 195,087 |
| Searchable GWAS association records                       | 132,765 |
| Publications represented by searchable records            |     152 |
| Study contexts                                            |     828 |
| Trait label/ontology combinations                         |     357 |
| Quarantined association records                           |   2,298 |

Reference assembly: **Os-Nipponbare-Reference-IRGSP-1.0**, nuclear chromosomes 1–12. All local coordinates are **1-based, inclusive**. Predicted-only genes and organelle genes are outside this first dataset.

Sources: [RAP-DB downloads](https://rapdb.dna.naro.go.jp/download/irgsp1.html), [RAP/MSU mapping explanation](https://rapdb.dna.naro.go.jp/faq.html), [GWAS Atlas rice download](https://ngdc.cncb.ac.cn/gwas/download), and [GWAS Atlas curation documentation](https://ngdc.cncb.ac.cn/gwas/documentation). Exact download URLs and SHA-256 hashes are in `data/provenance.json`. Scientific release dates are source versions, not authorship metadata.

The GWAS Atlas download is a curated association export, not complete genome-wide summary statistics. Its documentation describes harmonization of historical coordinates. This prototype uses exported `Chr`/`Pos` only when `Ref ver` explicitly identifies IRGSP-1.0; it preserves the original reported location, which may use an older assembly. Upstream coordinate conversions have not been independently revalidated against original studies. No local liftover was performed.

Of 135,063 source records, 2,291 have the unresolved reference label `1`; five additional records have p-values outside (0,1]. Two further positions exceed the reference chromosome lengths. These 2,298 records remain in `quarantine`, with original rows and reasons. Unicode minus signs in 139 valid p-values are normalized for numeric parsing only; source text is retained. All source rows are retained. The interface groups exact repeats by default and can show every source row; neither count implies independent replication.

## Try a real example

Search **Os03g0407400**, a locus with GS3 among its aliases (the annotation's preferred display symbol is LK3). At ±25 kb:

- Gene: chr3:16,729,501–16,735,109, negative strand.
- Query region: chr3:16,704,501–16,760,109.
- Three overlapping gene spans, 240 source rows and 237 displayed evidence records after grouping exact repeats.
- One record reports grain length at chr3:16,733,441, p = 7.82 × 10⁻⁵⁴, study ST415, [publication PMID 28385155](https://pubmed.ncbi.nlm.nih.gov/28385155/).

This example demonstrates retrieval of published evidence. It does not establish causality or independently reproduce that study's analysis. Record counts include correlated observations and potentially repeated source entries.

Try **LOC_Os01g01010** to check MSU-to-RAP mapping. Missing cross-references are not guessed: the RAP/MSU source gives no mapping for Os03g0407400 in this release.

## Rebuild and validate

```sh
python3 fetch_data.py
python3 build_database.py
python3 validate_data.py
python3 -m unittest discover -s tests -v
```

The fetcher verifies the included snapshot checksums and refuses silently changed downloads. The builder writes a replacement database before swapping it into place. Keep raw sources and `provenance.json` to reproduce this snapshot. `validate_data.py` compares every accepted coordinate and p-value against its source row, verifies that every source row is imported or quarantined, and checks database integrity and foreign keys. Eleven automated integration test groups cover identifier ambiguity/mappings, boundary/strand arithmetic, region/CSV consistency, eight-record pagination, window/trait filters, quarantine handling and invalid query input.

## Project documentation and next steps

- [Project brief](docs/project_brief.md): current objectives and first-stage requirements.
- [Implementation plan](docs/implementation_plan.md): response to the revised project objectives, architecture, initial sources, milestones and open decisions.
- [Database design](docs/database_design.md): table relationships, coordinate rules, missing values and planned expression extension.
- [Learning notes](docs/learning_notes.md): a short explanation of the Python/database concepts used here.
- `schema.sql`, `build_database.py`, `app.py`, `static/`: database definition, import, API and interface.
- `data/provenance.json`, `data/validation.json`: reproducibility and import validation.

Before treating this as a research database: expand the selected-paper review, resolve excluded reference labels with supporting evidence, assess cross-study independence and allele semantics, and agree on source coverage and update policies. Exact-row duplication and chromosome bounds are now audited. Public hosting and redistribution require review of source-specific reuse terms. This local pilot is not a complete rice GWAS catalog or a production service.

## Review package

Start with the [review walkthrough](docs/review_guide.md) and [requirements checklist](docs/requirements_checklist.md). The first-stage functional checklist has 19 of 20 checkpoints implemented (95%); supervisor acceptance remains open. This does not describe the proportion of scientifically validated records.

The source audit groups 2,803 exact repeated rows into 129,962 distinct evidence records. Direct comparisons with Table 1 in two original papers match coordinates/statistics for 93 source records, including 29 with trait or model-label conflicts. Checked adjusted p-values are identified and conflicting metadata is flagged in the table and exports. Unchecked records are explicitly labelled. See `data/publication_review.json`, `data/validation.json` and `data/review_examples.json`.

To repeat the optional primary-table review, run `python3 audit_publications.py --download`. It downloads bounded open-access XML, not full supplementary archives. This is selected-record curation, not validation of the entire literature collection.

## Code quality

Every Python and JavaScript function has a one-sentence comment immediately above it. The repository includes pinned development tools, formatting configuration and automated checks for that convention. GitHub Actions checks source formatting, function comments, downloads, import validation and application tests. See [CONTRIBUTING.md](CONTRIBUTING.md) for the commands.

No software license has been selected for the project. Upstream datasets retain their own reuse terms; downloading data for this local academic prototype does not grant permission to redistribute it. The original email correspondence and prior CATCH work are outside this repository.

# Implementation approach for the rice evidence platform

## Objective and first deliverable

Build a gene-centered evidence platform beginning with rice. A search connects an identified gene to its genomic neighborhood and existing published GWAS associations. Single-cell expression is the next stage. CATCH is optional future analysis and does not block database development.

The first prototype now runs locally with real public annotations and published association records. It supports identifier resolution, coordinate retrieval, window selection, neighboring genes, association strengths, publication links and evidence exports. No new GWAS or RNA-seq processing is triggered by search.

## Initial datasets

| Source                                                                  | Use                                                                            | Choice and limitation                                                                                        |
| ----------------------------------------------------------------------- | ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------ |
| RAP-DB IRGSP-1.0 supported gene loci and annotation, release 2026-02-05 | Coordinates, descriptions, transcripts and symbols                             | Versioned, compact files; excludes predicted-only and organelle genes in the initial import                  |
| RAP-DB RAP/MSU correspondence, same release                             | Recognized MSU locus and transcript identifiers                                | Preserve multiple matches and absent mappings; do not equate differently named annotations automatically     |
| CNCB GWAS Atlas rice association export                                 | Published trait associations, statistics, study metadata and publication links | Downloadable curated significant associations; not full summary statistics or exhaustive literature coverage |
| RiceVarMap                                                              | Candidate second evidence source                                               | Assess overlap, reference/cultivar, programmatic access and reuse terms before importing; not loaded         |
| Published processed rice single-cell atlas                              | Next-stage expression source                                                   | Select tissue, reference, matrix/metadata availability and reuse terms with the supervisor; not loaded       |

Primary source pages: [RAP-DB](https://rapdb.dna.naro.go.jp/download/irgsp1.html), [identifier semantics](https://rapdb.dna.naro.go.jp/faq.html), [GWAS Atlas download](https://ngdc.cncb.ac.cn/gwas/download), [curation and coordinate harmonization](https://ngdc.cncb.ac.cn/gwas/documentation), [RiceVarMap](https://ricevarmap.ncpgr.cn/v3/).

Genome compatibility is the first quality gate. The initial coordinate space is Os-Nipponbare-Reference-IRGSP-1.0. Explicit equivalent labels are accepted; ambiguous `1` labels are held out. Historical source positions are preserved separately from Atlas's harmonized positions. The current import relies on upstream harmonization; it does not claim independent confirmation of those mappings.

## Database and application

Use SQLite for the local prototype: a file-based relational database supports indexed region searches without a separate database server. Python handles controlled imports and a read-only API; an HTML/JavaScript interface displays results and a genomic neighborhood chart. This avoids the memory-heavy analytical steps of CATCH. The current dataset can be explored on the existing Mac; no HPC service is required for this prototype.

Separate assembly, source release, gene, aliases, publication, study context, trait and association records. A variant position is not automatically assigned to one causal gene. Region overlap is computed when searching; source-reported gene annotations remain source assertions. Every association keeps its original source row and original fields for audit.

For a shared service, migrate to PostgreSQL with indexed genomic intervals and an application framework, add deployment configuration, concurrency/load tests, import permissions, accessibility checks and backup/update procedures. Preserve the same evidence model; do not add analytical jobs to interactive gene searches.

## Research interpretation

- Window means distance from the entire gene span, not from its midpoint or transcription start site.
- Overlap includes introns and does not establish a functional effect.
- Upstream/downstream depend on strand; genomic left/right alone is insufficient.
- Association p-values are retained as reported. Different studies use different populations, models, thresholds and phenotypes; do not interpret a pooled list as a meta-analysis.
- R² is variance explained as reported, not a beta coefficient. The current source export does not provide beta or the other allele; those fields remain missing.
- Multiple records can share a SNP, publication, phenotype or sample. Counts are evidence records, not independent studies or independent causal signals.
- An empty result means no matching record in the imported collection, not evidence against a gene's involvement.

## Milestones

1. **Local prototype — implemented.** Rice gene/alias search, reference coordinates, published regional evidence, neighboring-gene navigation, selectable windows, trait filter, publication provenance and CSV export.
2. **Curation review — initial audit implemented.** Chromosome boundaries and exact repeated rows are checked across the import. Selected original-table comparisons cover 93 source records in two papers and flag 29 trait/model conflicts. Ambiguous builds, zero p-values and out-of-bounds coordinates remain quarantined. Extend original-study review, assess study independence and agree trait/tissue priorities with Lisa.
3. **Evidence expansion.** Add a second rice GWAS source after documenting source keys, licensing, duplicates and genome conversions. Keep both source assertions and cross-source links. Add explicit study-defined regions/credible sets only where the original source provides them; a display window is not an LD block.
4. **Single-cell expression.** Import a processed rice dataset with dataset/tissue/cultivar/condition labels, cell annotations and gene identifier mapping. Store precomputed cell-type summaries, number of cells, fraction expressing and normalized-expression summaries with their method. Distinguish unmeasured, filtered and genuinely zero observations. Add dataset-specific expression panels without implying that expression confirms causality.
5. **Shared platform.** Choose hosting with the lab, review redistribution terms, migrate database as needed, add monitoring, access controls for imports, backups, update/version policy and external user testing. Multi-species support requires assembly-scoped identifiers and explicit orthology relationships.

## What to review with Lisa

The prototype makes the proposed approach concrete. The next useful review is whether the gene detail, neighborhood and evidence table answer the intended questions, which rice trait/tissue should guide validation, and which processed expression dataset should be prioritized. The prototype does not require a paid service or a public deployment.

## Prototype validation

The source has 135,063 records: 132,765 are searchable and 2,298 quarantined. All accepted source coordinates and p-values were checked against the raw export; source hashes, row accounting, foreign keys and SQLite integrity passed. Eleven automated test groups cover gene resolution and core query behavior. Browser checks cover ambiguous symbols, window/trait changes, pagination and selecting a neighboring gene. Selected primary-table review now covers 93 coordinate/statistic matches in two papers; 29 trait/model metadata conflicts are flagged. Broader publication curation remains an ongoing milestone.

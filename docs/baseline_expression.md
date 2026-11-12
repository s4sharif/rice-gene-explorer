# Baseline leaf expression in RiceLENS

## Imported coverage

The pilot contains **36,455 annotated well-watered nuclei**, **14 original cell types**, **33,438 distinct processed-matrix features**, and **32,017 exact matches** to the supported RAP gene catalog. There are **1,421 retained features without a supported-catalog match**. The annotated object contains 126,139 nuclei in total; only the 36,455 control nuclei contribute expression. Of 41,083 cells in the three processed control matrices, 4,628 lack retained author annotations and are excluded with explicit counts.

## Scope

This module uses Robertson et al. (2026), [GSE313726](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313726), and the authors' original cell annotations. The only eligible samples are:

| GEO sample | Original library  | Condition                        |
| ---------- | ----------------- | -------------------------------- |
| GSM9373072 | Well.Watered.Rep1 | Well-watered; 3 g water / g soil |
| GSM9373073 | Well.Watered.Rep2 | Well-watered; 3 g water / g soil |
| GSM9373074 | Well.Watered.Rep3 | Well-watered; 3 g water / g soil |

These are developing third leaves (~6 cm), Nipponbare, assayed using single-nucleus RNA-seq against IRGSP1.0. Eight drought/salinity samples are excluded from expression calculations. The authors' cell labels were defined in their combined atlas; retaining those labels does not mean importing stress-cell expression or performing a stress comparison.

## Inputs and processing

- The published annotated RDS object supplies cell metadata only. A narrow streaming XDR reader skips assay payloads instead of loading the entire Seurat object into RAM. Unsupported serialization types fail explicitly.
- The three control MEX matrices supply full available processed RNA counts, feature IDs and barcodes. Counts are not taken from integrated/scaled assays, and variable-gene or marker lists do not determine feature coverage.
- Every included cell must match its original annotation exactly within its sample. Duplicate keys and unmatched annotation cells stop the import. Matrix cells without retained author annotations are excluded and counted in the join report.
- Each replicate keeps its original feature list. A feature absent from one replicate is unavailable there, not a zero measurement. The database retains unmapped features and maps only unambiguous recognized RAP loci for the search interface.
- Sparse matrices are processed one replicate at a time, avoiding a dense all-gene-by-all-cell matrix. Only precomputed summaries are served to the website.

## Metrics

For nucleus i and feature g, expression is `log1p(10000 × count[g,i] / library_size[i])`. Library size is the sum across all features available in that replicate's processed MEX matrix. The displayed mean includes zero-count annotated nuclei. Detection means `count > 0`; the percentage is detected nuclei divided by eligible nuclei in that cell type for replicates containing that feature.

Per-sample summaries preserve counts and summed log expression. The combined display pools nuclei across eligible replicates, so larger replicates receive more weight. It is a descriptive summary, not an inferential test. Cell types are ranked separately for each gene within this study; equal means share a rank. Small cell groups and groups represented in few replicates need caution.

No differential expression, stress contrast, GWAS enrichment or causal inference is computed. These values must not be compared numerically with expression from a different study or normalization pipeline.

## Display

A searched gene loads its expression automatically. Users can select up to twelve neighboring genes from the current GWAS window. Dot area represents detection percentage, and a shared sequential blue scale represents mean expression for the currently selected genes within this dataset. Numeric values, cell counts, ranks, tissue-system labels and replicate coverage appear in an accessible table. Missing feature/group combinations are labelled NA, distinct from measured zero detection. CSV export includes the publication and metric definition.

## Reproduction

With the main rice annotation database already built:

```sh
python3 -m pipelines.fetch_expression
python3 scripts/extract_rds_metadata.py data/raw/expression/GSE313726_OsLeafStressIntegrated.rds data/raw/expression/cell_metadata.tsv
Rscript scripts/summarize_expression.R
python3 -m pipelines.build_expression
python3 -m pipelines.validate_expression
```

The extraction configuration is versioned in `data/expression_selection.json`; the fetch script copies it to `data/raw/expression/selection.json`. R needs Matrix and jsonlite; Python steps use the standard library. Downloads total approximately 1 GB. Raw inputs and SQLite databases are ignored by Git. `data/expression_inputs.json` pins filenames, public URLs and SHA-256 checksums; a changed download is rejected.

## Outputs and validation

- `data/raw/expression/cell_metadata.tsv`: original per-cell annotations; local only.
- `data/raw/expression/features.tsv`: union of available control features.
- `data/raw/expression/summaries.tsv.gz`: per-gene, per-sample, per-cell-type sufficient summaries.
- `data/raw/expression/extraction.json`: eligibility, join counts, group counts and per-sample coverage.
- `data/databases/expression.sqlite`: indexed features, groups, sample summaries and pooled summaries.
- `data/expression_manifest.json`: source provenance, metric and observed import counts.
- `data/expression_validation.json`: independent count-based validation results.

The independent Python validator rereads the three original MEX matrices and recalculates four example genes across sample/cell-type groups. It checks denominator counts, detection counts and normalized expression against the R summaries. Unit tests separately check tied ranks, measured zero, missing features and comparison limits. The metadata reader is checked against an R-created serialized fixture.

Observed validation: all **168** independent gene-by-replicate-by-cell-type checks passed (four genes × three replicates × fourteen cell types). For example, Os07g0141400 has mean log expression 1.4807 and 63.05% detection in mesophyll nuclei; these are descriptive values from this baseline dataset, not a cross-study or differential-expression result.

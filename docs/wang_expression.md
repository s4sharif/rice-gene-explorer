# Wang et al. (2025): normal-root expression subset

## Imported scope

This second expression dataset uses [GSE273875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273875), a separately deposited root experiment associated with [Wang et al., Nature (2025)](https://doi.org/10.1038/s41586-025-09251-0). It is **not the complete eight-organ atlas** in GSE232863. The main eight-organ object is about 5.1 GB; this annotated root object is about 390 MB compressed, and approximately 2.2 GB in R memory.

| GEO sample | Original sample ID | Source condition       | Inclusion             |
| ---------- | ------------------ | ---------------------- | --------------------- |
| GSM8438488 | LTP_G_1            | osltpl120 mutant / LTP | Excluded              |
| GSM8438489 | LTP_G_2            | osltpl120 mutant / LTP | Excluded              |
| GSM8438490 | NIP_G_1            | Nipponbare / WT        | Included: 8,427 cells |
| GSM8438491 | NIP_G_2            | Nipponbare / WT        | Included: 9,753 cells |

The source object contains 36,160 cells. Only the **18,180 WT cells** contribute expression. Both source fields must agree: `orig.ident` is NIP_G_1 or NIP_G_2, and `condition` is WT. No mutant comparisons or differential expression are run. The source sample records identify roots and genotype; developmental age is not independently established here and is not borrowed from other experiments in the paper.

## Cell labels and available genes

The deposited object contains a `cell_type` column with nine original labels: cortex, endodermis, epidermis, epidermis (neat root hair), exodermis, meristem, root cap, sclerenchyma and vascular. The unusual spelling “epidermis (neat root hair)” is retained as supplied. The authors' [analysis archive](https://doi.org/10.5281/zenodo.15501781), script 15, documents the cluster-to-cell-type assignment. Those annotations were derived from the combined WT/mutant object, but only WT counts enter RiceLENS summaries.

The RNA count layer has **30,206 features**, whereas the scaled layer has only 2,000 variable features. The import uses all 30,206 RNA-count features, including measured zeros. It uses neither the scaled layer nor a marker list. The retained count layer reflects the authors' feature/cell QC, so this does not imply every rice gene was measured or retained.

The GEO methods name MSU7 as the alignment reference, but the deposited nuclear feature IDs are already RAP-style IDs. **28,636 features match the supported RiceLENS RAP catalog exactly**. Other feature labels are retained as unavailable to RAP search; the import does not guess their identity or transfer MSU coordinates into the GWAS browser. The versioned RAP-MSU alias resource is available for explicitly named MSU identifiers, but the observed mapped features in this subset use direct RAP matches.

## Expression and study separation

Metrics use the same explicitly stated log1p CP10K definition as the leaf pilot, recalculated from uncorrected integer RNA counts. Detection is count > 0. Means include zero-count WT cells. Sample-level summaries are retained; the displayed within-study summary weights replicates by their cell counts. Cell types are ranked separately for each gene. Normalization similarity does not make values comparable between studies, organs, or single-cell and single-nucleus assays.

The new database is `data/databases/expression_wang.sqlite`, separate from Robertson's `data/databases/expression.sqlite`. The website's expression-study selector loads one dataset at a time; it never pools the two studies. The CSV identifies the selected study and its own publication/metric. Neighboring genes can be compared within either dataset.

## Reproduction

```sh
python3 -m pipelines.fetch_expression --study wang
R_MAX_VSIZE=3G Rscript scripts/summarize_wang.R
python3 -m pipelines.build_expression --study wang
python3 -m pipelines.validate_expression --study wang
```

The source file is gzip-compressed twice. The R reader explicitly opens both compression layers. R needs Matrix and jsonlite; Python uses the standard library. Source objects, exported counts and local databases are excluded from Git. `data/expression_wang_inputs.json` pins the object, original sample metadata and analysis archive by checksum.

The R import verifies unique genes/barcodes, count-layer dimensions, integer nonnegative counts, exact metadata-to-count-cell matching and agreement of normal-sample/WT labels. It exports unmodified WT counts for an independent Python calculation of selected genes across all available sample/cell-type groups. This checks the calculations using the deposited count layer; it is not validation against independently generated sequencing data.

Outputs include extraction counts, features and per-sample summaries under `data/raw/wang`, the separate SQLite database, `data/expression_wang_manifest.json` and `data/expression_wang_validation.json`. The full eight-organ atlas and the other papers remain future imports.

Validation passed **54** independent gene–replicate–cell-type calculations (three available example genes × two WT replicates × nine cell types). The requested fourth example, Os07g0141400, is absent from this deposited count assay and is not treated as zero. All 19 software tests passed, including explicit study-separation checks.

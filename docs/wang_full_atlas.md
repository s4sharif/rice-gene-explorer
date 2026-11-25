# Wang et al. (2025): full deposited eight-organ RNA atlas

The full processed Seurat object from [GSE232863](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232863), associated with [Wang et al. (2025)](https://doi.org/10.1038/s41586-025-09251-0), is imported separately from the normal-root experiment in GSE273875. This adds baseline RNA expression; it does not reproduce the paper’s ATAC, regulatory-network, trajectory or perturbation analyses.

## Actual deposited coverage

| Source organ code | Organ label           | Deposited cells |
| ----------------- | --------------------- | --------------: |
| Root              | Crown root            |          31,108 |
| ST                | Culm                  |          19,740 |
| SP                | Young spike           |          17,457 |
| SAM               | Shoot apical meristem |          13,791 |
| leaf              | Third leaf            |          13,771 |
| Bud               | Tiller bud            |           8,858 |
| Seed              | Spikelet              |           6,276 |
| Flag              | Flag leaf             |           4,394 |
| **Total**         | **8 organs**          |     **115,395** |

The paper and GEO summary report **116,564 cells**, but the deposited RNA count assay and cell metadata both contain **115,395**. All deposited cells are retained. The reason for the 1,169-cell discrepancy has not been established; absent cells are not invented or counted. The `Seed` code is labelled Spikelet according to the matched GEO library titles; original embryo-related cell labels in this group are retained without reinterpretation.

The count matrix contains **37,960 features** and **144,433,235 stored entries**. All features are retained, including measured zeros. **35,797** map uniquely to the supported RAP catalog; the other **2,163** remain explicitly unmapped/ambiguous rather than being guessed or labelled zero. There are **32 original cell-type labels**, **56 distinct organ–cell-type contexts**, **16 original samples**, **112 sample–context groups**, and **4,251,520 feature–sample–context summaries**.

Original `tissue`, `sample` and `cluster_names` fields are preserved. The same cell-type name in different organs remains a separate context. Source labels such as `SP-unknown` and `Seed-unknown` are not removed or replaced. Each of the 16 sample names is matched to its original GEO RNA library; those records state Nipponbare genotype and normal treatment. No stress selection or differential expression is performed.

## Counts, normalization and validation

Only `RNA@counts` is read. RNA normalized/scaled values, ATAC, `atacRNA`, `chromvar`, marker lists and variable-feature subsets are not used. Library sizes are calculated over every retained RNA-count feature and agree exactly with source `nCount_RNA` for all deposited cells.

For each gene and sample–organ–cell-type group, the import stores the number of cells, number with RNA count > 0, and sum of `log1p(10000 × count / RNA library size)`. Study-level means include all zero-count cells and use cell-weighted pooling of the two source samples per organ context. Detection percentage is `100 × n_expressing / n_cells`. Ranks are within this atlas; values are not pooled with Robertson, Li or the separate Wang root experiment.

An independent validation reads the original sparse counts column by column, recalculates four genes with scalar `math.log1p`, and compares every sample context: **448 independent checks** passed for Os01g0100100, Os01g0883800, Os03g0407400 and Os12g0528400. Database integrity, foreign keys, all-feature coverage, all mapped zero-expression contexts, bounds, source checksums and sample-to-study aggregation also passed. This validates calculations against the deposited counts; it is not independent sequencing or biological validation.

The new dataset is enabled only when a matching source validation report is present. An independent R-generated fixture verifies that the streaming reader extracts the count slot rather than a differently scaled data slot, rejects an ATAC substitute, and preserves original metadata.

## Memory-bounded extraction

The source is **5,521,719,712 bytes**, SHA-256 `7a308241e07575770f227057b1cd716b4cb9c255200c4ff4189c567bc2fbc958`. Loading the whole multiomics object is unnecessary. The streaming XDR reader writes the RNA matrix’s big-endian `i`, `p` and `x` arrays to disk, skips other assays, and retains the original feature names, barcodes and metadata. It verifies the assay name and sparse-matrix class. Unsupported serialization is rejected. Pointer/environment records are skipped without evaluating code, following [R’s serialization format](https://github.com/wch/r-source/blob/trunk/src/main/serialize.c).

NumPy memory maps those arrays and summarizes blocks of 512 cells. The web application reads only compact precomputed SQLite results and does not load the source matrix. This import completed locally on the 8 GB Mac; no Fir job or Render deployment was used.

## Reproduction

Use Python 3.11 or newer in an environment with the pinned offline dependency:

```sh
python3 -m pip install -r requirements/expression.txt
python3 -m pipelines.fetch_expression --study wang_atlas
python3 -m scripts.stream_wang_atlas data/raw/wang_atlas/GSE232863_scRNA_omics.Rds data/raw/wang_atlas
python3 -m pipelines.prepare_wang_atlas
python3 -m pipelines.build_expression --study wang_atlas
python3 -m pipelines.validate_wang_atlas
python3 scripts/organize_datasets.py
```

The source download, disk-backed arrays and generated SQLite database are excluded from Git. Allow approximately 10 GB of additional space for this study and temporary database construction. On Fir, run extraction/import through a compute allocation rather than the shared login node.

## Files and calculations

- `data/expression_wang_atlas_inputs.json`: pinned RDS and GEO sample-metadata sources.
- `data/raw/wang_atlas/matrix.json`, `i.bin`, `p.bin`, `x.bin`: original RNA sparse counts, dimensions, vector lengths and byte order.
- `data/raw/wang_atlas/features.tsv`, `barcodes.tsv`, `cell_metadata.tsv`: original features, cells and all deposited annotations.
- `data/raw/wang_atlas/library_sizes.tsv`: RNA library sizes used as normalization denominators.
- `data/raw/wang_atlas/summaries.tsv.gz`: all feature–sample–context summaries.
- `data/raw/wang_atlas/extraction.json`: groups, organs, original sample evidence, discrepancy note and extracted-file hashes.
- `data/databases/expression_wang_atlas.sqlite`: feature mapping, cell groups, sample summaries and gene/context summaries (about 731 MiB locally).
- `data/expression_wang_atlas_manifest.json`: dataset identity, methods, source hash and original context annotations.
- `data/expression_wang_atlas_validation.json`: independent-calculation results and checks.
- `datasets/wang_2025_eight_organ_atlas/`: descriptive relative links to all inputs and outputs, with calculations in the dataset inventory.

## Local website

Cell Expression uses the existing shared gene/window context. The [study selector](expression_study_selector.md) also exposes the separately imported Wang root experiment, Robertson controls and both Li annotation versions without pooling their measurements. It provides an eight-organ selector, a dot plot, exact values, publication/source links, and comparison of up to six neighbouring genes within this study. Dot size encodes detection percentage; a sequential blue scale encodes mean expression. Unknown labels and missing features remain explicit. The CSV exports all 56 contexts for each selected gene even when the figure shows one organ.

This is a local review feature. Existing expression datasets, GWAS functionality and source coordinates remain separate. Neither main nor the live deployment is changed automatically.

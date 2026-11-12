# Jing lifecycle atlas: supplied metadata integration

The supplied article is _A spatiotemporal lifecycle atlas decodes spatial coordination in rice_, DOI [10.1016/j.cell.2026.09.023](https://doi.org/10.1016/j.cell.2026.09.023). Lisa referred to its 2026 online release; the supplied proof cites **Cell 190, January 7, 2027**. RiceLENS records both rather than silently changing the study identity.

## What the supplied workbook contains

`1-s2.0-S0092867426011256-mmc2.xlsx` is Table S1. It contains genome statistics (S1a), spatial/single-nucleus sample and QC statistics (S1b), cell-type names and abbreviations (S1c), marker references (S1d), enrichment results (S1e), top-marker statistics and some MSU cross-references (S1f), and model training/evaluation metadata (S1g–h). It is **not an all-gene expression matrix** and has no per-nucleus barcode-to-cell-type table.

RiceLENS imports the original sample-table values, cell-type definitions and gene cross-references as supporting metadata. It does not import fold changes, `pct.1`, `pct.2` or adjusted marker p-values as baseline expression. In particular, a marker fold change is not a cell-type mean and a gene absent from a top-marker table is not unexpressed.

## Gene cross-references

The importer reads the explicit author ZH11/MSU pairs in Table S1f (H/K and Z/AA). It excludes empty/`NA` placeholders, retains multi-target mappings as unresolved, joins MSU loci to the frozen local RAP crosswalk and rejects collisions where multiple ZH11 genes lead to the same RAP locus. It makes no symbol-only or coordinate-based transfers.

The supplement provides **1,955 ZH11 genes with cross-references**, of which **1,639 have unique accepted RAP mappings** under these rules. This is a marker-table-limited mapping resource, not a complete genome conversion. All source targets and candidate RAP matches remain in the metadata output, including ambiguous mappings. The website shows matching source identifiers in “Additional studies and data coverage” for a searched gene. No cross-reference means unavailable mapping evidence, not zero expression.

```sh
python3 -m pipelines.prepare_jing_metadata /path/to/1-s2.0-S0092867426011256-mmc2.xlsx
```

The original workbook is unchanged. `data/expression_jing_metadata.json` records its checksum, extracted metadata and mapping outcomes. `/api/expression/sources?gene=Os02g0662700` exposes study coverage separately from quantitative expression. Tests protect multi-target and many-to-one cases.

## Full-text data routes checked

The article's resource table identifies CNSA **CNP0008591**, the [RSTA atlas](https://db.genomics.cn/stomics/rice/) and [author code](https://github.com/RiceSpatioTemporalAtlas/RSTA). The code tree was reviewed at commit `d6b3c29497d5afbcab01e447fc4b14e656135233`; it contains analysis/model scripts and a Google Drive training-resource link. The Drive link redirects to Google sign-in in this environment. Local paths inside example notebooks are not public downloads. The public CNSA FTP project includes sequencing folders and ZH11 genome FASTA/GFF/GTF; those are not an annotated expression object.

[STOmicsDB STDS0000421](https://db.cngb.org/stomics/datasets/STDS0000421/data) lists `ZH11_SC.h5ad` at 79.85 GB. Its page and RSTA download metadata currently expose no direct download path for that file. The public Cirro visualization schema exposes 851,725 nuclei, 39,914 features, cell-type/sample/organ labels and per-feature values, but those values are rounded. Their relationship to the original raw RNA assay and exact detection denominators has not been verified. They are therefore not silently substituted for count-derived RiceLENS expression.

The supplied paper describes Seurat normalization, while the spatial analysis also uses SCTransform. This does not establish which assay a visualization export contains. Spatial bins must remain separate from nuclei. Organ/stage context must be retained rather than pooling identical cell labels across organs or stages.

## Still required for quantitative integration

A downloadable annotated single-nucleus RNA object or raw-count matrix with final per-nucleus metadata is required, preferably in organ-sized subsets. It must identify sample, cultivar/genotype, organ, stage, original barcode, cell type, feature IDs and the raw RNA count layer. A comprehensive author ZH11-to-MSU/RAP crosswalk would extend coverage beyond the partial supplement mapping. Alternatively, a documented all-gene group-summary export with exact detection counts, denominators and normalization could be imported as a distinct source metric after validation.

**Current status: supporting metadata integrated; Jing expression measurements not imported.** It is not selectable as a measured expression dataset, and pending data do not appear as zero-valued dots.

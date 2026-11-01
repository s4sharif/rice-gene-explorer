# Remaining expression papers: file-level findings

Hosting is deferred at Sharif's request. The focus is additional baseline expression evidence.

## Liu et al. (2026), early seed atlas

The [paper's data-availability statement](https://onlinelibrary.wiley.com/doi/full/10.1111/jipb.70373) now identifies processed snRNA-seq deposits:

- [OMIX010779](https://ngdc.cncb.ac.cn/omix/release/OMIX010779): wild-type SJ2, 3 and 5 days after pollination, two replicates per stage.
- [OMIX010934](https://ngdc.cncb.ac.cn/omix/release/OMIX010934): SJ2 at 7 days after pollination, two replicates.
- OMIX013730 is the osbzr4 mutant experiment and is outside the baseline import scope.
- CRA027252 is bulk RNA-seq and must not be substituted for cell-type expression.

The inspected smallest D7 archive contains barcodes, features and a count matrix, but no barcode-to-cell-type annotation table. The supplementary workbook descriptions identify markers and other analyses; marker lists do not supply the all-gene per-cell annotation join required by RiceLENS. A verified metadata table or annotated object remains necessary. No Liu expression values have been imported.

The deposit pages also cite an earlier OsBZR4 paper, but the JIPB statement explicitly names these accessions. Any reused samples must retain their original identifiers and cannot be counted as independent replication across publications. SJ2 must not be labelled Nipponbare.

## Jing et al. (online 2026; Cell 2027), lifecycle atlas

The supplied full text and Table S1 are now inspected. Sample/cell-type definitions and a conservative partial gene crosswalk are integrated as supporting metadata. There are 1,639 unique accepted RAP mappings; marker statistics are not used as baseline expression. See [Jing integration and remaining inputs](jing_metadata.md).

The original 79.85 GB H5AD is listed, but the inspected download interfaces expose no direct download path. The public visualization values are rounded and have not been verified against the raw RNA assay. No quantitative Jing expression is imported.

## Li et al. (2023)

The supplied Table S2 provides the original author annotations and has now been joined exactly to the original diploid GEO matrices. The author version retains 8,077 pistil nuclei after excluding 96 stamen residues. See [author import and validation](li_author_expression.md). The earlier [scPlantDB-curated version](li_expression.md) is preserved separately; these are two versions of the same libraries, not independent replication.

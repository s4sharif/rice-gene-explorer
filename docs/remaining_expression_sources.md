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

## Jing et al. (2026), lifecycle atlas

The [paper](https://doi.org/10.1016/j.cell.2026.09.023) identifies the authors' [Rice Spatio Temporal Atlas](https://db.genomics.cn/stomics/rice/) and CNSA CNP0008591. [STOmicsDB dataset STDS0000421](https://db.cngb.org/stomics/datasets/STDS0000421/data) lists an original `ZH11_SC.h5ad` object of approximately **79.85 GB**. This confirms a processed object exists; it has not been downloaded or inspected locally.

A suitable subset or documented all-gene cell-type summary export would be needed for this machine. The T2T reference and ZH11 cultivar also require a verified mapping to RiceLENS gene identifiers; no name-only coordinate transfer should be attempted. Spatial bins must not be treated as individual nuclei. No Jing expression values have been imported.

## Li et al. (2023)

A usable curated annotation source was found in scPlantDB and the diploid subset is now implemented. See [the Li import](li_expression.md). Original-author final cell annotations remain unverified, and the 1N library remains separate and unimported.

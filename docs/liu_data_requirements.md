# Liu early-seed atlas: full-text assessment

The supplied paper is _A single-nucleus transcriptomic atlas reveals distinct cell identities and key regulators of cellular differentiation during early rice seed development_, DOI [10.1111/jipb.70373](https://doi.org/10.1111/jipb.70373). Its methods describe alignment to IRGSP-1.0.57, Seurat processing, filtering, integration and author cell-type annotation. Its data-availability statement identifies:

| Source     | Baseline scope                                                      |
| ---------- | ------------------------------------------------------------------- |
| OMIX010779 | SJ2, 3 and 5 days after pollination, two biological replicates each |
| OMIX010934 | SJ2, 7 days after pollination, two biological replicates            |
| OMIX013730 | osbzr4 mutant; excluded                                             |
| CRA027252  | Bulk RNA-seq; excluded from single-nucleus expression               |

The full text does not provide a final barcode-to-cell-type table or a direct link to an annotated object. The previously inspected D7 count archive contains a standard sparse matrix, feature IDs and barcodes. A matrix barcode identifies a nucleus, but does not supply its biological cell identity. The described supporting datasets contain marker lists and downstream analyses; these cannot replace final per-nucleus annotations or all-gene expression.

The newly supplied `mmc3.xlsx` is **Li's 2023 pistil Table S2**, not Liu's seed data. The other workbook is **Jing's Table S1**. Their labels cannot be applied to Liu's nuclei. RiceLENS therefore displays Liu as a reviewed, pending source and supplies the exact deposit links, without fabricating measurements or assigning cell identities from marker lists.

## Exact input needed

A final metadata table must contain the original count-matrix barcode, library/replicate identifier and author cell-type label for each retained nucleus, with organ, developmental stage, genotype and condition supplied either per nucleus or in a sample table. All-gene raw RNA counts must be linked to those rows through exact sample/barcode keys. An annotated Seurat/AnnData object retaining that raw assay and metadata is also suitable.

After receipt, the import must verify the six wild-type library identities, separate stages, reject duplicate/missing joins, preserve original labels and feature identifiers, document gene mapping, calculate detection from raw RNA counts, and validate selected gene/group summaries independently. Shared accessions also appear in an earlier OsBZR4 publication and must not be counted as independent data merely because two papers cite them.

## Request text for the authors or advisor

> We are integrating baseline rice expression into RiceLENS using the SJ2 wild-type snRNA-seq libraries deposited under OMIX010779 and OMIX010934. Could you share the final retained-cell annotation table (original barcode, library/replicate, stage and cell-type label), or the annotated Seurat/AnnData object with its all-gene raw RNA count assay? We need the final barcode labels used in the paper, rather than marker lists, so that gene expression and exact expressing-cell percentages can be summarized reproducibly. We will retain your original labels and cite the paper and deposits.

This is a draft request only; no message has been sent.

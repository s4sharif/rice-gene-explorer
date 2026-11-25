# Selecting baseline expression studies

Cell Expression has a Study / dataset selector above the organ controls. It lists installed datasets from the existing read-only expression catalog rather than hardcoding a single Wang endpoint. Five dataset/annotation versions currently represent three publications:

| Choice                                       | Cells or nuclei | Contexts | Scope                                                            |
| -------------------------------------------- | --------------: | -------: | ---------------------------------------------------------------- |
| Wang (2025), eight-organ atlas               |         115,395 |       56 | All deposited RNA-count cells from eight organs                  |
| Wang (2025), separate normal-root experiment |          18,180 |        9 | Wild-type normal root replicates, mutants excluded               |
| Robertson (2026), well-watered leaf          |          36,455 |       14 | Control leaf nuclei only                                         |
| Li (2023), author-annotated pistil           |           8,077 |        8 | Original author labels, diploid nuclei, stamen residues excluded |
| Li (2023), alternate scPlantDB annotations   |           9,769 |       11 | A separate curation of overlapping Li libraries                  |

The Li choices reuse source libraries and must not be counted as independent replication. Each has its own labels, retained cells, methods and dataset identifier. The warning is visible whenever either version is selected.

A selection replaces the plotted data, numeric table, percentage denominators, publication and repository links, normalization notes, annotation source and CSV URL. Means and colour scales are never pooled across studies. Neighbouring genes can still be compared within the selected dataset. The gene and genomic window remain unchanged when switching studies.

The study and organ choices are stored in URL parameters `expression` and `expression_organ`, so they survive navigation and reloads. The organ selector remembers a different choice for each study during the session. Single-organ studies show their one available organ as a disabled control. Robertson's tissue-system families (Epidermis, Meristem, Mesophyll, Vasculature) are preserved as tissue systems; they are not relabelled as separate organs. Original cell-type labels are preserved in every dataset and in CSV exports.

## Papers pending complete inputs

Jing and Liu are explicitly listed under Not yet integrated, with disabled expression options and an expandable requirements panel. Their existing supporting metadata are retained. This interface change does not import missing expression or convert marker statistics into baseline measurements.

- Jing: the installation lacks an accessible all-gene expression object/matrix with its matching final annotations and sufficient gene-ID mapping. Its limited published marker/cross-reference metadata are not quantitative expression.
- Liu: counts have been located, but the final per-nucleus barcode/library-to-cell-type annotation join is still missing. Bulk RNA-seq, mutant libraries and marker lists are not substituted.

These are data/annotation requirements, not a need for more CPU alone. No new Fir computation, Render deployment or source-matrix reanalysis was performed for the selector. See [remaining sources](remaining_expression_sources.md) for the recorded file-level findings.

## Validation

The existing expression separation and numeric tests remain in place. A regression test verifies that organ labels are independent of tissue-system families and that unknown organs are explicitly unspecified. Local browser checks exercise all five choices, their context counts, publication text and dataset-specific export URLs. The source matrices, SQLite measurements, input manifests and original validation reports are unchanged by this interface update.

# First-stage acceptance checklist

This checklist measures the rice gene-search and genomic-neighborhood prototype. It does not measure completeness of the rice literature, correctness of every upstream record, single-cell integration or production readiness. Nineteen of twenty review checkpoints are implemented and checked: **95% on this functional checklist**. The remaining checkpoint requires Lisa's review; this percentage is not a scientific validation rate.

| #   | Checkpoint                                                                         | Status                | Evidence                                                              |
| --- | ---------------------------------------------------------------------------------- | --------------------- | --------------------------------------------------------------------- |
| 1   | Search a rice RAP gene ID                                                          | Complete              | Identifier integration test and review examples                       |
| 2   | Search recognized gene symbols                                                     | Complete              | GS3 and SD1 searches; ambiguity retained                              |
| 3   | Resolve alternate identifiers                                                      | Complete              | RAP transcript and official MSU cross-references                      |
| 4   | Show gene name, description and identifier                                         | Complete              | Gene detail panel                                                     |
| 5   | Show chromosome, start, end and strand                                             | Complete              | Annotation import and boundary tests                                  |
| 6   | Preserve assembly and annotation versions                                          | Complete              | Gene details, provenance and exports                                  |
| 7   | Retrieve existing regional GWAS evidence                                           | Complete              | Indexed queries over imported published evidence                      |
| 8   | Show associated traits                                                             | Complete              | Evidence table, trait filter and conflict notes                       |
| 9   | Show SNP positions                                                                 | Complete              | Table, chart and source-row audit                                     |
| 10  | Show gene-relative distances                                                       | Complete              | Inclusive-boundary and strand tests                                   |
| 11  | Distinguish overlap, upstream and downstream                                       | Complete              | Table labels and position filters                                     |
| 12  | Show reported p-values and available statistics                                    | Complete              | Missing values explicit; checked FDR statistics labelled              |
| 13  | Identify studies and publications                                                  | Complete              | Study/publication filters and links                                   |
| 14  | Link evidence to original sources                                                  | Complete              | Original-row endpoint, download provenance and paper links            |
| 15  | Support 10/25/50/100 kb and larger windows                                         | Complete              | Window selector; chromosome bounds enforced                           |
| 16  | List and visualize neighboring genes                                               | Complete              | Reference-coordinate arrows and gene table                            |
| 17  | Select a neighboring gene and investigate it                                       | Complete              | Browser navigation and API checks                                     |
| 18  | Query prepared data without rerunning GWAS                                         | Complete              | Read-only API; imports run separately                                 |
| 19  | Supply schema, initial datasets, reproducible prototype and validation examples    | Complete              | Repository, import audits, primary-table comparisons and review guide |
| 20  | Confirm the prototype answers the lab's questions in a supervised usability review | Pending Lisa's review | Walkthrough in review_guide.md                                        |

## Scientific limitations that remain

Only selected records in two original papers received a direct coordinate/statistic comparison. Ninety-three source rows matched within the documented rounding tolerance; 29 of these had trait or model metadata conflicts. A numeric match is not confirmation of sample size, allele orientation, causal interpretation or every source field. Other records are labelled as not independently reviewed.

The excluded records have a documented disposition: 2,291 unresolved reference labels, five zero p-values with unknown underflow/censoring semantics, and two out-of-bounds positions. They remain in quarantine. They have not been guessed into the searchable collection.

A shared public service, comprehensive original-study curation, additional data sources, additional expression atlases and optional CATCH analysis remain separate stages. A first Robertson baseline leaf-expression module is now implemented; see baseline_expression.md. No percentage in this checklist implies those are complete.

# RiceLENS: advisor feedback review

## Search behavior

The frozen RAP-DB annotation assigns the exact alias `SD1` to both `Os01g0883800` (display symbol SD1; gibberellin 20 oxidase 2, plant height) and `Os03g0130100` (display symbol OSD1; cinnamate:CoA ligase, submerged seed germination). This is a shared source alias, not substring matching or automatic replacement of SD1 with OSD1. The separate meiotic-process gene `Os02g0591500` also has OSD1 as a symbol.

Search now prioritizes an exact display-symbol match, retains all ambiguous candidates, and shows the matching identifier, gene ID, chromosome and source description. GS3 still offers both the grain-size locus and glutathione synthetase. Users must choose; RiceLENS does not silently infer which biological function they intended. Source aliases are retained rather than deleted to make a familiar symbol look unique. Regression tests cover SD1 and GS3 ambiguity.

## GWAS interpretation and coordinates

The 132,765 accepted rows are **curated published association records**, not complete SNP-level GWAS summary statistics. They should not be described as all tested SNPs, or as a uniformly significant collection: source studies use different reporting and significance thresholds, and RiceLENS has not independently verified every threshold. The interface and CSV preserve the reported values. There is no shared significance line or platform-wide claim of significance.

The primary-table review covers 93 coordinate/statistic matches in two publications, with 29 trait/model conflicts flagged. Verified FDR-adjusted values appear as orange diamonds; values with unverified adjustment status appear as dark grey circles. Reconstructing a p-value from a rounded minus-log value does not establish whether it was adjusted, so those records stay in the unverified category. Text and shapes accompany color, including plot detail labels and table statistic types. No record is newly called unadjusted without source evidence.

Gene coordinates use Os-Nipponbare-Reference-IRGSP-1.0. Regional queries require the same assembly ID for SNPs and genes, enforce chromosome boundaries, and use 1-based inclusive positions. SNP compatibility relies on the Atlas's explicit reference labels. Original reported coordinates remain available, but original-study remapping is not independently validated. There are 2,298 excluded rows: 2,291 unresolved reference labels, five zero p-values, and two out-of-range SNPs. This is a compatibility check against the imported reference, not proof of every source liftover.

## Website and accessibility

The website title, banner and README now use RiceLENS. Plot categories use shape as well as color, gene arrows include names and strand, and plot points include statistic details for keyboard focus and hover. Primary buttons and the lab badge use a darker blue with white text. This is an accessibility improvement, not a claim of a complete WCAG audit.

## Baseline expression pilot

The Robertson controls now supply 36,455 nuclei across three replicates and 14 original cell types. RiceLENS supports within-study expression ranks, detection percentages, neighboring-gene dot plots, numeric tables and expression CSV export. All available processed-matrix features are retained; 32,017 match the search annotation. See [baseline expression](baseline_expression.md).

## Remaining work

- Confirm private hosting account, access method and budget; no hosted deployment has been created.
- Expand beyond the three imported subsets (Robertson controls, Wang normal roots, and Li diploid pistils with scPlantDB labels); seed and lifecycle imports still need suitable files and mapping checks.
- Continue source-specific GWAS statistic and coordinate curation.
- Obtain advisor feedback on the revised search and plots.

See [expression sources](expression_sources.md). The earlier 19/20 first-stage checklist does not measure completion of this expanded expression scope.

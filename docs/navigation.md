# RiceLENS pages and navigation

RiceLENS provides four directly addressable pages behind the same hosted reviewer authentication:

- `/`: Home introduces the platform and links to each tool.
- `/variation`: Genetic Variation displays association-linked SNP observations, reported allele fields, coordinate provenance and the genomic neighborhood. This is not a complete variant catalog, and effect alleles are not relabelled as reference/alternate alleles.
- `/gwas`: GWAS Associations displays published traits, statistics, publications, neighboring genes and CSV export.
- `/expression`: Cell Type Expression displays baseline expression, study-specific organ and cell-type choices, dot plots, numeric tables, publication evidence and neighboring-gene comparisons.

Navigation preserves the selected gene, window and expression study. In-page tool changes also retain active filters and comparison selections. Direct reloads restore the bookmarked gene, window and expression study; other filters and comparison selections are session state. Legacy root URLs containing a gene open GWAS. Home bookmarks use `page=home` to preserve the landing page on reload. The browser Back and Forward controls switch tool pages.

Advanced association filters and neighboring-gene comparison controls start collapsed. Expression choices come only from the selected study response and reset when that response changes. The current imported subsets each cover one organ; these controls do not imply that the full Wang atlas is already integrated. Existing Li versions remain available while Wang and Robertson are the current validation priority.

Expression context filters affect displayed plots and tables, not the exported CSV. The export is explicitly labelled as the full study comparison. Displayed ranks remain the original within-study ranks; selecting a cell type does not recalculate its rank. Organ labels use explicit group organ metadata when supplied, otherwise the study's organ. Tissue-system labels such as Vasculature are not substituted for organs.

Missing features are distinguished from measured zero detection. Plots retain sequential blue expression colors, shaped association markers and numeric tables. Navigation includes current-page indicators, a skip link, keyboard focus styles, labelled controls and responsive wrapping.

The footer includes the known Wilkins Lab / University of Manitoba affiliation and city. The approved logo and full institutional address remain pending; no substitute logo or street address is invented. Tutorial remains deferred. No new datasets or scientific calculations were introduced by this layout update.

## Verification

The automated suite covers direct tool routes and authentication, alongside the existing gene, GWAS, expression and export tests. Browser checks cover gene-preserving navigation, direct expression reload, cell-type filtering, switching from Robertson to Wang with dependent-choice reset, and the separate variant table. Hosting auto-deployment is disabled, so pushing these changes does not update the live Render website.

## Wide genomic windows

Regions with more than 20 genes or a flank greater than 100 kb use a four-row compact gene track with only the selected gene labelled. All returned genes retain focusable/clickable marks and full descriptions in the readable hover/focus detail area. Closely spaced marks can overlap at overview scale; zooming narrows the supported flanking window around the selected gene. Smaller, sparse regions retain full gene labels. The plot still shows at most the 500 strongest matching association records, as before; its regional count and export describe the complete matching evidence.

Zoom controls step through the same supported windows as the flanking selector, preserving filters. The neighboring-gene search filters the complete returned table by identifier, symbol or description without changing the plotted genes. These display changes do not alter genomic coordinates, expression values or association statistics.

## Recorded-variant search

The Variation table has an explicitly submitted search plus an autocomplete list. Suggestions query all matching association records in the selected gene window and existing scientific filters, with up to eight suggestions showing position, record ID and trait. Commas, whitespace and case are normalized for matching, while typed input remains unchanged. Position text supports substring matching; a complete numeric `study:source-row` record ID selects that exact association. The IDs remain the source's association identifiers, not unique SNP identifiers.

Typing updates suggestions without changing the table. Search/Enter filters results across all pages and the evidence export; selecting a suggestion explicitly applies its exact record ID. Clear restores regional results, pagination restarts at page one, and no-match guidance suggests checking the input or widening the window. Arrow keys select suggestions, Enter applies the highlighted option and Escape closes the list. No silent spelling or coordinate correction occurs. The search is scoped to Variation; switching to GWAS removes its effect there while preserving the typed search for a return to Variation. Variant search text is session state and is not restored by a bookmark reload.

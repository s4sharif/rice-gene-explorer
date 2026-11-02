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

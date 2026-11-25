# GWAS-focused layout

This review branch implements Lisa’s revised Home, Locus and Genes, and GWAS workflow. Cell Expression and Tutorial intentionally have empty main content; existing expression databases and API methods are retained. Homolog Search has an independent input and explicit unavailable-data notice, not fabricated mappings.

## Navigation and selection

Navigation order: Home, Locus and Genes, GWAS, Cell Expression, Homolog Search, Download, Tutorial, Contact Us. `/variation` remains a compatible alias for `/locus`. The sidebar is sticky on desktop and becomes a full-width search block on small screens. It is shared across Locus and Genes, GWAS and Cell Expression. Home has no search field.

Gene searches retain ambiguous matches and descriptions. Coordinate searches accept `chr3:16729501`, `3:16729501-16735109` and comma-formatted equivalents, with inclusive IRGSP-1.0 bounds and a maximum 2 Mb input interval. Flanking windows are applied to the entire selected gene or interval. Coordinate searches never silently select the nearest gene. Upstream/downstream are strand-aware for genes; interval searches use lower/higher coordinate labels.

The URL preserves the selected ID or interval, window, trait, maximum P-value, position filter, sort direction, results page and duplicate-row choice. Clicking a neighbour updates the selected gene, sidebar, tables, plot, downloads and statistics; it retains the window and evidence filters. Neighbour choices and locus tables remain independent of GWAS trait filters.

## Counts and source limitations

Homepage counts are calculated by `/api/overview` from the serving databases. The present snapshot contains 828 study contexts in 152 publications, 357 trait/ontology entries and 132,765 source association rows. Study contexts are not presented as unique publications.

The displayed atlas count is 62,712 cells/nuclei in the selected annotation versions: Robertson controls (36,455), Wang normal-root subset (18,180), and Li author-labelled diploid pistil nuclei (8,077). The overlapping Li scPlantDB annotation version (9,769) is not added. The homepage disclosure lists the contributing versions. This does not claim full Wang eight-organ coverage, a union of all annotation versions, or resumed expression browsing.

The source `QTL` field is retained in results and CSV exports: 6,599 imported source rows contain a nonempty label. These rows remain GWAS evidence, with visible QTL names. This field does not establish an independently mapped linkage-QTL dataset. GWAS Atlas describes a curated significant-GWAS association resource: <https://ngdc.cncb.ac.cn/gwas/documentation>. No separate linkage-QTL data are inferred or invented.

The imported export does not verify reference/alternative allele orientation. Ref/Alt cells remain explicitly unreported; effect alleles, where provided, are labelled as effect alleles in the overlapping SNP details. Locus-table SNPs overlap each gene span, including introns; absence here is not absence of variation.

## Tables and downloads

Associations are sorted on the server before eight-row pagination, in either P-value direction. CSV includes every filtered record in that order. The plot always shows the strongest 500 matching records, independently of table sort. The searched span is highlighted with a blue band and dashed boundaries. FDR-reviewed evidence remains an orange diamond and adjustment-unverified evidence a grey circle; shape supplements colour.

The summary uses all filtered evidence, not the current table page: total records, the most represented trait (all ties), and minimum reported P-value. The minimum is not described as a common adjusted statistic. Exact duplicate grouping and original-row expansion remain available. Source details preserve coordinate-review qualifications and curation conflicts.

Download offers separate gene and association CSV exports for the retained selection. Gene exports are not restricted by GWAS filters. Contact details are copied from the public Wilkins Lab website, <https://www.oliviawilkinslab.com/>. The footer contains only the pending citation section; no citation has been invented.

## Review

Run `python3 -m unittest discover -s tests -v` and the standard formatting/comment checks. New regression cases cover coordinate bounds and empty intervals, sort order across pages and exports, complete-result summaries, alias/variant overlap, unverified allele orientation, QTL-label preservation, live counts, new routes and authentication. Existing expression tests remain in place although expression UI is paused.

No source data, normalized expression values, SQL schema or scientific source checksums are changed by this layout revision. Review locally before merging to main or manually deploying to Render.

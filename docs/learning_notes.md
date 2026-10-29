# Concepts used in this prototype

A database is a collection of organized tables. Each row describes one item; a key identifies it. A foreign key connects it to another table. Keeping publications in one table and associations in another avoids copying publication metadata into every analysis result and helps us check consistency.

A gene identifier behaves like a lookup key, but a gene symbol behaves more like a name: it can have synonyms or refer to several records. The alias table translates recognized identifiers into candidate genes. We ask the researcher to choose if a symbol is ambiguous.

An index helps the database locate relevant rows without scanning the whole table. Our important index is on assembly, chromosome and position, because this is what a regional search needs.

The import script performs the work once: read downloaded files, standardize allowed formats, validate fields, and insert records. The API handles searches later. This is the difference between data preparation and retrieval; clicking Search does not redo statistical analysis.

The browser sends a small request such as `/api/search?q=GS3`. Python returns JSON: structured text containing rows and labels. JavaScript turns that response into tables and a plot. SQLite stores the underlying evidence locally.

A regional join is an interval comparison. If a gene spans positions 100–200, a SNP at 150 overlaps it. A SNP at 205 is five bases beyond the gene. Whether it is upstream or downstream also depends on the gene's strand. These are geometric relationships, not causal conclusions.

The vertical plot axis is −log10(p). A p-value of 0.000001 becomes 6; smaller p-values appear higher. This improves visibility across orders of magnitude. It does not make results from different study designs directly comparable.

The tests check problems that could change an interpretation: wrong strand direction, off-by-one boundaries, mistaken symbol resolution, incompatible references or lost source rows. They complement source curation; passing software tests does not prove that every published association or upstream coordinate conversion is correct.

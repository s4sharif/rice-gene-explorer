# Project brief

## Current scope

Develop a searchable, gene-centered platform that connects existing plant GWAS evidence to genomic regions and, in a later phase, single-cell expression. Begin with rice. CATCH may become an optional analytical component but is not required for the first platform.

## Questions the platform should answer

1. What published GWAS evidence exists within or around the searched gene?
2. Which other genes are in the surrounding region?
3. In which cell types is the gene expressed once processed single-cell datasets are integrated?

## First-stage requirements

Accept a rice gene ID, symbol or recognized alternative identifier. Retrieve the gene's coordinates, strand and annotation version. Show regional trait associations, variant positions, distance and orientation relative to the gene, reported p-values, available allele/effect statistics, study provenance and publication links.

Allow 10, 25, 50 and 100 kb flanking windows and larger regions when useful. Display neighboring genes and GWAS positions together and let researchers select a neighboring gene to investigate it. Query processed, organized evidence rather than running a new GWAS during a search.

Explain that an association near a gene does not establish causality. Keep overlap, upstream and downstream distinct. Preserve assembly and annotation compatibility throughout the database.

## Review materials

The repository contains the working local prototype, initial source imports, schema, implementation approach, validation scripts and learning notes. The original project correspondence is not included; this brief records the functional requirements.

The next review should assess the first-stage interface and its source coverage, followed by study-level curation and selection of a processed rice single-cell dataset. Expression integration, public hosting and CATCH analysis are outside the implemented prototype.

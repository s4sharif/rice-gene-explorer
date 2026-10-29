# Database design

## Implemented relationships

```text
assembly ──< gene ──< alias >── source
    │          │
    │          └─ annotation source release
    └──< association >── source
              │
              ├── study >── publication
              └── trait
source ──< quarantine
```

| Table       | Purpose and important fields                                                                                                                                                                          |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| assembly    | Species, cultivar, assembly ID, coordinate convention                                                                                                                                                 |
| source      | Source URL, version/release and file SHA-256                                                                                                                                                          |
| gene        | Canonical RAP ID, assembly, annotation source, chromosome, inclusive start/end, strand, display symbol, description                                                                                   |
| alias       | Identifier or symbol, type, source and gene; one alias may point to several genes                                                                                                                     |
| publication | PMID, title, journal, author, publication label and original publication link                                                                                                                         |
| study       | Source StudyId plus publication/model/population/sample-size/environment/tissue context                                                                                                               |
| trait       | Source label and ontology accession; distinct label/accession pairs are not silently collapsed                                                                                                        |
| association | Source row, study, trait, assembly, chromosome, position, original/numeric p-value, available alleles/frequency/R², missing effect-size fields, original reference/location and complete original row |
| quarantine  | Source row, reason and complete original row for records excluded from coordinate searches                                                                                                            |

`schema.sql` is the executable specification. Foreign keys enforce relationships. Region indexes include assembly and chromosome, so matching a chromosome number alone is insufficient. Aliases use case-insensitive matching and return all candidate loci.

Association primary keys combine source GaP ID and physical source line because the source reuses some IDs. A study key combines its source ID and context hash to avoid overwriting records when a source ID has several contexts. This preserves evidence; it does not solve study independence or cross-source deduplication. The full raw source row preserves sampling conditions, year, QTL, reported gene and other fields not yet normalized for filtering.

The first prototype contains one species and one annotation release. Before adding assemblies or alternative annotations, change the gene key to an internal key with uniqueness over species/assembly/annotation/accession, and scope all aliases accordingly. Add source snapshot keys to association uniqueness before retaining multiple refreshes of the same dataset. Do not append a new version into this single-snapshot schema without that migration.

## Coordinate query

For gene span [s,e] and flanking window w:

```text
query_start = max(1, s - w)
query_end   = e + w
```

A neighboring gene overlaps the display region when `gene.start <= query_end AND gene.end >= query_start`. A SNP is included when its harmonized position lies within the inclusive display region on the same assembly and chromosome. The import checks every position against the official RAP-DB chromosome lengths; display regions clamp to those lengths. Out-of-bounds source positions remain quarantined.

Distance is zero for a SNP in [s,e]; otherwise it is the distance to the nearest gene boundary. For a positive-strand gene, left is upstream and right downstream; reverse those labels for negative-strand genes. Gene-body overlap is not exon/CDS overlap. No LD regions, credible sets, regulatory contacts, nearest-gene assignment or gene-based p-values are inferred.

## Missing values and evidence

Store absent effect size, other allele and absent source statistics as NULL. Display missing values explicitly. Do not replace missing effects with zero, infer alleles from position or substitute R² for effect size. Keep original p-value text; parse a separate numeric value for sorting/plotting. Unicode U+2212 is replaced with ASCII minus only during parsing. Zero/invalid p-values remain quarantined; plotting does not clip them to an invented number.

The current import accepts explicit IRGSP-1.0 `Ref ver` labels and uses exported Chr/Pos, following Atlas's documented harmonization. `Reported location` is retained even when its older assembly differs. The numeric label `1` is unresolved and excluded. No coordinates are guessed or lifted over locally. Curated published associations are not equivalent to full GWAS summary statistics.

## Next-stage extension (not implemented)

Add `expression_dataset`, `cell_type`, `expression_summary`, `gene_mapping` and `processing_method`. Each expression summary should key to a specific dataset, gene mapping and cell-type label and retain number of cells, number expressing, normalization method and summary unit. Tissue, cultivar, condition, assembly and annotation belong to the dataset. Preserve ambiguous/unmapped IDs with a mapping status. Later add sample/donor information if needed to interpret variability.

Add `genomic_region` and `region_evidence` for source-defined loci, credible sets or QTL intervals, preserving how boundaries were established. Keep a user-selected display window separate from a source-defined association region. Optional CATCH outputs should be stored as versioned analysis results linked to input datasets and methods, separately from published SNP association evidence.

## Evidence grouping and review overlays

`evidence_summary` stores the representative source row and repeat count for a fingerprint over every source field except GaP ID. `association_display` exposes one representative per exact evidence group. The original `association` and `quarantine` tables retain source rows unchanged. The API offers `include_repeats=true` to retrieve all accepted source rows.

A versioned `publication_review.json` overlay records selected original-paper comparisons by association ID. It identifies checked statistic types and metadata conflicts without overwriting original traits or models. `validate_data.py` checks that the overlay refers to the same GWAS source checksum and existing accepted records. These selected checks do not validate unmatched records.

The original-row endpoint `/api/association?id=...` returns the complete imported source assertion and any primary-table review. Study/publication/trait/statistic/position filters apply consistently to the table, chart and CSV; gene neighbors continue to describe the selected reference window.

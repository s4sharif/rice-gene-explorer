# Walkthrough for Lisa

## Start

Follow the README's clone/download/build instructions, then open the local app. The first download uses about 21 MB of compressed annotation and GWAS data. Single-cell analysis is deliberately outside this first-stage review.

## Five-minute functional review

1. Search `GS3`. Two annotation loci match this symbol; choose `Os03g0407400` rather than assuming which match was intended.
2. At ±25 kb, inspect the three neighboring gene spans and the evidence table. The current snapshot has 240 source rows represented by 237 unique evidence records. Eight records appear per page. Try Next and Go to page.
3. Select publication `33367255`. Inspect chr3:16,733,441 at 3.40E-43. The primary paper identifies this statistic as FDR-adjusted and the trait as grain length; the imported source says grain width. The warning preserves and explains the discrepancy.
4. Reset filters, select upstream or overlap, then a p-value threshold. These filter reported records; a threshold is not a universal genome-wide significance criterion. Use the study selector to examine a specific study context.
5. Increase the window to ±50 kb. Click a neighboring gene in the table or chart and inspect its own coordinates and regional evidence.
6. Search `Os05g0187500` (GW5) and select publication `33367255`. The paper's Table 1 reports chr5:5,371,686 with FDR-adjusted p = 9.3E-34 for grain width; the import retains that coordinate and statistic. This is regional evidence around GW5, not proof that it is causal.
7. Search `Os06g0133000` (Waxy/WX1) and select the same publication. Table 1 reports chr6:1,765,761 with FDR-adjusted p = 6.4E-08 for endosperm type. Inspect its location relative to the gene rather than assigning it to the nearest gene automatically.
8. Try `LOC_Os01g01010` to test an official MSU cross-reference. Try a nonexistent identifier and an impossible trait filter to check empty results.
9. Export CSV with a filter active. The export contains all matching records across pages and preserves review notes, statistic type, original-source row and coordinates. Enable repeated source rows if every original row is needed.

## What the primary-source check establishes

`audit_publications.py` compares imported records with Table 1 of [Kumar et al., PMID 33367255](https://pmc.ncbi.nlm.nih.gov/articles/PMC7748012/) and [Zhao et al., PMID 29616055](https://pmc.ncbi.nlm.nih.gov/articles/PMC5864933/). Article identities should be read from their linked publications; matching is based on coordinate and statistic, with 2% relative tolerance for rounded values. For the latter table, p-values are reconstructed from the reported rounded −log10(p), and table footnotes distinguish model columns.

Ninety-three source records match the compared coordinates/statistics. Twenty-nine have conflicting trait or model labels in the export. The application retains the source assertion and flags the conflict. No unverified correction is written into the source data. Unmatched records, including most of the database, remain outside this primary-table validation scope.

Repeat the audit with:

```sh
python3 audit_publications.py --download
```

The article XML downloads remain local. The review report records source URLs, file hashes, comparison tolerance and per-record outcomes. Updating a source snapshot requires repeating the audit; the validation script rejects a report tied to a different association-file checksum.

## Data-quality disposition

- All 36,061 gene spans and 132,765 accepted association positions pass reference chromosome bounds, using the [official RAP-DB chromosome list](https://rapdb.dna.naro.go.jp/genome-wide_variations/genome-wide_var_workflow_v2.html).
- The 2,298 excluded source rows retain their full original fields and a reason. No implicit reference conversion, replacement p-value or guessed allele is applied.
- There are 129,962 distinct evidence records after grouping 2,803 exact repeated rows across 2,702 groups. Grouping compares every original field except the source GaP identifier; different models, conditions or other source fields remain separate. This is not an independence assessment or cross-study meta-analysis.
- Only 56 source records contain an unambiguous A/C/G/T effect-allele label; another 22 say ambiguous. No export rows provide effect-allele frequency, beta or the other allele. Missing values remain explicit.

## Review decision

Confirm whether gene lookup, regional navigation, evidence filters and exports answer the first-stage questions. Record any corrections to identifier handling, wording or desired default windows. Then agree on a rice tissue/trait and a processed single-cell dataset for the next stage. The functional checklist records 19 completed checkpoints and leaves supervisor acceptance open.

# Dataset and calculation directory

Open this folder locally to browse descriptive dataset names. The naming convention is `author_year_parameter`, with annotation versions included when necessary. Shared reference resources use their provider name.

| Folder                                         | Scope                                                    | Calculation status                                       |
| ---------------------------------------------- | -------------------------------------------------------- | -------------------------------------------------------- |
| `robertson_2026_well_watered_leaf`             | Three well-watered leaf controls                         | Baseline expression integrated                           |
| `wang_2025_normal_root`                        | Two normal Nipponbare root replicates                    | Baseline expression integrated                           |
| `li_2023_diploid_pistil_author_annotations`    | Diploid pistil, original author cell labels              | Baseline expression integrated                           |
| `li_2023_diploid_pistil_scplantdb_annotations` | Same diploid libraries, curated cell labels              | Alternative annotation version, not independent evidence |
| `jing_2026_gene_mapping_metadata`              | Supporting annotations and conservative gene mapping     | Metadata only; expression pending                        |
| `liu_2026_seed_data_review`                    | Developing-seed source availability review               | Expression pending                                       |
| `rapdb_gene_coordinates`                       | Shared RAP-DB gene reference and GWAS Atlas associations | Gene search, windows and association evidence            |

## Find the files behind a result

Open [calculation_index.csv](calculation_index.csv) in a spreadsheet, or search [inventory.json](inventory.json). Each row identifies a descriptive path, the original path, its role and its calculation or usage. These are file inventories: downloaded archives, supporting code and excluded-sample metadata can be present without contributing cells to the calculations. Inclusion is determined by the input, selection and extraction manifests, not by file presence.

For each integrated expression dataset:

1. The `inputs` JSON identifies downloaded sources, accessions and SHA-256 checksums.
2. Selection, cell metadata and extraction files identify retained samples and cell types. Original download objects may contain excluded treatments.
3. `features.tsv` retains the available feature universe; `summaries.tsv.gz` contains per-sample/cell-type cell counts, expressing-cell counts and sums of normalized log expression.
4. `pipelines/build_expression.py` maps unambiguous gene IDs and aggregates the summaries into SQLite. Mean expression is the cell-weighted mean of `log1p(10000 * count / library_size)`; percent expressing is `100 * n_expressing / n_cells`. Ranks are calculated within the selected study.
5. Validation reports document independent recomputation checks. They are not input expression data.

The SQLite `sample_summary` table preserves sample-level summaries, `summary` preserves gene/cell-type aggregates, `feature` records ID mapping and `metadata` records study provenance. Hosted snapshots omit `sample_summary` to reduce storage; local databases retain it. Neither SQLite version stores the complete original single-cell matrix.

See [Robertson methods](../docs/baseline_expression.md), [Wang methods](../docs/wang_expression.md), [Li curated methods](../docs/li_expression.md), [Li author methods](../docs/li_author_expression.md), [Jing metadata](../docs/jing_metadata.md) and [Liu requirements](../docs/liu_data_requirements.md) for exact source selection and limitations. Shared GWAS/RAP inputs and checksums are recorded in [provenance](../data/provenance.json); their import is implemented in `pipelines/build_database.py`.

## Storage and regeneration

The descriptive files in study folders are **relative symbolic links**, not additional copies. Original source filenames, runtime paths, database IDs and website bookmarks remain intact. This organization therefore adds negligible storage and does not change calculations. Local links and their large targets are not uploaded to GitHub; the index, documentation and generator are versioned.

Recreate or refresh the folders after fetching/building data:

```sh
python3 scripts/organize_datasets.py
```

The index reflects files present when generated. A fresh clone must download/build its data before all links can be recreated. Missing Jing/Liu expression data are explicitly pending; their review files are not expression datasets.

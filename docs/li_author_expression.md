# Li pistil expression with original author annotations

The supplied `mmc3.xlsx` is **Table S2 from Li et al. (2023)**, _Single-nucleus sequencing deciphers developmental trajectories in rice pistils_, DOI [10.1016/j.devcel.2023.03.004](https://doi.org/10.1016/j.devcel.2023.03.004). It is not Liu's seed atlas. The supplied `mmc7.pdf` contains that Li article and its supplementary figures.

## Included data

The `Rice_Pistil_index` worksheet contains all 8,173 nuclei in the authors' combined diploid atlas. The article explains its selection of 6,004 nuclei from replicate 1 and 2,169 from replicate 2 (article page 2). Every published sample/barcode key matches the original GEO count matrices exactly.

The import excludes 96 nuclei explicitly labelled `sta` (stamen residue), retaining **8,077 nuclei**: 5,939 from GSM6346349 and 2,138 from GSM6346350. It retains eight author cell groups. Codes expand only as defined in the article's Figure 3 legend: chalaza, epidermis, nucellus, ovule, ovary wall and stigma. C12 and C13 remain numbered author cluster labels. Original cluster numbers and temporal-map labels remain in the intermediate annotation table, without treating a pseudotime label as a physical cell type.

All **38,726 original count-matrix features** are processed, with **35,797** unambiguous supported RAP mappings. This includes all available features, not just marker genes. Mean expression is mean(log1p(10,000 × RNA count / per-nucleus RNA library size)), including zeros. Detection is original count > 0. Replicates are pooled with nuclei-weighted means, as in the other local expression imports. No stress or differential-expression analysis is performed.

## Reproduce

Place the supplied publisher supplement at `data/raw/li_author/mmc3.xlsx`; its checksum is recorded in `data/expression_li_author_inputs.json`.

```sh
python3 -m pipelines.fetch_expression --study li_author
python3 -m pipelines.prepare_li_author
Rscript scripts/summarize_li.R author
python3 -m pipelines.build_expression --study li_author
python3 -m pipelines.validate_expression --study li_author
```

The fetch command retrieves the original GEO inputs and verifies the supplied workbook. The parser reads literal XLSX values without Excel or added Python packages. The join fails on duplicate, missing or inconsistent sample/barcode keys.

Outputs are `data/databases/expression_li_author.sqlite`, the tracked manifest and join/independent-validation reports, plus intermediate matched-cell and sparse summary tables under the ignored raw directory. The author version appears separately in the website study selector and CSV exports.

## Validation and interpretation

All 8,173 source barcodes matched. Four genes were independently recalculated from the original MEX matrices across both samples and all eight retained groups: **64 gene/sample/group checks passed**. SQLite integrity and numeric bounds also passed.

The existing scPlantDB-curated import is preserved. These two annotation versions use the same original libraries and **are not independent evidence or biological replication**. Their cell selections and labels differ; RiceLENS keeps them separate. The author dataset remains enriched for diploid nuclei by the experiment's sorting design and does not estimate natural cell-type proportions. Missing/unmapped features are unavailable, not measured zero.

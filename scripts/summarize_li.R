# Summarize original diploid pistil counts using exactly matched curated cell
# annotations.
library(Matrix)
library(jsonlite)
root <- "data/raw/li"
metadata <- read.delim(file.path(root, "matched_cells.tsv"), stringsAsFactors = FALSE)
controls <- c("GSM6346349", "GSM6346350")
stopifnot(all(metadata$sample %in% controls), !anyDuplicated(paste(metadata$sample, metadata$count_barcode)),
  !anyNA(metadata$cell_type), all(nchar(metadata$cell_type) > 0))
config <- fromJSON(file.path(root, "selection.json"))
output <- gzfile(file.path(root, "summaries.tsv.gz"), "wt")
writeLines("feature_id\tsample\tcell_type\tn_cells\tn_expressing\tsum_log", output)
features_all <- character()
features_per_sample <- list()
sample_groups <- list()
join_report <- list()
for (sample in controls) {
  prefix <- file.path(root, config$prefixes[[sample]])
  features <- read.delim(gzfile(paste0(prefix, "_features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)[[1]]
  barcodes <- readLines(gzfile(paste0(prefix, "_barcodes.tsv.gz")))
  stopifnot(!anyDuplicated(features), !anyDuplicated(barcodes))
  mat <- as(readMM(gzfile(paste0(prefix, "_matrix.mtx.gz"))), "CsparseMatrix")
  stopifnot(nrow(mat) == length(features), ncol(mat) == length(barcodes), all(mat@x > 0),
    all(mat@x == floor(mat@x)))
  cells <- metadata[metadata$sample == sample, , drop = FALSE]
  indices <- match(cells$count_barcode, barcodes)
  stopifnot(!anyNA(indices), nrow(cells) > 0)
  mat <- mat[, indices, drop = FALSE]
  sizes <- colSums(mat)
  stopifnot(all(sizes > 0))
  features_all <- union(features_all, features)
  features_per_sample[[sample]] <- length(features)
  join_report[[sample]] <- list(original_matrix_nuclei = length(barcodes), matched_curated_nuclei = nrow(cells),
    unmatched_annotation_nuclei = sum(is.na(indices)), excluded_without_curated_label = length(barcodes) -
      nrow(cells))
  for (label in sort(unique(cells$cell_type))) {
    ix <- which(cells$cell_type == label)
    part <- mat[, ix, drop = FALSE]
    detected <- rowSums(part > 0)
    part@x <- log1p(part@x * rep(10000/sizes[ix], diff(part@p)))
    sums <- rowSums(part)
    write.table(data.frame(feature_id = features, sample = sample, cell_type = label,
      n_cells = length(ix), n_expressing = detected, sum_log = sums), output, sep = "\t",
      row.names = FALSE, col.names = FALSE, quote = FALSE)
    sample_groups[[length(sample_groups) + 1]] <- list(sample = sample, cell_type = label,
      tissue_system = "Pistil (2N)", n_cells = length(ix))
  }
  cat(sample, nrow(cells), "matched nuclei", length(features), "original features\n")
  rm(mat, part)
  gc()
}
close(output)
groups <- list()
for (label in sort(unique(metadata$cell_type))) groups[[length(groups) + 1]] <- list(cell_type = label,
  tissue_system = "Pistil (2N)", n_cells = sum(metadata$cell_type == label))
write.table(data.frame(feature_id = sort(features_all)), file.path(root, "features.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write_json(list(eligible_samples = controls, n_cells = nrow(metadata), groups = groups, sample_groups = sample_groups,
  features_per_sample = features_per_sample, join_report = join_report, annotation_source = "scPlantDB SRP386976",
  selection = "Two wild-type diploid libraries only; haploid library excluded", normalization = "log1p CP10K from all original GEO matrix features"),
  file.path(root, "extraction.json"), pretty = TRUE, auto_unbox = TRUE, digits = 15)

# Summarize all available processed RNA-count features in annotated well-watered nuclei.
library(Matrix)
library(jsonlite)

root <- "data/raw/expression"
metadata <- read.delim(file.path(root, "cell_metadata.tsv"), check.names = FALSE,
                       stringsAsFactors = FALSE, row.names = 1)
stopifnot(!anyDuplicated(rownames(metadata)))
print(names(metadata))
# The sample and barcode columns are selected explicitly after source metadata inspection.
config <- fromJSON(file.path(root, "selection.json"))
stopifnot(all(c(config$sample_column, "CellAnnotation", "TissueSystem") %in% names(metadata)))
controls <- c("Well.Watered.Rep1", "Well.Watered.Rep2", "Well.Watered.Rep3")
metadata <- metadata[metadata[[config$sample_column]] %in% unlist(config$sample_values), , drop = FALSE]
stopifnot(all(metadata$Condition == config$required_condition))
stopifnot(nrow(metadata) > 0, !anyNA(metadata$CellAnnotation), !anyNA(metadata$TissueSystem))
stopifnot(all(nchar(metadata$CellAnnotation) > 0))
output <- gzfile(file.path(root, "summaries.tsv.gz"), "wt")
writeLines("feature_id\tsample\tcell_type\tn_cells\tn_expressing\tsum_log", output)
feature_union <- character()
feature_by_sample <- list()
sample_groups <- list()
join_report <- list()
retained_metadata <- list()

for (index in seq_along(controls)) {
  sample <- controls[index]
  acc <- paste0("GSM", 9373071 + index)
  prefix <- file.path(root, paste0(acc, "_", sample))
  features <- read.delim(gzfile(paste0(prefix, ".genes.tsv.gz")), header = FALSE,
                         stringsAsFactors = FALSE)[[1]]
  barcodes <- readLines(gzfile(paste0(prefix, ".barcodes.tsv.gz")))
  stopifnot(!anyDuplicated(features), !anyDuplicated(barcodes))
  mat <- as(readMM(gzfile(paste0(prefix, ".matrix.mtx.gz"))), "CsparseMatrix")
  stopifnot(nrow(mat) == length(features), ncol(mat) == length(barcodes),
            all(mat@x > 0), all(mat@x == floor(mat@x)))
  cells <- metadata[metadata[[config$sample_column]] == config$sample_values[[sample]], , drop = FALSE]
  cell_barcodes <- if (config$barcode_column == "rownames") rownames(cells) else cells[[config$barcode_column]]
  if (nchar(config$barcode_remove_prefix) > 0) {
    cell_barcodes <- sub(config$barcode_remove_prefix, "", cell_barcodes)
  }
  stopifnot(!anyDuplicated(cell_barcodes))
  column_indices <- match(cell_barcodes, barcodes)
  stopifnot(!anyNA(column_indices))
  mat <- mat[, column_indices, drop = FALSE]
  rownames(mat) <- features
  colnames(mat) <- rownames(cells)
  library_sizes <- colSums(mat)
  stopifnot(all(library_sizes > 0))
  feature_union <- union(feature_union, features)
  feature_by_sample[[sample]] <- features
  join_report[[sample]] <- list(matrix_cells = length(barcodes), annotated_baseline_cells = nrow(cells),
                                unmatched_annotation_cells = sum(is.na(column_indices)),
                                excluded_unannotated_matrix_cells = length(barcodes) - nrow(cells))
  retained_metadata[[sample]] <- cells
  for (cell_type in sort(unique(cells$CellAnnotation))) {
    indices <- which(cells$CellAnnotation == cell_type)
    tissues <- unique(cells$TissueSystem[indices])
    stopifnot(length(tissues) == 1)
    counts <- mat[, indices, drop = FALSE]
    detected <- rowSums(counts > 0)
    counts@x <- log1p(counts@x * rep(10000 / library_sizes[indices], diff(counts@p)))
    totals <- rowSums(counts)
    summary <- data.frame(feature_id = features, sample = sample, cell_type = cell_type,
                          n_cells = length(indices), n_expressing = detected, sum_log = totals)
    write.table(summary, output, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    sample_groups[[length(sample_groups) + 1]] <- list(sample = sample, cell_type = cell_type,
                                                      tissue_system = tissues, n_cells = length(indices))
  }
  cat(sample, "annotated nuclei", nrow(cells), "features", length(features), "\n")
  rm(mat, counts)
  gc()
}
close(output)
retained <- do.call(rbind, retained_metadata)
groups <- list()
for (label in sort(unique(retained$CellAnnotation))) {
  ix <- which(retained$CellAnnotation == label)
  tissue <- unique(retained$TissueSystem[ix])
  stopifnot(length(tissue) == 1)
  groups[[length(groups) + 1]] <- list(cell_type = label, tissue_system = tissue, n_cells = length(ix))
}
write.table(data.frame(feature_id = sort(feature_union)), file.path(root, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write_json(list(eligible_samples = controls, n_cells = nrow(retained), groups = groups,
                sample_groups = sample_groups, join_report = join_report,
                features_per_sample = lapply(feature_by_sample, length),
                normalization = "log1p CP10K; library sizes use all features available in each processed MEX matrix",
                feature_missingness = "Features absent from a sample matrix are unavailable for that sample, not imputed zero"),
           file.path(root, "extraction.json"), pretty = TRUE, auto_unbox = TRUE, digits = 15)

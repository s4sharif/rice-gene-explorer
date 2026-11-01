# Import only the published normal Nipponbare roots from the annotated Wang object.
library(Matrix)
library(jsonlite)
root <- "data/raw/wang"
x <- readRDS(gzcon(gzfile(file.path(root, "GSE273875_LTP_NIP_label_240622.rds.gz"), "rb")))
metadata <- attr(x, "meta.data")
stopifnot(!is.null(metadata), !anyDuplicated(rownames(metadata)), all(c("orig.ident", "condition",
  "cell_type") %in% names(metadata)))
print(table(metadata$orig.ident, metadata$condition))
assay <- attr(x, "assays")[["RNA"]]
counts <- attr(assay, "counts")
if (is.null(counts)) {
  layers <- attr(assay, "layers")
  stopifnot("counts" %in% names(layers))
  counts <- layers[["counts"]]
  features <- attr(attr(assay, "features"), "dimnames")[[1]]
  cells <- attr(attr(assay, "cells"), "dimnames")[[1]]
  stopifnot(length(features) == nrow(counts), length(cells) == ncol(counts))
  rownames(counts) <- features
  colnames(counts) <- cells
}
stopifnot(inherits(counts, "sparseMatrix"), all(counts@x >= 0), all(counts@x == floor(counts@x)),
  !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)), all(rownames(metadata) %in%
    colnames(counts)))
controls <- c("NIP_G_1", "NIP_G_2")
keep <- metadata$orig.ident %in% controls
stopifnot(all(metadata$condition[keep] == "WT"), !anyNA(metadata$cell_type[keep]), all(nchar(as.character(metadata$cell_type[keep])) >
  0))
metadata <- metadata[keep, , drop = FALSE]
counts <- counts[, rownames(metadata), drop = FALSE]
rm(x, assay)
if (exists("layers")) rm(layers)
gc()
write.table(data.frame(feature_id = rownames(counts)), file.path(root, "features.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE)
output <- gzfile(file.path(root, "summaries.tsv.gz"), "wt")
writeLines("feature_id\tsample\tcell_type\tn_cells\tn_expressing\tsum_log", output)
sample_groups <- list()
features_per_sample <- list()
for (sample in controls) {
  md <- metadata[metadata$orig.ident == sample, , drop = FALSE]
  mat <- counts[, rownames(md), drop = FALSE]
  sizes <- colSums(mat)
  stopifnot(all(sizes > 0))
  features_per_sample[[sample]] <- nrow(mat)
  writeMM(mat, file.path(root, paste0(sample, ".matrix.mtx")))
  status <- system2("gzip", c("-f", file.path(root, paste0(sample, ".matrix.mtx"))))
  stopifnot(status == 0)
  feature_connection <- gzfile(file.path(root, paste0(sample, ".features.tsv.gz")), "wt")
  writeLines(rownames(mat), feature_connection)
  close(feature_connection)
  write.table(data.frame(barcode = rownames(md), cell_type = as.character(md$cell_type)),
    file.path(root, paste0(sample, ".cells.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  for (label in sort(unique(as.character(md$cell_type)))) {
    cols <- which(md$cell_type == label)
    part <- mat[, cols, drop = FALSE]
    detected <- rowSums(part > 0)
    part@x <- log1p(part@x * rep(10000/sizes[cols], diff(part@p)))
    totals <- rowSums(part)
    write.table(data.frame(feature_id = rownames(mat), sample = sample, cell_type = label,
      n_cells = length(cols), n_expressing = detected, sum_log = totals), output, sep = "\t",
      row.names = FALSE, col.names = FALSE, quote = FALSE)
    sample_groups[[length(sample_groups) + 1]] <- list(sample = sample, cell_type = label,
      tissue_system = "Root", n_cells = length(cols))
  }
  cat(sample, ncol(mat), "normal cells", nrow(mat), "RNA-count features\n", sep = " ")
  rm(mat, part)
  gc()
}
close(output)
groups <- list()
for (label in sort(unique(as.character(metadata$cell_type)))) {
  groups[[length(groups) + 1]] <- list(cell_type = label, tissue_system = "Root", n_cells = sum(metadata$cell_type ==
    label))
}
write_json(list(eligible_samples = controls, n_cells = nrow(metadata), groups = groups, sample_groups = sample_groups,
  features_per_sample = features_per_sample, normalization = "log1p CP10K from original RNA counts",
  selection = "orig.ident NIP_G_1/NIP_G_2 and condition WT; both mutant replicates excluded",
  feature_missingness = "All features in the joined original RNA count assay; no marker/HVG selection"),
  file.path(root, "extraction.json"), pretty = TRUE, auto_unbox = TRUE, digits = 15)

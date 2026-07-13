#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(limma)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

abundance_file <- args[1]
metadata_file <- args[2]
out_raw <- args[3]
out_batch_corrected <- args[4]

read_table_auto <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    read_xlsx(path)
  } else {
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  }
}

clr_transform <- function(mat, pseudocount = 1) {
  log_mat <- log(mat + pseudocount)
  sweep(log_mat, 2, colMeans(log_mat, na.rm = TRUE), "-")
}

plot_pca <- function(mat, meta, outfile, title) {
  pca <- prcomp(t(mat), scale. = TRUE)

  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

  df <- as.data.frame(pca$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(meta, by = "sample")

  p <- ggplot(df, aes(PC1, PC2, color = condition, shape = batch)) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = title,
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(outfile, p, width = 7, height = 5)
}

abund <- read_table_auto(abundance_file)
meta <- read_table_auto(metadata_file)

required_meta <- c("sample", "condition", "batch")
missing_meta <- setdiff(required_meta, colnames(meta))

if (length(missing_meta) > 0) {
  stop("metadata is missing required columns: ", paste(missing_meta, collapse = ", "))
}

taxon_col <- colnames(abund)[1]

mat <- abund
rownames(mat) <- mat[[taxon_col]]
mat[[taxon_col]] <- NULL
mat <- as.matrix(mat)
mode(mat) <- "numeric"

samples <- intersect(colnames(mat), meta$sample)

if (length(samples) < 3) {
  stop("Fewer than 3 matching samples between abundance and metadata")
}

mat <- mat[, samples, drop = FALSE]
meta <- meta[match(samples, meta$sample), , drop = FALSE]

# Remove taxa with zero variance, otherwise prcomp can fail.
keep <- apply(mat, 1, var, na.rm = TRUE) > 0
mat <- mat[keep, , drop = FALSE]

mat_clr <- clr_transform(mat)

plot_pca(
  mat_clr,
  meta,
  out_raw,
  "CLR PCA, KIS excluded"
)

design <- model.matrix(~ condition, data = meta)

batch_corrected <- removeBatchEffect(
  mat_clr,
  batch = meta$batch,
  design = design
)

plot_pca(
  batch_corrected,
  meta,
  out_batch_corrected,
  "CLR PCA after batch correction, KIS excluded"
)

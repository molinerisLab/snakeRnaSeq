#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

abundance_file <- args[1]
metadata_file <- args[2]
top_table_file <- args[3]
out_heatmap <- args[4]
top_n <- as.integer(args[5])

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

abund <- read_table_auto(abundance_file)
meta <- read_table_auto(metadata_file)
top <- read_table_auto(top_table_file)

if ("Geneid" %in% colnames(abund) && !"GeneID" %in% colnames(abund)) {
  colnames(abund)[colnames(abund) == "Geneid"] <- "GeneID"
}

if ("Geneid" %in% colnames(top) && !"GeneID" %in% colnames(top)) {
  colnames(top)[colnames(top) == "Geneid"] <- "GeneID"
}

required_meta <- c("sample", "condition", "batch")
missing_meta <- setdiff(required_meta, colnames(meta))
if (length(missing_meta) > 0) {
  stop("metadata is missing required columns: ", paste(missing_meta, collapse = ", "))
}

required_top <- c("GeneID", "Pvalue")
missing_top <- setdiff(required_top, colnames(top))
if (length(missing_top) > 0) {
  stop("top_table is missing required columns: ", paste(missing_top, collapse = ", "))
}

taxon_col <- colnames(abund)[1]

mat_df <- abund
rownames(mat_df) <- mat_df[[taxon_col]]
mat_df[[taxon_col]] <- NULL

mat <- as.matrix(mat_df)
mode(mat) <- "numeric"

samples <- intersect(colnames(mat), meta$sample)

if (length(samples) < 2) {
  stop("Fewer than 2 matching samples between abundance and metadata")
}

mat <- mat[, samples, drop = FALSE]
meta <- meta[match(samples, meta$sample), , drop = FALSE]

top$Pvalue <- as.numeric(top$Pvalue)

top_taxa <- top %>%
  arrange(Pvalue) %>%
  head(top_n) %>%
  pull(GeneID)

top_taxa <- intersect(top_taxa, rownames(mat))

if (length(top_taxa) < 2) {
  stop("Fewer than 2 top taxa from top_table were found in abundance matrix")
}

mat_clr <- clr_transform(mat)
mat_top <- mat_clr[top_taxa, , drop = FALSE]

# Clean taxa labels
clean_taxa <- rownames(mat_top)

clean_taxa <- gsub("_\\d+_S$", "", clean_taxa)
clean_taxa <- gsub("_", " ", clean_taxa)

rownames(mat_top) <- clean_taxa

mat_scaled <- t(scale(t(mat_top)))
mat_scaled[is.na(mat_scaled)] <- 0

annotation_col <- meta[, c("sample", "condition", "batch"), drop = FALSE]
rownames(annotation_col) <- annotation_col$sample
annotation_col$sample <- NULL
annotation_col <- annotation_col[colnames(mat_scaled), , drop = FALSE]

pdf(out_heatmap, width = 10, height = 9)
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 7,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  main = paste0("Top ", length(top_taxa), " taxa by nominal P value")
)
dev.off()
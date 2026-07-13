#!/usr/bin/env Rscript

library(tidyverse)

# Input arguments from Snakemake
clr_file <- snakemake@input[["clr"]]
meta_file <- snakemake@input[["metadata"]]
out_pca <- snakemake@output[["clr_pca"]]

meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  mutate(
    group = factor(group),
    condition = factor(condition),
    batch = factor(batch)
  )

load_matrix <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)

  if (!"Geneid" %in% colnames(df)) {
    stop(
      "Expected a 'Geneid' column. Found: ",
      paste(colnames(df), collapse = ", ")
    )
  }

  mat <- df %>%
    column_to_rownames("Geneid") %>%
    as.matrix()

  storage.mode(mat) <- "numeric"

  if (any(!is.finite(mat))) {
    stop("The abundance matrix contains NA, NaN, or infinite values.")
  }

  mat
}

# Load CLR abundance matrix
clr_mat <- load_matrix(clr_file)

# Check that sample names are unique
if (anyDuplicated(colnames(clr_mat))) {
  stop("The CLR matrix contains duplicated sample names.")
}

if (anyDuplicated(meta$sample)) {
  stop("The metadata contains duplicated sample names.")
}

# Check correspondence between matrix and metadata
missing_metadata <- setdiff(colnames(clr_mat), meta$sample)

if (length(missing_metadata) > 0) {
  stop(
    "Samples missing from metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# Retain and order metadata to match the abundance matrix
meta <- meta %>%
  filter(sample %in% colnames(clr_mat)) %>%
  slice(match(colnames(clr_mat), sample))

if (!identical(meta$sample, colnames(clr_mat))) {
  stop("Metadata could not be aligned with the CLR matrix.")
}

# Remove invariant taxa, because scale.=TRUE cannot scale zero-variance rows
taxon_sd <- apply(clr_mat, 1, sd)

if (any(!is.finite(taxon_sd))) {
  stop("Non-finite taxon standard deviations encountered.")
}

clr_mat <- clr_mat[taxon_sd > 0, , drop = FALSE]

if (nrow(clr_mat) < 2) {
  stop("Fewer than two variable taxa remain for PCA.")
}

# Samples are observations; taxa are variables
pca <- prcomp(
  t(clr_mat),
  center = TRUE,
  scale. = TRUE
)

variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta, by = "sample")

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    shape = batch
  )
) +
  geom_point(size = 4) +
  theme_bw() +
  labs(
    title = "PCA of CLR-transformed abundances",
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2]),
    color = "Condition",
    shape = "Batch"
  )

ggsave(
  filename = out_pca,
  plot = p_pca,
  width = 7,
  height = 5,
  units = "in"
)
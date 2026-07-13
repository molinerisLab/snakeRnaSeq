#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# -----------------------------
# Read Snakemake inputs/outputs
# -----------------------------

input_matrix <- snakemake@input[["matrix"]]

output_relative <- snakemake@output[["relative"]]
output_clr <- snakemake@output[["clr"]]

log_file <- snakemake@log[[1]]

pseudocount <- snakemake@params[["pseudocount"]]

if (is.null(pseudocount) || is.na(pseudocount) || pseudocount <= 0) {
  stop("Parameter 'pseudocount' must be a positive numeric value.")
}

# -----------------------------
# Set up logging
# -----------------------------

dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")

sink(log_con, split = TRUE)
sink(log_con, type = "message")

on.exit({
  while (sink.number(type = "message") > 0) {
    sink(type = "message")
  }
  while (sink.number(type = "output") > 0) {
    sink(type = "output")
  }
  close(log_con)
}, add = TRUE)

message("Starting abundance normalization")
message("Input matrix: ", input_matrix)
message("Relative abundance output: ", output_relative)
message("CLR output: ", output_clr)
message("CLR pseudocount: ", pseudocount)

# -----------------------------
# Create output directories
# -----------------------------

dir.create(dirname(output_relative), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_clr), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read count table
# -----------------------------

counts_raw <- read_tsv(input_matrix, show_col_types = FALSE)

if (nrow(counts_raw) == 0) {
  stop("Input table has zero rows.")
}

if (ncol(counts_raw) < 2) {
  stop("Input table must contain one feature column and at least one sample column.")
}

# Accept common feature/taxon column names.
feature_candidates <- c("Geneid", "GeneID", "geneid", "FeatureID", "feature_id", "taxon", "Taxon", "name")

feature_col <- intersect(feature_candidates, names(counts_raw))[1]

if (is.na(feature_col)) {
  feature_col <- names(counts_raw)[1]
  warning("No standard feature column found. Using first column as feature column: ", feature_col)
}

message("Feature column: ", feature_col)

# Preserve optional taxonomy annotation columns if present.
annotation_cols <- intersect(
  c("taxonomy_id", "taxonomy_lvl", "taxid", "rank"),
  names(counts_raw)
)

id_cols <- unique(c(feature_col, annotation_cols))

sample_cols <- setdiff(names(counts_raw), id_cols)

if (length(sample_cols) == 0) {
  stop("No sample columns found after excluding feature and annotation columns.")
}

message("Number of sample columns detected: ", length(sample_cols))
message("Number of taxa/features before duplicate collapsing: ", nrow(counts_raw))

# -----------------------------
# Convert sample columns to numeric
# -----------------------------

counts_raw <- counts_raw %>%
  mutate(across(all_of(sample_cols), as.numeric))

if (any(is.na(as.matrix(counts_raw[, sample_cols])))) {
  warning("NA values found in count columns. Replacing NA values with 0.")
  counts_raw[sample_cols] <- lapply(counts_raw[sample_cols], function(x) {
    x[is.na(x)] <- 0
    x
  })
}

if (any(as.matrix(counts_raw[, sample_cols]) < 0)) {
  stop("Negative values found in count table. Counts must be non-negative.")
}

# -----------------------------
# Collapse duplicate features, if present
# -----------------------------

annotation_df <- counts_raw %>%
  select(all_of(id_cols)) %>%
  distinct(.data[[feature_col]], .keep_all = TRUE)

counts_summed <- counts_raw %>%
  group_by(.data[[feature_col]]) %>%
  summarise(
    across(all_of(sample_cols), sum),
    .groups = "drop"
  )

message("Number of taxa/features after duplicate collapsing: ", nrow(counts_summed))

count_matrix <- counts_summed %>%
  column_to_rownames(feature_col) %>%
  as.matrix()

# -----------------------------
# Check sample totals
# -----------------------------

sample_totals <- colSums(count_matrix)

zero_total_samples <- names(sample_totals)[sample_totals == 0]

if (length(zero_total_samples) > 0) {
  stop(
    "These samples have total count of zero and cannot be normalized: ",
    paste(zero_total_samples, collapse = ", ")
  )
}

message("Minimum sample total: ", min(sample_totals))
message("Maximum sample total: ", max(sample_totals))

# -----------------------------
# Relative abundance normalization
# -----------------------------
# Formula:
# relative abundance = count / total counts in that sample

relative_matrix <- sweep(
  count_matrix,
  2,
  sample_totals,
  FUN = "/"
)

relative_df <- relative_matrix %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column(var = feature_col)

# Reattach annotation columns if present.
if (length(annotation_cols) > 0) {
  relative_df <- annotation_df %>%
    select(all_of(id_cols)) %>%
    right_join(relative_df, by = feature_col)
}

write_tsv(relative_df, output_relative)

message("Relative abundance normalization complete.")
message("Relative abundance table written to: ", output_relative)

# -----------------------------
# CLR normalization
# -----------------------------
# Formula per sample:
#
# x_pseudo = count + pseudocount
# clr(x) = log(x_pseudo) - mean(log(x_pseudo))
#
# This centers each sample around its geometric mean.

pseudo_matrix <- count_matrix + pseudocount

log_matrix <- log(pseudo_matrix)

clr_matrix <- sweep(
  log_matrix,
  2,
  colMeans(log_matrix),
  FUN = "-"
)

clr_df <- clr_matrix %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column(var = feature_col)

# Reattach annotation columns if present.
if (length(annotation_cols) > 0) {
  clr_df <- annotation_df %>%
    select(all_of(id_cols)) %>%
    right_join(clr_df, by = feature_col)
}

write_tsv(clr_df, output_clr)

message("CLR normalization complete.")
message("CLR table written to: ", output_clr)

# -----------------------------
# Final checks
# -----------------------------

relative_sample_sums <- colSums(relative_matrix)

message("Relative abundance sample sums:")
print(relative_sample_sums)

clr_sample_means <- colMeans(clr_matrix)

message("CLR sample means, should be approximately zero:")
print(clr_sample_means)

message("Normalization finished successfully.")
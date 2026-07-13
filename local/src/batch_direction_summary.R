#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

abundance_file <- args[1]
metadata_file <- args[2]
top_table_file <- args[3]
out_summary <- args[4]
out_plot <- args[5]
top_n <- as.integer(args[6])

dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

read_table_auto <- function(path) {
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

required_top <- c("GeneID", "logFC", "Pvalue", "Pvalue_adj")
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
top$Pvalue_adj <- as.numeric(top$Pvalue_adj)
top$logFC <- as.numeric(top$logFC)

top_taxa <- top %>%
  arrange(Pvalue) %>%
  head(top_n) %>%
  pull(GeneID)

top_taxa <- intersect(top_taxa, rownames(mat))

if (length(top_taxa) < 1) {
  stop("No top taxa from top_table were found in abundance matrix")
}

mat_clr <- clr_transform(mat)

long_df <- as.data.frame(mat_clr[top_taxa, , drop = FALSE]) %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "sample",
    values_to = "clr_abundance"
  ) %>%
  left_join(meta, by = "sample")

batch_fc <- long_df %>%
  group_by(GeneID, batch, condition) %>%
  summarise(mean_clr = mean(clr_abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_clr)

if (!all(c("Control", "Resistant") %in% colnames(batch_fc))) {
  stop("Both Control and Resistant must be present within the retained batches")
}

batch_fc <- batch_fc %>%
  mutate(batch_logFC = Resistant - Control) %>%
  select(GeneID, batch, batch_logFC)

wide_fc <- batch_fc %>%
  pivot_wider(names_from = batch, values_from = batch_logFC)

batch_cols <- setdiff(colnames(wide_fc), "GeneID")

summary_df <- wide_fc %>%
  left_join(
    top %>% select(GeneID, pooled_logFC = logFC, Pvalue, Pvalue_adj),
    by = "GeneID"
  )

summary_df$n_batches <- rowSums(!is.na(summary_df[, batch_cols, drop = FALSE]))
summary_df$n_positive <- rowSums(summary_df[, batch_cols, drop = FALSE] > 0, na.rm = TRUE)
summary_df$n_negative <- rowSums(summary_df[, batch_cols, drop = FALSE] < 0, na.rm = TRUE)

summary_df$direction_consistent <- with(
  summary_df,
  n_batches > 0 & (n_positive == n_batches | n_negative == n_batches)
)

summary_df <- summary_df %>%
  arrange(Pvalue)

write.table(
  summary_df,
  out_summary,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (length(batch_cols) >= 2) {
  plot_df <- summary_df

  x_batch <- batch_cols[1]
  y_batch <- batch_cols[2]

  p <- ggplot(plot_df, aes_string(x = x_batch, y = y_batch)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(aes(shape = direction_consistent), size = 2.5, alpha = 0.8) +
    labs(
      title = "Within-batch logFC consistency",
      x = paste0("CLR mean difference in ", x_batch, " (Resistant - Control)"),
      y = paste0("CLR mean difference in ", y_batch, " (Resistant - Control)")
    ) +
    theme_bw()

  ggsave(out_plot, p, width = 6, height = 6)
} else {
  p <- ggplot(summary_df, aes(pooled_logFC)) +
    geom_histogram(bins = 30) +
    labs(
      title = "Only one batch available; consistency plot not applicable",
      x = "Pooled logFC",
      y = "Number of taxa"
    ) +
    theme_bw()

  ggsave(out_plot, p, width = 6, height = 5)
}
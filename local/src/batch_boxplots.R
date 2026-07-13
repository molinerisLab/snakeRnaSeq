#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

abundance_file <- args[1]
metadata_file <- args[2]
top_table_file <- args[3]
out_boxplots <- args[4]
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

if (length(top_taxa) < 1) {
  stop("No top taxa from top_table were found in abundance matrix")
}

mat_clr <- clr_transform(mat)

plot_df <- as.data.frame(mat_clr[top_taxa, , drop = FALSE]) %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(
    cols = -GeneID,
    names_to = "sample",
    values_to = "clr_abundance"
  ) %>%
  left_join(meta, by = "sample")

plot_df$GeneID <- gsub("_\\d+_S$", "", plot_df$GeneID)
plot_df$GeneID <- gsub("_", " ", plot_df$GeneID)
# plot_df$GeneID <- sub("^([^ ]+).*$", "\\1", plot_df$GeneID) # Remove species-level resolution entirely and keeps only the genus (or first token)
plot_df$condition <- factor(plot_df$condition, levels = c("Control", "Resistant"))

p <- ggplot(plot_df, aes(condition, clr_abundance, fill = condition)) +
  geom_boxplot(width = 0.65, alpha = 0.75, outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = condition), width = 0.15, alpha = 0.8, size = 2) +
  facet_grid(GeneID ~ batch, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#56B4E9", "Resistant" = "#D55E00")) +
  scale_color_manual(values = c("Control" = "#0072B2", "Resistant" = "#A63603")) +
  labs(
    title = "Top taxa abundance by condition within batch",
    x = NULL,
    y = "CLR abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text.y = element_text(size = 7, face = "italic"),
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

ggsave(out_boxplots, p, width = 12, height = max(6, top_n * 1.2))

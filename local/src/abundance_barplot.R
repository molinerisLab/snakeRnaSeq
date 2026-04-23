#!/usr/bin/env Rscript

library(tidyverse)
library(pheatmap)
library(vegan)
library(RColorBrewer)

# -----------------------------
# Input arguments from Snakemake
# -----------------------------
rel_file <- snakemake@input[["relative"]]
meta_file <- snakemake@input[["metadata"]]

outdir <- snakemake@output[["outdir"]]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load metadata
# -----------------------------
meta <- read_tsv(meta_file) %>%
  mutate(
    condition = factor(condition),
    batch = factor(batch)
  )

# -----------------------------
# Function to load matrix
# -----------------------------
load_matrix <- function(file) {
  df <- read_tsv(file)
  
  df <- df %>%
    select(name, ends_with(".bracken_num"))
  
  colnames(df) <- colnames(df) %>%
    str_replace("\\.bracken_num$", "")
  
  mat <- df %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  return(mat)
}

# -----------------------------
# Load data
# -----------------------------
rel_mat <- load_matrix(rel_file)

# Ensure sample order consistency
meta <- meta %>%
  filter(sample %in% colnames(rel_mat)) %>%
  arrange(match(sample, colnames(rel_mat)))

# -----------------------------
# Stacked Bar Plot (Relative)
# -----------------------------
rel_long <- as.data.frame(rel_mat) %>%
  rownames_to_column("taxon") %>%
  pivot_longer(-taxon, names_to = "sample", values_to = "abundance") %>%
  left_join(meta, by = "sample")

top_taxa <- rel_long %>%
  group_by(taxon) %>%
  summarise(total = sum(abundance)) %>%
  slice_max(total, n = 20) %>%
  pull(taxon)

rel_long <- rel_long %>%
  mutate(taxon = ifelse(taxon %in% top_taxa, taxon, "Other"))

p_bar <- ggplot(rel_long, aes(sample, abundance, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_wrap(~condition, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  legend.text = element_text(face = "italic")) +
  labs(title = "Relative Abundance (Top Taxa)", fill = "Species", x = "Sample", y = "Relative Abundance")

ggsave(file.path(outdir, "stacked_barplot.pdf"), p_bar, width = 12, height = 6)

# -----------------------------
# Done
# -----------------------------
cat("Plots generated in:", outdir, "\n")
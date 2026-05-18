#!/usr/bin/env Rscript

library(tidyverse)
library(pheatmap)
library(vegan)
library(RColorBrewer)


# Input arguments from Snakemake
rel_file <- snakemake@input[["relative"]]
meta_file <- snakemake@input[["metadata"]]
out_bar <- snakemake@output[["barplot"]]


# Load metadata
meta <- read_tsv(meta_file) %>%
  mutate(
    group=factor(group),
    condition = factor(condition),
    batch = factor(batch)
  )

# Function to load matrix
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

# Load data
rel_mat <- load_matrix(rel_file)

# Ensure sample order consistency
meta <- meta %>%
  filter(sample %in% colnames(rel_mat)) %>%
  arrange(match(sample, colnames(rel_mat)))


# Stacked Bar Plot (Relative)
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
  facet_wrap(~group, scales = "free_x") + # You can adjust this based on your metadata. For example, you might want to facet by condition or batch instead.
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
  legend.text = element_text(face = "italic")) +
  labs(title = "Relative Abundance (Top Taxa)", fill = "Species", x = "Sample", y = "Relative Abundance")

ggsave(out_bar, p_bar, width = 12, height = 6)

cat("Plots generated in:", out_bar, "\n")
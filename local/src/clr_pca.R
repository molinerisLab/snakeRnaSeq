#!/usr/bin/env Rscript

library(tidyverse)
library(pheatmap)
library(vegan)
library(RColorBrewer)


# Input arguments from Snakemake

clr_file <- snakemake@input[["clr"]]
meta_file <- snakemake@input[["metadata"]]
out_pca <- snakemake@output[["clr_pca"]]


# Load metadata
meta <- read_tsv(meta_file) %>%
  mutate(
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
clr_mat <- load_matrix(clr_file)

# Align metadata with matrix columns
meta <- meta %>%
  filter(sample %in% colnames(clr_mat)) %>%
  arrange(match(sample, colnames(clr_mat)))


# PCA
pca <- prcomp(t(clr_mat), scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  left_join(meta, by = "sample")

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "PCA (CLR-transformed)")

ggsave(out_pca, p_pca, width = 7, height = 5)
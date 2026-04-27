#!/usr/bin/env Rscript

library(tidyverse)
library(pheatmap)
library(vegan)
library(RColorBrewer)


# Input arguments from Snakemake
rel_file <- snakemake@input[["relative"]]
meta_file <- snakemake@input[["metadata"]]
out_heatmap <- snakemake@output[["rel_heatmap"]]


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
rel_mat <- load_matrix(rel_file)

# Ensure sample order consistency
meta <- meta %>%
  filter(sample %in% colnames(rel_mat)) %>%
  arrange(match(sample, colnames(rel_mat)))


# Relative Heatmap
top_taxa_rel <- names(sort(rowSums(rel_mat), decreasing = TRUE))[1:30]
rel_sub <- rel_mat[top_taxa_rel, ]

annotation <- meta %>%
  column_to_rownames("sample")

pheatmap(rel_sub,
         scale = "row",
         annotation_col = annotation,
         show_rownames = FALSE,
         filename = out_heatmap,
         width = 8, height = 10)
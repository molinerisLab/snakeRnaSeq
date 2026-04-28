#!/usr/bin/env Rscript

library(tidyverse)
library(pheatmap)
library(vegan)
library(RColorBrewer)


# Input arguments from Snakemake

clr_file <- snakemake@input[["clr"]]
meta_file <- snakemake@input[["metadata"]]
out_heatmap <- snakemake@output[["clr_heatmap"]]


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

# Create annotation
annotation <- meta %>%
  column_to_rownames("sample")

# CLR Heatmap
pheatmap(clr_mat,
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = annotation,
         show_rownames = FALSE,
         filename = out_heatmap,
         width = 8, height = 10)
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

abundance_file <- args[1]
metadata_file  <- args[2]
top_table_file <- args[3]
out_abundance  <- args[4]
out_metadata   <- args[5]
out_top_table  <- args[6]

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

abund <- read_table_auto(abundance_file)
meta  <- read_table_auto(metadata_file)

if (!"sample" %in% colnames(meta)) {
  stop("metadata.txt must contain a 'sample' column")
}

if (!"batch" %in% colnames(meta)) {
  stop("metadata.txt must contain a 'batch' column")
}

if (!"condition" %in% colnames(meta)) {
  stop("metadata.txt must contain a 'condition' column")
}

# Exclude KIS/Kisumu if present. This is safe even when KIS is absent.
meta_no_kis <- meta %>%
  filter(!grepl("KIS|Kisumu", batch, ignore.case = TRUE),
         !grepl("^KIS", sample, ignore.case = TRUE))

sample_cols <- intersect(meta_no_kis$sample, colnames(abund))

if (length(sample_cols) == 0) {
  stop("No abundance sample columns matched metadata sample names")
}

taxon_col <- colnames(abund)[1]

abund_no_kis <- abund %>%
  select(all_of(taxon_col), all_of(sample_cols))

# Normalize first column name for downstream scripts
colnames(abund_no_kis)[1] <- "GeneID"

if (grepl("\\.xlsx$", top_table_file, ignore.case = TRUE)) {
  top_table <- as.data.frame(read_excel(top_table_file))
} else {
  top_table <- read_table_auto(top_table_file)
}

if ("Geneid" %in% colnames(top_table) && !"GeneID" %in% colnames(top_table)) {
  colnames(top_table)[colnames(top_table) == "Geneid"] <- "GeneID"
}

write.table(
  abund_no_kis,
  out_abundance,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  meta_no_kis,
  out_metadata,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  top_table,
  out_top_table,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
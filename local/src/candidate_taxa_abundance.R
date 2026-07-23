#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop(
    "Usage: Rscript plot_candidate_taxa_abundance.R ",
    "<diff_table.xlsx/tsv> <relative_abundance.tsv> <metadata.tsv> ",
    "<out_plot.pdf> <out_selected_taxa.tsv> <out_long_abundance.tsv> ",
    "<top_n> <padj_cutoff> [group_col]"
  )
}

diff_file <- args[1]
abundance_file <- args[2]
metadata_file <- args[3]
out_pdf <- args[4]
out_selected <- args[5]
out_long <- args[6]
top_n <- as.integer(args[7])
padj_cutoff <- as.numeric(args[8])
group_col <- ifelse(length(args) >= 9, args[9], "condition")

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_selected), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_long), recursive = TRUE, showWarnings = FALSE)

message("Differential table: ", diff_file)
message("Abundance table: ", abundance_file)
message("Metadata table: ", metadata_file)
message("Grouping column: ", group_col)

# -----------------------------
# Helper functions
# -----------------------------

read_table_auto <- function(file) {
  ext <- tolower(tools::file_ext(file))

  if (ext %in% c("xlsx", "xls")) {
    message("Reading Excel file: ", file)
    read_excel(file, sheet = 1)
  } else {
    message("Reading delimited text file: ", file)
    read_tsv(file, show_col_types = FALSE)
  }
}

make_taxon_label <- function(x) {
  x %>%
    str_remove("_S$") %>%
    str_remove("_[0-9]+$") %>%
    str_replace_all("_", " ")
}

# -----------------------------
# Read differential abundance table
# -----------------------------

diff <- read_table_auto(diff_file)

message("Columns in differential table:")
message(paste(names(diff), collapse = ", "))

required_diff_cols <- c("GeneID", "logFC", "Pvalue", "Pvalue_adj")
missing_diff_cols <- setdiff(required_diff_cols, names(diff))

if (length(missing_diff_cols) > 0) {
  stop(
    "Differential abundance table is missing required columns: ",
    paste(missing_diff_cols, collapse = ", "),
    "\nDetected columns were: ",
    paste(names(diff), collapse = ", ")
  )
}

if (!"significance" %in% names(diff)) {
  diff$significance <- as.numeric(diff$Pvalue_adj < padj_cutoff)
}

diff <- diff %>%
  mutate(
    GeneID = as.character(GeneID),
    logFC = as.numeric(logFC),
    Pvalue = as.numeric(Pvalue),
    Pvalue_adj = as.numeric(Pvalue_adj),
    significance = suppressWarnings(as.numeric(significance)),
    abs_logFC = abs(logFC)
  ) %>%
  filter(!is.na(GeneID))


# -----------------------------
# Select candidate taxa
# -----------------------------
# Use FDR-adjusted results only when at least top_n taxa are significant.
# Otherwise, fall back to the top_n taxa ranked by raw P value.

sig_taxa <- diff %>%
  filter(
    (!is.na(Pvalue_adj) & Pvalue_adj < padj_cutoff) |
      (!is.na(significance) & significance != 0)
  ) %>%
  distinct(GeneID, .keep_all = TRUE)

n_sig <- nrow(sig_taxa)

message(
  "Number of significant taxa at cutoff ",
  padj_cutoff,
  ": ",
  n_sig
)

if (n_sig >= top_n) {
  selected <- sig_taxa %>%
    arrange(Pvalue_adj, Pvalue, desc(abs_logFC)) %>%
    slice_head(n = top_n) %>%
    mutate(
      selection_rule = paste0(
        "top_",
        top_n,
        "_FDR_significant_by_adjusted_Pvalue_less_than_",
        padj_cutoff
      )
    )
} else {
  selected <- diff %>%
    arrange(Pvalue, desc(abs_logFC)) %>%
    distinct(GeneID, .keep_all = TRUE) %>%
    slice_head(n = top_n) %>%
    mutate(
      selection_rule = paste0(
        "top_",
        top_n,
        "_nominal_by_raw_Pvalue_due_to_only_",
        n_sig,
        "_FDR_significant_taxa"
      )
    )
}

# -----------------------------
# Read relative abundance table
# -----------------------------

abundance <- read_tsv(abundance_file, show_col_types = FALSE)

message("Columns in abundance table:")
message(paste(names(abundance), collapse = ", "))

if ("Geneid" %in% names(abundance)) {
  abundance_gene_col <- "Geneid"
} else if ("GeneID" %in% names(abundance)) {
  abundance_gene_col <- "GeneID"
} else {
  abundance_gene_col <- names(abundance)[1]
  warning("No Geneid/GeneID column found in abundance table. Using first column: ", abundance_gene_col)
}

abundance <- abundance %>%
  rename(GeneID = all_of(abundance_gene_col)) %>%
  mutate(GeneID = as.character(GeneID))

sample_cols <- setdiff(names(abundance), "GeneID")

if (length(sample_cols) == 0) {
  stop("No sample columns found in abundance table.")
}

abundance <- abundance %>%
  mutate(across(all_of(sample_cols), as.numeric))

# Keep only selected taxa
abundance_selected <- abundance %>%
  semi_join(selected, by = "GeneID")

missing_taxa <- setdiff(selected$GeneID, abundance_selected$GeneID)

if (length(missing_taxa) > 0) {
  warning(
    "These selected taxa were not found in the abundance table: ",
    paste(missing_taxa, collapse = ", ")
  )
}

if (nrow(abundance_selected) == 0) {
  stop("None of the selected taxa were found in the abundance table.")
}

# -----------------------------
# Convert abundance table to long format
# -----------------------------

abundance_long <- abundance_selected %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  mutate(
    sample = as.character(sample),
    abundance = as.numeric(abundance)
  )

# -----------------------------
# Read metadata
# -----------------------------

metadata <- read_tsv(metadata_file, show_col_types = FALSE)

message("Columns in metadata:")
message(paste(names(metadata), collapse = ", "))

if (!"sample" %in% names(metadata)) {
  stop("Metadata file must contain a column named 'sample'.")
}

if (!group_col %in% names(metadata)) {
  stop(
    "Metadata file does not contain requested grouping column: ",
    group_col,
    "\nAvailable metadata columns are: ",
    paste(names(metadata), collapse = ", ")
  )
}

metadata <- metadata %>%
  transmute(
    sample = as.character(sample),
    group = as.character(.data[[group_col]])
  )

# Drop samples that are not in metadata
# This prevents failure if the abundance table has extra samples.
metadata_samples <- metadata$sample

extra_abundance_samples <- setdiff(unique(abundance_long$sample), metadata_samples)

if (length(extra_abundance_samples) > 0) {
  warning(
    "These abundance samples are not in metadata and will be dropped: ",
    paste(extra_abundance_samples, collapse = ", ")
  )
}

abundance_long <- abundance_long %>%
  filter(sample %in% metadata_samples)

missing_abundance_samples <- setdiff(metadata_samples, unique(abundance_long$sample))

if (length(missing_abundance_samples) > 0) {
  warning(
    "These metadata samples were not found in the abundance table: ",
    paste(missing_abundance_samples, collapse = ", ")
  )
}

# -----------------------------
# Join abundance, metadata, and differential results
# -----------------------------

selected_for_join <- selected %>%
  mutate(
    taxon_label = make_taxon_label(GeneID),
    plot_label = paste0(
      taxon_label,
      "\nlogFC=", round(logFC, 2),
      "; P=", signif(Pvalue, 2),
      "; FDR=", signif(Pvalue_adj, 2)
    )
  ) %>%
  select(
    GeneID,
    taxon_label,
    plot_label,
    logFC,
    Pvalue,
    Pvalue_adj,
    significance,
    selection_rule
  )

plot_df <- abundance_long %>%
  left_join(metadata, by = "sample") %>%
  left_join(selected_for_join, by = "GeneID")

if (any(is.na(plot_df$group))) {
  bad_samples <- plot_df %>%
    filter(is.na(group)) %>%
    distinct(sample) %>%
    pull(sample)

  stop(
    "Some samples have no group assignment: ",
    paste(bad_samples, collapse = ", ")
  )
}

# Put Control before Resistant if both exist
if (all(c("Control", "Resistant") %in% unique(plot_df$group))) {
  plot_df <- plot_df %>%
    mutate(group = factor(group, levels = c("Control", "Resistant")))
} else {
  plot_df <- plot_df %>%
    mutate(group = factor(group))
}

plot_df <- plot_df %>%
  mutate(
    plot_label = factor(plot_label, levels = selected_for_join$plot_label)
  )

# -----------------------------
# Write selected taxa and long table
# -----------------------------

write_tsv(
  selected_for_join %>%
    select(
      GeneID,
      taxon_label,
      logFC,
      Pvalue,
      Pvalue_adj,
      significance,
      selection_rule
    ),
  out_selected
)

write_tsv(
  plot_df %>%
    select(
      GeneID,
      taxon_label,
      sample,
      group,
      abundance,
      logFC,
      Pvalue,
      Pvalue_adj,
      significance,
      selection_rule
    ),
  out_long
)

# -----------------------------
# Plot
# -----------------------------

nonzero_abundances <- plot_df$abundance[plot_df$abundance > 0]

if (length(nonzero_abundances) > 0) {
  pseudocount <- min(nonzero_abundances, na.rm = TRUE) / 2
} else {
  pseudocount <- 1e-12
}

plot_df <- plot_df %>%
  mutate(abundance_for_plot = abundance + pseudocount)

plot_title <- paste0(
  "Relative abundance of candidate species\n",
  "Candidate selection: ",
  unique(selected$selection_rule)[1]
)

p <- ggplot(plot_df, aes(x = group, y = abundance_for_plot)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    alpha = 0.35
  ) +
  geom_point(
    aes(shape = group),
    position = position_jitter(width = 0.12, height = 0),
    size = 2,
    alpha = 0.85
  ) +
  facet_wrap(~ plot_label, scales = "free_y") +
  scale_y_log10(labels = label_scientific()) +
  labs(
    title = plot_title,
    x = NULL,
    y = paste0(
      "Relative abundance + pseudocount",
      "\nlog10 scale; pseudocount = ",
      signif(pseudocount, 3)
    ),
    shape = group_col
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 8)
  )

plot_height <- max(6, ceiling(length(unique(plot_df$GeneID)) / 3) * 3.2)

ggsave(
  filename = out_pdf,
  plot = p,
  width = 12,
  height = plot_height,
  limitsize = FALSE
)

message("Saved plot to: ", out_pdf)
message("Saved selected taxa to: ", out_selected)
message("Saved long-format abundance table to: ", out_long)
message("Done.")
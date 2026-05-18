#!/usr/bin/env Rscript

# Plot order-level taxa from a folder of *.barplot.txt / *.txt files
#
# Input format expected per file:
#   o__Flavobacteriales    217537.0
#   o__Acetobacterales      69887.0
#
# Usage:
#   Rscript plot_order_barplots.R /path/to/barplot_folder output_order_barplot.pdf
#
# Optional:
#   Rscript plot_order_barplots.R /path/to/barplot_folder output_order_barplot.pdf 20
#
# The optional third argument is the number of top orders to show before grouping
# the rest as "Other". Default = 20.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript plot_order_barplots.R <barplot_folder> <output_pdf> [top_n]\n",
    "Example: Rscript plot_order_barplots.R ./barplots order_barplot.pdf 20"
  )
}

input_dir <- args[1]
output_pdf <- args[2]
top_n <- ifelse(length(args) >= 3, as.integer(args[3]), 20)

if (!dir.exists(input_dir)) {
  stop("Input folder does not exist: ", input_dir)
}

files <- list.files(
  input_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No .txt files found in folder: ", input_dir)
}

read_barplot_file <- function(file) {
  sample_name <- basename(file)
  sample_name <- str_remove(sample_name, "\\.txt$")
  sample_name <- str_remove(sample_name, "\\.barplot$")

  df <- read.table(
    file,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    comment.char = "#",
    fill = TRUE
  )

  if (ncol(df) < 2) {
    warning("Skipping malformed file: ", file)
    return(data.frame(sample = character(), order = character(), count = numeric()))
  }

  df <- df[, 1:2]
  colnames(df) <- c("taxon", "count")

  df$count <- suppressWarnings(as.numeric(df$count))

  df %>%
    filter(!is.na(taxon), !is.na(count)) %>%
    filter(str_starts(taxon, "o__")) %>%
    mutate(
      sample = sample_name,
      order = str_remove(taxon, "^o__")
    ) %>%
    select(sample, order, count)
}

order_counts <- bind_rows(lapply(files, read_barplot_file))

if (nrow(order_counts) == 0) {
  stop("No order-level taxa found. Expected taxa beginning with 'o__'.")
}

# Sum duplicated order entries within each sample, if any
order_counts <- order_counts %>%
  group_by(sample, order) %>%
  summarise(count = sum(count), .groups = "drop")

# Convert to relative abundance per sample
order_rel <- order_counts %>%
  group_by(sample) %>%
  mutate(
    total_count = sum(count),
    relative_abundance = count / total_count
  ) %>%
  ungroup()

# Keep top N orders overall; group the rest as Other
top_orders <- order_rel %>%
  group_by(order) %>%
  summarise(total_abundance = sum(relative_abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = top_n) %>%
  pull(order)

plot_df <- order_rel %>%
  mutate(order_plot = ifelse(order %in% top_orders, order, "Other")) %>%
  group_by(sample, order_plot) %>%
  summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
  mutate(percent_abundance = relative_abundance * 100)

# Put Other at the bottom/top consistently
order_levels <- plot_df %>%
  group_by(order_plot) %>%
  summarise(total = sum(percent_abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(order_plot)

order_levels <- c(setdiff(order_levels, "Other"), "Other")

plot_df$order_plot <- factor(plot_df$order_plot, levels = rev(order_levels))

p <- ggplot(plot_df, aes(x = sample, y = percent_abundance, fill = order_plot)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(
    x = "Sample",
    y = "Relative abundance (%)",
    fill = "Order",
    title = "Order-level taxonomic composition"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = output_pdf,
  plot = p,
  width = max(8, length(unique(plot_df$sample)) * 0.6),
  height = 6
)

message("Saved order-level barplot to: ", output_pdf)
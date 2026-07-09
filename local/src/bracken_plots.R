#!/usr/bin/env Rscript
# ---------- 1. Data Preparation ----------
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggrepel)
library(scales)
library(ggExtra)

# ---------- CLI Arguments ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript bracken_plots.R <bracken_merged_file> [species_name] [top_n] [output_dir]")
}

input_file    <- args[1]
species_focus <- ifelse(length(args) >= 2, args[2], "Staphylococcus aureus")
top_n_species <- ifelse(length(args) >= 3, as.integer(args[3]), 30)
output_dir    <- ifelse(length(args) >= 4, args[4], ".")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read the merged Bracken output
# check.names = FALSE prevents R from adding 'X' in front of sample names that start with numbers
raw_df <- read.delim(input_file, header = TRUE, sep = "\t", check.names = FALSE)

# Auto-detect ID columns (handles both taxonomy_lv and taxonomy_lvl)
id_cols <- intersect(names(raw_df), c("taxonomy_id", "taxonomy_lv", "taxonomy_lvl"))

# If file has both _num and _frac columns, keep only _num
has_num_frac <- any(grepl("_num$", names(raw_df))) && any(grepl("_frac$", names(raw_df)))
if (has_num_frac) {
  sample_cols <- grep("_num$", names(raw_df), value = TRUE)
  raw_df <- raw_df[, c("name", id_cols, sample_cols)]
  # Clean sample names: remove the _num suffix
  names(raw_df) <- gsub("\\.braken_num$", "", names(raw_df))
}

# Pivot the wide matrix into a long format
df <- raw_df %>%
  select(-all_of(id_cols)) %>%
  # Rename 'name' to 'species' to match our plotting variables
  rename(species = name) %>%
  # Pivot all columns EXCEPT 'species' into our long format
  pivot_longer(
    cols = -species, 
    names_to = "sample",
    values_to = "abundance"
  )


# Now generate the species_summary exactly as before...
species_summary <- df %>%
  group_by(species) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    prevalence = sum(abundance > 0) / n_distinct(df$sample) * 100
  ) %>%
  arrange(desc(mean_abundance))

# Use top_n_species from CLI args (default: 30)
top_species_list <- head(species_summary$species, top_n_species)


# ---------- Plot 1: Clustered Heatmap ----------
# Filter for top species and reshape to wide format (samples as columns, species as rows)
df_wide <- df %>%
  filter(species %in% top_species_list) %>%
  select(sample, species, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0) %>%
  as.data.frame()

# Set row names to species and remove the species column
rownames(df_wide) <- df_wide$species
df_wide <- df_wide[, -1]

# Log-transform the data (log10(x + 1e-5)) to handle zeros and massive abundance skews
log_matrix <- log10(as.matrix(df_wide) + 1e-5)

# Generate Heatmap (saves automatically if filename is provided)
pheatmap(log_matrix, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         show_colnames = FALSE, # Hide 270 sample names for clean UI
         fontsize_row = 10,
         main = paste0("Log10 Abundance of Top ", top_n_species, " Species Across All Samples"),
         filename = file.path(output_dir, paste0("heatmap_top", top_n_species, "_species.png")),
         width = 12, height = 8)


sample_data_cols <- setdiff(names(raw_df), c("name", id_cols))
total_reads <- colSums(raw_df[, sample_data_cols, drop = FALSE], na.rm = TRUE)

species_df <- raw_df %>%
  filter(name == species_focus) %>%
  select(-all_of(id_cols)) %>%
  pivot_longer(
    cols = -name,
    names_to = "sample",
    values_to = "raw_counts"
  ) %>%
  mutate(
    total_sample_reads = total_reads[sample],
    fraction = raw_counts / total_sample_reads
  ) %>%
  mutate(sample = gsub("_001$", "", sample))

# ---------- 2. The Base Scatter Plot ----------
p_base <- ggplot(species_df, aes(x = raw_counts + 1, y = fraction)) +
  # Using alpha=0.6 means overlapping points will appear darker red automatically
  geom_point(alpha = 0.6, color = "firebrick3", size = 3) +
  
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
  
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000, 300000),
    labels = c("0", "10", "100", "1,000", "10,000", "100,000", "300,000")
  ) +
  scale_y_continuous(labels = scales::percent) + 
  
  labs(
    title = paste(species_focus, "Profiling Across Cohort"),
    subtitle = "Marginal histograms highlight cohort density at zero vs. extreme outliers",
    x = paste0("Total ", species_focus, " Reads (Log10 Scale)"),
    y = "Relative Abundance within Unmapped Fraction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank()
  )

# ---------- 3. Adding the Marginals ----------
# This wraps the base plot with histograms on the top and right
p_marginal <- ggMarginal(p_base, type = "histogram", fill = "steelblue", color = "white", bins = 40)

# Save the marginal plot
species_tag <- gsub(" ", "_", tolower(species_focus))
ggsave(file.path(output_dir, paste0(species_tag, "_marginal_scatter.png")),
       p_marginal, width = 12, height = 7, dpi = 300)
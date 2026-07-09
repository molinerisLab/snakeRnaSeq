#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
# ---------- CLI Args ----------
args <- commandArgs(trailingOnly = TRUE)
unmapped_file <- ifelse(length(args) >= 1, args[1], "plots/unmapped_depth.csv")
total_file    <- ifelse(length(args) >= 2, args[2], "plots/sequencing_depth.csv")
output_dir    <- ifelse(length(args) >= 3, args[3], "plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- Load ----------
clean_sample <- function(x) {
  gsub("_001(_unmapped)?_R1\\.fastq$", "", x)
}

df_unmapped <- read.csv(unmapped_file, header = TRUE)
df_unmapped$name <- clean_sample(df_unmapped$name)

# Rank by unmapped depth
df <- df_unmapped %>%
  arrange(unmapped_depth)

df$name <- factor(df$name, levels = df$name)

# Thresholds from UNMAPPED depth
median_depth <- median(df$unmapped_depth)
low_threshold <- median_depth * 0.75
high_threshold <- median_depth * 1.25

# Color groups for palette/legend
df <- df %>%
  mutate(group = case_when(
    unmapped_depth < low_threshold  ~ "Low",
    unmapped_depth > high_threshold ~ "High",
    TRUE                            ~ "Median range"
  ))

df$group <- factor(df$group, levels = c("Low", "Median range", "High"))

palette_cov <- c(
  "Low" = "#d95f02",
  "Median range" = "#bdbdbd",
  "High" = "#1b9e77"
)

# ---------- Plot 1: all samples, unmapped only ----------
p_unmapped_all <- ggplot(df, aes(x = name, y = unmapped_depth, fill = group)) +
  geom_col() +
  geom_hline(yintercept = median_depth, linetype = "dashed", color = "red") +
  geom_hline(yintercept = low_threshold, linetype = "dotted", color = "orange") +
  geom_hline(yintercept = high_threshold, linetype = "dotted", color = "orange") +
  annotate(
    "text",
    x = nrow(df) * 0.02,
    y = median_depth,
    label = paste0("Median Unmapped: ", format(round(median_depth), big.mark = ",")),
    vjust = -0.8, hjust = 0, color = "red", fontface = "bold"
  ) +
  scale_fill_manual(values = palette_cov) +
  labs(
    title = "Unmapped Reads Across All Samples",
    subtitle = "Samples ranked by unmapped depth",
    x = "Samples ranked by unmapped reads",
    y = "Number of unmapped reads",
    fill = "Coverage group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "unmapped_reads_overview.jpg"), p_unmapped_all, width = 14, height = 6, dpi = 300)

# ---------- Plot 2: outliers only, 16:9 ----------
outliers <- df %>%
  filter(group != "Median range") %>%
  arrange(unmapped_depth)

outliers$name <- factor(outliers$name, levels = outliers$name)

p_unmapped_out <- ggplot(outliers, aes(x = name, y = unmapped_depth, fill = group)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  geom_col(width = 0.85) +
  geom_hline(yintercept = median_depth, linetype = "dashed", color = "red") +
  geom_hline(yintercept = low_threshold, linetype = "dotted", color = "orange") +
  geom_hline(yintercept = high_threshold, linetype = "dotted", color = "orange") +
  annotate(
    "text",
    x = 1,
    y = median_depth,
    label = paste0("Median Unmapped: ", format(round(median_depth), big.mark = ",")),
    vjust = -0.5,
    hjust = 0,
    color = "red",
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(
values = c("Low" = "#d95f02", "High" = "#1b9e77"),
breaks = c("Low", "High"),
drop = TRUE
) +


  labs(
    title = "Outlier Samples by Unmapped Reads",
    subtitle = "Samples outside +/-25% of median unmapped depth",
    x = "Sample",
    y = "Number of unmapped reads",
    fill = "Coverage group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "unmapped_reads_outliers_16x9.jpg"), p_unmapped_out, width = 22, height = 10, dpi = 300)

# ---------- Plot 3: Percentage Mapped vs Unmapped ----------
# Load the total sequencing depth for the baseline 100%
df_total <- read.csv(total_file, header = TRUE)
df_total$name <- clean_sample(df_total$name)
df_total <- df_total %>% rename(total_depth = sequencing_depth)

# Merge datasets
df_merged <- merge(df, df_total, by = "name", all.x = TRUE)

# Calculate percentages for mapped and unmapped reads
df_merged <- df_merged %>%
  mutate(
    mapped_depth = total_depth - unmapped_depth,
    mapped_pct = (mapped_depth / total_depth) * 100,
    unmapped_pct = (unmapped_depth / total_depth) * 100
  )

# Reshape the data for a stacked barplot using tidyr::pivot_longer
df_stacked <- df_merged %>%
  select(name, mapped_pct, unmapped_pct) %>%
  pivot_longer(
    cols = c(mapped_pct, unmapped_pct),
    names_to = "read_status",
    values_to = "percentage"
  ) %>%
  mutate(
    # Set factors to control stacking order (Unmapped on top, Mapped on bottom)
    read_status = factor(
      read_status, 
      levels = c("unmapped_pct", "mapped_pct"), 
      labels = c("Unmapped", "Mapped")
    )
  )

# Keep the same sample ordering established in earlier plots
df_stacked$name <- factor(df_stacked$name, levels = levels(df$name))

# Create the percentage stacked bar plot (zoomed to 90-100%)
p_pct <- ggplot(df_stacked, aes(x = name, y = percentage, fill = read_status)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  coord_cartesian(ylim = c(80, 100)) + # Zoom in on the 90-100% range
  labs(
    title = "Mapped vs Unmapped Reads (Zoomed >90%)",
    subtitle = "100% represents the initial sequencing depth",
    x = "Samples",
    y = "Percentage of reads (%)",
    fill = "Read Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# Save the plot
ggsave(file.path(output_dir, "mapped_unmapped_percentage.jpg"), p_pct, width = 14, height = 6, dpi = 300)

# ---------- Plot 4: Top 25% by unmapped percentage ----------
q75_unmapped <- quantile(df_merged$unmapped_pct, 0.75, na.rm = TRUE)

top25_unmapped <- df_merged %>%
  filter(unmapped_pct >= q75_unmapped) %>%
  arrange(unmapped_pct) %>%
  mutate(
    name = factor(name, levels = name),
    read_status = "Unmapped"
  )

# Use the same mapped/unmapped palette as Plot 3
map_unmap_colors <- c("Unmapped" = "#F8766D", "Mapped" = "#00BFC4")

p_top25_unmapped <- ggplot(top25_unmapped, aes(x = name, y = unmapped_pct, fill = read_status)) +
  geom_col() +
  scale_fill_manual(values = map_unmap_colors, breaks = c("Unmapped", "Mapped"), drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Top 25% Samples by Unmapped Percentage",
    subtitle = paste0("Samples with unmapped % >= 75th percentile (", round(q75_unmapped, 2), "%)"),
    x = "Sample",
    y = "Unmapped reads (%)",
    fill = "Read Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "top25_unmapped_percentage.jpg"), p_top25_unmapped, width = 14, height = 7, dpi = 300)
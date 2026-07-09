#!/usr/bin/env Rscript
library(ggplot2)

# ---------- CLI Args ----------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[1], "plots/sequencing_depth.csv")
output_dir <- ifelse(length(args) >= 2, args[2], "plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- Load ----------
df <- read.csv(input_file, header = TRUE)
df$name <- gsub("_001_R1.*", "", df$name)
n_samples <- nrow(df)

# ---------- Rank all samples ----------
df <- df[order(df$sequencing_depth), ]
df$name <- factor(df$name, levels = df$name)

median_depth <- median(df$sequencing_depth)
low_threshold <- median_depth * 0.75
high_threshold <- median_depth * 1.25

# ---------- Plot 1: full distribution (clean) ----------
p_all <- ggplot(df, aes(x = name, y = sequencing_depth)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = median_depth, linetype = "dashed", color = "red") +
  geom_hline(yintercept = low_threshold, linetype = "dotted", color = "orange") +
  geom_hline(yintercept = high_threshold, linetype = "dotted", color = "orange") +
  
  # Add the median annotation back, positioned on the left side
  annotate("text", x = nrow(df) * 0.02, y = median_depth, 
           label = paste("Median:", round(median_depth)), 
           vjust = -0.8, hjust = 0, color = "red", fontface = "bold") +
           
  labs(
    title = "Sequencing Depth Across All Samples",
    subtitle = "Names hidden for readability. See outlier plot for sample IDs.",
    x = "Samples ranked by depth",
    y = "Number of reads"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "sequencing_depth_overview.jpg"), p_all, width = 14, height = 6, dpi = 300)

# ---------- Plot 2: only outliers with readable names ----------
outliers <- df[df$sequencing_depth < low_threshold | df$sequencing_depth > high_threshold, ]
outliers$group <- ifelse(outliers$sequencing_depth < low_threshold, "Low", "High")

# sort for clean reading
outliers <- outliers[order(outliers$sequencing_depth), ]
outliers$name <- factor(outliers$name, levels = outliers$name)

p_out <- ggplot(outliers, aes(x = name, y = sequencing_depth, fill = group)) +
  geom_col() +
  coord_flip() +
  
  # Add the EXACT same reference lines to keep consistency across plots
  geom_hline(yintercept = median_depth, linetype = "dashed", color = "red") +
  geom_hline(yintercept = low_threshold, linetype = "dotted", color = "orange") +
  geom_hline(yintercept = high_threshold, linetype = "dotted", color = "orange") +
  
  # Add the median annotation (placed at the very top of the categorical axis)
  annotate("text", x = nrow(outliers), y = median_depth, 
           label = paste("Median:", round(median_depth)), 
           vjust = -0.5, hjust = -0.1, color = "red", fontface = "bold", size = 4) +
           
  scale_fill_manual(values = c("Low" = "#d95f02", "High" = "#1b9e77")) +
  labs(
    title = "Outlier Samples (Outside +/-25% of Median)",
    subtitle = paste0("Median = ", round(median_depth), " reads"),
    x = "Sample",
    y = "Number of reads"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "top"
  )

# height scales with number of outliers so labels stay readable
h <- max(6, 0.28 * nrow(outliers))
ggsave(file.path(output_dir, "sequencing_depth_outliers_readable.jpg"), p_out, width = 11, height = h, dpi = 300)

# ---------- Plot 2: outliers on X axis (16:9 friendly) ----------
outliers <- df[df$sequencing_depth < low_threshold | df$sequencing_depth > high_threshold, ]
outliers$group <- ifelse(outliers$sequencing_depth < low_threshold, "Low", "High")
outliers <- outliers[order(outliers$sequencing_depth), ]
outliers$name <- factor(outliers$name, levels = outliers$name)

p_out <- ggplot(outliers, aes(x = name, y = sequencing_depth, fill = group)) +
  geom_col(width = 0.85) +
  geom_hline(yintercept = median_depth, linetype = "dashed", color = "red") +
  geom_hline(yintercept = low_threshold, linetype = "dotted", color = "orange") +
  geom_hline(yintercept = high_threshold, linetype = "dotted", color = "orange") +
  annotate(
    "text",
    x = 1,
    y = median_depth,
    label = paste0("Median: ", format(round(median_depth), big.mark = ",")),
    vjust = -0.5,
    hjust = 0,
    color = "red",
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(values = c("Low" = "#d95f02", "High" = "#1b9e77")) +
  labs(
    title = "Outlier Samples (Outside +/-25% of Median)",
    subtitle = "Samples ranked by sequencing depth",
    x = "Sample",
    y = "Number of reads"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# 16:9 export for slides
ggsave(file.path(output_dir, "sequencing_depth_outliers_readable_16x9.jpg"), p_out, width = 16, height = 9, dpi = 300)

# ---------- Load additional requirement for axis formatting ----------
# If not installed, run: install.packages("scales")
library(scales) 

# ---------- Plot 3: Density Plot ----------
# Visualizes the continuous distribution of the cohort's read depth
p_density <- ggplot(df, aes(x = sequencing_depth)) +
  geom_density(fill = "steelblue", alpha = 0.5, color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = median_depth, linetype = "dashed", color = "red", linewidth = 1) +
  
  # Annotate the median line
  annotate("text", x = median_depth, y = Inf,
           label = paste("Median:", format(round(median_depth), big.mark = ",")),
           vjust = 2, hjust = -0.1, color = "red", fontface = "bold", size = 5) +
           
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "Distribution of Sequencing Depth",
    subtitle = paste0("Density estimation across all ", n_samples, " samples"),
    x = "Number of reads",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "sequencing_depth_density_16x9.jpg"), p_density, width = 16, height = 9, dpi = 300)


# ---------- Plot 4: Violin + Boxplot Hybrid ----------
# Combines the density shape with clear quartile boundaries and outlier dots
# We create a dummy categorical variable to plot a single aggregated violin
df$cohort <- "Full Cohort"

p_violin <- ggplot(df, aes(x = cohort, y = sequencing_depth)) +
  # The violin layer (density shape)
  geom_violin(fill = "steelblue", alpha = 0.3, color = NA) +
  
  # The embedded boxplot layer (quartiles and outliers)
  geom_boxplot(width = 0.15, fill = "white", color = "black", 
               outlier.color = "#d95f02", outlier.size = 2.5, outlier.alpha = 0.7) +
               
  geom_hline(yintercept = median_depth, linetype = "dashed", color = "red") +
  
  annotate("text", x = 1.15, y = median_depth,
           label = paste("Median:", format(round(median_depth), big.mark = ",")),
           vjust = -0.5, color = "red", fontface = "bold", size = 5) +
           
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Sequencing Depth Spread",
    subtitle = "Violin plot showing data density with embedded quartiles",
    x = "",
    y = "Number of reads"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(file.path(output_dir, "sequencing_depth_violin_16x9.jpg"), p_violin, width = 10, height = 9, dpi = 300)
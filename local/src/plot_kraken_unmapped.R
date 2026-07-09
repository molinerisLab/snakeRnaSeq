#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)

# ---------- CLI Args ----------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[1], "plots/sequencing_unmapped_kraken.csv")
output_dir <- ifelse(length(args) >= 2, args[2], "plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- Load ----------
clean_sample <- function(x) {
  gsub("_001(_unmapped)?_R1\\.fastq$", "", x)
}

df <- read.csv(input_file, header = TRUE)
required_cols <- c("name", "unmapped_depth", "kraken_classified_reads")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

df$name <- clean_sample(df$name)
df$unmapped_depth        <- as.numeric(df$unmapped_depth)
df$kraken_classified_reads <- as.numeric(df$kraken_classified_reads)

# Guard against kraken > unmapped
n_bad <- sum(df$kraken_classified_reads > df$unmapped_depth, na.rm = TRUE)
if (n_bad > 0) {
  warning(paste(n_bad, "samples have kraken_classified_reads > unmapped_depth; clamping."))
  df$kraken_classified_reads <- pmin(df$kraken_classified_reads, df$unmapped_depth)
}

# Rank by kraken classified reads
df <- df %>% arrange(kraken_classified_reads)
df$name <- factor(df$name, levels = df$name)

# Thresholds (raw values for grouping)
median_kraken  <- median(df$kraken_classified_reads, na.rm = TRUE)
low_threshold  <- median_kraken * 0.75
high_threshold <- median_kraken * 1.25

# Values used for plotting on log scale
df <- df %>%
  mutate(
    unmapped_plot = log10(unmapped_depth + 1),
    kraken_plot   = log10(kraken_classified_reads + 1)
  )

median_kraken_plot <- log10(median_kraken + 1)
low_threshold_plot <- log10(low_threshold + 1)
high_threshold_plot <- log10(high_threshold + 1)

df <- df %>%
  mutate(group = case_when(
    kraken_classified_reads < low_threshold  ~ "Low",
    kraken_classified_reads > high_threshold ~ "High",
    TRUE                                     ~ "Median range"
  ))
df$group <- factor(df$group, levels = c("Low", "Median range", "High"))

# Palette: grey for total unmapped, blue for kraken subset
metric_palette <- c(
  "Unmapped reads"          = "#bdbdbd",  # grey
  "Kraken classified reads" = "#0577b4"   # light blue
)

# ---------- Plot helper (two-layer identity bars) ----------
make_plot <- function(data, hide_x = FALSE, x_angle = 60) {
  data <- data %>%
    arrange(kraken_classified_reads) %>%
    mutate(name = factor(name, levels = name))

  ggplot(data, aes(x = name)) +
    # Full unmapped bar (grey) — drawn first so blue sits on top
    geom_col(aes(y = unmapped_plot, fill = "Unmapped reads"),
             position = "identity", width = 0.85) +
    # Kraken subset bar (blue) — drawn second, starts from 0
    geom_col(aes(y = kraken_plot, fill = "Kraken classified reads"),
             position = "identity", width = 0.85) +
    geom_hline(yintercept = median_kraken_plot, linetype = "dashed", color = "#b71c1c") +
    geom_hline(yintercept = low_threshold_plot, linetype = "dotted", color = "#ef6c00") +
    geom_hline(yintercept = high_threshold_plot, linetype = "dotted", color = "#ef6c00") +
    scale_fill_manual(
      values = metric_palette,
      breaks = c("Unmapped reads", "Kraken classified reads")   # legend order
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position    = "top",
      axis.text.x = if (hide_x)
        element_blank()
      else
        element_text(angle = x_angle, hjust = 1, vjust = 1, size = 7),
      axis.ticks.x = if (hide_x) element_blank() else element_line()
    )
}

# ---------- Plot helper (linear bars) ----------
make_plot_linear <- function(data, hide_x = FALSE, x_angle = 60) {
  data <- data %>%
    arrange(kraken_classified_reads) %>%
    mutate(name = factor(name, levels = name))

  ggplot(data, aes(x = name)) +
    # Full unmapped bar (grey) — drawn first so blue sits on top
    geom_col(aes(y = unmapped_depth, fill = "Unmapped reads"),
             position = "identity", width = 0.85) +
    # Kraken subset bar (blue) — drawn second, starts from 0
    geom_col(aes(y = kraken_classified_reads, fill = "Kraken classified reads"),
             position = "identity", width = 0.85) +
    geom_hline(yintercept = median_kraken, linetype = "dashed", color = "#b71c1c") +
    geom_hline(yintercept = low_threshold, linetype = "dotted", color = "#ef6c00") +
    geom_hline(yintercept = high_threshold, linetype = "dotted", color = "#ef6c00") +
    scale_fill_manual(
      values = metric_palette,
      breaks = c("Unmapped reads", "Kraken classified reads")   # legend order
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position    = "top",
      axis.text.x = if (hide_x)
        element_blank()
      else
        element_text(angle = x_angle, hjust = 1, vjust = 1, size = 7),
      axis.ticks.x = if (hide_x) element_blank() else element_line()
    )
}

# ---------- Plot 1: all samples ----------
p_all <- make_plot(df, hide_x = TRUE) +
  annotate("text",
           x     = nrow(df) * 0.02,
           y     = median_kraken_plot,   # use log-scale y position
           label = paste0("Median Kraken: ", format(round(median_kraken), big.mark = ",")),
           vjust = -0.8, hjust = 0, color = "#b71c1c", fontface = "bold") +
  labs(
    title    = "Kraken Classified Reads vs Total Unmapped Reads Across All Samples",
    subtitle = "Grey = total unmapped reads; blue = Kraken-classified reads ",
    x        = "Samples ranked by Kraken classified reads",
    y        = "Read count",
    fill     = NULL
  )

ggsave(file.path(output_dir, "stacked_kraken_unmapped_overview.jpg"), p_all, width = 14, height = 6, dpi = 300)

p_all_lin <- make_plot_linear(df, hide_x = TRUE) +
  annotate("text",
           x     = nrow(df) * 0.02,
           y     = median_kraken,
           label = paste0("Median Kraken: ", format(round(median_kraken), big.mark = ",")),
           vjust = -0.8, hjust = 0, color = "#b71c1c", fontface = "bold") +
  labs(
    title    = "Kraken Classified Reads vs Total Unmapped Reads Across All Samples (Linear Scale)",
    x        = "Samples ranked by Kraken classified reads",
    y        = "Read count",
    fill     = NULL
  )

ggsave(file.path(output_dir, "stacked_kraken_unmapped_overview_linear.jpg"), p_all_lin, width = 14, height = 6, dpi = 300)

# ---------- Plot 2: outliers only ----------
df_out <- df %>%
  filter(group != "Median range") %>%
  arrange(kraken_classified_reads)

df_out$name <- factor(df_out$name, levels = df_out$name)

coverage_palette <- c(
  "Low"  = "#d95f02",
  "High" = "#1b9e77"
)

p_out <- ggplot(df_out, aes(x = name)) +
  # total unmapped as grey background column
  geom_col(aes(y = unmapped_plot), fill = "#bdbdbd", width = 0.85) +
  # kraken subset filled by Low/High group (no blue)
  geom_col(aes(y = kraken_plot, fill = group), width = 0.85) +
  geom_hline(yintercept = median_kraken_plot, linetype = "dashed", color = "#b71c1c") +
  geom_hline(yintercept = low_threshold_plot, linetype = "dotted", color = "#ef6c00") +
  geom_hline(yintercept = high_threshold_plot, linetype = "dotted", color = "#ef6c00") +
  annotate(
    "text",
    x = 1,
    y = median_kraken_plot,   # use log-scale y position
    label = paste0("Median Kraken: ", format(round(median_kraken), big.mark = ",")),
    vjust = -0.5, hjust = 0, color = "#b71c1c", fontface = "bold", size = 4
  ) +
  scale_fill_manual(
    values = coverage_palette,
    breaks = c("Low", "High"),
    drop = TRUE,
    name = "Coverage group"
  ) +
  labs(
    title = "Outlier Samples by Kraken Classified Reads",
    subtitle = "Grey = total unmapped; Low/High colors fill Kraken subset",
    x = "Sample",
    y = "log10(read count + 1)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "stacked_kraken_unmapped_outliers_16x9.jpg"), p_out, width = 22, height = 10, dpi = 600)

p_out_lin <- ggplot(df_out, aes(x = name)) +
  # total unmapped as grey background column
  geom_col(aes(y = unmapped_depth), fill = "#bdbdbd", width = 0.85) +
  # kraken subset filled by Low/High group (no blue)
  geom_col(aes(y = kraken_classified_reads, fill = group), width = 0.85) +
  geom_hline(yintercept = median_kraken, linetype = "dashed", color = "#b71c1c") +
  geom_hline(yintercept = low_threshold, linetype = "dotted", color = "#ef6c00") +
  geom_hline(yintercept = high_threshold, linetype = "dotted", color = "#ef6c00") +
  annotate(
    "text",
    x = 1,
    y = median_kraken,
    label = paste0("Median Kraken: ", format(round(median_kraken), big.mark = ",")),
    vjust = -0.5, hjust = 0, color = "#b71c1c", fontface = "bold", size = 4
  ) +
  scale_fill_manual(
    values = coverage_palette,
    breaks = c("Low", "High"),
    drop = TRUE,
    name = "Coverage group"
  ) +
  labs(
    title = "Outlier Samples by Kraken Classified Reads (Linear Scale)",
    subtitle = "Grey = total unmapped; Low/High colors fill Kraken subset",
    x = "Sample",
    y = "Read count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "stacked_kraken_unmapped_outliers_16x9_linear.jpg"), p_out_lin, width = 22, height = 10, dpi = 600)
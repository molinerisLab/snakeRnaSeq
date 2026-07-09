#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------- CLI Arguments ----------
args       <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot.R <frac_file> <num_file> [top_n] [threshold_pct] [output_dir] [names_file]\n",
       "  threshold_pct: prevalence threshold as a fraction, e.g. 0.9 for 90%")
}
frac_file  <- args[1]
num_file   <- args[2]
top_n      <- ifelse(length(args) >= 3, as.integer(args[3]), 10)
threshold  <- ifelse(length(args) >= 4, as.numeric(args[4]),  0.9)
output_dir <- ifelse(length(args) >= 5, args[5], ".")
names_file <- ifelse(length(args) >= 6, args[6], NA)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

read_tsv_compat <- function(path) {
  if ("show_col_types" %in% names(formals(readr::read_tsv))) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    suppressMessages(readr::read_tsv(path, progress = FALSE))
  }
}

# Load data
frac <- read_tsv_compat(frac_file)
num  <- read_tsv_compat(num_file)

id_cols <- c("name", "taxonomy_id", "taxonomy_lvl")
sample_cols <- setdiff(names(frac), id_cols)

# Ensure same samples
if (!setequal(sample_cols, setdiff(names(num), id_cols))) {
  stop("Sample columns differ between fraction and count files.")
}

# Convert to numeric
frac <- frac %>% mutate(across(all_of(sample_cols), as.numeric))
num  <- num  %>% mutate(across(all_of(sample_cols), as.numeric))

n_samples <- length(sample_cols)

# ---- Species-level stats ----
species_frac <- frac %>%
  filter(taxonomy_lvl == "S") %>%
  mutate(
    n_present_frac = rowSums(across(all_of(sample_cols), ~ !is.na(.) & . > 0)),
    mean_fraction = rowMeans(across(all_of(sample_cols)), na.rm = TRUE),
    total_fraction = rowSums(across(all_of(sample_cols)), na.rm = TRUE)
  )

species_num <- num %>%
  filter(taxonomy_lvl == "S") %>%
  mutate(
    n_present_num = rowSums(across(all_of(sample_cols), ~ !is.na(.) & . > 0)),
    total_counts = rowSums(across(all_of(sample_cols)), na.rm = TRUE)
  )

# Merge stats
species_stats <- species_frac %>%
  inner_join(
    species_num %>% select(all_of(id_cols), n_present_num, total_counts),
    by = id_cols
  )


# ---- Top N species OVERALL (by mean fraction) ----
top10_overall <- species_stats %>%
  arrange(desc(mean_fraction)) %>%
  slice_head(n = top_n)

# Save summary
top10_overall_summary <- top10_overall %>%
  select(
    all_of(id_cols),
    n_present_frac,
    n_present_num,
    mean_fraction,
    total_fraction,
    total_counts
  )

write_tsv(top10_overall_summary, file.path(output_dir, paste0("top", top_n, "_species_overall.tsv")))

# FRACTION
plot_frac_overall <- frac %>%
  inner_join(top10_overall %>% select(all_of(id_cols)), by = id_cols) %>%
  select(name, all_of(sample_cols)) %>%
  pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
  mutate(type = "Fraction")

# COUNTS
plot_num_overall <- num %>%
  inner_join(top10_overall %>% select(all_of(id_cols)), by = id_cols) %>%
  select(name, all_of(sample_cols)) %>%
  pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
  mutate(type = "Counts")

# Combine
plot_df_overall <- bind_rows(plot_frac_overall, plot_num_overall)

# Order samples for overall plot
sample_order_overall <- plot_frac_overall %>%
  group_by(sample) %>%
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(total) %>%
  pull(sample)

plot_df_overall <- plot_df_overall %>%
  mutate(sample = factor(sample, levels = sample_order_overall))

# ---- Plot ----
p_overall <- ggplot(plot_df_overall, aes(x = sample, y = value, fill = name)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Top ", top_n, " species overall (by mean fraction)"),
    subtitle = "Stacked abundance: fraction and raw counts",
    x = "Samples",
    y = "Abundance",
    fill = "Species"
  )

ggsave(
  file.path(output_dir, paste0("top", top_n, "_species_overall_stacked.png")),
  p_overall,
  width = 16,
  height = 8,
  dpi = 300
)


# ---- Strict filter: present in ALL samples in BOTH tables ----
strict_species <- species_stats %>%
  filter(n_present_frac == n_samples, n_present_num == n_samples)

if (nrow(strict_species) == 0) {
  warning("No species present in all samples in both tables — skipping the 'all samples' plot.")
} else {

  if (nrow(strict_species) < top_n) {
    warning("Only ", nrow(strict_species), " species meet strict criterion (requested top_n = ", top_n, ").")
  }

  # ---- Top N by mean fraction (present in ALL samples) ----
  top10 <- strict_species %>%
    arrange(desc(mean_fraction)) %>%
    slice_head(n = top_n)

  # Save summary
  top10_summary <- top10 %>%
    select(
      all_of(id_cols),
      n_present_frac,
      n_present_num,
      mean_fraction,
      total_fraction,
      total_counts
    )

  write_tsv(top10_summary, file.path(output_dir, paste0("top", top_n, "_species_present_all", n_samples, ".tsv")))

  # ---- Prepare plotting data ----

  # FRACTION
  plot_frac <- frac %>%
    inner_join(top10 %>% select(all_of(id_cols)), by = id_cols) %>%
    select(name, all_of(sample_cols)) %>%
    pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
    mutate(type = "Fraction")

  # COUNTS
  plot_num <- num %>%
    inner_join(top10 %>% select(all_of(id_cols)), by = id_cols) %>%
    select(name, all_of(sample_cols)) %>%
    pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
    mutate(type = "Counts")

  # Combine
  plot_df <- bind_rows(plot_frac, plot_num)

  # Order samples consistently (by total fraction)
  sample_order <- plot_frac %>%
    group_by(sample) %>%
    summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(total) %>%
    pull(sample)

  plot_df <- plot_df %>%
    mutate(sample = factor(sample, levels = sample_order))

  # ---- Plot ----
  p <- ggplot(plot_df, aes(x = sample, y = value, fill = name)) +
    geom_col(width = 1, color = "black", linewidth = 0.2) +
    facet_wrap(~type, scales = "free_y", ncol = 1) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste0("Top ", top_n, " species present in all ", n_samples, " samples"),
      subtitle = "Stacked abundance: fraction and raw counts",
      x = "Samples",
      y = "Abundance",
      fill = "Species"
    )

  ggsave(
    file.path(output_dir, paste0("top", top_n, "_species_present_all", n_samples, "_stacked.png")),
    p,
    width = 16,
    height = 8,
    dpi = 300
  )
}



# ---- Top N species present in >= threshold% of samples ----

threshold_count <- ceiling(threshold * n_samples)
threshold_pct   <- round(threshold * 100)

species_90 <- species_stats %>%
  filter(n_present_frac >= threshold_count, n_present_num >= threshold_count)

if (nrow(species_90) == 0) {
  stop("No species found in >=", threshold_pct, "% of samples.")
}

if (nrow(species_90) < top_n) {
  warning("Only ", nrow(species_90), " species meet the ", threshold_pct, "% criterion (requested top_n = ", top_n, ").")
}

top10_90 <- species_90 %>%
  arrange(desc(mean_fraction)) %>%
  slice_head(n = top_n)

# Save summary
top10_90_summary <- top10_90 %>%
  select(
    all_of(id_cols),
    n_present_frac,
    n_present_num,
    mean_fraction,
    total_fraction,
    total_counts
  )

write_tsv(top10_90_summary, file.path(output_dir, paste0("top", top_n, "_species_present_", threshold_pct, "perc.tsv")))

# ---- Prepare plotting data ----

# FRACTION
plot_frac_90 <- frac %>%
  inner_join(top10_90 %>% select(all_of(id_cols)), by = id_cols) %>%
  select(name, all_of(sample_cols)) %>%
  pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
  mutate(type = "Fraction")

# COUNTS
plot_num_90 <- num %>%
  inner_join(top10_90 %>% select(all_of(id_cols)), by = id_cols) %>%
  select(name, all_of(sample_cols)) %>%
  pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
  mutate(type = "Counts")

# Combine
plot_df_90 <- bind_rows(plot_frac_90, plot_num_90)

# Order samples for the 90% plot
sample_order_90 <- plot_frac_90 %>%
  group_by(sample) %>%
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(total) %>%
  pull(sample)

plot_df_90 <- plot_df_90 %>%
  mutate(sample = factor(sample, levels = sample_order_90))

# ---- Plot ----
p_90 <- ggplot(plot_df_90, aes(x = sample, y = value, fill = name)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Top ", top_n, " species present in ≥", threshold_pct, "% of samples"),
    subtitle = paste0("Threshold: ≥", threshold_count, " / ", n_samples, " samples"),
    x = "Samples",
    y = "Abundance",
    fill = "Species"
  )

ggsave(
  file.path(output_dir, paste0("top", top_n, "_species_present_", threshold_pct, "perc_stacked.png")),
  p_90,
  width = 16,
  height = 8,
  dpi = 300
)


# ---- Optional: subset plot from names_file ----
if (!is.na(names_file) && file.exists(names_file)) {

  target_samples <- readLines(names_file, warn = FALSE)
  target_samples <- trimws(target_samples)
  target_samples <- sub("_R1\\.fastq$", "", target_samples)
  target_samples <- target_samples[nzchar(target_samples)]
  target_samples <- target_samples[target_samples != "name"]  # safe if header exists
  target_samples <- unique(target_samples)

  if (length(target_samples) == 0) {
    warning(names_file, " is empty or malformed — skipping subset plot.")
  } else {

    missing_samples <- setdiff(target_samples, sample_cols)
    if (length(missing_samples) > 0) {
      stop("These samples are missing in your data:\n", paste(missing_samples, collapse = ", "))
    }

    sample_cols_subset <- target_samples
    n_subset <- length(sample_cols_subset)

    species_frac_subset <- frac %>%
      filter(taxonomy_lvl == "S") %>%
      mutate(
        n_present_frac = rowSums(across(all_of(sample_cols_subset), ~ !is.na(.) & . > 0)),
        mean_fraction   = rowMeans(across(all_of(sample_cols_subset)), na.rm = TRUE),
        total_fraction  = rowSums(across(all_of(sample_cols_subset)), na.rm = TRUE)
      )

    species_num_subset <- num %>%
      filter(taxonomy_lvl == "S") %>%
      mutate(
        n_present_num = rowSums(across(all_of(sample_cols_subset), ~ !is.na(.) & . > 0)),
        total_counts  = rowSums(across(all_of(sample_cols_subset)), na.rm = TRUE)
      )

    species_stats_subset <- species_frac_subset %>%
      inner_join(
        species_num_subset %>% select(all_of(id_cols), n_present_num, total_counts),
        by = id_cols
      )

    strict_subset <- species_stats_subset %>%
      filter(n_present_frac == n_subset, n_present_num == n_subset)

    if (nrow(strict_subset) == 0) {
      warning("No species present in all selected samples — skipping subset plot.")
    } else {

      top_subset <- strict_subset %>%
        arrange(desc(mean_fraction)) %>%
        slice_head(n = top_n)

      plot_frac_subset <- frac %>%
        inner_join(top_subset %>% select(all_of(id_cols)), by = id_cols) %>%
        select(name, all_of(sample_cols_subset)) %>%
        pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
        mutate(type = "Fraction")

      plot_num_subset <- num %>%
        inner_join(top_subset %>% select(all_of(id_cols)), by = id_cols) %>%
        select(name, all_of(sample_cols_subset)) %>%
        pivot_longer(cols = -name, names_to = "sample", values_to = "value") %>%
        mutate(type = "Counts")

      plot_df_subset <- bind_rows(plot_frac_subset, plot_num_subset)

      sample_order_subset <- plot_frac_subset %>%
        group_by(sample) %>%
        summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
        arrange(total) %>%
        pull(sample)

      plot_df_subset <- plot_df_subset %>%
        mutate(sample = factor(sample, levels = sample_order_subset))

      p_subset <- ggplot(plot_df_subset, aes(x = sample, y = value, fill = name)) +
        geom_col(width = 1, color = "black", linewidth = 0.2) +
        facet_wrap(~type, scales = "free_y", ncol = 1) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
          axis.ticks.x = element_line(),
          panel.grid.major.x = element_blank(),
          legend.position = "right"
        ) +
        labs(
          title    = paste0("Top ", top_n, " species present in selected samples"),
          subtitle = paste("Strict: present in all", n_subset, "selected samples"),
          x = "Samples",
          y = "Abundance",
          fill = "Species"
        )

      ggsave(
        file.path(output_dir, paste0("top", top_n, "_species_selected_samples.png")),
        p_subset,
        width = 16,
        height = 8,
        dpi = 300
      )
    }
  }

} else if (!is.na(names_file)) {
  warning("names_file '", names_file, "' not found — skipping subset plot.")
}
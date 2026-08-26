#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})


# ============================================================
# Arguments
# ============================================================

option_list <- list(

  make_option(
    "--results_dir",
    type = "character",
    help = "Directory containing per-sample Kraken comparison CSV files"
  ),

  make_option(
    "--pdf",
    type = "character",
    help = "Output PDF report"
  ),

  make_option(
    "--summary",
    type = "character",
    help = "Output combined summary CSV"
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)


# ============================================================
# Validate arguments
# ============================================================

required_args <- c(
  "results_dir",
  "pdf",
  "summary"
)

missing_args <- required_args[
  vapply(
    required_args,
    function(x) {
      is.null(opt[[x]]) || !nzchar(opt[[x]])
    },
    logical(1)
  )
]

if (length(missing_args) > 0) {
  stop(
    "Missing required argument(s): ",
    paste(missing_args, collapse = ", ")
  )
}

if (!dir.exists(opt$results_dir)) {
  stop(
    "Results directory does not exist: ",
    opt$results_dir
  )
}

dir.create(
  dirname(opt$pdf),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(opt$summary),
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Safe correlation
# ============================================================

safe_cor <- function(x, y, method = "spearman") {

  keep <- complete.cases(x, y)

  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < 3 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(NA_real_)
  }

  cor(
    x,
    y,
    method = method
  )
}


# ============================================================
# Locate per-sample files
# ============================================================

summary_files <- list.files(
  opt$results_dir,
  pattern = "_summary\\.csv$",
  full.names = TRUE
)

taxon_files <- list.files(
  opt$results_dir,
  pattern = "_taxon_comparison\\.csv$",
  full.names = TRUE
)


if (length(summary_files) == 0) {
  stop("No per-sample summary files found.")
}

if (length(taxon_files) == 0) {
  stop("No per-sample taxon comparison files found.")
}


# ============================================================
# Read per-sample summaries
# ============================================================

sample_summary <- map_dfr(
  summary_files,
  function(file) {

    x <- read_csv(
      file,
      show_col_types = FALSE
    )

    # Keep sample from file if present;
    # otherwise infer it from filename.
    if (!"sample" %in% names(x)) {
      x$sample <- str_remove(
        basename(file),
        "_summary\\.csv$"
      )
    }

    x
  }
) %>%
  distinct(
    sample,
    .keep_all = TRUE
  ) %>%
  arrange(sample)


# ============================================================
# Additional summary metrics
# ============================================================

sample_summary <- sample_summary %>%
  mutate(

    classified_percent_change =
      new_classified_percent -
      old_classified_percent,

    classified_reads_change =
      new_classified_reads -
      old_classified_reads,

    detected_taxid_change =
      new_taxa -
      old_taxa,

    new_to_old_unique_ratio =
      if_else(
        old_only_taxa > 0,
        new_only_taxa / old_only_taxa,
        NA_real_
      )
  )


# ============================================================
# Save combined overview
# ============================================================

write_csv(
  sample_summary,
  opt$summary
)


# ============================================================
# Read all taxon comparison tables
# ============================================================

taxon_all <- map_dfr(
  taxon_files,
  function(file) {

    sample_name <- str_remove(
      basename(file),
      "_taxon_comparison\\.csv$"
    )

    read_csv(
      file,
      show_col_types = FALSE
    ) %>%
      mutate(
        sample = sample_name
      )
  }
)


# ============================================================
# Sample ordering
# ============================================================

sample_levels <- sample_summary$sample

sample_summary <- sample_summary %>%
  mutate(
    sample = factor(
      sample,
      levels = sample_levels
    )
  )

taxon_all <- taxon_all %>%
  mutate(
    sample = factor(
      sample,
      levels = sample_levels
    )
  )


# ============================================================
# Plot 1
# Classification percentage: old versus new
# ============================================================

classification_long <- sample_summary %>%
  select(
    sample,
    old_classified_percent,
    new_classified_percent
  ) %>%
  pivot_longer(
    cols = c(
      old_classified_percent,
      new_classified_percent
    ),
    names_to = "database",
    values_to = "classified_percent"
  ) %>%
  mutate(
    database = recode(
      database,
      old_classified_percent = "Old database",
      new_classified_percent = "New database"
    )
  )


p1 <- ggplot(
  classification_long,
  aes(
    x = sample,
    y = classified_percent,
    group = database,
    shape = database
  )
) +
  geom_point(
    size = 2
  ) +
  geom_line(
    aes(
      linetype = database
    )
  ) +
  labs(
    title = "Classification rate across samples",
    x = "Sample",
    y = "Classified reads (%)",
    shape = NULL,
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# ============================================================
# Plot 2
# Change in classification percentage
# ============================================================

p2 <- ggplot(
  sample_summary,
  aes(
    x = classified_percent_change,
    y = reorder(
      sample,
      classified_percent_change
    )
  )
) +
  geom_col() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Change in classification rate",
    subtitle = "Positive values indicate higher classification with the new database",
    x = "New database − Old database (percentage points)",
    y = "Sample"
  ) +
  theme_bw()


# ============================================================
# Plot 3
# Overall abundance agreement
# ============================================================

p3 <- ggplot(
  sample_summary,
  aes(
    x = sample,
    y = spearman_rho
  )
) +
  geom_col() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_cartesian(
    ylim = c(-1, 1)
  ) +
  labs(
    title = "Abundance agreement between databases",
    subtitle = "Spearman correlation among TaxIDs detected by both databases",
    x = "Sample",
    y = "Spearman rho"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# ============================================================
# Plot 4
# Number of detected TaxIDs
# ============================================================

taxa_long <- sample_summary %>%
  select(
    sample,
    old_taxa,
    new_taxa
  ) %>%
  pivot_longer(
    cols = c(
      old_taxa,
      new_taxa
    ),
    names_to = "database",
    values_to = "detected_taxids"
  ) %>%
  mutate(
    database = recode(
      database,
      old_taxa = "Old database",
      new_taxa = "New database"
    )
  )


p4 <- ggplot(
  taxa_long,
  aes(
    x = sample,
    y = detected_taxids,
    group = database,
    shape = database
  )
) +
  geom_point(
    size = 2
  ) +
  geom_line(
    aes(
      linetype = database
    )
  ) +
  labs(
    title = "Detected TaxIDs across samples",
    subtitle = "Counts include all classified taxonomic ranks",
    x = "Sample",
    y = "Number of detected TaxIDs",
    shape = NULL,
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# ============================================================
# Plot 5
# Database-specific detections
# ============================================================

unique_taxa_long <- sample_summary %>%
  select(
    sample,
    old_only_taxa,
    new_only_taxa
  ) %>%
  pivot_longer(
    cols = c(
      old_only_taxa,
      new_only_taxa
    ),
    names_to = "database",
    values_to = "taxa"
  ) %>%
  mutate(
    database = recode(
      database,
      old_only_taxa = "Old database only",
      new_only_taxa = "New database only"
    )
  )


p5 <- ggplot(
  unique_taxa_long,
  aes(
    x = sample,
    y = taxa,
    fill = database
  )
) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title = "Database-specific taxonomic detections",
    x = "Sample",
    y = "Number of TaxIDs",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# ============================================================
# Plot 6
# Rank-specific abundance agreement
# ============================================================

rank_order <- c(
  "D", "P", "C", "O", "F", "G", "S"
)


rank_agreement <- taxon_all %>%
  filter(
    rank %in% rank_order
  ) %>%
  group_by(
    sample,
    rank
  ) %>%
  summarise(

    shared_taxa = sum(
      clade_reads_old > 0 &
        clade_reads_new > 0
    ),

    spearman = {

      keep <-
        clade_reads_old > 0 &
        clade_reads_new > 0

      safe_cor(
        abundance_old[keep],
        abundance_new[keep],
        method = "spearman"
      )
    },

    .groups = "drop"
  ) %>%
  mutate(
    rank = factor(
      rank,
      levels = rank_order
    )
  )


p6 <- ggplot(
  rank_agreement,
  aes(
    x = rank,
    y = sample,
    fill = spearman
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0,
    na.value = "grey90"
  ) +
  labs(
    title = "Database agreement by taxonomic rank",
    subtitle = "Spearman abundance correlation for each sample",
    x = "Kraken rank",
    y = "Sample",
    fill = "Spearman\nrho"
  ) +
  theme_bw()


# ============================================================
# Plot 7
# Largest genus/species abundance changes across samples
#
# Restricting this plot to genus/species avoids root,
# domain, phylum, etc. dominating the result.
# ============================================================

largest_changes <- taxon_all %>%
  filter(
    rank %in% c("G", "S"),
    clade_reads_old > 0 |
      clade_reads_new > 0
  ) %>%
  mutate(
    label = paste0(
      as.character(sample),
      " | ",
      taxon
    )
  ) %>%
  slice_max(
    order_by = abs(abundance_difference),
    n = 30,
    with_ties = FALSE
  )


p7 <- ggplot(
  largest_changes,
  aes(
    x = abundance_difference,
    y = reorder(
      label,
      abundance_difference
    )
  )
) +
  geom_col() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Largest genus/species abundance changes",
    subtitle = "Across all samples",
    x = "New database − Old database (percentage points)",
    y = NULL
  ) +
  theme_bw()


# ============================================================
# Plot 8
# Classification change versus abundance agreement
# ============================================================

p8 <- ggplot(
  sample_summary,
  aes(
    x = classified_percent_change,
    y = spearman_rho,
    label = sample
  )
) +
  geom_point(
    size = 3
  ) +
  geom_text(
    nudge_y = 0.015,
    check_overlap = TRUE
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Overall impact of database change",
    x = "Change in classification rate (percentage points)",
    y = "Abundance agreement (Spearman rho)"
  ) +
  theme_bw()


# ============================================================
# Save combined PDF
# ============================================================

pdf(
  opt$pdf,
  width = 12,
  height = 8
)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)

dev.off()


# ============================================================
# Console summary
# ============================================================

cat("\n==========================================\n")
cat("Kraken Database Comparison: All Samples\n")
cat("==========================================\n")

cat(
  "Samples:",
  nrow(sample_summary),
  "\n"
)

cat(
  "Median change in classification rate:",
  round(
    median(
      sample_summary$classified_percent_change,
      na.rm = TRUE
    ),
    4
  ),
  "percentage points\n"
)

cat(
  "Median Spearman correlation:",
  round(
    median(
      sample_summary$spearman_rho,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Mean new-only TaxIDs:",
  round(
    mean(
      sample_summary$new_only_taxa,
      na.rm = TRUE
    ),
    1
  ),
  "\n"
)

cat(
  "Mean old-only TaxIDs:",
  round(
    mean(
      sample_summary$old_only_taxa,
      na.rm = TRUE
    ),
    1
  ),
  "\n"
)

cat("==========================================\n")
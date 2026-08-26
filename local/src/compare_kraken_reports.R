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
    "--old_db",
    type = "character",
    help = "Kraken 2 report generated using the old database"
  ),
  make_option(
    "--new_db",
    type = "character",
    help = "Kraken 2 report generated using the new database"
  ),
  make_option(
    "--sample",
    type = "character",
    help = "Sample name"
  ),
  make_option(
    "--pdf",
    type = "character",
    help = "Output PDF containing comparison plots"
  ),
  make_option(
    "--table",
    type = "character",
    help = "Output CSV containing taxon-level comparison"
  ),
  make_option(
    "--summary",
    type = "character",
    help = "Output CSV containing database-level summary"
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)


# ============================================================
# Validate arguments
# ============================================================

required_args <- c("old_db", "new_db", "sample", "pdf", "table", "summary")

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

if (!file.exists(opt$old_db)) {
  stop("Old Kraken report does not exist: ", opt$old_db)
}

if (!file.exists(opt$new_db)) {
  stop("New Kraken report does not exist: ", opt$new_db)
}


# ============================================================
# Create output directories
# ============================================================

for (path in c(opt$pdf, opt$table, opt$summary)) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# Function to read Kraken 2 report
# ============================================================

read_k2report <- function(file) {
  x <- read_tsv(
    file,
    col_names = FALSE,
    col_types = "dddccc",
    show_col_types = FALSE,
    trim_ws = TRUE,
    progress = FALSE
  )
  if (ncol(x) != 6) {
    stop(
      "Expected 6 columns in Kraken report, but found ",
      ncol(x),
      " in: ",
      file
    )
  }
  colnames(x) <- c("percentage", "clade_reads", "direct_reads", "rank", "taxid", "taxon")
  x %>%
    mutate(
      percentage = as.numeric(percentage),
      clade_reads = as.numeric(clade_reads),
      direct_reads = as.numeric(direct_reads),
      rank = as.character(rank),
      taxid = as.character(taxid),
      taxon = str_trim(taxon)
    )
}


# ============================================================
# Safe correlation function
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
  cor(x, y, method = method)
}


# ============================================================
# Read reports
# ============================================================

old <- read_k2report(opt$old_db)
new <- read_k2report(opt$new_db)


# ============================================================
# Calculate read totals
#
# IMPORTANT:
# Do not sum clade_reads across the report because Kraken
# clade counts are hierarchical and would count reads multiple
# times.
#
# direct_reads represents reads assigned directly to each node,
# so summing direct_reads gives the total number of reads.
# ============================================================

old_total_reads <- sum(old$direct_reads, na.rm = TRUE)

new_total_reads <- sum(new$direct_reads, na.rm = TRUE)

old_unclassified_reads <- sum(old$direct_reads[old$rank == "U"], na.rm = TRUE)

new_unclassified_reads <- sum(new$direct_reads[new$rank == "U"], na.rm = TRUE)

old_classified_reads <- old_total_reads - old_unclassified_reads

new_classified_reads <- new_total_reads - new_unclassified_reads


# Same sample should normally contain the same number of reads
if (old_total_reads != new_total_reads) {

  warning(
    "Total read counts differ between reports for sample ",
    opt$sample,
    ": old = ",
    old_total_reads,
    ", new = ",
    new_total_reads
  )
}


# ============================================================
# Keep detected classified taxa
#
# clade_reads > 0 prevents zero-count taxa from being treated
# as detected if Kraken was run with --report-zero-counts.
# ============================================================

old_classified <- old %>%
  filter(
    rank != "U",
    clade_reads > 0
  )

new_classified <- new %>%
  filter(
    rank != "U",
    clade_reads > 0
  )


# ============================================================
# Check TaxID uniqueness
# ============================================================

if (anyDuplicated(old_classified$taxid)) {
  stop("Duplicate TaxIDs detected in old Kraken report.")
}

if (anyDuplicated(new_classified$taxid)) {
  stop("Duplicate TaxIDs detected in new Kraken report.")
}


# ============================================================
# Join reports using NCBI TaxID
# ============================================================

comparison <- full_join(

  old_classified %>%
    select(
      taxid,
      taxon_old = taxon,
      rank_old = rank,
      percentage_reported_old = percentage,
      clade_reads_old = clade_reads,
      direct_reads_old = direct_reads
    ),

  new_classified %>%
    select(
      taxid,
      taxon_new = taxon,
      rank_new = rank,
      percentage_reported_new = percentage,
      clade_reads_new = clade_reads,
      direct_reads_new = direct_reads
    ),

  by = "taxid"
)


# ============================================================
# Identify differences
# ============================================================

comparison <- comparison %>%
  mutate(

    # Prefer current/new taxonomy when available
    taxon = coalesce(taxon_new, taxon_old),

    rank = coalesce(rank_new, rank_old),

    # Record taxonomy changes between databases
    taxon_name_changed =
      !is.na(taxon_old) &
      !is.na(taxon_new) &
      taxon_old != taxon_new,

    rank_changed =
      !is.na(rank_old) &
      !is.na(rank_new) &
      rank_old != rank_new,

    percentage_reported_old =
      replace_na(percentage_reported_old, 0),

    percentage_reported_new =
      replace_na(percentage_reported_new, 0),

    clade_reads_old =
      replace_na(clade_reads_old, 0),

    clade_reads_new =
      replace_na(clade_reads_new, 0),

    direct_reads_old =
      replace_na(direct_reads_old, 0),

    direct_reads_new =
      replace_na(direct_reads_new, 0),

    # Recalculate exact abundance from counts rather than using
    # Kraken's rounded report percentage
    abundance_old =
      if (old_total_reads > 0) {
        100 * clade_reads_old / old_total_reads
      } else {
        0
      },

    abundance_new =
      if (new_total_reads > 0) {
        100 * clade_reads_new / new_total_reads
      } else {
        0
      },

    abundance_difference =
      abundance_new - abundance_old,

    detection = case_when(

      clade_reads_old > 0 &
        clade_reads_new > 0 ~
        "Both",

      clade_reads_old > 0 &
        clade_reads_new == 0 ~
        "Old_database_only",

      clade_reads_old == 0 &
        clade_reads_new > 0 ~
        "New_database_only",

      TRUE ~
        "Neither"
    )
  )


# ============================================================
# Save taxon comparison table
# ============================================================

write_csv(comparison, opt$table)


# ============================================================
# Database-level summary
# ============================================================

summary <- tibble(
  sample = opt$sample,
  old_total_reads = old_total_reads,
  new_total_reads = new_total_reads,
  old_unclassified_reads = old_unclassified_reads,
  new_unclassified_reads = new_unclassified_reads,
  old_classified_reads = old_classified_reads,
  new_classified_reads = new_classified_reads,
  old_classified_percent = 
    if_else(
      old_total_reads > 0,
      100 * old_classified_reads / old_total_reads,
      NA_real_
    ),
  new_classified_percent =
    if_else(
      new_total_reads > 0,
      100 * new_classified_reads / new_total_reads,
      NA_real_
    ),
  old_taxa = n_distinct(old_classified$taxid),
  new_taxa = n_distinct(new_classified$taxid),
  shared_taxa = sum(comparison$detection == "Both"),
  old_only_taxa = sum(comparison$detection == "Old_database_only"),
  new_only_taxa = sum(comparison$detection == "New_database_only"),
  taxon_name_changes = sum(comparison$taxon_name_changed, na.rm = TRUE),
  rank_changes = sum(comparison$rank_changed, na.rm = TRUE)
)


# ============================================================
# Overall abundance similarity
# ============================================================

shared <- comparison %>%
  filter(clade_reads_old > 0, clade_reads_new > 0)

spearman_rho <- safe_cor(
  shared$abundance_old,
  shared$abundance_new,
  method = "spearman"
)
pearson_r <- safe_cor(
  shared$abundance_old,
  shared$abundance_new,
  method = "pearson"
)

summary <- summary %>%
  mutate(
    spearman_rho = spearman_rho,
    pearson_r = pearson_r
  )


write_csv(summary, opt$summary)

# ============================================================
# Plot 1: Top taxa
# ============================================================

top_taxa <- comparison %>%
  filter(rank %in% c("P", "C", "O", "F", "G", "S")) %>%
  mutate(max_abundance = pmax(abundance_old, abundance_new)) %>%
  slice_max(
    order_by = max_abundance,
    n = 20,
    with_ties = FALSE
  )

top_long <- comparison %>%
  filter(taxid %in% top_taxa$taxid) %>%
  select(taxid, taxon, abundance_old, abundance_new) %>%
  pivot_longer(
    cols = c(abundance_old, abundance_new),
    names_to = "database",
    values_to = "abundance"
  ) %>%
  mutate(database = recode(database, abundance_old = "Old database", abundance_new = "New database"))

p1 <- ggplot(top_long, aes(x = abundance, y = reorder(taxon, abundance))) +
  geom_col() +
  facet_wrap(~database) +
  labs(
    title = paste("Top taxa -", opt$sample),
    x = "Percentage of reads",
    y = NULL
  ) +
  theme_bw()

# ============================================================
# Plot 2: Database correlation
# ============================================================

p2 <- ggplot(shared, aes(x = abundance_old, y = abundance_new)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Abundance agreement between databases",
    subtitle = paste(
      "Spearman rho =",
      ifelse(
        is.na(spearman_rho),
        "NA",
        round(spearman_rho, 3)
      )
    ),
    x = "Old database (%)",
    y = "New database (%)"
  ) +
  theme_bw()

# ============================================================
# Plot 3: Largest abundance differences
# ============================================================

difference_data <- comparison %>%
  filter(clade_reads_old > 0 | clade_reads_new > 0) %>%
  slice_max(order_by = abs(abundance_difference), n = 20, with_ties = FALSE)

p3 <- ggplot(difference_data, aes(x = abundance_difference, y = reorder(taxon, abundance_difference))) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Largest abundance differences",
    x = "New database − Old database (%)",
    y = NULL
  ) +
  theme_bw()

# ============================================================
# Plot 4: Detection overlap
# ============================================================

overlap <- comparison %>%
  filter(detection != "Neither") %>%
  count(detection)

p4 <- ggplot(overlap, aes(x = detection, y = n)) +
  geom_col() +
  labs(
    title = "Taxonomic detection overlap",
    x = NULL,
    y = "Number of TaxIDs"
  ) +
  theme_bw()

# ============================================================
# Plot 5: Rank-specific agreement
# ============================================================

rank_comparison <- comparison %>%
  filter(rank %in% c("D", "P", "C", "O", "F", "G", "S")) %>%
  group_by(rank) %>%
  summarise(
    shared_taxa =
      sum(clade_reads_old > 0 & clade_reads_new > 0),
    spearman = {
      keep <- clade_reads_old > 0 & clade_reads_new > 0
      safe_cor(
        abundance_old[keep],
        abundance_new[keep],
        method = "spearman"
      )
    },
    .groups = "drop"
  )

p5 <- ggplot(rank_comparison, aes(x = rank, y = spearman)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Database agreement by taxonomic rank",
    x = "Kraken rank",
    y = "Spearman correlation"
  ) +
  theme_bw()

# ============================================================
# Save PDF
# ============================================================

pdf(opt$pdf, width = 11, height = 8)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

dev.off()


# ============================================================
# Console summary
# ============================================================

cat("\n==========================================\n")
cat("Kraken Database Comparison\n")
cat("==========================================\n")

cat("Sample:", opt$sample, "\n\n")

cat("Old database total reads:", summary$old_total_reads, "\n")

cat("New database total reads:", summary$new_total_reads, "\n\n")

cat("Old database classified reads:", summary$old_classified_reads, "\n")

cat("New database classified reads:", summary$new_classified_reads, "\n\n")

cat("Old database classified (%):", round(summary$old_classified_percent, 2), "\n")

cat("New database classified (%):", round(summary$new_classified_percent, 2), "\n\n")

cat("Old database taxa:", summary$old_taxa, "\n")

cat("New database taxa:", summary$new_taxa, "\n")

cat("Shared taxa:", summary$shared_taxa, "\n")

cat("Old database only:", summary$old_only_taxa, "\n")

cat("New database only:", summary$new_only_taxa, "\n\n")

cat("Taxon name changes:", summary$taxon_name_changes, "\n")

cat("Rank changes:", summary$rank_changes, "\n\n")

cat("Spearman correlation:", round(spearman_rho, 4), "\n")

cat("Pearson correlation:", round(pearson_r, 4), "\n")

cat("\n==========================================\n")
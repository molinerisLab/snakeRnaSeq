#!/usr/bin/env Rscript

library(tidyverse)

rel_file <- snakemake@input[["relative"]]
meta_file <- snakemake@input[["metadata"]]
out_bar <- snakemake@output[["barplot"]]

meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  mutate(
    group = factor(group),
    condition = factor(condition),
    batch = factor(batch)
  )

# Clean taxon names by removing trailing taxonomy IDs and rank suffixes, replacing underscores with spaces, and trimming whitespace.:
clean_taxon_name <- function(x) {
  x %>%
    # Remove trailing taxonomy ID and rank suffix.
    str_remove("_[0-9]+_[A-Za-z]+$") %>%

    # Replace underscores with spaces.
    str_replace_all("_", " ") %>%

    # Remove repeated whitespace.
    str_squish()
}

load_matrix <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)

  if (!"Geneid" %in% colnames(df)) {
    stop(
      "Expected a 'Geneid' column. Found: ",
      paste(colnames(df), collapse = ", ")
    )
  }

  mat <- df %>%
    column_to_rownames("Geneid") %>%
    as.matrix()

  storage.mode(mat) <- "numeric"

  if (any(!is.finite(mat))) {
    stop("The abundance matrix contains NA, NaN, or infinite values.")
  }

  if (any(mat < 0)) {
    stop("The abundance matrix contains negative values.")
  }

  mat
}

rel_mat <- load_matrix(rel_file)

missing_metadata <- setdiff(colnames(rel_mat), meta$sample)

if (length(missing_metadata) > 0) {
  stop(
    "Samples missing from metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

sample_sums <- colSums(rel_mat)

cat("Relative-abundance sample sums:\n")
print(sample_sums)

if (any(sample_sums <= 0)) {
  stop("At least one sample has total abundance equal to zero.")
}

# Normalize each sample to sum to 1.
rel_mat <- sweep(rel_mat, 2, sample_sums, "/")

meta <- meta %>%
  filter(sample %in% colnames(rel_mat)) %>%
  arrange(match(sample, colnames(rel_mat)))

rel_long <- as.data.frame(
  rel_mat,
  check.names = FALSE
) %>%
  rownames_to_column("taxon_raw") %>%
  pivot_longer(
    cols = -taxon_raw,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  left_join(meta, by = "sample")

# Select the top taxa using the original identifiers.
top_taxa <- rel_long %>%
  group_by(taxon_raw) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  slice_max(
    order_by = mean_abundance,
    n = 20,
    with_ties = FALSE
  ) %>%
  pull(taxon_raw)

rel_long <- rel_long %>%
  mutate(
    taxon = if_else(
      taxon_raw %in% top_taxa,
      clean_taxon_name(taxon_raw),
      "Other"
    )
  ) %>%
  group_by(
    sample,
    group,
    condition,
    batch,
    taxon
  ) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample = factor(
      sample,
      levels = colnames(rel_mat)
    )
  )

# Order taxa by mean abundance and place "Other" last.
taxon_levels <- rel_long %>%
  filter(taxon != "Other") %>%
  group_by(taxon) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  pull(taxon)

rel_long <- rel_long %>%
  mutate(
    taxon = factor(
      taxon,
      levels = c(taxon_levels, "Other")
    )
  )

dir.create(
  dirname(out_bar),
  recursive = TRUE,
  showWarnings = FALSE
)

p_bar <- ggplot(
  rel_long,
  aes(
    x = sample,
    y = abundance,
    fill = taxon
  )
) +
  geom_col(width = 0.9) +
  facet_wrap(
    ~condition,           # You can adjust this based on your metadata. For example, group, condition, or batch.
    scales = "free_x"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Relative Abundance — Top 20 Taxa",
    fill = "Species",
    x = "Sample",
    y = "Relative abundance"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.text = element_text(face = "italic"),
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = out_bar,
  plot = p_bar,
  width = 12,
  height = 6
)

cat("Plot generated:", out_bar, "\n")
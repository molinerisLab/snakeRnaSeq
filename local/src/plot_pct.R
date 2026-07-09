library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------- CLI Args ----------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[1], "plots/kraken_composition.csv")
output_dir <- ifelse(length(args) >= 2, args[2], "plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

df <- suppressMessages(readr::read_csv(input_file, show_col_types = FALSE)) %>%
  mutate(
    kraken_input_reads = as.numeric(kraken_input_reads),
    kraken_unclassified_reads = as.numeric(kraken_unclassified_reads),
    kraken_classified_reads = as.numeric(kraken_classified_reads),
    kraken_human_9606_reads = as.numeric(kraken_human_9606_reads),

    # Safety clamps
    kraken_human_9606_reads = pmin(kraken_human_9606_reads, kraken_classified_reads),
    kraken_classified_nonhuman_reads = pmax(kraken_classified_reads - kraken_human_9606_reads, 0),

    denom = pmax(kraken_input_reads, 1),

    # The three percentages you requested
    pct_human_of_input = 100 * kraken_human_9606_reads / denom,
    pct_classified_of_input = 100 * kraken_classified_reads / denom,
    pct_unclassified_of_input = 100 * kraken_unclassified_reads / denom,

    # Disjoint parts for a valid 100% stack
    pct_classified_nonhuman_of_input = 100 * kraken_classified_nonhuman_reads / denom
  )

sample_order <- df %>% arrange(pct_human_of_input) %>% pull(name)

plot_df <- df %>%
  select(
    name,
    pct_unclassified_of_input,
    pct_classified_nonhuman_of_input,
    pct_human_of_input
  ) %>%
  pivot_longer(-name, names_to = "segment", values_to = "percent") %>%
  mutate(
    name = factor(name, levels = sample_order),
    segment = factor(
      segment,
      levels = c(
        "pct_unclassified_of_input",
        "pct_classified_nonhuman_of_input",
        "pct_human_of_input"
      ),
      labels = c(
        "Kraken unclassified",
        "Kraken classified non-human",
        "Kraken 9606 Homo sapiens"
      )
    )
  )

p <- ggplot(plot_df, aes(x = name, y = percent, fill = segment)) +
  geom_col(width = 0.9) +
  scale_y_continuous(
  limits  = c(0, 100),
  expand  = c(0, 0),
  oob     = scales::squish,          # <-- questa riga
  labels  = function(x) paste0(x, "%")
) +
  labs(
    title = "Read composition per sample",
    subtitle = "100% stacked unmapped reads; Kraken input is split into unclassified, classified non-human, and classified 9606",
    x = "Sample",
    y = "Percent of reads",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "stacked_unmapped_kraken_percent_with_9606.png"), p, width = 16, height = 6, dpi = 300)
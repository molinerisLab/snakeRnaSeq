#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
})

# ---------------------------
# Parse args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

n <- length(args)

if (n < 6) {
  stop("Need at least: <inputs...> <pdf> <png> <rds> <lfc> <padj>")
}

# Last five arguments are outputs + thresholds
pdf_file    <- args[n-4]
png_file    <- args[n-3]
rds_file    <- args[n-2]
lfc_cutoff  <- as.numeric(args[n-1])
padj_cutoff <- as.numeric(args[n])

# All previous args are input files
input_files <- args[1:(n-5)]

# ---------------------------
# Load all contrasts
# ---------------------------
load_contrast <- function(f) {
  df <- read.table(f, header = FALSE)
  colnames(df) <- c("gene", "log2FoldChange", "pvalue", "padj")

  # extract contrast name from filename
  # e.g. deseq2.toptable_clean.contrast_resistant_vs_unexposed.gz
  contrast <- sub(".*contrast_(.*)\\.gz$", "\\1", f)

  df$contrast <- contrast
  df
}

combined <- map_df(input_files, load_contrast)

# ---------------------------
# Prepare data
# ---------------------------
combined <- combined %>%
  mutate(
    negLog10Padj = -log10(padj),
    significance = case_when(
      padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# ---------------------------
# Plot (facet per contrast)
# ---------------------------
p <- ggplot(combined, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = significance), alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("gray70", "red")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  facet_wrap(~ contrast, scales = "free") +
  labs(
    title = "Combined Volcano Plots",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal(base_size = 14)

# ---------------------------
# Save outputs
# ---------------------------
ggsave(pdf_file, p, width = 12, height = 8)
ggsave(png_file, p, width = 12, height = 5.5, dpi = 300)
saveRDS(p, file = rds_file)


#===========================================================================================
# Modified version
#===========================================================================================

#!/usr/bin/env Rscript

# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(purrr)
#   library(readr)
# })

# # ---------------------------
# # Parse args
# # ---------------------------
# args <- commandArgs(trailingOnly = TRUE)
# n <- length(args)

# if (n < 6) {
#   stop("Usage: <inputs...> <pdf> <png> <rds> <lfc> <padj>")
# }

# pdf_file    <- args[n-4]
# png_file    <- args[n-3]
# rds_file    <- args[n-2]
# lfc_cutoff  <- as.numeric(args[n-1])
# padj_cutoff <- as.numeric(args[n])

# input_files <- args[1:(n-5)]

# # ---------------------------
# # Load contrasts
# # ---------------------------
# load_contrast <- function(f) {

#   df <- read_tsv(
#     f,
#     col_names = c("gene", "log2FoldChange", "pvalue", "padj"),
#     show_col_types = FALSE
#   )

#   contrast <- sub(".*contrast_(.*)\\.gz$", "\\1", f)

#   df %>%
#     mutate(
#       contrast = contrast,
#       padj = ifelse(is.na(padj), 1, padj),
#       negLog10Padj = -log10(padj),
#       significance = ifelse(
#         padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff,
#         "Significant",
#         "Not Significant"
#       )
#     )
# }

# combined <- map_df(input_files, load_contrast)

# # ---------------------------
# # Preserve contrast order
# # ---------------------------
# combined$contrast <- factor(
#   combined$contrast,
#   levels = unique(combined$contrast)
# )

# # ---------------------------
# # Plot
# # ---------------------------
# p <- ggplot(
#   combined,
#   aes(x = log2FoldChange, y = negLog10Padj)
# ) +
#   geom_point(
#     aes(color = significance),
#     alpha = 0.6,
#     size = 1.1
#   ) +
#   scale_color_manual(
#     values = c("Not Significant" = "gray70",
#                "Significant"     = "red")
#   ) +
#   geom_vline(
#     xintercept = c(-lfc_cutoff, lfc_cutoff),
#     linetype = "dashed",
#     linewidth = 0.4
#   ) +
#   geom_hline(
#     yintercept = -log10(padj_cutoff),
#     linetype = "dashed",
#     linewidth = 0.4
#   ) +
#   facet_wrap(~ contrast, nrow = 2) +
#   labs(
#     x = expression(log[2]~Fold~Change),
#     y = expression(-log[10]~FDR)
#   ) +
#   theme_classic(base_size = 13) +
#   theme(
#     legend.position = "top",
#     strip.background = element_blank(),
#     strip.text = element_text(face = "bold")
#   )

# # ---------------------------
# # Save outputs
# # ---------------------------
# ggsave(pdf_file, p, width = 14, height = 9)
# ggsave(png_file, p, width = 14, height = 9, dpi = 300)
# saveRDS(p, file = rds_file)


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)

top_table_file <- args[1]
out_volcano <- args[2]
out_ma <- args[3]
out_pvalue_hist <- args[4]

read_table_auto <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    read_xlsx(path)
  } else {
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  }
}

df <- read_table_auto(top_table_file)

if ("Geneid" %in% colnames(df) && !"GeneID" %in% colnames(df)) {
  colnames(df)[colnames(df) == "Geneid"] <- "GeneID"
}

required <- c("GeneID", "logFC", "Pvalue", "Pvalue_adj")
missing <- setdiff(required, colnames(df))

if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

df <- df %>%
  mutate(
    logFC = as.numeric(logFC),
    Pvalue = as.numeric(Pvalue),
    Pvalue_adj = as.numeric(Pvalue_adj),
    neg_log10_p = -log10(Pvalue),
    status = case_when(
      Pvalue_adj < 0.05 ~ "FDR < 0.05",
      Pvalue < 0.05 ~ "nominal P < 0.05",
      TRUE ~ "not nominal"
    )
  ) %>%
  arrange(Pvalue)

df$GeneID <- gsub("_\\d+_S$", "", df$GeneID)
df$GeneID <- gsub("_", " ", df$GeneID)

df$label <- NA_character_
df$label[seq_len(min(10, nrow(df)))] <- df$GeneID[seq_len(min(10, nrow(df)))]


volcano <- ggplot(df, aes(logFC, neg_log10_p)) +
  geom_point(aes(color = status), alpha = 0.8, size = 2.2) +
  scale_color_manual(values = c("not nominal" = "grey70", "nominal P < 0.05" = "#E69F00", "FDR < 0.05" = "#D55E00")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  geom_text_repel(aes(label = label), max.overlaps = 20, size = 2.5) +
  labs(
    title = "Resistant vs Control: Volcano Plot",
    x = "logFC: Resistant vs Control",
    y = "-log10(P value)",
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(out_volcano, volcano, width = 8, height = 6)

if (all(c("Control", "Resistant") %in% colnames(df))) {
  df <- df %>%
    mutate(
      Control = as.numeric(Control),
      Resistant = as.numeric(Resistant),
      mean_abundance = rowMeans(cbind(Control, Resistant), na.rm = TRUE)
    )

  ma <- ggplot(df, aes(rank_mean, logFC)) +
    geom_point(aes(color = status), alpha = 0.8, size = 2.2) +
    scale_color_manual(values = c("not nominal" = "grey70", "nominal P < 0.05" = "#56B4E9", "FDR < 0.05" = "#0072B2")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = "Resistant vs Control: Effect-size Rank Plot",
      x = "Taxa ranked by top-table order",
      y = "logFC: Resistant vs Control",
      color = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

} else {
  df$rank_mean <- seq_len(nrow(df))

  ma <- ggplot(df, aes(rank_mean, logFC)) +
    geom_point(aes(color = status), alpha = 0.8, size = 2.2) +
    scale_color_manual(values = c("not nominal" = "grey70", "nominal P < 0.05" = "#56B4E9", "FDR < 0.05" = "#0072B2")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = "Resistant vs Control: Effect-size Rank Plot",
      x = "Taxa ranked by top-table order",
      y = "logFC: Resistant vs Control",
      color = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

ggsave(out_ma, ma, width = 8, height = 6)

p_hist <- ggplot(df, aes(Pvalue)) +
  geom_histogram(bins = 30, boundary = 0, fill = "#4DBBD5", color = "black", alpha = 0.85) +
  labs(
    title = "P-value Distribution",
    x = "P value",
    y = "Number of taxa"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(out_pvalue_hist, p_hist, width = 6, height = 5)
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Usage: pca_plot.R <counts_file> <metadata_file> <png_out> <pdf_out>")
}

counts_file  <- args[1]
metadata_file <- args[2]
png_out <- args[3]
pdf_out <- args[4]

# Load metadata
sample_info <- read.table(metadata_file, header = TRUE, sep="\t")
sample_info$condition <- factor(sample_info$condition,
                                levels = c("susceptible", "resistant", "unexposed"))
rownames(sample_info) <- sample_info$sample

# Load count matrix
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep="\t", check.names = FALSE)

# Ensure sample ordering matches metadata
counts <- counts[, rownames(sample_info)]

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info,
    design    = ~ condition
)

# Variance Stabilizing Transformation
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA Plot
p <- ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of RNA-seq Gene Expression Profiles") +
    scale_color_manual(
        name = " ",
        breaks = c("resistant", "susceptible", "unexposed"),
        labels = c("Resistant", "Susceptible", "Unexposed"),
        values = c("resistant" = "#E41A1C",
                   "susceptible" = "#377EB8",
                   "unexposed" = "#4DAF4A")
    ) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "bottom",
        legend.title = element_blank()
    )

# Save outputs
ggsave(png_out, p, width = 6, height = 4.5, dpi = 300)
ggsave(pdf_out, p, width = 6, height = 4.5)


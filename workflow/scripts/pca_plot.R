#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(yaml)
    library(dplyr)
    library(RColorBrewer)
})

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_flag_value <- function(flag) {
    idx <- match(flag, args)
    if (is.na(idx) || idx == length(args)) {
        return(NULL)
    }
    return(args[[idx + 1]])
}

counts_file  <- get_flag_value("--counts")
metadata_file <- get_flag_value("--metadata")
config_file  <- get_flag_value("--config")
png_out      <- get_flag_value("--png")
pdf_out      <- get_flag_value("--pdf")

if (is.null(counts_file) ||
    is.null(metadata_file) ||
    is.null(config_file) ||
    is.null(png_out) ||
    is.null(pdf_out)) {

    stop("Usage: pca_plot1.R --counts <counts_file> --metadata <metadata_file> --config <config_yaml> --png <png_out> --pdf <pdf_out>")
}

# -----------------------------
# Load config
# -----------------------------
config <- yaml::read_yaml(config_file)
pc     <- config$PCA

color_var  <- pc$COLOR_BY
shape_var  <- pc$SHAPE_BY
ntop       <- ifelse(is.null(pc$NTOP), 500, pc$NTOP)
add_labels <- ifelse(is.null(pc$LABEL_POINTS), FALSE, pc$LABEL_POINTS)

# -----------------------------
# Load metadata
# -----------------------------
sample_info <- read.table(metadata_file,
                          header = TRUE,
                          sep = "\t",
                          check.names = FALSE)

rownames(sample_info) <- sample_info$sample

# Ensure grouping variables exist
if (!color_var %in% colnames(sample_info)) {
    stop(paste("COLOR_BY variable not found in metadata:", color_var))
}

if (!is.null(shape_var) && !is.na(shape_var) &&
    shape_var != "null" &&
    !shape_var %in% colnames(sample_info)) {
    stop(paste("SHAPE_BY variable not found in metadata:", shape_var))
}

# Convert grouping columns to factors
sample_info[[color_var]] <- as.factor(sample_info[[color_var]])

if (!is.null(shape_var) && !is.na(shape_var) && shape_var != "null") {
    sample_info[[shape_var]] <- as.factor(sample_info[[shape_var]])
} else {
    shape_var <- NULL
}

# -----------------------------
# Load count matrix
# -----------------------------
counts <- read.table(counts_file,
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)

counts <- as.matrix(counts)

# Ensure counts are integers and non-negative
if (any(counts < 0)) {
    stop("Count matrix contains negative values. Use raw featureCounts output.")
}

mode(counts) <- "integer"

# Match sample order
counts <- counts[, rownames(sample_info)]

# -----------------------------
# Create DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info,
    design    = as.formula(paste("~", color_var))
)

# Variance stabilizing transform
vsd <- vst(dds, blind = TRUE)

# PCA using ntop most variable genes
pcaData <- plotPCA(
    vsd,
    intgroup = color_var,
    ntop = ntop,
    returnData = TRUE
)

percentVar <- round(100 * attr(pcaData, "percentVar"))

# If shape factor exists, merge it
if (!is.null(shape_var)) {
    pcaData[[shape_var]] <- colData(vsd)[[shape_var]]
}

# -----------------------------
# Build aesthetics (modern tidy eval)
# -----------------------------
if (is.null(shape_var)) {
    aes_mapping <- aes(
        x = PC1,
        y = PC2,
        color = .data[[color_var]]
    )
} else {
    aes_mapping <- aes(
        x = PC1,
        y = PC2,
        color = .data[[color_var]],
        shape = .data[[shape_var]]
    )
}

# -----------------------------
# Base PCA plot
# -----------------------------
p <- ggplot(pcaData, aes_mapping) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(pc$TITLE) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "bottom",
        legend.title = element_blank()
    )

# -----------------------------
# Optional ellipses (only if >=3 samples per group)
# -----------------------------
group_counts <- pcaData %>%
    group_by(.data[[color_var]]) %>%
    tally()

if (all(group_counts$n >= 3)) {
    p <- p +
        stat_ellipse(
            aes(group = .data[[color_var]]),
            linetype = 2,
            linewidth = 0.8,
            show.legend = FALSE
        )
}

# -----------------------------
# Optional sample labels
# -----------------------------
if (add_labels) {
    p <- p +
        ggrepel::geom_text_repel(
            aes(label = name),
            size = 3,
            show.legend = FALSE
        )
}

# -----------------------------
# Custom legend labels (if provided)
# -----------------------------
if (!is.null(pc$LABELS) &&
    !is.null(pc$LABELS[[color_var]])) {

    label_map <- unlist(pc$LABELS[[color_var]])

    p <- p +
        scale_color_discrete(labels = label_map)
}

# -----------------------------
# Save outputs
# -----------------------------
ggsave(png_out, p, width = 6, height = 4.5, dpi = 300)
ggsave(pdf_out, p, width = 6, height = 4.5)

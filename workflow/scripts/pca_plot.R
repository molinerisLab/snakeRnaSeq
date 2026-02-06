#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(yaml)
})

# -----------------------------
# CLI parsing
# -----------------------------
get_flag_value <- function(flag, default = NA_character_) {
    args <- commandArgs(trailingOnly = TRUE)
    idx <- match(flag, args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[[idx + 1]]
}

counts_file   <- get_flag_value("--counts")
metadata_file <- get_flag_value("--metadata")
config_file   <- get_flag_value("--config")
png_out       <- get_flag_value("--png")
pdf_out       <- get_flag_value("--pdf")

if (any(is.na(c(counts_file, metadata_file, config_file, png_out, pdf_out)))) {
    stop("Usage: pca_plot.R --counts <file> --metadata <file> --config <config.yaml> --png <out.png> --pdf <out.pdf>")
}

# -----------------------------
# Load config
# -----------------------------
cfg <- yaml::read_yaml(config_file)
if (is.null(cfg$PCA)) stop("Missing PCA block in config.yaml")
pc <- cfg$PCA

# -----------------------------
# Load metadata
# -----------------------------
sample_info <- read.table(
    metadata_file,
    header = TRUE,
    sep = pc$METADATA_SEP,
    check.names = FALSE
)

rownames(sample_info) <- sample_info[[pc$SAMPLE_COLUMN]]

# Apply factor levels if provided
if (!is.null(pc$LEVELS)) {
    for (var in names(pc$LEVELS)) {
        if (var %in% colnames(sample_info)) {
            sample_info[[var]] <- factor(
                sample_info[[var]],
                levels = pc$LEVELS[[var]]
            )
        }
    }
}

# -----------------------------
# Load counts
# -----------------------------
counts <- read.table(
    counts_file,
    header = TRUE,
    row.names = 1,
    sep = pc$COUNTS_SEP,
    check.names = FALSE
)

counts <- counts[, rownames(sample_info)]

# -----------------------------
# DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = as.formula(pc$DESIGN)
)

vsd <- vst(dds, blind = TRUE)
color_var <- pc$COLOR_BY
shape_var <- pc$SHAPE_BY

intgroup_vars <- color_var
if (!is.null(shape_var) && shape_var != "null") {
    intgroup_vars <- c(color_var, shape_var)
}

pcaData <- plotPCA(vsd, intgroup = intgroup_vars, returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

# -----------------------------
# Base ggplot
# -----------------------------
aes_mapping <- aes_string(
    x = "PC1",
    y = "PC2",
    color = color_var,
    label = "name"
)

if (!is.null(shape_var) && shape_var != "null") {
    aes_mapping <- aes_string(
        x = "PC1",
        y = "PC2",
        color = color_var,
        shape = shape_var,
        label = "name"
    )
}

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
# Ellipses (optional)
# -----------------------------
if (!is.null(pc$ADD_ELLIPSE) && pc$ADD_ELLIPSE) {
    p <- p +
        stat_ellipse(
            aes_string(group = color_var),
            type = pc$ELLIPSE_TYPE,
            linetype = 2
        )
}

# -----------------------------
# Apply manual colors
# -----------------------------
if (!is.null(pc$COLORS)) {
    p <- p +
        scale_color_manual(
            values = unlist(pc$COLORS),
            drop = FALSE
        )
}

# -----------------------------
# Apply manual shapes (if used)
# -----------------------------
if (!is.null(shape_var) && shape_var != "null" && !is.null(pc$SHAPES)) {
    p <- p +
        scale_shape_manual(
            values = unlist(pc$SHAPES),
            drop = FALSE
        )
}

# -----------------------------
# Save outputs
# -----------------------------
ggsave(
    png_out,
    p,
    width = pc$WIDTH,
    height = pc$HEIGHT,
    dpi = pc$DPI
)

ggsave(
    pdf_out,
    p,
    width = pc$WIDTH,
    height = pc$HEIGHT
)

#!/usr/bin/env R स्क्रिप्ट

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript extract_genes_from_heatmap_rds.R <input.rds> <output_genes.txt>")
}

rds_file <- args[1]
out_file <- args[2]

cat("Loading RDS file:", rds_file, "\n")
p <- readRDS(rds_file)

if (!is.null(p$tree_row)) {
    genes <- p$tree_row$labels[p$tree_row$order]
    cat("Found", length(genes), "clustered genes.\n")
} else {
    if (!is.null(p$tree_row$labels)) {
        genes <- p$tree_row$labels
        cat("Found", length(genes), "genes (not clustered).\n")
    } else {
        # Fallback if tree_row$labels isn't populated but rownames are in the data
        if (!is.null(rownames(p$gtable$grobs[[1]]$label))) {
           genes <- rownames(p$gtable$grobs[[1]]$label)
        } else {
            stop("Could not extract gene labels from the pheatmap object.")
        }
    }
}

write.table(genes, out_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("Gene list saved to:", out_file, "\n")

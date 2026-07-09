#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: extract_fasta_from_rds.R <input.rds> [output.fasta] [transcriptome.fa]")
}

rds_file <- args[1]
out_fasta <- ifelse(length(args) >= 2, args[2], "extracted_transcripts.fasta")
ref_fasta <- ifelse(length(args) >= 3, args[3], "kallisto_output/combined_transcriptome.fa")

if (!file.exists(rds_file)) stop(paste("File not found:", rds_file))
if (!file.exists(ref_fasta)) stop(paste("Reference FASTA not found:", ref_fasta))

# Read RDS
cat("Reading", rds_file, "...\n")
x <- readRDS(rds_file)

# Extract full names based on object type
if (inherits(x, "pheatmap")) {
  full_names <- x$tree_row$labels[x$tree_row$order]
} else if (is.matrix(x) || is.data.frame(x)) {
  full_names <- rownames(x)
} else {
  stop("Unknown RDS object format.")
}

# Keep the full transcript names since samtools faidx index requires exact matches
ids <- full_names

# Write to temp file
tmp_list <- tempfile("transcripts_", fileext = ".txt")
writeLines(ids, tmp_list)

# Remove the output file if it exists to append fresh
if (file.exists(out_fasta)) file.remove(out_fasta)

# Run samtools faidx
cat(sprintf("Extracting %d sequences from %s to %s...\n", length(ids), ref_fasta, out_fasta))
cmd <- paste0("while read -r id; do samtools faidx ", shQuote(ref_fasta), " \"$id\" >> ", shQuote(out_fasta), "; done < ", shQuote(tmp_list))
system(cmd)

unlink(tmp_list)
cat("Done!\n")

#!/usr/bin/env Rscript

library(rhdf5)
library(tximport)
library(data.table)

cat("Finding kallisto files...\n")
base_dir <- "kallisto_PsiCLASS_idx_combined"
files <- list.files(base_dir, pattern="abundance.h5", full.names=TRUE, recursive=TRUE)

names(files) <- basename(dirname(files))

if (length(files) == 0) {
  stop("No abundance.h5 files found in ", base_dir)
}
cat(sprintf("Found %d abundance.h5 files.\n", length(files)))

cat("Extracting transcript IDs and building tx2gene mapping...\n")
ids <- as.character(h5read(files[1], "aux/ids"))

# Vectorized splitting using data.table
split_dt <- as.data.table(tstrsplit(ids, "|", fixed=TRUE))

col1 <- as.character(split_dt[[1]])
col2 <- if (ncol(split_dt) >= 2) as.character(split_dt[[2]]) else rep(NA_character_, nrow(split_dt))
col6 <- if (ncol(split_dt) >= 6) as.character(split_dt[[6]]) else rep(NA_character_, nrow(split_dt))

gene_ids <- fifelse(!is.na(col2) & col2 != "-", col2, col1)
symbol_ids <- fifelse(!is.na(col6) & col6 != "-", col6, gene_ids)

rm(split_dt, col1, col2, col6)

tx2gene <- data.frame(TXNAME = ids, GENEID = gene_ids)

# Strict validation
stopifnot(sum(is.na(tx2gene$GENEID)) == 0)
stopifnot(sum(duplicated(tx2gene$TXNAME)) == 0)
stopifnot(sum(tx2gene$GENEID == "-") == 0)
stopifnot(sum(tx2gene$GENEID == "") == 0)
cat(sprintf("tx2gene built: %d transcripts -> %d genes\n",
            nrow(tx2gene), length(unique(tx2gene$GENEID))))

cat("Running tximport at gene level...\n")
# "no" = raw estimated counts; correct for DESeq2, wrong for edgeR/voom without length offsets
txi <- tximport(files, type="kallisto", txOut=FALSE, tx2gene=tx2gene,
                countsFromAbundance="no")  

cat("Writing outputs to tximport_genes_counts.tsv and tximport_genes_tpm.tsv...\n")
out_counts <- data.table(Geneid = rownames(txi$counts))
out_counts <- cbind(out_counts, txi$counts)
fwrite(out_counts, "tximport_genes_counts.tsv", sep="\t")

out_tpm <- data.table(Geneid = rownames(txi$abundance))
out_tpm <- cbind(out_tpm, txi$abundance)
fwrite(out_tpm, "tximport_genes_tpm.tsv", sep="\t")

# Free massive memory footprint
rm(txi, out_counts, out_tpm)

# STRICT 1-to-1 annotation map
anno_dt <- data.table(Geneid = gene_ids, Symbol = symbol_ids)

conflicts <- anno_dt[, .N, by = Geneid][N > 1]
if (nrow(conflicts) > 0) {
  message(nrow(conflicts), " Geneids had conflicting symbols \u2014 kept first match deterministically")
}

# Sort by Symbol to make deduplication stable and deterministic across runs
setorder(anno_dt, Geneid, Symbol)
gene_annotation_map <- anno_dt[, .(Symbol = Symbol[1]), by = Geneid]

fwrite(gene_annotation_map, "gene_id_to_symbol.tsv", sep="\t")
cat("Annotation map written to gene_id_to_symbol.tsv\n")

cat("Done.\n")
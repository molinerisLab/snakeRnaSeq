library(data.table)

# 0. Input validation
stopifnot(file.exists("transcripts_counts.tsv"))

# 1. Load the transcript matrix
tx_data <- fread("transcripts_counts.tsv")

# Explicitly pull the first column as a 1D character array
tx_names <- tx_data[[1]]

# Drop the first column and convert explicitly to a numeric matrix
tx_data[, names(tx_data)[1] := NULL]
tx_mat <- as.matrix(tx_data)
storage.mode(tx_mat) <- "numeric"
rm(tx_data)

## Use as.data.table() to correctly unpack the list into N rows x K columns
split_dt <- as.data.table(tstrsplit(tx_names, "|", fixed=TRUE))

# Safely extract columns
col1 <- as.character(split_dt[[1]])
col2 <- if (ncol(split_dt) >= 2) as.character(split_dt[[2]]) else rep(NA_character_, nrow(split_dt))
col6 <- if (ncol(split_dt) >= 6) as.character(split_dt[[6]]) else rep(NA_character_, nrow(split_dt))

# Group by the ENSG ID (2nd element) with fallback for unannotated loci (1st element)
gene_key <- fifelse(!is.na(col2) & col2 != "-", col2, col1)

# Extract the Gene Symbol (6th element) with fallback to gene_key
symbol_key <- fifelse(!is.na(col6) & col6 != "-", col6, gene_key)
# Clean up split_dt to free more memory
rm(split_dt)

# 3. Sum transcript expression by the gene keys
stopifnot(
  "row mismatch between tx_mat and gene_key" =
    nrow(tx_mat) == length(gene_key)
)
gene_mat <- rowsum(tx_mat, group=gene_key)
rm(tx_mat)

# 4. Format and write the primary count output
out_dt <- data.table(Geneid = rownames(gene_mat))
out_dt <- cbind(out_dt, gene_mat)
fwrite(out_dt, "updated_genes_counts.tsv", sep="\t")

# 5. Create a STRICT 1-to-1 annotation map
anno_dt <- data.table(Geneid = gene_key, Symbol = symbol_key)

# Check for and warn about conflicts before deduplicating
conflicts <- anno_dt[, .N, by = Geneid][N > 1]
if (nrow(conflicts) > 0) {
  message(nrow(conflicts), " Geneids had conflicting symbols \u2014 kept first match")
}

# Deduplicate: if there are conflicts, just take the first symbol
gene_annotation_map <- anno_dt[, .(Symbol = Symbol[1]), by = Geneid]

fwrite(gene_annotation_map, "gene_id_to_symbol.tsv", sep="\t")
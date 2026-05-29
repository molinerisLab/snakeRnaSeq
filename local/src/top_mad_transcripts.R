library(data.table)
library(matrixStats)
cat("Loading TMM normalized data...\n")
genes <- fread("tximport_genes_counts.tmm.tsv")
tx <- fread("transcripts_counts.tmm.tsv")
setnames(genes, 1, "Gene_ID")
setnames(tx, 1, "Transcript_ID")
cat("Mapping Transcripts to Genes...\n")
split_names <- strsplit(tx$Transcript_ID, "\\|")
tx[, Gene_ID := sapply(split_names, function(x) {
  if (length(x) >= 2) return(x[2]) # Extracts Ensembl Gene ID
  else return(x[1])
})]
cat("Calculating log2(TMM + 1) and MAD...\n")
cat("Filtering lowly expressed genes and transcripts...\n")
sample_cols <- setdiff(names(tx), c("Transcript_ID", "Gene_ID"))
genes_mat <- as.matrix(genes[, ..sample_cols])
tx_mat <- as.matrix(tx[, ..sample_cols])
# Keep genes/transcripts with an average TMM > 1 across samples
genes_keep <- rowMeans(genes_mat) > 1
tx_keep <- rowMeans(tx_mat) > 1
genes_mat <- genes_mat[genes_keep, ]
tx_mat <- tx_mat[tx_keep, ]
# Subset the original data.tables to match
genes <- genes[genes_keep, ]
tx <- tx[tx_keep, ]
cat("Calculating log2(TMM + 1) and MAD...\n")
genes_log <- log2(genes_mat + 1)
tx_log <- log2(tx_mat + 1)
gene_mad_values <- rowMads(genes_log)
tx_mad_values <- rowMads(tx_log)
genes_mad_dt <- data.table(Gene_ID = genes$Gene_ID, Gene_MAD = gene_mad_values)
tx_mad_dt <- data.table(Transcript_ID = tx$Transcript_ID, Gene_ID = tx$Gene_ID, Transcript_MAD = tx_mad_values)
merged_dt <- merge(tx_mad_dt, genes_mad_dt, by = "Gene_ID", all.x = TRUE)
merged_dt <- merged_dt[Transcript_MAD > 0]
cat("\n=========================================\n")
cat("TOP 10 TRANSCRIPTS WITH HIGHEST MAD:\n")
cat("=========================================\n")
top10_high_tx <- merged_dt[order(-Transcript_MAD)][1:10]
print(top10_high_tx)
fwrite(top10_high_tx, "plots/top10_highest_transcript_MAD.csv")
cat("\n=========================================\n")
cat("TOP 10 HIGH VAR TRANSCRIPTS + LOW VAR GENES:\n")
cat("=========================================\n")
# To find high var transcript + low var gene, we can look for the largest difference
# between Transcript_MAD and Gene_MAD, or the highest ratio.
# Since Transcript_MAD can be higher than Gene_MAD if the gene has multiple transcripts
# that are anti-correlated (so the total gene sum is stable but transcripts switch).
# This is classic "Isoform Switching".
merged_dt[, MAD_diff := Transcript_MAD - Gene_MAD]
top10_switchers <- merged_dt[order(-MAD_diff)][1:1000]
print(top10_switchers)
fwrite(top10_switchers, "plots/top1000_isoform_switching.csv")
fwrite(merged_dt[order(-MAD_diff)], "plots/all_MAD_diff.csv")
cat("\nData saved to plots/top10_highest_transcript_MAD.csv and plots/top10_isoform_switching.csv\n")
cat("\n=========================================\n")
cat("EXPORTING TOP N VARIABLE GENES FOR ENRICHR:\n")
cat("=========================================\n")
# 1. Define N and extract the top variable genes by MAD
N_genes <- 500
top_genes_dt <- genes_mad_dt[order(-Gene_MAD)][1:N_genes]
# 2. Load the dictionary we made earlier
anno_map <- fread("gene_id_to_symbol.tsv")
# 3. Merge to get the official Gene Symbols
top_genes_annotated <- merge(top_genes_dt, anno_map, by.x="Gene_ID", by.y="Geneid", all.x=TRUE)
# 4. Filter for EnrichR
# We only want genes that actually have a recognized symbol. 
enrichr_targets <- top_genes_annotated[Symbol != Gene_ID & !is.na(Symbol)]$Symbol
cat(sprintf("Extracted %d valid Gene Symbols out of the top %d variable genes.\n", 
            length(enrichr_targets), N_genes))
# 5. Export as a clean text file for easy copy-pasting
write.table(enrichr_targets, "plots/EnrichR_input_list.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("List exported to plots/EnrichR_input_list.txt. Ready for the web interface!\n")
# Extract all valid gene symbols from your merged dataset (where MAD was calculated)
# We filter out the unannotated loci just like we did for the targets
background_targets <- anno_map[Symbol != Geneid & !is.na(Symbol)]$Symbol
write.table(unique(background_targets), "plots/EnrichR_background_list.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("Background list exported: %d genes.\n", length(unique(background_targets))))
cat("\n=========================================\n")
cat("PREPARING .rnk FILE FOR SnakeGSEA:\n")
cat("=========================================\n")
# 1. Merge your MAD data with the Gene Symbol dictionary
# We will use the merged_dt from your earlier script
gsea_data <- merge(merged_dt, anno_map, by.x="Gene_ID", by.y="Geneid", all.x=TRUE)
# 2. Filter out unannotated genes (GSEA requires official symbols)
gsea_data <- gsea_data[Symbol != Gene_ID & !is.na(Symbol)]
# 3. Create the two-column format required by SnakeGSEA
# We are using MAD_diff to look for pathways enriched for isoform switching
rnk_df <- gsea_data[, .(Symbol, MAD_diff)]
# 4. Sort by the ranking metric in descending order (CRITICAL for GSEA)
# 4. Collapse duplicate Gene Symbols by taking the maximum MAD_diff
# This is CRITICAL for GSEA, which expects exactly one value per gene symbol
rnk_df <- rnk_df[, .(MAD_diff = max(MAD_diff, na.rm = TRUE)), by = Symbol]
# 5. Sort by the ranking metric in descending order
rnk_df <- rnk_df[order(-MAD_diff)]
# 5. Export as a tab-separated file WITHOUT column names
# 6. Export as a tab-separated file WITHOUT column names
# SnakeGSEA will pick up any .rnk file in the directory
fwrite(rnk_df, "isoform_switching.rnk", sep="\t", col.names=FALSE)
cat("Exported isoform_switching.rnk successfully.\n")
cat("Exported isoform_switching.rnk successfully (Duplicates collapsed).\n")
cat("Ready to execute: snakemake -p -j N_CORES all\n")

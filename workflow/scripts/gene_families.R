#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

read_any <- function(path, header = TRUE) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    fread(cmd = sprintf("zcat %s", shQuote(path)),
          header = header, data.table = FALSE, sep = "\t")
  } else {
    fread(path, header = header, data.table = FALSE, sep = "\t")
  }
}

extract_attr <- function(x, key) {
  pat <- paste0("(^|;)", key, "=([^;]+)")
  out <- sub(paste0(".", pat, "."), "\\2", x, perl = TRUE)
  out[!grepl(pat, x, perl = TRUE)] <- NA_character_
  out
}

first_non_na <- function(...) {
  xs <- list(...)
  out <- xs[[1]]
  for (i in 2:length(xs)) out[is.na(out)] <- xs[[i]][is.na(out)]
  out
}

extract_agap <- function(x) {
  m <- regmatches(x, regexpr("AGAP[0-9]+", x, perl = TRUE))
  m[m == ""] <- NA_character_
  m
}

# Clean transcript-style names like CYP6Z3-RA -> CYP6Z3
strip_isoform <- function(x) {
  x <- gsub("-R[A-Z]$", "", x)
  x
}

classify_family <- function(symbol, desc) {
  s <- tolower(ifelse(is.na(symbol), "", symbol))
  d <- tolower(ifelse(is.na(desc), "", desc))

  # Detox
  if (grepl("^cyp", s) || grepl("cytochrome p450|\\bp450\\b", d)) return("CYP")
  if (grepl("^gst", s) || grepl("glutathione.*transferase|glutathione s-transferase", d)) return("GST")
  # COE / esterase family (also catch "CCE" style)
  if (grepl("^coe", s) || grepl("^cce", s) || grepl("carboxylesterase|cholinesterase|juvenile hormone esterase|\\besterase\\b", d)) return("COE")

  # Cuticular proteins
  if (grepl("^cpr", s) || grepl("^cpap", s) || grepl("^cplcp", s) || grepl("^cpcfc", s) ||
      grepl("^twdl", s) || grepl("cuticular protein|\\bcuticle\\b", d)) return("CP")

  # Salivary gland proteins (SG*, D7*)
  if (grepl("^sg", s) || grepl("^d7", s) || grepl("\\bsalivary\\b|\\bd7\\b", d)) return("SGP")

  "Other"
}

option_list <- list(
  make_option(c("--all_contrasts"), type="character", help="ALL_contrast header-added .gz", metavar="FILE"),
  make_option(c("--gff"), type="character", help="VectorBase GFF (.gff or .gff.gz)", metavar="FILE"),
  make_option(c("--out"), type="character", help="Output TSV: GeneID<tab>family", metavar="FILE"),
  make_option(c("--gene_col"), type="character", default="GeneID", help="Gene column in ALL-contrast file")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$all_contrasts) || is.null(opt$gff) || is.null(opt$out)) {
  stop("Required: --all_contrasts, --gff, --out")
}

# Gene set from your ALL-contrast DE file
deg <- read_any(opt$all_contrasts, header = TRUE)
if (!(opt$gene_col %in% colnames(deg))) stop("Gene column not found: ", opt$gene_col)

gene_ids <- unique(deg[[opt$gene_col]])
gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]

# Read GFF
gff <- read_any(opt$gff, header = FALSE)
if (ncol(gff) < 9) stop("GFF appears to have <9 tab-separated columns: ", opt$gff)
gff <- gff[, 1:9]
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
gff <- gff[!grepl("^#", gff$seqid), , drop = FALSE]

# --- Gene-like rows in this VectorBase file ---
genes <- gff[gff$type %in% c("protein_coding_gene", "ncRNA_gene"), , drop = FALSE]

# Extract GeneID and descriptions from gene rows
g_attrs <- genes$attributes
g_id_raw <- extract_attr(g_attrs, "ID")
g_desc <- first_non_na(extract_attr(g_attrs, "description"),
                       extract_attr(g_attrs, "product"),
                       extract_attr(g_attrs, "Note"))

GeneID <- extract_agap(g_id_raw)

gene_tbl <- data.frame(GeneID = GeneID, gene_desc = g_desc, stringsAsFactors = FALSE)
gene_tbl <- gene_tbl[!is.na(gene_tbl$GeneID), , drop = FALSE]
gene_tbl <- gene_tbl[!duplicated(gene_tbl$GeneID), , drop = FALSE]

# --- mRNA rows may carry Name/Alias/product fields (optional enrichment) ---
tx <- gff[gff$type == "mRNA", , drop = FALSE]
t_attrs <- tx$attributes
t_parent <- extract_attr(t_attrs, "Parent")
t_name   <- extract_attr(t_attrs, "Name")
t_alias  <- extract_attr(t_attrs, "Alias")
t_desc   <- first_non_na(extract_attr(t_attrs, "description"),
                         extract_attr(t_attrs, "product"),
                         extract_attr(t_attrs, "Note"))

tx_gene <- extract_agap(t_parent)
tx_symbol <- first_non_na(t_name, t_alias)
tx_symbol <- strip_isoform(tx_symbol)

tx_tbl <- data.frame(GeneID = tx_gene, symbol = tx_symbol, tx_desc = t_desc, stringsAsFactors = FALSE)
tx_tbl <- tx_tbl[!is.na(tx_tbl$GeneID), , drop = FALSE]

# Collapse multiple transcripts per gene: take first non-NA symbol/desc
# (data.table makes this easy)
dt <- as.data.table(tx_tbl)
tx_collapsed <- dt[, .(
  symbol = symbol[which(!is.na(symbol) & symbol != "")][1],
  tx_desc = tx_desc[which(!is.na(tx_desc) & tx_desc != "")][1]
), by = GeneID]
tx_collapsed <- as.data.frame(tx_collapsed)

# Join gene + transcript info
ann <- merge(gene_tbl, tx_collapsed, by = "GeneID", all.x = TRUE)

# Combine descriptions (prefer informative transcript desc if gene desc is "unspecified product")
ann$desc <- ann$gene_desc
is_unspec <- !is.na(ann$desc) & tolower(ann$desc) == "unspecified product"
ann$desc[is_unspec] <- ann$tx_desc[is_unspec]
ann$desc <- first_non_na(ann$desc, ann$tx_desc)

# Keep only genes that appear in DE
ann <- ann[ann$GeneID %in% gene_ids, , drop = FALSE]

missing_n <- length(setdiff(gene_ids, ann$GeneID))
message("Genes in DE results: ", length(gene_ids))
message("Genes matched in GFF: ", nrow(ann))
message("Genes missing from GFF match: ", missing_n)

# Classify
ann$family <- mapply(classify_family, ann$symbol, ann$desc, USE.NAMES = FALSE)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
write.table(ann[, c("GeneID","family")], file = opt$out, sep = "\t", row.names = FALSE, quote = FALSE)

print(sort(table(ann$family), decreasing = TRUE))

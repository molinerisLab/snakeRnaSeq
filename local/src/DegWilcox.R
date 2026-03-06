#!/usr/bin/env Rscript
sink(stderr())#the stdout has to be kept clean form messages

'
Usage:
  DEGwilcox <gene_expression_matrix> <sample_metadata> [--condition=<condition>] [--g1=<g1>] [--g2=<g2>] [--min_exp=<min_exp>] [--min_samples_ratio=<min_sampes_ratio>]

The file <gene_expression_matrix> shoud have gene names on the first column and samples as column names.
The expression data are expected as raw (not normalized) data.

The file <sample_metadata> shoud have gene sample names on the first column.

Options:
  --condition <condition>        Name of dichotomous variable column in the metadata file to discriminate the two group of samples [default: condition].
  --g1 <g1>                      Level 1 of the dichotomous variable, indentifing group1_samples samples [default: case].
  --g2 <g2>                      Level 2 of the dichotomous variable, indentifing group2_samples samples [default: control].
  --min_exp <min_exp>            Minimum expression levels (cpm) to consider a gene as expressed [default: 1].
  --min_samples_ratio <min_sampes_ratio>            If less that <min_exp_sampes_ratio> (rounded up) samples show a gene expressed in at least one group, then the gene if filtered out as not expressed [default: 0.5].
' -> doc

# Required libraries
suppressMessages(library(edgeR))
suppressMessages(library(readr))
suppressMessages(library(docopt))

# Get command-line arguments
opt <- docopt(doc)

opt$min_exp <- as.numeric(opt$min_exp)
opt$min_samples_ratio <- as.numeric(opt$min_samples_ratio)

warning(paste("min_samples_ratio: ", opt$min_samples_ratio))

# to debug
#opt<-list(condition="SampleType", g1="D",g2="H",min_exp=1,min_samples_ratio=0.5)

test<-FALSE
if(test){
    opt<-list()
    setwd("/home/epigen/Proserpio_Vulvodinia/dataset/analisi_first_draft")
    opt$min_exp <- 1
    opt$min_samples_ratio <- 0.5
    use_raw_counts <- FALSE
    opt$gene_expression_matrix <- "GEP.count.gz"
    opt$sample_metadata<-"metadata.txt"
    opt$condition<-"disease"
    opt$g1<-"D"
    opt$g2<-"H"
}

###########################################
#
# Read in data
#
counts <- read.table(opt$gene_expression_matrix, header=T, sep="\t", quote="",  row.names = 1)
metadata <- read.table(opt$sample_metadata, header=T, sep="\t", quote="", row.names=1) #row.names=1 check for uniquenes of sample names
metadata$sample <- rownames(metadata)  # ensure there's a sample colum$sample

###########################
#
# Check input data
# 

#remove columns form counts if the sample is not present in metadata
initial_samples <- colnames(counts)
common_samples <- intersect(colnames(counts), metadata$sample)

counts = counts[,common_samples]
removed_colnames <- setdiff(initial_samples, common_samples)
if (length(removed_colnames) > 0) {
  warning("The following columns in the <gene_expression_matrix> were removed because the corresponding sample is not present in <sample_metadata>: ", 
          paste(removed_colnames, collapse = ", "))
}

#remove rows from metadata if the sample is not present in counts
initial_samples <- metadata$sample
metadata <- metadata[common_samples,]
removed_rownames <- setdiff(initial_samples, common_samples)
if (length(removed_rownames) > 0) {
  warning("The following samples in <sample_metadata> were removed because not preset in the ccommon_samples)",
          paste(removed_rownames, collapse = ", "))
}


###########################################
#
# Gene filtering
#

# CPM conversion
counts_cpm <- cpm(counts, normalized.lib.sizes=FALSE)

# Identify samples belonging to each group
group1_samples <- metadata$sample[metadata[[opt$condition]] == opt$g1]
group2_samples <- metadata$sample[metadata[[opt$condition]] == opt$g2]
group1 <- counts_cpm[, metadata[[opt$condition]] == opt$g1]
group2 <- counts_cpm[, metadata[[opt$condition]] == opt$g2]


min_samples_g1 <- ceiling(length(group1_samples) * opt$min_samples_ratio)
min_samples_g2 <- ceiling(length(group2_samples) * opt$min_samples_ratio)

gene_exp_g1 <- group1[rowSums(group1 >= opt$min_exp) >= min_samples_g1, ]
gene_exp_g2 <- group2[rowSums(group2 >= opt$min_exp) >= min_samples_g2, ]

# Union of expressed gene in the two groups
expressed_genes <- unique(c(rownames(gene_exp_g1), rownames(gene_exp_g2)))

warning(paste("min_samples_g1: ", min_samples_g1))
warning(paste("min_samples_g2: ", min_samples_g2))
warning(paste("expressed_genes: ", length(expressed_genes)))



# Filter original gene expression matrix
counts <- counts[rownames(counts) %in% expressed_genes,]

y <- DGEList(counts=counts);
y <- calcNormFactors(y,method="TMM");
counts_norm <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE)

group1 <- counts_norm[, metadata[[opt$condition]] == opt$g1]
group2 <- counts_norm[, metadata[[opt$condition]] == opt$g2]




###########################################
#
# Compute statistics
#

if(!identical(rownames(group1), rownames(group2))){
  stop("Dataframes should have the same row names")
}

# Function to compute Mann-Whitney test
mann_whitney <- function(gene) {
  test <- wilcox.test(gene[group1_samples], gene[group2_samples])
  return(test$p.value)
}

median_change <- apply(group1, 1, median) - apply(group2, 1, median)

# Compute p-values
p_values <- apply(counts_norm, 1, mann_whitney) 

# Adjust p-values using Benjamini-Hochberg method
p_values_adjusted <- p.adjust(p_values, method = "BH")


###########################################
#
# Write output
#


# Prepare results data frame
results <- data.frame(
  Geneid = rownames(counts),
  log_median_change = median_change,
  p_value = p_values,
  p_value_adjusted = p_values_adjusted
)
results <- results[order(results$p_value), ]

# Output

sink()

write.table(results, "/dev/stdout", sep="\t", quote=F, row.names = F)

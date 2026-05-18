#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

input_file <- args[1]
output_file <- args[2]
cols_to_remove <- args[3:length(args)]

# Read the input file
df <- read.delim(input_file)

# Remove specified columns
df <- df[, !(names(df) %in% cols_to_remove)]

# Write the modified table
write.table(
  df,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Removed columns: ", paste(cols_to_remove, collapse=", "), "\n")
cat("Output saved to: ", output_file, "\n")
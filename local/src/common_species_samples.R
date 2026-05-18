# # Load required library
# library(readxl)
# library(yaml)

# # Import and process xlsx files
# process_bracken_files <- function(dir_path = "Reports") {
#   # Get all xlsx files
#   xlsx_files <- list.files(path = dir_path, pattern = "\\.frac.txt.xlsx$", full.names = TRUE)
  
#   all_data <- list()
  
#   # Process each file
#   for (file_path in xlsx_files) { #allow for more than one file to be processed
#     file_name <- tools::file_path_sans_ext(basename(file_path))
#     df <- as.data.frame(read_excel(file_path))
#     all_data[[file_name]] <- df
#   }
  
#   # Get taxonomy IDs from column 2 of each dataframe
#   taxonomy_ids <- lapply(all_data, function(df) df[, 2])
  
#   # Find common taxonomy IDs across all files
#   common_ids <- Reduce(intersect, taxonomy_ids)
  
#   # Subset each dataframe to only include rows with common taxonomy IDs
#   subsetted_data <- lapply(all_data, function(df) {
#     df[df[, 2] %in% common_ids, ]
#   })
  
#   return(list(
#     original = all_data,
#     common_taxonomy_ids = common_ids,
#     subsets = subsetted_data
#   ))
# }

# # Run the function and store results
# results <- process_bracken_files()

# # Access results:
# # - Original dataframes: results$original
# # - Common taxonomy IDs: results$common_taxonomy_ids
# # - Subsetted dataframes: results$subsets

# # Convert columns 4 to 9 to numeric for all subset matrices
# subset_matrix_names <- names(results$subsets)

# for (subset_name in subset_matrix_names) {
#   # Get the current subset matrix
#   current_matrix <- results$subsets[[subset_name]]

#   # Check if the object is a matrix or data frame before proceeding
#   if (!is.matrix(current_matrix) && !is.data.frame(current_matrix)) {
#     warning(paste0("Object ", subset_name, " is not a matrix or data frame. Type: ", class(current_matrix)))
#     next
#   }

#   # Now safely check the number of columns
#   n_cols <- ncol(current_matrix)
#   if (is.null(n_cols)) {
#     warning(paste0("Could not determine number of columns for ", subset_name))
#     next
#   }

#   if (n_cols >= 9) {
#     # Convert columns 4 to 9 to numeric
#     for (col in 4:9) {
#       current_matrix[, col] <- as.numeric(current_matrix[, col])
#     }

#     # Reassign the modified matrix back to the results list
#     results$subsets[[subset_name]] <- current_matrix
#     cat("Converted columns 4-9 to numeric in", subset_name, "\n")
#   } else {
#     warning(paste0("Matrix ", subset_name, " has only ", n_cols, " columns (fewer than 9). Skipping conversion."))
#   }
# }

# # Assuming results$subsets contains your subsetted dataframes
# # Extract sample-specific taxonomy IDs that are non-zero across all datasets

# # Define the sample columns to extract
# samples <- snakemake@params[["samples"]]
# cat("Using samples:", paste(samples, collapse = ", "), "\n")

# # Create a list to store results for each sample
# sample_results <- list()

# # Create a list to store the sample dataframes
# sample_dataframes <- list()

# # Process each sample
# for (sample in samples) {
#   # Find all dataframes that contain columns matching this sample
#   matching_dfs <- list()
#   matching_cols <- list()
  
#   for (dataset_name in names(results$subsets)) {
#     df <- results$subsets[[dataset_name]]
#     # Find columns containing the sample name (case insensitive)
#     col_idx <- grep(paste0("_", sample, "$"), colnames(df), ignore.case = TRUE)
    
#     if (length(col_idx) > 0) {
#       matching_dfs[[dataset_name]] <- df
#       matching_cols[[dataset_name]] <- col_idx
#     }
#   }
  
#   # Skip if sample not found in any dataset
#   if (length(matching_dfs) == 0) {
#     warning(paste("Sample", sample, "not found in any dataset. Skipping."))
#     next
#   }
  
#   # Get taxonomy IDs with non-zero values for this sample across all datasets
#   non_zero_ids <- lapply(names(matching_dfs), function(dataset_name) {
#     df <- matching_dfs[[dataset_name]]
#     cols <- matching_cols[[dataset_name]]
    
#     # Check if any sample column has value > 0
#     rows_with_values <- rowSums(df[, cols, drop=FALSE] > 0) > 0
#     df[rows_with_values, 2]  # Column 2 is taxonomy ID
#   })
  
#   # Find common non-zero taxonomy IDs across all datasets where this sample appears
#   common_non_zero <- Reduce(intersect, non_zero_ids)
  
#   # Create a combined dataframe for this sample
#   if(length(common_non_zero) > 0) {
#     # Get all columns for this sample from all datasets
#     all_sample_cols <- list()
    
#     for (dataset_name in names(matching_dfs)) {
#       df <- matching_dfs[[dataset_name]]
#       sample_cols <- matching_cols[[dataset_name]]
      
#       # Extract only rows with the common non-zero IDs and columns we want
#       subset_df <- df[df[, 2] %in% common_non_zero, c(1, 2, sample_cols)]
      
#       # If this is the first dataset, use it as our base
#       if (is.null(all_sample_cols[[sample]])) {
#         all_sample_cols[[sample]] <- subset_df
#       } else {
#         # Merge with existing data by name and taxonomy_id
#         existing_df <- all_sample_cols[[sample]]
#         merged_df <- merge(existing_df, subset_df, by=c(1,2))
#         all_sample_cols[[sample]] <- merged_df
#       }
#     }
    
#     # Store the final result
#     sample_dataframes[[sample]] <- all_sample_cols[[sample]]
#   }
# }

# # After the sample_dataframes have been created, add this code to save the results

# # Create output directory if it doesn't exist
# output_dir <- "Reports"
# if (!dir.exists(output_dir)) {
#   dir.create(output_dir, recursive = TRUE)
# }

# # Save common species list (taxonomy IDs common across all files)
# common_species_file <- file.path(output_dir, "common_taxa.csv")
# common_species_df <- data.frame(taxonomy_id = results$common_taxonomy_ids)
# write.csv(common_species_df, common_species_file, row.names = FALSE)
# cat("Saved common taxa to", common_species_file, "\n")

# # Save each sample dataframe as CSV
# for (sample in names(sample_dataframes)) {
#   if (!is.null(sample_dataframes[[sample]])) {
#     sample_file <- file.path(output_dir, paste0("sample_", sample, ".csv"))
#     write.csv(sample_dataframes[[sample]], sample_file, row.names = FALSE)
#     cat("Saved sample", sample, "data to", sample_file, "\n")
#   } else {
#     warning(paste("No data available for sample", sample))
#   }
# }

# # Print summary
# cat("\nSummary of saved files:\n")
# cat("- Common species across all files:", nrow(common_species_df), "species\n")
# for (sample in names(sample_dataframes)) {
#   if (!is.null(sample_dataframes[[sample]])) {
#     cat("- Sample", sample, ":", nrow(sample_dataframes[[sample]]), "species\n")
#   }
# }


#======================================================================================

# Load required library
library(readxl)
library(yaml)

# Import and process tsv files
process_bracken_files <- function(dir_path = "Reports") {
  # Get all tsv files
  tsv_files <- list.files(path = dir_path, pattern = "\\.frac.tsv$", full.names = TRUE)
  
  all_data <- list()
  
  # Process each file
  for (file_path in tsv_files) { #allow for more than one file to be processed
    file_name <- tools::file_path_sans_ext(basename(file_path))
    df <- as.data.frame(read_csv(file_path))
    all_data[[file_name]] <- df
  }
  
  # Get taxonomy IDs from column 2 of each dataframe
  taxonomy_ids <- lapply(all_data, function(df) df[, 2])
  
  # Find common taxonomy IDs across all files
  common_ids <- Reduce(intersect, taxonomy_ids)
  
  # Subset each dataframe to only include rows with common taxonomy IDs
  subsetted_data <- lapply(all_data, function(df) {
    df[df[, 2] %in% common_ids, ]
  })
  
  return(list(
    original = all_data,
    common_taxonomy_ids = common_ids,
    subsets = subsetted_data
  ))
}

# Run the function and store results
results <- process_bracken_files()

# Access results:
# - Original dataframes: results$original
# - Common taxonomy IDs: results$common_taxonomy_ids
# - Subsetted dataframes: results$subsets

# Convert columns 4 to 9 to numeric for all subset matrices
subset_matrix_names <- names(results$subsets)

for (subset_name in subset_matrix_names) {
  # Get the current subset matrix
  current_matrix <- results$subsets[[subset_name]]

  # Check if the object is a matrix or data frame before proceeding
  if (!is.matrix(current_matrix) && !is.data.frame(current_matrix)) {
    warning(paste0("Object ", subset_name, " is not a matrix or data frame. Type: ", class(current_matrix)))
    next
  }

  # Now safely check the number of columns
  n_cols <- ncol(current_matrix)
  if (is.null(n_cols)) {
    warning(paste0("Could not determine number of columns for ", subset_name))
    next
  }

  if (n_cols >= 9) {
    # Convert columns 4 to 9 to numeric
    for (col in 4:9) {
      current_matrix[, col] <- as.numeric(current_matrix[, col])
    }

    # Reassign the modified matrix back to the results list
    results$subsets[[subset_name]] <- current_matrix
    cat("Converted columns 4-9 to numeric in", subset_name, "\n")
  } else {
    warning(paste0("Matrix ", subset_name, " has only ", n_cols, " columns (fewer than 9). Skipping conversion."))
  }
}

# Assuming results$subsets contains your subsetted dataframes
# Extract sample-specific taxonomy IDs that are non-zero across all datasets

# Define the sample columns to extract
samples <- snakemake@params[["samples"]]
cat("Using samples:", paste(samples, collapse = ", "), "\n")

# Create a list to store results for each sample
sample_results <- list()

# Create a list to store the sample dataframes
sample_dataframes <- list()

# Process each sample
for (sample in samples) {
  # Find all dataframes that contain columns matching this sample
  matching_dfs <- list()
  matching_cols <- list()
  
  for (dataset_name in names(results$subsets)) {
    df <- results$subsets[[dataset_name]]
    # Find columns containing the sample name (case insensitive)
    col_idx <- grep(paste0("_", sample, "$"), colnames(df), ignore.case = TRUE)
    
    if (length(col_idx) > 0) {
      matching_dfs[[dataset_name]] <- df
      matching_cols[[dataset_name]] <- col_idx
    }
  }
  
  # Skip if sample not found in any dataset
  if (length(matching_dfs) == 0) {
    warning(paste("Sample", sample, "not found in any dataset. Skipping."))
    next
  }
  
  # Get taxonomy IDs with non-zero values for this sample across all datasets
  non_zero_ids <- lapply(names(matching_dfs), function(dataset_name) {
    df <- matching_dfs[[dataset_name]]
    cols <- matching_cols[[dataset_name]]
    
    # Check if any sample column has value > 0
    rows_with_values <- rowSums(df[, cols, drop=FALSE] > 0) > 0
    df[rows_with_values, 2]  # Column 2 is taxonomy ID
  })
  
  # Find common non-zero taxonomy IDs across all datasets where this sample appears
  common_non_zero <- Reduce(intersect, non_zero_ids)
  
  # Create a combined dataframe for this sample
  if(length(common_non_zero) > 0) {
    # Get all columns for this sample from all datasets
    all_sample_cols <- list()
    
    for (dataset_name in names(matching_dfs)) {
      df <- matching_dfs[[dataset_name]]
      sample_cols <- matching_cols[[dataset_name]]
      
      # Extract only rows with the common non-zero IDs and columns we want
      subset_df <- df[df[, 2] %in% common_non_zero, c(1, 2, sample_cols)]
      
      # If this is the first dataset, use it as our base
      if (is.null(all_sample_cols[[sample]])) {
        all_sample_cols[[sample]] <- subset_df
      } else {
        # Merge with existing data by name and taxonomy_id
        existing_df <- all_sample_cols[[sample]]
        merged_df <- merge(existing_df, subset_df, by=c(1,2))
        all_sample_cols[[sample]] <- merged_df
      }
    }
    
    # Store the final result
    sample_dataframes[[sample]] <- all_sample_cols[[sample]]
  }
}

# After the sample_dataframes have been created, add this code to save the results

# Create output directory if it doesn't exist
output_dir <- "Reports"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save common species list (taxonomy IDs common across all files)
common_species_file <- file.path(output_dir, "common_taxa.csv")
common_species_df <- data.frame(taxonomy_id = results$common_taxonomy_ids)
write.csv(common_species_df, common_species_file, row.names = FALSE)
cat("Saved common taxa to", common_species_file, "\n")

# Save each sample dataframe as CSV
for (sample in names(sample_dataframes)) {
  if (!is.null(sample_dataframes[[sample]])) {
    sample_file <- file.path(output_dir, paste0("sample_", sample, ".csv"))
    write.csv(sample_dataframes[[sample]], sample_file, row.names = FALSE)
    cat("Saved sample", sample, "data to", sample_file, "\n")
  } else {
    warning(paste("No data available for sample", sample))
  }
}

# Print summary
cat("\nSummary of saved files:\n")
cat("- Common species across all files:", nrow(common_species_df), "species\n")
for (sample in names(sample_dataframes)) {
  if (!is.null(sample_dataframes[[sample]])) {
    cat("- Sample", sample, ":", nrow(sample_dataframes[[sample]]), "species\n")
  }
}

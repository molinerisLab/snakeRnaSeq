library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(gridExtra)
library(readxl)
library(purrr)

##########################################
# Part 1: Top Species Analysis (from top_species.R)
##########################################

message("*** PART 1: Processing top species analysis ***")

# List all braken files
files <- list.files("boutputs_filtered", pattern = "\\.braken$", full.names = TRUE)

all_data <- lapply(files, function(f) {
  df <- read_tsv(f)
  sample <- str_replace(basename(f), "\\.braken$", "")
  if (nrow(df) > 0 && all(c("name", "fraction_total_reads") %in% colnames(df))) {
    df %>%
      select(species = name, fraction = fraction_total_reads) %>%
      mutate(sample = sample)
  } else {
    NULL
  }
})

all_df <- bind_rows(Filter(Negate(is.null), all_data))

if (nrow(all_df) > 0) {
  # Top 10 per sample
  top10_per_sample <- all_df %>%
    group_by(sample) %>%
    arrange(desc(fraction), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  write_csv(top10_per_sample, "top10_taxa_per_sample.csv")
  
  # Find common taxa among top 10 of all samples
  top10_lists <- split(top10_per_sample$species, top10_per_sample$sample)
  common_species <- Reduce(intersect, top10_lists)
  
  # Save common species to a CSV
  write_csv(tibble(common_species = common_species), "common_top10_taxa_all_samples.csv")
  
  # Overall top 10 species (sum fractions across all samples)
  top10_overall_species <- all_df %>%
    group_by(species) %>%
    summarise(total_fraction = sum(fraction, na.rm = TRUE)) %>%
    arrange(desc(total_fraction)) %>%
    slice_head(n = 10) %>%
    pull(species)
  
  # Get individual fractions for these species in each sample
  top10_overall <- all_df %>%
    filter(species %in% top10_overall_species) %>%
    select(species, sample, fraction) %>%
    pivot_wider(names_from = sample, values_from = fraction, values_fill = 0)
  
  write_csv(top10_overall, "top10_taxa_overall.csv")
} else {
  write_csv(tibble(), "top10_taxa_per_sample.csv")
  write_csv(tibble(), "top10_taxa_overall.csv")
}

# Read the top 10 species per sample
df <- tryCatch({
  data <- read_csv("top10_taxa_per_sample.csv")
  # Check which columns are available
  message("Columns in CSV: ", paste(colnames(data), collapse=", "))
  
  # Rename columns if they don't match expectations
  if ("species" %in% colnames(data) && "sample" %in% colnames(data) && "fraction" %in% colnames(data)) {
    data
  } else if (ncol(data) == 3) {
    # Assume the columns are in the right order but named incorrectly
    message("Renaming columns to expected format")
    colnames(data) <- c("species", "fraction", "sample")
    data
  } else {
    message("CSV file has unexpected format")
    NULL
  }
}, error = function(e) {
  message("Error reading top10_taxa_per_sample.csv: ", e$message)
  NULL
})

if (!is.null(df) && "species" %in% colnames(df)) {
  # Replace the current factor assignment with alphabetical sorting
  df <- df %>%
    mutate(species = factor(species, levels = sort(unique(species))))
  
  # Stacked bar plot
  p <- ggplot(df, aes(x = sample, y = fraction, fill = species)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
    labs(x = "Sample", y = "Fraction", fill = "Species",
         title = "Top 10 taxa per Sample (Stacked Barplot)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("top10_taxa_per_sample_barplot.png", p, width = 10, height = 6)
  message("Created top10_taxa_per_sample_barplot.png")
} else {
  message("Skipping top10_taxa_per_sample_barplot.png creation - missing required data")
}

# Create a faceted plot for all patients
create_patients_faceted_plot <- function() {
  # Find all patient directories
  patient_dirs <- list.dirs("../", recursive = FALSE)
  patient_dirs <- grep("P_", patient_dirs, value = TRUE)
  message(paste("Found", length(patient_dirs), "patient directories"))
  
  # Create a list to store data from all patients
  all_patients_data <- list()
  
  # Process each patient directory
  for (patient_dir in patient_dirs) {
    patient_id <- basename(patient_dir)
    message(paste("Processing", patient_id))
    
    # Path to braken files
    braken_path <- file.path(patient_dir, "boutputs_filtered")
    
    # Check if the directory exists
    if (!dir.exists(braken_path)) {
      message(paste("Directory not found:", braken_path))
      next
    }
    
    # List all braken files for this patient
    files <- list.files(braken_path, pattern = "\\.braken$", full.names = TRUE)
    
    if (length(files) == 0) {
      message(paste("No .braken files found in", braken_path))
      next
    }
    
    # Process each braken file
    patient_data <- lapply(files, function(f) {
      df <- read_tsv(f)
      sample <- str_replace(basename(f), "\\.braken$", "")
      
      if (nrow(df) > 0 && all(c("name", "fraction_total_reads") %in% colnames(df))) {
        df %>%
          select(species = name, fraction = fraction_total_reads) %>%
          mutate(
            sample = sample,
            patient_id = patient_id  # Add patient ID
          )
      } else {
        NULL
      }
    })
    
    # Combine data for this patient
    patient_df <- bind_rows(Filter(Negate(is.null), patient_data))
    
    if (nrow(patient_df) > 0) {
      # Get top 10 species for this patient (across all samples)
      top_species <- patient_df %>%
        group_by(species) %>%
        summarise(total_fraction = sum(fraction, na.rm = TRUE)) %>%
        arrange(desc(total_fraction)) %>%
        slice_head(n = 10) %>%
        pull(species)
      
      # Filter for just the top 10 species
      filtered_df <- patient_df %>%
        filter(species %in% top_species)
      
      # Add to all patients data
      all_patients_data[[patient_id]] <- filtered_df
    }
  }
  
  # Combine all patient data
  combined_df <- bind_rows(all_patients_data)
  
  if (nrow(combined_df) > 0) {
    # Create output directory if it doesn't exist
    if (!dir.exists("plots")) {
      dir.create("plots", recursive = TRUE)
    }
    
    # Create a faceted plot with one facet per patient
    p <- ggplot(combined_df, aes(x = sample, y = fraction, fill = species)) +
      geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
      facet_wrap(~ patient_id, scales = "free_x") +
      labs(
        x = "Sample", 
        y = "Fraction", 
        fill = "Species",
        title = "Top 10 taxa per Patient (Faceted Plot)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        panel.spacing = unit(1, "lines")
      )
    
    # Save the combined faceted plot WITH THE EXACT FILENAME EXPECTED BY SNAKEMAKE
    ggsave("plots/all_patients_top_taxa_faceted.png", p, width = 16, height = 10, dpi = 300)
    message("Created faceted plot with all patients' data")
  } else {
    message("No valid data found across patient directories")
  }
}

# Call the function to create the faceted plot
create_patients_faceted_plot()

###########################################
# Generate top species barplot for all samples normalized to 1 million reads (CPM)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

message("\n*** Generating top species barplot normalized to 1 million reads (CPM) ***")
bmerged = read_excel("bracken_merged_abbundances.num.txt.xlsx")
bmerged = bmerged[,-c(2,3)]
counts = read.table("Results/reads_summary.csv", header = TRUE, sep = "\t")
raw_counts = data.frame(name=readnum$sample,initial_reads= readnum$initial_reads)
bmerged_long <- bmerged %>%
pivot_longer(cols = -name, names_to = "sample", values_to = "fraction") %>%
left_join(raw_counts, by = c("sample" = "name")) %>%
mutate(CPM = (fraction / initial_reads) * 1000000)

sored_bmerged <- bmerged_long %>%
arrange(desc(CPM))

sored_bmerged <- bmerged_long %>%
  group_by(sample) %>%
  arrange(desc(CPM), .by_group = TRUE) %>%  # Sort by CPM within each sample
  mutate(rank = row_number()) %>%
  filter(rank <= 10) %>%
  ungroup() %>%
  select(-rank)

p <- ggplot(sored_bmerged, aes(x = sample, y = CPM, fill = name)) +
geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
labs(x = "Sample", y = "CPM", fill = "Species",
title = "Top 10 Species per Sample (Stacked Barplot)") +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("barplot_CPM.png", p, width = 10, height = 6)

Now i want to plot the same plot but stackbar to 100(so i see the percentage of each plotted)

# Create 100% stacked barplot (percentage view)
p_percent <- ggplot(sored_bmerged, aes(x = sample, y = CPM, fill = name)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
    labs(x = "Sample", y = "Percentage", fill = "Taxa",
         title = "Top 10 Taxa per Sample (100% Stacked Barplot)") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave("barplot_CPM_percentage.png", p_percent, width = 10, height = 6)

#################Facet barplot for top 10 taxa per sample with CPM for all patients ###############

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Directory containing all patient files
data_dir <- "../analysis_P/bracken_num"
xlsx_files <- list.files(data_dir, pattern = "bracken_merged_abbundances.num.txt.xlsx$", full.names = TRUE)
reads_files <- list.files(data_dir, pattern = "reads_summary.csv$", full.names = TRUE)

# Helper: get patient/sample prefix from filename
get_patient_id <- function(path) str_split(basename(path), "_", simplify = TRUE)[1]

all_data <- list()

for (xlsx in xlsx_files) {
  patient_id <- get_patient_id(xlsx)
  # Find corresponding reads_summary.csv
  reads_file <- reads_files[grepl(patient_id, reads_files)]
  if (length(reads_file) == 0) next
  
  bmerged <- read_excel(xlsx)
  if (ncol(bmerged) > 2) bmerged <- bmerged[,-c(2,3)]
  reads_summary <- read.table(reads_file, header = TRUE, sep = "\t")
  raw_counts <- data.frame(name = reads_summary$sample, initial_reads = reads_summary$initial_reads)
  
  bmerged_long <- bmerged %>%
    pivot_longer(cols = -name, names_to = "sample", values_to = "fraction") %>%
    left_join(raw_counts, by = c("sample" = "name")) %>%
    mutate(CPM = (fraction / initial_reads) * 1e6, patient = patient_id)
  
  # Top 10 taxa per sample
  top10 <- bmerged_long %>%
    group_by(patient, sample) %>%
    arrange(desc(CPM), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  all_data[[patient_id]] <- top10
}

combined <- bind_rows(all_data)

p <- ggplot(combined, aes(x = sample, y = CPM, fill = name)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
  facet_wrap(~ patient, scales = "free_x") +
  labs(x = "Sample", y = "CPM", fill = "Taxa",
       title = "Top 10 Taxa per Sample (CPM, Faceted by Patient)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom")

ggsave("facet_barplot_CPM_patients.png", p, width = 16, height = 10, dpi = 300)

p_percent <- ggplot(combined, aes(x = sample, y = CPM, fill = name)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  facet_wrap(~ patient, scales = "free_x") +
  labs(x = "Sample", y = "Percentage", fill = "Taxa",
       title = "Top 10 Taxa per Sample (100% Stacked Barplot, Faceted by Patient)") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom")

ggsave("facet_barplot_CPM_percentage_patients.png", p_percent, width = 16, height = 10, dpi = 300)

##########################################
# Part 2: Plot Reads Analysis (from plot_reads.R)
##########################################

message("\n*** PART 2: Processing reads data ***")

if (file.exists("reads_summary.csv")) {
  df <- read_tsv("reads_summary.csv")
  
  df$classified_kraken <- as.numeric(df$classified_kraken)
  df$Kraken_human_assigned_reads <- as.numeric(df$Kraken_human_assigned_reads)
  df$contaminant_reads <- as.numeric(df$contaminant_reads)
  
  plot_df <- df %>%
    mutate(
      non_human = classified_kraken - Kraken_human_assigned_reads - contaminant_reads,
      # Ensure non_human doesn't go negative
      non_human = pmax(non_human, 0)
    ) %>%
    select(sample, Kraken_human_assigned_reads, contaminant_reads, non_human) %>%
    pivot_longer(cols = c("Kraken_human_assigned_reads", "contaminant_reads", "non_human"),
                 names_to = "type", values_to = "count")
  
  p <- ggplot(plot_df, aes(x = sample, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = format(count, big.mark = ",", scientific = FALSE)),
              position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
    labs(x = "Sample", y = "Classified Reads", fill = "Type",
         title = "Classified Kraken Reads: Human, Contaminant, and Non-human per Sample") +
    scale_fill_manual(values = c("Kraken_human_assigned_reads" = "red", 
                                "contaminant_reads" = "orange",
                                "non_human" = "steelblue"),
                      labels = c("Human", "Contaminant", "Non-human")) +
    theme_minimal()
  
  ggsave("classified_vs_human_contaminant_barplot.png", p, width = 14, height = 10)
  message("Created classified vs. human vs. contaminant reads barplot")
  
  # Create a 100% stacked bar chart with CPM values including contaminants
  cpm_df <- df %>%
    mutate(
      non_human = pmax(classified_kraken - Kraken_human_assigned_reads - contaminant_reads, 0),
      total_reads = Kraken_human_assigned_reads + contaminant_reads + non_human,
      # Calculate CPM for human, contaminant, and non-human reads
      human_cpm = (Kraken_human_assigned_reads / total_reads) * 1000000,
      contaminant_cpm = (contaminant_reads / total_reads) * 1000000,
      non_human_cpm = (non_human / total_reads) * 1000000
    ) %>%
    select(sample, human_cpm, contaminant_cpm, non_human_cpm) %>%
    pivot_longer(cols = c("human_cpm", "contaminant_cpm", "non_human_cpm"),
                names_to = "type", values_to = "cpm")
  
  # Convert type labels to be more readable
  cpm_df$type <- case_when(
    cpm_df$type == "human_cpm" ~ "Human",
    cpm_df$type == "contaminant_cpm" ~ "Contaminant",
    cpm_df$type == "non_human_cpm" ~ "Non-human",
    TRUE ~ cpm_df$type
  )
  
  # Create the 100% stacked bar chart
  p_percent <- ggplot(cpm_df, aes(x = sample, y = cpm, fill = type)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f CPM", cpm/1000)),
              position = position_fill(vjust = 0.5), color = "black", fontface = "bold", size = 4.5) +
    labs(x = "Sample", y = "Proportion", fill = "Type",
         title = "Proportion of Human vs. Contaminant vs. Non-human Reads (CPM)") +
    scale_fill_manual(values = c("Human" = "red", "Contaminant" = "orange", "Non-human" = "steelblue")) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_blank()
    )
  
  ggsave("human_vs_contaminant_vs_nonhuman_percent_cpm.png", p_percent, width = 14, height = 10)
  message("Created 100% stacked bar chart with CPM values including contaminants")
} else {
  warning("reads_summary.csv not found - skipping classified vs. human reads plots")
}

# Extract read counts from top10 species
message("Processing read counts for top 10 species...")

# Function to extract read count from k2report files
extract_read_count <- function(sample_name, species_name) {
  # Build the path to the breport_file
  k2report_path <- file.path("breports_filtered", paste0(sample_name, ".breport"))
  
  # Check if the file exists
  if (!file.exists(k2report_path)) {
    warning(paste("File not found:", k2report_path))
    return(NA)
  }
  
  # Read the k2report file
  k2report <- read.delim(k2report_path, header = FALSE, sep = "\t", 
                        stringsAsFactors = FALSE)
  
  # Extract the required columns (read count is in column 2, name in column 6)
  k2report <- data.frame(
    reads = k2report[, 2],
    name = trimws(k2report[, 6]),
    stringsAsFactors = FALSE
  )
  
  # Clean the species names by removing rank prefixes (e.g., "S 9606 Homo sapiens" becomes "Homo sapiens")
  k2report$clean_name <- str_replace(k2report$name, "^[A-Z][0-9]?\\s+\\d+\\s+", "")
  
  # Find the row matching the species name
  matches <- grep(paste0("^", species_name, "$"), k2report$clean_name, fixed = TRUE)
  
  if (length(matches) > 0) {
    return(k2report$reads[matches[1]])
  } else {
    # Try a more flexible match if exact match fails
    matches <- grep(species_name, k2report$clean_name, fixed = TRUE)
    if (length(matches) > 0) {
      return(k2report$reads[matches[1]])
    } else {
      warning(paste("Species not found:", species_name, "in file:", k2report_path))
      return(NA)
    }
  }
}

# Read the top10 species data generated in Part 1
if (file.exists("top10_species_per_sample.csv")) {
  top10 <- read.csv("top10_species_per_sample.csv", stringsAsFactors = FALSE)
  
  # Ensure we have the expected three columns
  if (ncol(top10) != 3) {
    warning("Expected top10_species_per_sample.csv to have 3 columns: species,fraction,sample")
  } else {
    colnames(top10) <- c("species", "fraction", "sample")
    
    # Initialize the total_reads column
    top10$total_reads <- NA
    
    # Process each row
    for (i in 1:nrow(top10)) {
      sample_name <- top10$sample[i]
      species_name <- top10$species[i]
      
      # Extract read count and add to dataframe
      read_count <- extract_read_count(sample_name, species_name)
      top10$total_reads[i] <- read_count
      
      # Print progress
      cat(sprintf("Processing %s - %s: %s reads\n", 
                  sample_name, species_name, 
                  ifelse(is.na(read_count), "NA", read_count)))
    }
    
    # Write the updated dataframe to CSV
    write.csv(top10, "top10_species_with_reads.csv", row.names = FALSE)
    
    # Create a barplot showing species composition per sample using total reads
    top10_with_reads <- top10  # We already have the data
    
    # Replace the current factor assignment with alphabetical sorting
    top10_with_reads$species <- factor(top10_with_reads$species, 
                                      levels = sort(unique(top10_with_reads$species)))
    
    # Create a stacked bar plot showing species composition for each sample
    p_species_per_sample <- ggplot(top10_with_reads, aes(x = sample, y = total_reads, fill = species)) +
      geom_bar(stat = "identity", position = "stack", color = "black") +  # Added black border
      theme_minimal() +
      labs(title = "Species Composition per Sample",
           x = "Sample", 
           y = "Total Reads",
           fill = "Species") +
      scale_y_continuous(labels = scales::comma) +
      theme(legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5)) +
      guides(fill = guide_legend(ncol = 1))
    
    # Save the plot
    ggsave("species_composition_per_sample.png", p_species_per_sample, width = 10, height = 7, dpi = 300)
    
    # Create a modified version with percentages instead of absolute values
    p_species_percentage <- top10_with_reads %>%
      group_by(sample) %>%
      mutate(percentage = total_reads / sum(total_reads, na.rm = TRUE) * 100) %>%
      ungroup() %>%
      ggplot(aes(x = sample, y = percentage, fill = species)) +
      geom_bar(stat = "identity", position = "stack", color = "black") +  # Added black border
      theme_minimal() +
      labs(title = "Species Composition per Sample (Percentage)",
           x = "Sample", 
           y = "Percentage (%)",
           fill = "Species") +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      theme(legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5)) +
      guides(fill = guide_legend(ncol = 1))
    
    # Save the percentage plot
    ggsave("species_percentage_per_sample.png", p_species_percentage, width = 10, height = 7, dpi = 300)
    
    message("Analysis complete. Created species composition plots")
  }
} else {
  warning("top10_species_per_sample.csv not found - skipping species composition plots")
}

##########################################
# Add faceted barplot for top10_species_with_reads.csv
##########################################

message("\n*** Creating faceted barplot from top10_species_with_reads.csv ***")

if (file.exists("top10_species_with_reads.csv")) {
  # Read the data
  top10_reads <- read_csv("top10_species_with_reads.csv")
  
  # Verify the data has the expected columns
  if (!all(c("species", "sample", "total_reads") %in% colnames(top10_reads))) {
    warning("top10_species_with_reads.csv is missing expected columns")
  } else {
    # For each sample, reorder species by total_reads to ensure consistent ordering
    top10_reads <- top10_reads %>%
      group_by(sample) %>%
      mutate(
        # Create a factor of species ordered by read count within each sample
        species_ordered = factor(species, 
                               levels = species[order(total_reads, decreasing = TRUE)])
      ) %>%
      ungroup()
    
    # Create the faceted barplot
    p_faceted_reads <- ggplot(top10_reads, aes(x = species_ordered, y = total_reads, fill = species_ordered)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
      facet_wrap(~ sample, scales = "free_y", ncol = 2) +  # Put samples in facets
      labs(
        title = "Top 10 Species by Sample",
        x = "Species",
        y = "Total Reads"
      ) +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        legend.position = "none",  # Remove the legend as it's redundant
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(1, "lines")
      )
    
    # Create output directory if it doesn't exist
    if (!dir.exists("plots")) {
      dir.create("plots", recursive = TRUE)
    }
    
    # Save the plot
    output_file_png <- file.path("plots", "top10_species_reads_faceted.png")
    ggsave(output_file_png, p_faceted_reads, width = 14, height = 10, dpi = 300)
    
    message(paste("Created faceted barplot at", output_file_png))
  }
} else {
  warning("top10_species_with_reads.csv not found - skipping faceted barplot")
}

##########################################
# Part 3: Plot Sample Species Analysis (from plot_sample_species.R)
##########################################

message("\n*** PART 3: Processing sample species analysis ***")
# Add code to process bracken files and create results_csv directory


# Function to process all CSV files in the results_csv folder and create a single faceted plot
process_all_csv_files <- function() {
  # Find all directories named "results_csv"
  base_dir <- "."
  csv_dirs <- list.dirs(base_dir, recursive = TRUE) %>%
    grep("results_csv$", ., value = TRUE)
  
  if (length(csv_dirs) == 0) {
    message("No 'results_csv' directories found - skipping CSV processing")
    return(NULL)
  }
  
  # Get all CSV files from these directories
  csv_files <- lapply(csv_dirs, function(dir) {
    list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  }) %>% unlist()
  
  if (length(csv_files) == 0) {
    message("No CSV files found in results_csv directories")
    return(NULL)
  }
  
  message(paste("Found", length(csv_files), "CSV files to process"))
  
  # Create empty list to store combined data from all CSV files
  all_data <- list()
  
  # Process each CSV file
  for (i in seq_along(csv_files)) {
    csv_path <- csv_files[i]
    # Extract filename without extension for output naming
    filename <- tools::file_path_sans_ext(basename(csv_path))
    
    # Read the CSV file
    data <- tryCatch(
      read_csv(csv_path),
      error = function(e) {
        warning(paste("Error reading", csv_path, ":", e$message))
        return(NULL)
      }
    )
    
    if (is.null(data)) next
    
    # Check if the file has the expected structure
    if (!("name" %in% colnames(data))) {
      warning(paste("Skipping", csv_path, "- 'name' column not found"))
      next
    }
    
    # Get sample columns (all except name and taxonomy_id)
    sample_cols <- setdiff(colnames(data), c("name", "taxonomy_id"))
    
    # Select top 10 species across all samples
    top_species <- data %>%
      mutate(total = rowSums(across(all_of(sample_cols)), na.rm = TRUE)) %>%
      arrange(desc(total)) %>%
      slice_head(n = 10) %>%
      pull(name)
    
    # Filter data to include only top species
    plot_data <- data %>%
      filter(name %in% top_species) %>%
      select(name, all_of(sample_cols)) %>%
      pivot_longer(cols = all_of(sample_cols), 
                  names_to = "Sample", 
                  values_to = "Reads") %>%
      mutate(source_file = filename)  # Add source file information
    
    # Add to the combined dataset
    all_data[[i]] <- plot_data
  }
  
  if (length(all_data) == 0) {
    message("No valid data found in CSV files")
    return(NULL)
  }
  
  # Combine all datasets
  combined_data <- bind_rows(all_data)
  
  # Create output directory if it doesn't exist
  if (!dir.exists("plots")) {
    dir.create("plots", recursive = TRUE)
  }
  
  # Create the faceted plot with better layout for 6 samples
  p_faceted <- ggplot(combined_data, aes(x = Sample, y = Reads, fill = name)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    facet_wrap(~ source_file, scales = "free_x", ncol = 3) +  # 3 columns for better layout
    labs(
      title = "Species Composition per Sample (All Datasets)",
      x = "Sample",
      y = "Total Reads",
      fill = "Species"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom", 
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_line(color = "grey98"),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    guides(fill = guide_legend(ncol = 3))  # Arrange legend in 3 columns
  
  # Save the faceted plot
  output_file_png <- file.path("plots", "all_species_composition_faceted.png")
  ggsave(output_file_png, p_faceted, width = 18, height = 12, dpi = 300)
  message(paste("Created combined faceted plot at", output_file_png))
  
  # Create individual plots for each sample type (AD, AP, M, T, OD, OP)
  sample_types <- c("AD", "AP", "M", "T", "OD", "OP")
  
  for (sample_type in sample_types) {
    # Filter data for the specific sample type
    sample_data <- combined_data %>%
      filter(str_detect(Sample, sample_type))
    
    if (nrow(sample_data) == 0) {
      message(paste("No data found for sample type:", sample_type))
      next
    }
    
    # Create individual plot for this sample type
    p_individual <- ggplot(sample_data, aes(x = Sample, y = Reads, fill = name)) +
      geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
      facet_wrap(~ source_file, scales = "free_x") +
      labs(
        title = paste("Species Composition -", sample_type, "Samples"),
        x = "Sample",
        y = "Total Reads",
        fill = "Species"
      ) +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightblue", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right", 
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_line(color = "grey98"),
        panel.spacing = unit(1.5, "lines"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      ) +
      guides(fill = guide_legend(ncol = 1))
    
    # Save individual plot
    individual_output <- file.path("plots", paste0("species_composition_", sample_type, ".png"))
    ggsave(individual_output, p_individual, width = 14, height = 8, dpi = 300)
    message(paste("Created individual plot for", sample_type, "at", individual_output))
  }
  
  return(p_faceted)
}
# Function to create faceted barplots of top 10 species per sample
plot_top_species_per_sample <- function() {
  # Find all directories named "results_csv"
  base_dir <- "."
  csv_dirs <- list.dirs(base_dir, recursive = TRUE) %>%
    grep("results_csv$", ., value = TRUE)
  
  if (length(csv_dirs) == 0) {
    message("No 'results_csv' directories found - skipping top species per sample plot")
    return(NULL)
  }
  
  # Get all CSV files from these directories
  csv_files <- lapply(csv_dirs, function(dir) {
    list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  }) %>% unlist()
  
  if (length(csv_files) == 0) {
    message("No CSV files found in results_csv directories")
    return(NULL)
  }
  
  # Create empty list to store sample data
  samples_data <- list()
  
  # Process each CSV file
  for (csv_path in csv_files) {
    # Read the CSV file
    data <- tryCatch(
      read_csv(csv_path),
      error = function(e) {
        warning(paste("Error reading", csv_path, ":", e$message))
        return(NULL)
      }
    )
    
    if (is.null(data) || !("name" %in% colnames(data))) {
      next
    }
    
    # Get sample columns (all except name and taxonomy_id)
    sample_cols <- setdiff(colnames(data), c("name", "taxonomy_id"))
    
    # Process each sample column
    for (sample_col in sample_cols) {
      # Get simple sample name without dataset prefix
      sample_name <- str_extract(sample_col, "[^_]+$")
      
      if (is.na(sample_name)) sample_name <- sample_col
      
      # Select top 10 species for this sample
      top_species_data <- data %>%
        arrange(desc(!!sym(sample_col))) %>%
        slice_head(n = 10) %>%
        select(name, !!sample_col) %>%
        mutate(
          Sample = sample_name,
          Species = name,
          Reads = !!sym(sample_col)
        ) %>%
        select(Sample, Species, Reads)
      
      # Add to the list
      samples_data[[length(samples_data) + 1]] <- top_species_data
    }
  }
  
  # Combine all sample data
  if (length(samples_data) == 0) {
    message("No valid data found for plotting top species per sample")
    return(NULL)
  }
  
  all_samples_data <- bind_rows(samples_data)
  
  # Reorder species by read count within each sample
  all_samples_data <- all_samples_data %>%
    group_by(Sample) %>%
    mutate(
      # Add this line to handle duplicate species names
      Species_unique = make.unique(as.character(Species)),
      # Use the unique species names for creating the factor
      Species = factor(Species_unique, 
                      levels = Species_unique[order(Reads, decreasing = TRUE)])
    ) %>%
    ungroup() %>%
    # Remove the temporary column after factoring
    select(-Species_unique)
  
  # Create faceted barplot
  p <- ggplot(all_samples_data, aes(x = Species, y = Reads, fill = Species)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    facet_wrap(~ Sample, scales = "free_y", ncol = 2) +
    labs(
      title = "Top 10 Species per Sample",
      x = "Species",
      y = "Total Reads"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines")
    )
  
  # Create output directory if it doesn't exist
  if (!dir.exists("plots")) {
    dir.create("plots", recursive = TRUE)
  }
  
  # Save the plot as a single PNG file
  output_file_png <- file.path("plots", "top_species_per_sample.png")
  ggsave(output_file_png, p, width = 16, height = 12, dpi = 300)
  
  message(paste("Created top species per sample plot at", output_file_png))
  
  return(p)
}

# Run the main functions from Part 3
process_all_csv_files()
plot_top_species_per_sample()

message("\n*** Processing complete! ***")
message("All plots saved to the 'plots' directory")
message("CSV files saved to the current working directory")


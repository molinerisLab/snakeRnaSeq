library(readxl)
library(tidyr)
library(yaml)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)

#-----------------------#
#--- Helper Functions --#
#-----------------------#

read_counts <- function(filepath) {
  counts <- read.table(filepath, header = TRUE, sep = "\t")
  data.frame(name = counts$sample, initial_reads = counts$initial_reads)
}

prepare_long_data <- function(bmerged, raw_counts, patient_id = NULL) {
  bmerged %>%
    pivot_longer(cols = -name, names_to = "sample", values_to = "fraction") %>%
    left_join(raw_counts, by = c("sample" = "name")) %>%
    mutate(
      CPM = (fraction / initial_reads) * 1e6,
      patient = patient_id
    )
}

get_top10 <- function(df) {
  df %>%
    group_by(across(any_of(c("sample", "patient")))) %>%
    arrange(desc(CPM), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
}

plot_stacked <- function(data, file, title, percent = FALSE, facet = FALSE) {
  p <- ggplot(data, aes(x = sample, y = CPM, fill = name)) +
    geom_bar(stat = "identity", position = ifelse(percent, "fill", "stack"),
             color = "black", size = 0.3) +
    labs(
      x = "Sample",
      y = ifelse(percent, "Percentage", "CPM"),
      fill = "Taxa",
      title = title
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (percent) {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  if (facet) {
    p <- p + facet_wrap(~patient, scales = "free_x")
  }

  ggsave(file, p, width = ifelse(facet, 16, 30), height = 10, dpi = 300)
}

# Create Plots directory if it doesn't exist
if (!dir.exists("Plots")) dir.create("Plots")

#------------------------------#
#--- Main: Single Dataset ----#
#------------------------------#

bmerged <- read_excel("Reports/bracken_merged_abbundances.num.txt.xlsx")[,-c(2,3)]
raw_counts <- read_counts("Reports/reads_summary.csv")
bmerged_long <- prepare_long_data(bmerged, raw_counts)
top10_data <- get_top10(bmerged_long)

plot_stacked(top10_data, "Plots/barplot_CPM.png", "Top 10 Species per Sample")
plot_stacked(top10_data, "Plots/barplot_CPM_percentage.png", "Top 10 Taxa per Sample (100% Stacked)", percent = TRUE)

#---------------------------------------------#
#--- Optional: Multiple Patients via Config --#
#---------------------------------------------#

config <- yaml.load_file("config.yaml")
if (isTRUE(config$multiple_files)) {
  xlsx_files <- list.files(".", pattern = "bracken_merged_abbundances.num.txt.xlsx$", full.names = TRUE)
  reads_files <- list.files(".", pattern = "reads_summary.csv$", full.names = TRUE)
  
  get_patient_id <- function(path) str_split(basename(path), "_", simplify = TRUE)[1]
  
  all_data <- list()
  for (xlsx in xlsx_files) {
    patient_id <- get_patient_id(xlsx)
    reads_file <- reads_files[grepl(patient_id, reads_files)]
    if (length(reads_file) == 0) next

    bmerged <- read_excel(xlsx)
    if (ncol(bmerged) > 2) bmerged <- bmerged[,-c(2,3)]
    raw_counts <- read_counts(reads_file)
    
    long_data <- prepare_long_data(bmerged, raw_counts, patient_id)
    all_data[[patient_id]] <- get_top10(long_data)
  }

  combined <- bind_rows(all_data)
  plot_stacked(combined, "Plots/facet_barplot_CPM_patients.png", "Top 10 Taxa per Sample (Faceted)", facet = TRUE)
  plot_stacked(combined, "Plots/facet_barplot_CPM_percentage_patients.png", "Top 10 Taxa per Sample (100% Stacked, Faceted)", percent = TRUE, facet = TRUE)
}

#------------------------------------------------#
#--- Optional: Metadata-Split Visualization -----#
#------------------------------------------------#

metadata <- read.csv2("metadata.txt", sep = "\t")
split_by_col <- as.character(config$split_by)

if (!(split_by_col %in% colnames(metadata))) {
  stop(paste("Column", split_by_col, "not found in metadata"))
}

metadata_split_list <- split(metadata, metadata[[split_by_col]])

for (split_name in names(metadata_split_list)) {
  samples_sub <- metadata_split_list[[split_name]]$sample
  bmerged_sub <- bmerged[, c("name", intersect(colnames(bmerged), samples_sub)), drop = FALSE]
  
  if (ncol(bmerged_sub) <= 1) {
    message(paste("Skipping", split_name, "- no matching sample columns"))
    next
  }

  long_data <- prepare_long_data(bmerged_sub, raw_counts)
  top10_split <- get_top10(long_data)

  plot_stacked(top10_split, paste0("Plots/barplot_CPM_", split_name, ".png"), 
               paste("Top 10 Species per Sample -", split_name))
  
  plot_stacked(top10_split, paste0("Plots/barplot_CPM_percentage_", split_name, ".png"),
               paste("Top 10 Taxa per Sample (100% Stacked) -", split_name), percent = TRUE)
}

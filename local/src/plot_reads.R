library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)   
library(yaml)

if (file.exists("Reports/reads_summary.csv")) {
  df <- read_tsv("Reports/reads_summary.csv")
  
  df$classified_kraken <- as.numeric(df$classified_kraken)
  df$Kraken_human_assigned_reads <- as.numeric(df$Kraken_human_assigned_reads)
  df$contaminant_reads <- as.numeric(df$contaminant_reads)
  df$initial_reads <- as.numeric(df$initial_reads)
  df$unclassified <- as.numeric(df$unmapped_reads - df$classified_kraken)
  
  plot_df <- df %>%
    mutate(
      non_human = classified_kraken - Kraken_human_assigned_reads - contaminant_reads,
      # Ensure non_human doesn't go negative
      non_human = pmax(non_human, 0),
      Kraken_human_assigned_reads = Kraken_human_assigned_reads / initial_reads * 1e6,
      contaminant_reads = contaminant_reads / initial_reads * 1e6,
      non_human = non_human / initial_reads * 1e6,
      non_classified = unclassified / initial_reads * 1e6
    ) %>%
    select(sample, Kraken_human_assigned_reads, contaminant_reads, non_human, non_classified) %>%
    pivot_longer(cols = c("Kraken_human_assigned_reads", "contaminant_reads", "non_human", "non_classified"),
                 names_to = "type", values_to = "count") %>%
    filter(count > 0) %>%
    mutate(type = factor(type, levels = c(
      "contaminant_reads",
      "non_classified",                
      "Kraken_human_assigned_reads",   
      "non_human"                    
                    
    )))
  
  p <- ggplot(plot_df, aes(x = sample, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = format(round(count), big.mark = ",", scientific = FALSE)),
              position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3) + # smaller text
    labs(x = "Sample", y = "Normalized Reads (per million initial reads)", fill = "Type",
         title = "Normalized Classified and Unclassified Kraken Reads: Human, Contaminant, and Assigned to taxa") +
    scale_fill_manual(
      values = c(
        "Kraken_human_assigned_reads" = "red", 
        "contaminant_reads" = "orange",
        "non_human" = "steelblue", 
        "non_classified" = "grey"
      ),
      labels = c(
        "Kraken_human_assigned_reads" = "Human",
        "contaminant_reads" = "Contaminant",
        "non_human" = "Assigned to taxa",
        "non_classified" = "Unclassified"
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1) # tilt x axis labels
    )
  
  ggsave(filename="Plots/classified_vs_human_contaminant_barplot_normalized.png", p, width = 25, height = 10)
  message("Created normalized classified vs. human vs. contaminant reads barplot")
}

######Plot_stacked_barplot_percentage######

plot_dfper <- df %>%
    mutate(
      non_human = classified_kraken - Kraken_human_assigned_reads - contaminant_reads,
      non_human = pmax(non_human, 0),
      Kraken_human_assigned_reads = Kraken_human_assigned_reads / initial_reads * 1e6,
      contaminant_reads = contaminant_reads / initial_reads * 1e6,
      non_human = non_human / initial_reads * 1e6, 
      non_classified = unclassified / initial_reads * 1e6
    ) %>%
    select(sample, Kraken_human_assigned_reads, contaminant_reads, non_human,non_classified) %>%
    pivot_longer(cols = c("Kraken_human_assigned_reads", "contaminant_reads", "non_human", "non_classified"),
                 names_to = "type", values_to = "count") %>%
    filter(count > 0) %>%
    mutate(type = factor(type, levels = c(
    "non_classified",                # Unclassified (bottom)
    "Kraken_human_assigned_reads",   # Human
    "non_human",                     # Assigned to taxa
    "contaminant_reads"              # Contaminant (top)
  ))) %>%
    group_by(sample) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  pp <- ggplot(plot_dfper, aes(x = sample, y = percentage, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 2.5) + # smaller text
    labs(x = "Sample", y = "Percentage of Normalized Reads", fill = "Type",
         title = "Percentage of Normalized Classified and Unclassified Kraken Reads per Sample") +
    scale_fill_manual(
      values = c(
        "Kraken_human_assigned_reads" = "red", 
        "contaminant_reads" = "orange",
        "non_human" = "steelblue", 
        "non_classified" = "grey"
      ),
      labels = c(
        "Kraken_human_assigned_reads" = "Human",
        "contaminant_reads" = "Contaminant",
        "non_human" = "Assigned to taxa",
        "non_classified" = "Unclassified"
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1) # tilt x axis labels
    )
  
  ggsave(filename="Plots/classified_vs_human_contaminant_barplot_percentage.png", pp, width = 25, height = 10)
  message("Created 100% stacked barplot of normalized classified vs. human vs. contaminant reads")






##################################### Split plot by metadata #####################################
config <- yaml.load_file("config.yaml")
split_by = config$split_by

# Read metadata and join with main dataframe
metadata <- read.csv2("metadata.txt", sep="\t")
  df <- read_tsv("Reports/reads_summary.csv")
  
  df$classified_kraken <- as.numeric(df$classified_kraken)
  df$Kraken_human_assigned_reads <- as.numeric(df$Kraken_human_assigned_reads)
  df$contaminant_reads <- as.numeric(df$contaminant_reads)
  df$initial_reads <- as.numeric(df$initial_reads)
  df$unclassified <- as.numeric(df$unmapped_reads - df$classified_kraken)
if (!split_by %in% colnames(metadata)) {
  stop(paste("Column", split_by, "not found in metadata.txt"))
}
metadata_nodup <- metadata[, setdiff(colnames(metadata), colnames(df)), drop = FALSE]

metadata_nodup$sample <- metadata$sample 

df <- left_join(df, metadata_nodup, by = "sample")


plot_df <- df %>%
  mutate(
    non_human = classified_kraken - Kraken_human_assigned_reads - contaminant_reads,
    non_human = pmax(non_human, 0),
    Kraken_human_assigned_reads = Kraken_human_assigned_reads / initial_reads * 1e6,
    contaminant_reads = contaminant_reads / initial_reads * 1e6,
    non_human = non_human / initial_reads * 1e6,
    non_classified = unclassified / initial_reads * 1e6
  ) %>%
  select(sample, !!split_by, Kraken_human_assigned_reads, contaminant_reads, non_human, non_classified) %>%
  pivot_longer(cols = c("Kraken_human_assigned_reads", "contaminant_reads", "non_human", "non_classified"),
               names_to = "type", values_to = "count") %>%
  filter(count > 0) %>%
  mutate(type = factor(type, levels = c(
    "contaminant_reads",
    "non_classified",                
    "Kraken_human_assigned_reads",   
    "non_human"                    
  )))

plot_dfper <- plot_df %>%
  group_by(sample) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Split and plot by split_by variable
groups <- unique(df[[split_by]])
for (grp in groups) {
  grp_plot_df <- plot_df %>% filter(.data[[split_by]] == grp)
  grp_plot_dfper <- plot_dfper %>% filter(.data[[split_by]] == grp)
  
  # Normalized barplot
  p_grp <- ggplot(grp_plot_df, aes(x = sample, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = format(round(count), big.mark = ",", scientific = FALSE)),
              position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 3) +
    labs(x = "Sample", y = "Normalized Reads (per million initial reads)", fill = "Type",
         title = paste("Normalized Reads by", split_by, "=", grp)) +
    scale_fill_manual(
      values = c(
        "Kraken_human_assigned_reads" = "red", 
        "contaminant_reads" = "orange",
        "non_human" = "steelblue", 
        "non_classified" = "grey"
      ),
      labels = c(
        "Kraken_human_assigned_reads" = "Human",
        "contaminant_reads" = "Contaminant",
        "non_human" = "Assigned to taxa",
        "non_classified" = "Unclassified"
      )
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(filename = paste0("Plots/classified_vs_human_contaminant_barplot_normalized_", split_by, "_", grp, ".png"),
         p_grp, width = 25, height = 10)
  
  # Percentage barplot
  pp_grp <- ggplot(grp_plot_dfper, aes(x = sample, y = percentage, fill = type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 2.5) +
    labs(x = "Sample", y = "Percentage of Normalized Reads", fill = "Type",
         title = paste("Percentage of Normalized Reads by", split_by, "=", grp)) +
    scale_fill_manual(
      values = c(
        "Kraken_human_assigned_reads" = "red", 
        "contaminant_reads" = "orange",
        "non_human" = "steelblue", 
        "non_classified" = "grey"
      ),
      labels = c(
        "Kraken_human_assigned_reads" = "Human",
        "contaminant_reads" = "Contaminant",
        "non_human" = "Assigned to taxa",
        "non_classified" = "Unclassified"
      )
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(filename = paste0("Plots/classified_vs_human_contaminant_barplot_percentage_", split_by, "_", grp, ".png"),
         pp_grp, width = 25, height = 10)
}


#!/usr/bin/env Rscript
# Exploratory QC for mosquito Bracken abundance profiles.
#
# Designed for Snakemake's script: directive, with named inputs:
#   clr, relative, metadata, totals
# Optional input:
#   counts = "qc/counts_all.taxid_collapsed.txt"
#   Strongly recommended: use original, non-pseudocounted counts for richness,
#   detected-taxon counts, zero fractions and diversity estimates.
# and named outputs:
#   pca, heatmap, summary
# Recommended additional named outputs:
#   group_summary, pca_coordinates
# Optional Snakemake params:
#   condition_column = "condition"
#   group_column     = "site"  # choose an actual metadata column
#   batch_column     = "batch"
#   outlier_z        = 3.5
#
# No external R packages are required.
# This script suggests candidates for QC review; it DOES NOT remove samples.
#
# Standalone usage (optional):
# Rscript explore_abundance_qc.R CLR.tsv RELATIVE.tsv metadata.txt \
#   sample_totals.txt pca.pdf sample_distances.pdf sample_qc_summary.tsv \
#   [group_qc_summary.tsv] [pca_coordinates.tsv] [counts.tsv]

fail <- function(...) stop(paste0(...), call. = FALSE)
info <- function(...) message("[explore_abundance_qc] ", paste0(...))

option_value <- function(name, default = NULL) {
  if (exists("snakemake", inherits = TRUE)) {
    x <- snakemake@params[[name]]
    if (!is.null(x) && length(x) > 0 && !is.na(x[1])) return(x[[1]])
  }
  default
}

if (exists("snakemake", inherits = TRUE)) {
  path_in <- function(name) {
    val <- snakemake@input[[name]]
    if (is.null(val) || length(val) != 1L) fail("Missing Snakemake input: ", name)
    as.character(val)
  }
  path_out <- function(name, fallback = NULL) {
    val <- snakemake@output[[name]]
    if (is.null(val) || length(val) == 0L) return(fallback)
    as.character(val[[1]])
  }
  counts_file <- snakemake@input[["counts"]]
  inputs <- list(
    clr = path_in("clr"), relative = path_in("relative"),
    metadata = path_in("metadata"), totals = path_in("totals"),
    counts = if (is.null(counts_file)) NULL else as.character(counts_file[[1]])
  )
  outputs <- list(
    pca = path_out("pca"), heatmap = path_out("heatmap"),
    summary = path_out("summary"),
    group_summary = path_out("group_summary"),
    pca_coordinates = path_out("pca_coordinates")
  )
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 7L) {
    fail("Expected >=7 positional arguments (see script header), or use Snakemake.")
  }
  inputs <- setNames(as.list(args[1:4]), c("clr", "relative", "metadata", "totals"))
  inputs$counts <- if (length(args) >= 10L) args[10] else NULL
  outputs <- list(
    pca = args[5], heatmap = args[6], summary = args[7],
    group_summary = if (length(args) >= 8L) args[8] else NULL,
    pca_coordinates = if (length(args) >= 9L) args[9] else NULL
  )
}

# Record progress in the Snakemake log file, when one is provided.
log_file <- NULL
if (exists("snakemake", inherits = TRUE) && length(snakemake@log) > 0L) {
  log_file <- as.character(snakemake@log[[1]])
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  writeLines("Exploratory abundance QC started", con = log_file)
}
info <- function(...) {
  message_text <- paste0("[explore_abundance_qc] ", paste0(...))
  message(message_text)
  if (!is.null(log_file)) cat(message_text, "\n", file = log_file, append = TRUE)
}

if (is.null(outputs$group_summary)) {
  outputs$group_summary <- file.path(dirname(outputs$summary), "group_qc_summary.tsv")
}
if (is.null(outputs$pca_coordinates)) {
  outputs$pca_coordinates <- file.path(dirname(outputs$summary), "pca_coordinates.tsv")
}
for (filename in unlist(inputs, use.names = FALSE)) {
  if (!file.exists(filename)) fail("Input does not exist: ", filename)
}
for (filename in unlist(outputs, use.names = FALSE)) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
}

# Base R parser: preserve sample-column names exactly as they appear in the files.
read_tsv_matrix <- function(path, label) {
  tab <- tryCatch(
    read.delim(path, header = TRUE, sep = "\t", check.names = FALSE,
               stringsAsFactors = FALSE, quote = "\"", comment.char = ""),
    error = function(e) fail("Cannot read ", label, " matrix: ", conditionMessage(e))
  )
  if (ncol(tab) < 3L) {
    fail(label, " matrix must have a feature-ID column and >=2 sample columns: ", path)
  }
  ids <- trimws(as.character(tab[[1]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    fail(label, " matrix contains missing or duplicate feature identifiers: ", path)
  }
  sample_ids <- colnames(tab)[-1]
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    fail(label, " matrix contains missing or duplicate sample column names: ", path)
  }
  values <- as.matrix(tab[-1])
  suppressWarnings(storage.mode(values) <- "double")
  if (any(!is.finite(values))) {
    fail(label, " matrix contains nonnumeric, missing, or infinite values: ", path)
  }
  rownames(values) <- ids
  colnames(values) <- sample_ids
  values
}

# Metadata may be tab-, comma-, or whitespace-delimited. CSV/TSV are preferred.
read_metadata <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  if (length(first) == 0L) fail("Empty metadata file: ", path)
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else
    if (grepl(",", first, fixed = TRUE)) "," else ""
  dat <- tryCatch(
    read.table(path, header = TRUE, sep = sep, check.names = FALSE,
               stringsAsFactors = FALSE, quote = "\"", comment.char = "",
               fill = FALSE),
    error = function(e) fail("Cannot read sample metadata: ", conditionMessage(e))
  )
  if (ncol(dat) < 1L || nrow(dat) < 1L) fail("Metadata is empty: ", path)
  dat
}

column_matching <- function(df, candidates) {
  names_low <- tolower(trimws(names(df)))
  for (candidate in candidates) {
    pos <- match(tolower(candidate), names_low)
    if (!is.na(pos)) return(names(df)[pos])
  }
  NULL
}

locate_sample_column <- function(metadata, samples) {
  preferred <- option_value("sample_column", NULL)
  if (!is.null(preferred)) {
    col <- column_matching(metadata, as.character(preferred))
    if (is.null(col)) fail("Configured sample_column not found in metadata: ", preferred)
    if (!all(samples %in% as.character(metadata[[col]]))) {
      fail("Configured sample_column does not contain all matrix samples: ", col)
    }
    return(col)
  }
  candidates <- c("Sample", "sample", "sample_id", "sampleid", "sample_name", "ID", "id")
  candidates <- unique(c(column_matching(metadata, candidates), names(metadata)))
  for (col in candidates) {
    v <- trimws(as.character(metadata[[col]]))
    if (!anyDuplicated(v) && all(samples %in% v)) return(col)
  }
  fail("No unique metadata column contains all matrix sample IDs. ",
       "Check spelling/case or set params.sample_column.")
}

resolve_attribute <- function(metadata, explicit, candidates, mandatory = FALSE) {
  if (!is.null(explicit) && nzchar(as.character(explicit))) {
    col <- column_matching(metadata, as.character(explicit))
    if (is.null(col)) fail("Requested metadata column not found: ", explicit)
    return(col)
  }
  col <- column_matching(metadata, candidates)
  if (is.null(col) && mandatory) fail("Cannot find metadata column among: ", paste(candidates, collapse = ", "))
  col
}

safe_label <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values) | !nzchar(values)] <- "Not recorded"
  values
}

robust_z <- function(x, min_n = 5L) {
  x <- as.numeric(x)
  result <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < min_n) return(result)
  med <- median(x[ok])
  dispersion <- stats::mad(x[ok], center = med, constant = 1.4826)
  if (!is.finite(dispersion) || dispersion <= 1e-12) return(result)
  result[ok] <- (x[ok] - med) / dispersion
  result
}

save_table <- function(data, filename) {
  write.table(data, filename, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "NA")
}

clr <- read_tsv_matrix(inputs$clr, "CLR")
relative <- read_tsv_matrix(inputs$relative, "Relative abundance")
samples <- colnames(clr)
if (!setequal(samples, colnames(relative))) {
  fail("CLR and relative abundance matrices have different sample IDs. ",
       "Only in CLR: ", paste(setdiff(samples, colnames(relative)), collapse = ", "), 
       "; only in relative: ", paste(setdiff(colnames(relative), samples), collapse = ", "))
}
if (!setequal(rownames(clr), rownames(relative))) {
  fail("CLR and relative matrices contain different taxon/feature IDs. ",
       "Both inputs must be derived from the same count matrix.")
}
relative <- relative[rownames(clr), samples, drop = FALSE]
if (any(relative < 0)) fail("Relative abundance matrix contains negative values.")
if (ncol(clr) < 3L) fail("At least 3 samples are required for exploratory QC.")

# When available, use pre-normalization Bracken counts for ecological summaries.
# Normalization may add pseudocounts, making all taxa appear to be present.
composition <- relative
composition_source <- "normalized_relative_abundance"
if (!is.null(inputs$counts)) {
  counts <- read_tsv_matrix(inputs$counts, "Pre-normalization count")
  if (!setequal(rownames(counts), rownames(clr)) ||
      !setequal(colnames(counts), samples)) {
    fail("Pre-normalization counts must contain the same features and samples ",
         "as the CLR and relative-abundance matrices.")
  }
  counts <- counts[rownames(clr), samples, drop = FALSE]
  if (any(counts < 0)) fail("Pre-normalization count matrix contains negative counts.")
  composition <- counts
  composition_source <- "pre_normalization_bracken_counts"
} else {
  info("WARNING: No optional 'counts' input supplied. Richness, zeros, and ",
       "diversity use the relative matrix; if pseudocounts were added during ",
       "normalization, these values do not measure actual taxon detection.")
}

metadata <- read_metadata(inputs$metadata)
sample_col <- locate_sample_column(metadata, samples)
metadata[[sample_col]] <- trimws(as.character(metadata[[sample_col]]))
if (anyDuplicated(metadata[[sample_col]])) {
  fail("Metadata has duplicate sample IDs in column ", sample_col)
}
extra_meta <- setdiff(metadata[[sample_col]], samples)
if (length(extra_meta)) info("Ignoring ", length(extra_meta), " metadata samples not in matrices.")
metadata <- metadata[match(samples, metadata[[sample_col]]), , drop = FALSE]
if (anyNA(metadata[[sample_col]])) fail("Metadata and abundance matrix sample IDs do not match.")

condition_col <- resolve_attribute(
  metadata, option_value("condition_column", NULL),
  c("condition", "phenotype", "status", "resistance_status")
)
group_col <- resolve_attribute(
  metadata, option_value("group_column", NULL),
  c("collection_group", "site", "collection_site", "location", "group", "batch")
)
if (is.null(group_col)) group_col <- condition_col
batch_col <- resolve_attribute(metadata, option_value("batch_column", NULL), c("batch", "run", "sequencing_batch"))

condition <- if (!is.null(condition_col)) safe_label(metadata[[condition_col]]) else rep("Not supplied", length(samples))
group <- if (!is.null(group_col)) safe_label(metadata[[group_col]]) else rep("All samples", length(samples))
batch <- if (!is.null(batch_col)) safe_label(metadata[[batch_col]]) else rep("Not supplied", length(samples))
if (is.null(group_col)) info("No group/site column detected; group outlier assessment unavailable.")
info("Matrix: ", nrow(clr), " taxa x ", length(samples), " samples")
info("Metadata sample column: ", sample_col,
     "; condition: ", if (is.null(condition_col)) "none" else condition_col,
     "; group: ", if (is.null(group_col)) "none" else group_col)

# qc/sample_totals.txt is produced by summing Bracken 'new_est_reads'.
# It is NOT the total number of sequenced FASTQ reads.
totals <- tryCatch(
  read.table(inputs$totals, header = TRUE, sep = "", check.names = FALSE,
             stringsAsFactors = FALSE, comment.char = "", quote = "\""),
  error = function(e) fail("Cannot read Bracken sample totals: ", conditionMessage(e))
)
if (ncol(totals) < 2L) fail("Totals must contain sample ID and total reads columns.")
totals_col_sample <- column_matching(totals, c("Sample", "sample", "sample_id"))
if (is.null(totals_col_sample)) totals_col_sample <- names(totals)[1]
totals_col_reads <- column_matching(totals, c("Total_Reads", "total_reads", "total", "reads"))
if (is.null(totals_col_reads)) totals_col_reads <- names(totals)[2]
totals_ids <- trimws(as.character(totals[[totals_col_sample]]))
if (anyDuplicated(totals_ids)) fail("Duplicate sample identifiers in Bracken totals file.")
if (!all(samples %in% totals_ids)) {
  fail("Bracken totals missing these samples: ", paste(setdiff(samples, totals_ids), collapse = ", "))
}
bracken_reads <- suppressWarnings(as.numeric(totals[[totals_col_reads]][match(samples, totals_ids)]))
if (any(!is.finite(bracken_reads)) || any(bracken_reads < 0)) {
  fail("Bracken totals contain missing, nonnumeric, infinite, or negative counts.")
}

# Shannon/Simpson are calculated for taxa included in the input matrix.
# They are not estimates of diversity of every organism in the original sample.
compositional_sums <- colSums(composition)
proportions <- composition
nonzero <- compositional_sums > 0
if (any(nonzero)) {
  proportions[, nonzero] <- sweep(composition[, nonzero, drop = FALSE], 2L,
                                 compositional_sums[nonzero], "/")
}
if (any(!nonzero)) proportions[, !nonzero] <- 0
log_input <- proportions
log_input[log_input <= 0] <- 1
richness <- colSums(composition > 0)
shannon <- colSums(-proportions * log(log_input))
simpson <- 1 - colSums(proportions^2)
zero_fraction <- colMeans(composition == 0)

# Aitchison distance = Euclidean distance between CLR vectors. Do not scale
# features independently; doing so would change the Aitchison geometry.
aitchison <- as.matrix(stats::dist(t(clr), method = "euclidean"))
dimnames(aitchison) <- list(samples, samples)
median_distance <- vapply(seq_along(samples), function(i) {
  median(aitchison[i, -i])
}, numeric(1))

# PCA using taxa with nonzero variance (constant taxa add no PCA information).
feature_variance <- rowSums((clr - rowMeans(clr))^2)
varying <- feature_variance > 1e-12
if (!any(varying)) fail("No variation between samples in CLR matrix; PCA is undefined.")
pca <- stats::prcomp(t(clr[varying, , drop = FALSE]), center = TRUE,
                     scale. = FALSE)
variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pc1 <- pca$x[, 1L]
pc2 <- if (ncol(pca$x) >= 2L) pca$x[, 2L] else rep(0, length(samples))
names(pc2) <- samples

outlier_threshold <- suppressWarnings(as.numeric(option_value("outlier_z", 3.5)))
if (length(outlier_threshold) != 1L || !is.finite(outlier_threshold) || outlier_threshold <= 0) {
  fail("params.outlier_z must be a positive number.")
}
z_distance <- robust_z(median_distance)
z_depth <- robust_z(log10(bracken_reads + 1))
z_richness <- robust_z(richness)
flag_distance <- !is.na(z_distance) & z_distance > outlier_threshold
flag_depth <- !is.na(z_depth) & z_depth < -outlier_threshold
flag_richness <- !is.na(z_richness) & z_richness < -outlier_threshold
flag_zero <- compositional_sums <= 0
review <- flag_distance | flag_depth | flag_richness | flag_zero
reasons <- vapply(seq_along(samples), function(i) {
  r <- character(0)
  if (flag_distance[i]) r <- c(r, "unusual_composition")
  if (flag_depth[i]) r <- c(r, "low_bracken_estimated_reads")
  if (flag_richness[i]) r <- c(r, "low_detected_richness")
  if (flag_zero[i]) r <- c(r, "zero_relative_abundance_total")
  if (length(r)) paste(r, collapse = ";") else ""
}, character(1))

# Group-level exploration is only a screening tool, never grounds for automatic
# exclusion. Avoid a numerical outlier claim when fewer than FIVE groups exist.
unique_groups <- unique(group)
group_members <- lapply(unique_groups, function(g) which(group == g))
names(group_members) <- unique_groups
group_centroids <- vapply(group_members, function(idx) {
  rowMeans(clr[, idx, drop = FALSE])
}, numeric(nrow(clr)))
if (!is.matrix(group_centroids)) {
  group_centroids <- matrix(group_centroids, ncol = length(unique_groups))
}
colnames(group_centroids) <- unique_groups
centroid_d <- if (length(unique_groups) > 1L) {
  as.matrix(stats::dist(t(group_centroids)))
} else matrix(0, 1L, 1L, dimnames = list(unique_groups, unique_groups))

between_group_median <- vapply(seq_along(unique_groups), function(i) {
  if (length(unique_groups) <= 1L) return(NA_real_)
  median(centroid_d[i, -i])
}, numeric(1))
group_z <- if (length(unique_groups) >= 5L) robust_z(between_group_median) else rep(NA_real_, length(unique_groups))
within_group_median <- vapply(group_members, function(idx) {
  if (length(idx) < 2L) return(NA_real_)
  d <- aitchison[idx, idx, drop = FALSE]
  median(d[upper.tri(d)])
}, numeric(1))
group_sizes <- lengths(group_members)
group_review <- !is.na(group_z) & group_z > outlier_threshold & group_sizes >= 2L

group_reason <- vapply(seq_along(unique_groups), function(i) {
  if (length(unique_groups) < 5L) return("insufficient_number_of_groups_for_robust_flagging")
  if (group_sizes[i] < 2L) return("single_sample_group;cannot_assess_replication")
  if (is.na(group_z[i])) return("robust_group_score_unavailable")
  if (group_review[i]) return("unusual_group_centroid;investigate_batch_or_biology")
  ""
}, character(1))

group_table <- data.frame(
  group_column = if (is.null(group_col)) "Not defined" else group_col,
  group = unique_groups,
  n_samples = as.integer(group_sizes),
  median_bracken_estimated_reads = vapply(group_members, function(idx) median(bracken_reads[idx]), numeric(1)),
  median_taxa_detected = vapply(group_members, function(idx) median(richness[idx]), numeric(1)),
  median_within_group_aitchison_distance = as.numeric(within_group_median),
  median_centroid_distance_to_other_groups = between_group_median,
  robust_z_centroid_distance = group_z,
  group_review_flag = group_review,
  group_review_reason = group_reason,
  stringsAsFactors = FALSE, check.names = FALSE
)

sample_table <- data.frame(
  sample = samples,
  condition = condition,
  group = group,
  batch = batch,
  bracken_estimated_reads = bracken_reads,
  taxa_detected = as.integer(richness),
  shannon = shannon,
  simpson = simpson,
  fraction_zero_taxa = zero_fraction,
  compositional_summary_source = composition_source,
  relative_abundance_column_sum = colSums(relative),
  median_aitchison_distance_to_other_samples = median_distance,
  robust_z_aitchison_distance = z_distance,
  robust_z_log_bracken_reads = z_depth,
  robust_z_richness = z_richness,
  PC1 = as.numeric(pc1[samples]),
  PC2 = as.numeric(pc2[samples]),
  sample_review_flag = review,
  sample_review_reason = reasons,
  group_review_flag = group_review[match(group, unique_groups)],
  stringsAsFactors = FALSE, check.names = FALSE
)

pca_table <- data.frame(
  sample = samples, condition = condition, group = group, batch = batch,
  PC1 = as.numeric(pc1[samples]), PC2 = as.numeric(pc2[samples]),
  PC1_variance_percent = variance_explained[1],
  PC2_variance_percent = if (length(variance_explained) >= 2L) variance_explained[2] else 0,
  sample_review_flag = review,
  stringsAsFactors = FALSE, check.names = FALSE
)

save_table(sample_table, outputs$summary)
save_table(group_table, outputs$group_summary)
save_table(pca_table, outputs$pca_coordinates)

# Consistent colors per factor; base R only.
condition_levels <- unique(condition)
group_levels <- unique(group)
condition_colors <- setNames(grDevices::hcl.colors(length(condition_levels), "Dark 3"), condition_levels)
group_colors <- setNames(grDevices::hcl.colors(length(group_levels), "Set 2"), group_levels)

plot_pca_page <- function(coloring, palette, heading, legend_heading) {
  colors <- unname(palette[coloring])
  x_label <- sprintf("PC1 (%.1f%%)", variance_explained[1])
  y_label <- sprintf("PC2 (%.1f%%)", if (length(variance_explained) >= 2L) variance_explained[2] else 0)
  graphics::par(mar = c(5, 5, 4, 2) + 0.1)
  graphics::plot(pc1, pc2, type = "n", xlab = x_label, ylab = y_label,
                 main = heading, cex.main = 1.1)
  graphics::abline(h = 0, v = 0, col = "grey85", lty = 3)
  graphics::points(pc1, pc2, pch = ifelse(review, 17, 19),
                   col = colors, cex = ifelse(review, 1.7, 1.3))
  if (length(samples) <= 60L) {
    graphics::text(pc1, pc2, labels = samples, pos = 3, cex = 0.62, col = "grey20")
  } else if (any(review)) {
    graphics::text(pc1[review], pc2[review], labels = samples[review], pos = 3, cex = 0.65)
  }
  legend_labels <- names(palette)
  # Place a small legend outside the plotting area when there are many groups.
  graphics::legend("topright", legend = legend_labels, col = unname(palette),
                   pch = 19, cex = if (length(legend_labels) > 10) 0.55 else 0.75,
                   bty = "n", title = legend_heading)
  graphics::mtext("Triangles = samples flagged for QC review; not automatically excluded",
                  side = 1, line = 3.5, cex = 0.7)
}

render_pca <- function(filename) {
  grDevices::pdf(filename, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_pca_page(condition, condition_colors, "CLR PCA: coloured by condition",
                if (is.null(condition_col)) "Condition" else condition_col)
  if (!identical(group, condition)) {
    plot_pca_page(group, group_colors, "CLR PCA: coloured by QC group",
                  if (is.null(group_col)) "Group" else group_col)
  }
  graphics::par(mar = c(5, 5, 4, 2) + 0.1)
  ncomp <- min(15L, length(variance_explained))
  graphics::barplot(variance_explained[seq_len(ncomp)],
                    names.arg = paste0("PC", seq_len(ncomp)),
                    col = "steelblue", xlab = "Principal component",
                    ylab = "Variance explained (%)",
                    main = "CLR PCA: variance explained")
  graphics::plot(log10(bracken_reads + 1), richness,
                 pch = ifelse(review, 17, 19),
                 col = unname(condition_colors[condition]),
                 cex = ifelse(review, 1.7, 1.3),
                 xlab = "log10(Bracken estimated reads + 1)",
                 ylab = "Number of detected taxa",
                 main = "Bracken coverage versus observed richness")
  graphics::text(log10(bracken_reads + 1), richness, labels = samples,
                 pos = 3, cex = 0.65)
  graphics::legend("topleft", legend = names(condition_colors),
                   col = unname(condition_colors), pch = 19, bty = "n", cex = 0.75)
}

render_heatmap <- function(filename) {
  grDevices::pdf(filename, width = 11, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  dendro <- as.dendrogram(stats::hclust(stats::as.dist(aitchison), method = "average"))
  heat_colors <- grDevices::colorRampPalette(c("#ffffe5", "#fd8d3c", "#800026"))(100)
  colors <- unname(condition_colors[condition])
  stats::heatmap(
    aitchison, Rowv = dendro, Colv = dendro, symm = TRUE,
    scale = "none", col = heat_colors, margins = c(10, 10),
    RowSideColors = colors, ColSideColors = colors,
    main = "Sample distances (Euclidean on CLR = Aitchison)",
    xlab = "Sample", ylab = "Sample"
  )
  # Second page: compare QC-group centroid profiles when groups are defined.
  if (length(unique_groups) >= 3L) {
    gr_dendro <- as.dendrogram(stats::hclust(stats::as.dist(centroid_d),
                                             method = "average"))
    stats::heatmap(
      centroid_d, Rowv = gr_dendro, Colv = gr_dendro, symm = TRUE,
      scale = "none", col = heat_colors, margins = c(10, 10),
      main = "Distances between QC-group CLR centroids",
      xlab = "Group", ylab = "Group"
    )
  }
}

render_pca(outputs$pca)
render_heatmap(outputs$heatmap)

info("Generated PCA: ", outputs$pca)
info("Generated Aitchison-distance heatmap: ", outputs$heatmap)
info("Sample summary: ", outputs$summary)
info("Group summary: ", outputs$group_summary)
info("PCA scores: ", outputs$pca_coordinates)
info("Sample candidates for review: ", sum(review), "/", length(review))
info("Groups flagged for review: ", sum(group_review), "/", length(group_review),
     " (requires >=5 distinct groups and >=2 samples per flagged group)")
info("Review flags are descriptive screening indicators. No exclusions were performed.")

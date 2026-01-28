#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(ggplot2)
  library(grid)
})

# -----------------------------
# minimal CLI parsing
# -----------------------------
get_flag_value <- function(flag, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

cfg_path <- get_flag_value("--config", NA_character_)
out_pdf  <- get_flag_value("--out_pdf", NA_character_)
out_png  <- get_flag_value("--out_png", NA_character_)
if (is.na(cfg_path)) stop("Missing --config config.yaml")

cfg <- yaml::read_yaml(cfg_path)
if (is.null(cfg$VOLCANO)) stop("Missing VOLCANO block in config.yaml")
vc <- cfg$VOLCANO

if (is.na(out_pdf)) out_pdf <- vc$OUTPUT$pdf
if (is.na(out_png)) out_png <- vc$OUTPUT$png

dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# IO
# -----------------------------
read_any <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    fread(cmd = sprintf("zcat %s", shQuote(path)), sep = "\t", header = TRUE, data.table = FALSE)
  } else {
    fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  }
}

stop_missing_cols <- function(df, cols, where = "") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required column(s)%s: %s\nAvailable columns: %s",
      ifelse(where == "", "", paste0(" in ", where)),
      paste(missing, collapse = ", "),
      paste(colnames(df), collapse = ", ")
    ))
  }
}

# -----------------------------
# COLORS parsing (config-driven, no hardcoding)
# -----------------------------
colors_raw <- vc$COLORS
if (is.null(colors_raw) || is.null(names(colors_raw)) || any(names(colors_raw) == "")) {
  stop("VOLCANO.COLORS must be a YAML mapping with names (e.g. COE: 'red').")
}
colors <- setNames(
  vapply(names(colors_raw), function(k) as.character(colors_raw[[k]]), character(1)),
  names(colors_raw)
)

other_label <- if (!is.null(vc$OTHER_LABEL)) as.character(vc$OTHER_LABEL) else "Other"
if (!(other_label %in% names(colors))) {
  stop(sprintf("VOLCANO.COLORS must include OTHER_LABEL '%s' (e.g. Other: 'grey70').", other_label))
}

# Optional order (if you want legend/factor order controlled from config)
if (!is.null(vc$COLOR_ORDER)) {
  ord <- as.character(unlist(vc$COLOR_ORDER))
  ord <- ord[ord %in% names(colors)]
  ord <- c(ord, setdiff(names(colors), ord))
  colors <- colors[ord]
}

# -----------------------------
# Read data
# -----------------------------
all_df <- read_any(vc$ALL_CONTRASTS_FILE)
fam_df <- read_any(vc$FAMILY_MAP_FILE)

cols <- vc$COLUMNS
need_all <- c(cols$contrast, cols$gene, cols$log2fc, cols$fdr)
stop_missing_cols(all_df, need_all, where = vc$ALL_CONTRASTS_FILE)
stop_missing_cols(fam_df, c("GeneID", "family"), where = vc$FAMILY_MAP_FILE)

contrasts_keep <- as.character(unlist(vc$CONTRASTS))
if (length(contrasts_keep) != 5) stop("VOLCANO.CONTRASTS must contain exactly 5 contrasts for the 5-panel layout.")

layout_order <- as.character(unlist(vc$LAYOUT$order))
if (length(layout_order) != 5) stop("VOLCANO.LAYOUT.order must contain exactly 5 contrasts.")
if (!all(layout_order %in% contrasts_keep)) stop("VOLCANO.LAYOUT.order must be a subset of VOLCANO.CONTRASTS.")

titles <- NULL
if (!is.null(vc$TITLES) && !is.null(names(vc$TITLES))) {
  titles <- setNames(
    vapply(names(vc$TITLES), function(k) as.character(vc$TITLES[[k]]), character(1)),
    names(vc$TITLES)
  )
}

# Join and normalize
df <- data.frame(
  contrast = all_df[[cols$contrast]],
  GeneID   = all_df[[cols$gene]],
  log2fc   = as.numeric(all_df[[cols$log2fc]]),
  fdr      = as.numeric(all_df[[cols$fdr]]),
  stringsAsFactors = FALSE
)
df <- df[df$contrast %in% contrasts_keep, , drop = FALSE]
df$contrast <- factor(df$contrast, levels = layout_order)

fam_df <- fam_df[, c("GeneID", "family")]
df <- merge(df, fam_df, by = "GeneID", all.x = TRUE)

df$family <- ifelse(is.na(df$family) | df$family == "", other_label, df$family)
df$family <- ifelse(df$family %in% names(colors), df$family, other_label)
df$family <- factor(df$family, levels = names(colors))

min_fdr <- as.numeric(vc$OUTPUT$min_fdr)
y_cap   <- as.numeric(vc$OUTPUT$y_cap)
df$fdr  <- pmax(df$fdr, min_fdr)
df$neglog10_fdr <- pmin(-log10(df$fdr), y_cap)

# -----------------------------
# ggplot + legend extraction
# -----------------------------
extract_legend_grob <- function(p) {
  g <- ggplotGrob(p)
  # Be tolerant across ggplot2 versions: match anything containing "guide"
  idx <- which(grepl("guide", sapply(g$grobs, function(x) x$name), ignore.case = TRUE))
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

make_volcano <- function(d, title, show_legend = FALSE) {
  p <- ggplot(d, aes(x = log2fc, y = neglog10_fdr)) +
    geom_point(aes(color = family), size = as.numeric(vc$OUTPUT$point_size)) +
    geom_hline(yintercept = -log10(as.numeric(vc$PADJ_CUTOFF)), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "solid") +
    geom_vline(xintercept = c(-as.numeric(vc$LFC_CUTOFF), as.numeric(vc$LFC_CUTOFF)), linetype = "dashed") +
    scale_color_manual(values = colors, breaks = names(colors), limits = names(colors), drop = FALSE) +
    scale_y_continuous(limits = c(0, y_cap), breaks = seq(0, y_cap, by = 10)) +
    labs(
      title = title,
      x = expression(log[2](FC)),
      y = expression(-log[10](FDR)),
      color = NULL
    ) +
    theme_bw(base_size = as.numeric(vc$OUTPUT$base_font_size)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )
  p
}

# Build the 5 plots (in layout order)
plots <- list()
for (cn in layout_order) {
  title <- if (!is.null(titles) && !is.null(titles[[cn]])) titles[[cn]] else cn
  plots[[cn]] <- make_volcano(df[df$contrast == cn, , drop = FALSE], title, show_legend = FALSE)
}

# Build ONE plot with legend, just to extract it (same scales/colors)
p_for_legend <- make_volcano(df[df$contrast == layout_order[[1]], , drop = FALSE],
                             title = if (!is.null(titles) && !is.null(titles[[layout_order[[1]]]])) titles[[layout_order[[1]]]] else layout_order[[1]],
                             show_legend = TRUE)

p_for_legend <- p_for_legend +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    legend.text = element_text(size = as.numeric(vc$OUTPUT$base_font_size)),
    legend.key.height = unit(1.2, "lines"),
    legend.key.width  = unit(1.2, "lines")
  )
  
legend_grob <- extract_legend_grob(p_for_legend)
if (is.null(legend_grob)) {
  stop("Could not extract ggplot2 legend grob. (No 'guide' grob found).")
}

g1 <- ggplotGrob(plots[[layout_order[[1]]]])
g2 <- ggplotGrob(plots[[layout_order[[2]]]])
g3 <- ggplotGrob(plots[[layout_order[[3]]]])
g4 <- ggplotGrob(plots[[layout_order[[4]]]])
g5 <- ggplotGrob(plots[[layout_order[[5]]]])

rel_widths <- unlist(vc$LAYOUT$rel_widths)
if (length(rel_widths) != 3) stop("VOLCANO.LAYOUT.rel_widths must have length 3.")

lay <- grid.layout(
  nrow = 2, ncol = 3,
  widths = unit(rel_widths, "null"),
  heights = unit(c(1, 1), "null")
)

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(layout = lay))

  draw_cell <- function(grob, r, c) {
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))
    grid.draw(grob)
    popViewport()
  }

  draw_cell(g1, 1, 1)
  draw_cell(g2, 1, 2)
  draw_cell(g3, 1, 3)
  draw_cell(g4, 2, 1)
  draw_cell(g5, 2, 2)
  draw_cell(legend_grob, 2, 3)

  popViewport()
}

# Save PDF/PNG
pdf(out_pdf, width = as.numeric(vc$OUTPUT$width_in), height = as.numeric(vc$OUTPUT$height_in), onefile = TRUE)
draw_all()
dev.off()

png(out_png,
    width = as.numeric(vc$OUTPUT$width_in),
    height = as.numeric(vc$OUTPUT$height_in),
    units = "in",
    res = as.integer(vc$OUTPUT$dpi))
draw_all()
dev.off()

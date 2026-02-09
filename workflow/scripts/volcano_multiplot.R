#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(ggplot2)
  library(grid)
  library(RColorBrewer)
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
    fread(cmd = sprintf("zcat %s", shQuote(path)), sep = "\t",
          header = TRUE, data.table = FALSE)
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
# Automatic COLOR handling
# -----------------------------

colors_cfg <- vc$COLORS

if (is.null(colors_cfg$FAMILIES) || is.null(colors_cfg$PALETTE)) {
  stop("VOLCANO.COLORS must contain FAMILIES and PALETTE fields.")
}

families <- as.character(unlist(colors_cfg$FAMILIES))
palette_name <- as.character(colors_cfg$PALETTE)

other_label <- if (!is.null(colors_cfg$OTHER_LABEL)) {
  as.character(colors_cfg$OTHER_LABEL)
} else {
  "Other"
}

other_color <- if (!is.null(colors_cfg$OTHER_COLOR)) {
  as.character(colors_cfg$OTHER_COLOR)
} else {
  "grey70"
}

n_fam <- length(families)

# ---- Determine palette type ----
if (palette_name %in% rownames(brewer.pal.info)) {

  max_colors <- brewer.pal.info[palette_name, "maxcolors"]

  # Validation
  if (n_fam > max_colors) {
    stop(sprintf(
      "Palette '%s' supports maximum %d colors but %d families were requested.\nChoose a larger palette (e.g. Set3) or reduce FAMILIES.",
      palette_name, max_colors, n_fam
    ))
  }

  palette_vals <- brewer.pal(max(n_fam, 3), palette_name)[1:n_fam]

} else if (palette_name %in% c("viridis","magma","plasma","cividis","inferno")) {

  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("viridisLite package required for viridis palettes.")
  }

  palette_vals <- viridisLite::viridis(n_fam, option = palette_name)

} else {
  stop(sprintf("Unknown palette '%s'. Use RColorBrewer or viridis palettes.", palette_name))
}

# Build named color vector
colors <- setNames(palette_vals, families)
colors[[other_label]] <- other_color

# Optional manual order
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

layout_order <- as.character(unlist(vc$LAYOUT$order))

if (length(layout_order) == 0) {
  stop("VOLCANO.LAYOUT.order must contain at least one contrast.")
}

if (!all(layout_order %in% contrasts_keep)) {
  stop("VOLCANO.LAYOUT.order must be a subset of VOLCANO.CONTRASTS.")
}

titles <- NULL
if (!is.null(vc$TITLES) && !is.null(names(vc$TITLES))) {
  titles <- setNames(
    vapply(names(vc$TITLES), function(k) as.character(vc$TITLES[[k]]), character(1)),
    names(vc$TITLES)
  )
}

# -----------------------------
# Join and normalize
# -----------------------------
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
  idx <- which(grepl("guide", sapply(g$grobs, function(x) x$name), ignore.case = TRUE))
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

make_volcano <- function(d, title, show_legend = FALSE) {

  n_categories <- length(colors)
  legend_ncol <- if (!is.null(vc$LEGEND_NCOL)) {
  as.integer(vc$LEGEND_NCOL)
  } else {
  ceiling(n_categories / 4)   # auto wrap: ~4 rows max
  }
  # legend_ncol <- if (!is.null(vc$LEGEND_NCOL)) as.integer(vc$LEGEND_NCOL) else 1
  legend_pt   <- if (!is.null(vc$LEGEND_POINT_MM)) as.numeric(vc$LEGEND_POINT_MM) else 4

  d_other <- d[d$family == other_label, , drop = FALSE]
  d_high  <- d[d$family != other_label, , drop = FALSE]

  p <- ggplot() +
    # background genes first
    geom_point(data = d_other, aes(x = log2fc, y = neglog10_fdr, color = family),
               size = as.numeric(vc$OUTPUT$point_size), alpha = 0.6)+
    # highlighted families on top
    geom_point(data = d_high, aes(x = log2fc, y = neglog10_fdr, color = family),
               size = as.numeric(vc$OUTPUT$point_size)) +
    geom_hline(yintercept = -log10(as.numeric(vc$PADJ_CUTOFF)), linetype = "dashed") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(-as.numeric(vc$LFC_CUTOFF),
                              as.numeric(vc$LFC_CUTOFF)), linetype = "dashed") +
    scale_color_manual(
      values = colors,
      drop = FALSE
    ) +
    guides(
      color = guide_legend(
        ncol = legend_ncol,
        override.aes = list(size = legend_pt)
      )
    ) +
    labs(
      title = title,
      x = expression(log[2](FC)),
      y = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = as.numeric(vc$OUTPUT$base_font_size)) +
    theme(
      # plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend) "right" else "none"
    )

  p
}

# -----------------------------
# Build plots dynamically
# -----------------------------
plots <- list()
for (cn in layout_order) {
  title <- if (!is.null(titles) && !is.null(titles[[cn]])) titles[[cn]] else cn
  plots[[cn]] <- make_volcano(
    df[df$contrast == cn, , drop = FALSE],
    title,
    show_legend = FALSE
  )
}

plot_grobs <- lapply(layout_order, function(cn) {
  ggplotGrob(plots[[cn]])
})

# Legend from first plot
p_for_legend <- make_volcano(
  df[df$contrast == layout_order[[1]], , drop = FALSE],
  title = "",
  show_legend = TRUE
)

legend_grob <- extract_legend_grob(p_for_legend)
if (is.null(legend_grob)) {
  stop("Could not extract ggplot2 legend grob.")
}

# -----------------------------
# Dynamic grid layout
# -----------------------------
n_plots <- length(plot_grobs)
total_cells <- n_plots + 1   # include legend

rel_widths <- unlist(vc$LAYOUT$rel_widths)
ncol <- length(rel_widths)
nrow <- ceiling(total_cells / ncol)

lay <- grid.layout(
  nrow = nrow,
  ncol = ncol,
  widths = unit(rel_widths, "null"),
  heights = unit(rep(1, nrow), "null")
)

draw_all <- function() {

  grid.newpage()
  pushViewport(viewport(layout = lay))

  draw_cell <- function(grob, r, c) {
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))
    grid.draw(grob)
    popViewport()
  }

  # Draw plots
  for (i in seq_along(plot_grobs)) {
    row <- ceiling(i / ncol)
    col <- ((i - 1) %% ncol) + 1
    draw_cell(plot_grobs[[i]], row, col)
  }

  # Draw legend
  legend_index <- n_plots + 1
  row <- ceiling(legend_index / ncol)
  col <- ((legend_index - 1) %% ncol) + 1
  draw_cell(legend_grob, row, col)

  popViewport()
}

# -----------------------------
# Save outputs
# -----------------------------
pdf(out_pdf,
    width = as.numeric(vc$OUTPUT$width_in),
    height = as.numeric(vc$OUTPUT$height_in),
    onefile = TRUE)
draw_all()
dev.off()

png(out_png,
    width = as.numeric(vc$OUTPUT$width_in),
    height = as.numeric(vc$OUTPUT$height_in),
    units = "in",
    res = as.integer(vc$OUTPUT$dpi))
draw_all()
dev.off()

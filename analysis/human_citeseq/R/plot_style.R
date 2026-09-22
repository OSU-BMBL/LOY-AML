# R/plot_style.R
# Shared colour palettes, display labels, figure export and donor identifiers.

library(ggplot2)

# Color scheme (consistent across all figures)
COLORS <- list(
  group = c("LOY" = "#E41A1C", "CN-Male" = "#377EB8", "CN-Female" = "#4DAF4A"),
  cd8_state = c(Naive = "#66C2A5", Effector = "#FC8D62", Dysfunctional = "#8DA0CB"),
  # Broad cell-type palette for overview UMAPs (Brewer Set1 family, matching the
  # group colours). Leukemia = the project red; immune lineages get distinct
  # Set1 hues. "Unidentified" defined (grey) though it is dropped from F1A.
  cell_type = c("Leukemia" = "#E41A1C", "CD8 T" = "#377EB8", "CD4 T" = "#4DAF4A",
                "B cell" = "#984EA3", "NK" = "#FF7F00",
                "Unidentified" = "#999999", "Other" = "#BDBDBD"),
  diverging = c(high = "#B2182B", low = "#2166AC")
)

# ---------------------------------------------------------------------------
# Group labels use upright superscripts in plotmath expressions.
# Data retain LOY, CN-Male and CN-Female; only display labels change.
# Parse with label_parsed or parse(text = ...).
GROUP_LABELS <- c("LOY"       = "M^{plain(LOY)}",
                  "CN-Male"   = "M^{plain(CN)}",
                  "CN-Female" = "F^{plain(CN)}")

#' Map raw group values to plotmath display labels (M^LOY / M^CN / F^CN).
#' Unknown values pass through unchanged. Returns a character vector suitable
#' for parsing (label_parsed / parse(text=)).
group_label <- function(x) {
  x <- as.character(x)
  out <- GROUP_LABELS[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

# ---------------------------------------------------------------------------
# Display the Naive state as Naive/stem-like; retain the original
# state values in analysis objects. These labels use plain text.
CD8_STATE_LABELS <- c("Naive"         = "Naive/stem-like",
                      "Effector"      = "Effector",
                      "Dysfunctional" = "Dysfunctional")

#' Map raw CD8 state values to display labels. Unknown values pass through.
cd8_state_label <- function(x) {
  x <- as.character(x)
  out <- CD8_STATE_LABELS[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

# Publication-ready theme
theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold"),
      axis.title = element_text(size = base_size - 1),
      axis.text = element_text(size = base_size - 3),
      legend.text = element_text(size = base_size - 5),
      legend.title = element_text(size = base_size - 3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# Save a PNG figure with standard dimensions.
save_fig <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = "white")
}

# ---------------------------------------------------------------------------
# Export conventions: editable PDF, 600-dpi TIFF and 300-dpi PNG.
# Dense point layers are rasterized; text remains editable in PDF.
# ---------------------------------------------------------------------------
#' Nature-style theme for submission panels (small Arial, editable text)
theme_nature_contract <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line   = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.35, colour = "black"),
      axis.title  = element_text(size = base_size),
      axis.text   = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text  = element_text(size = base_size - 0.7),
      legend.key.size = unit(3, "mm"),
      strip.text  = element_text(size = base_size - 0.3, face = "bold"),
      strip.background = element_blank(),
      plot.title  = element_text(size = base_size + 0.5, face = "bold"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank()
    )
}

#' Export a single submission panel as PDF (primary) + 600-dpi TIFF (final)
#' + 300-dpi PNG (shareable).
#'
#' Writes `<filename>.pdf` (cairo_pdf, Arial, editable text/vector),
#' `<filename>.tiff` (ragg, 600 dpi flattened), and `<filename>.png`
#' (300 dpi, white background). Accepts ggplot, patchwork, or a
#' grid grob / ComplexHeatmap (drawn via grid::grid.draw / ComplexHeatmap::draw).
#' Width/height in mm; Nature columns are 89 mm (1-col) / 183 mm (2-col).
#'
#' For ggplot/patchwork inputs the PNG is written with ggplot2::ggsave;
#' for the function-form (ComplexHeatmap path) it is written with
#' ragg::agg_png at 300 dpi.
#'
#' @param plot ggplot/patchwork object OR a function() that draws to the device
#'   (use the function form for ComplexHeatmap: \code{function() draw(ht)}).
#' @param family PDF font family. Default "Arial" -> cairo_pdf with embedded
#'   Arial. Pass
#'   `family = NULL` to use the plain default device (grDevices::pdf, Helvetica)
#'   for Helvetica font metrics.
save_pub_r <- function(plot, filename, width_mm = 89, height_mm = 70, dpi = 600,
                       png_dpi = 300, family = "Arial") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  draw_it <- function() {
    if (is.function(plot)) plot() else print(plot)
  }
  # Editable PDF; family = NULL selects the default Helvetica device.
  if (is.null(family)) {
    grDevices::pdf(paste0(filename, ".pdf"), width = w, height = h)
  } else {
    grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = family)
  }
  draw_it()
  grDevices::dev.off()
  # TIFF — flattened 600-dpi final deliverable
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h,
                 units = "in", res = dpi, background = "white")
  draw_it()
  grDevices::dev.off()
  # PNG — 300-dpi preview.
  if (is.function(plot)) {
    ragg::agg_png(paste0(filename, ".png"), width = w, height = h,
                  units = "in", res = png_dpi, background = "white")
    draw_it()
    grDevices::dev.off()
  } else {
    ggplot2::ggsave(paste0(filename, ".png"), plot, width = w, height = h,
                    dpi = png_dpi, bg = "white")
  }
  invisible(c(paste0(filename, ".pdf"), paste0(filename, ".tiff"),
              paste0(filename, ".png")))
}

#' Unified donor identifier
#'
#' LOY samples use `sample_id` (S1-S6); CN-Male/CN-Female use `donor_id`
#' (AMLxxxx). Pseudobulk DESeq2 needs a single donor column across groups.
#' Falls back to `orig.ident` if neither is set.
#'
#' @param obj Seurat object OR data.frame meta.data
#' @return character vector, length = ncol(obj) or nrow(meta)
get_donor_id <- function(obj) {
  meta <- if (inherits(obj, "Seurat")) obj@meta.data else obj
  donor <- ifelse(!is.na(meta$donor_id) & nzchar(as.character(meta$donor_id)),
                  as.character(meta$donor_id),
                  ifelse(!is.na(meta$sample_id) & nzchar(as.character(meta$sample_id)),
                         as.character(meta$sample_id),
                         as.character(meta$orig.ident)))
  stopifnot(!any(is.na(donor)), !any(donor == ""))
  donor
}

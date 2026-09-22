# Shared plotting components for expression comparisons, UMAPs and
# directional GO enrichment. Grouped comparisons can display donor-level
# pseudobulk statistics; cell-level expression is used for the distribution.
# Depends on R/plot_style.R for palettes, labels and export functions.
library(ggplot2)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Format a padj for a panel subtitle. Never let a nonzero padj round to
#' "0.000" (reads as an impossible p = 0).
.ff_fmt_pad <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
  if (p < 1e-3) return(formatC(p, format = "e", digits = 1))
  formatC(p, format = "f", digits = 3)
}

#' Build the per-gene stat label from a DESeq2 data.frame (gene, log2FC, padj).
#' Returns NULL when no stat table is supplied (graceful "no stat line").
#' Returns "not tested" when the table exists but the gene is absent.
.ff_stat_label <- function(deg_deseq2, gene) {
  if (is.null(deg_deseq2)) return(NULL)
  if (!(gene %in% deg_deseq2$gene)) return("not tested (low pseudobulk counts)")
  row <- deg_deseq2[deg_deseq2$gene == gene, , drop = FALSE]
  sprintf("log2FC %+.2f, padj %s", row$log2FC[1], .ff_fmt_pad(row$padj[1]))
}

# ---------------------------------------------------------------------------
# ff_violin_compare — grouped gene-comparison panel
# ---------------------------------------------------------------------------

#' Grouped gene-comparison panel with optional pseudobulk annotations.
#'
#' One code path for every comparison panel. Produces a patchwork grid of
#' per-gene violin+box facets. Each facet is filled by group, carries the
#' donor-level pseudobulk DESeq2 log2FC + padj (LOY vs CN-Male) as its subtitle,
#' and shows the %-positive (fraction of cells with expr > 0) under each group.
#' NO cell-level Wilcoxon stars/brackets (n-inflated -> overclaim).
#'
#' Two modes:
#'   * facet_celltype = NULL  -> one violin+box per gene (genes laid out in a
#'       grid of `ncol` columns). Used by F4D, F1C, F5E1/F5E2 (single facet var).
#'   * facet_celltype = <col> -> each gene panel is itself faceted by that cell
#'       grouping (e.g. Leukemia vs CD8 T), giving a genes x cell-type grid.
#'       Used by F3supp. In this mode %-positive is shown in the x labels of the
#'       FIRST facet only (kept compact) and the stat subtitle is the blast/whole
#'       pseudobulk stat passed in `deg_deseq2`.
#'
#' @param seurat Seurat object (already subset to the groups/cell type you want,
#'   OR pass a pre-built data.frame via `df` — see Details). Expression is read
#'   from the "data" layer.
#' @param genes character vector of genes to plot (one facet each).
#' @param cell_type optional. If supplied AND `facet_celltype` is NULL, cells are
#'   restricted to this cell_type (meta column `cell_type`). Ignored when a
#'   pre-built `df` is supplied or when faceting by cell type.
#' @param deg_deseq2 data.frame with columns gene, log2FC, padj (pseudobulk
#'   DESeq2, LOY vs CN-Male). NULL -> stat line omitted gracefully.
#' @param groups character vector of groups to show, in display order. Default
#'   c("LOY","CN-Male"). Add "CN-Female" only for Y/X dosage panels.
#' @param facet_celltype optional name of a cell-type grouping column to facet
#'   each gene panel by (genes x cell_type grid). When set, `df` MUST carry it.
#' @param ncol columns in the gene grid (ignored when facet_celltype set and a
#'   single combined plot is returned).
#' @param title,subtitle plot_annotation title/subtitle. Defaults to NULL for
#'   both — when NULL no title/subtitle block is rendered (Nature-minimal style).
#'   Pass an explicit string to add a title or subtitle.
#' @param df optional pre-built long data.frame with columns group, expr, gene
#'   (and facet_celltype col if faceting). Bypasses Seurat extraction; lets
#'   callers reuse cached dosage_df without re-touching the Seurat object.
#' @param y_lab y-axis label (default "log-norm expr").
#' @param base_size Base font size in points. Increase for panels that will be
#'   reduced substantially during final figure assembly.
#' @param cnf_note logical; append the CN-Female caveat to the subtitle when
#'   "CN-Female" is among `groups` AND subtitle is non-NULL. Default TRUE.
#' @return a patchwork object (facet_celltype = NULL) or a single faceted ggplot
#'   (facet_celltype set), ready for save_pub_r().
ff_violin_compare <- function(seurat = NULL, genes, cell_type = NULL,
                              deg_deseq2 = NULL,
                              groups = c("LOY", "CN-Male"),
                              facet_celltype = NULL, ncol = 6,
                              title = NULL, subtitle = NULL, df = NULL,
                              y_lab = "log-norm expr", cnf_note = TRUE,
                              group_labels = NULL,
                              sig_brackets = FALSE, sig_method = "wilcox.test",
                              sig_comparisons = NULL, title_colors = NULL,
                              base_size = 6.5) {
  group_cols <- COLORS$group[groups]
  # Opt-in pairwise significance brackets (mode 1). When TRUE, each panel gets
  # ggpubr::stat_compare_means brackets (p.signif stars) for every group pair
  # (or `sig_comparisons` if supplied) instead of the log2FC/padj subtitle.
  # Significance annotations use cells as the statistical unit.
  if (sig_brackets && is.null(sig_comparisons))
    sig_comparisons <- utils::combn(groups, 2, simplify = FALSE)
  # Optional superscript x-axis labels (mode 1 only). Named vector mapping each
  # group level -> a Markdown/HTML string (e.g. "M<sup>LOY</sup>"), rendered via
  # ggtext::element_markdown. NULL (default) keeps the plain group names, so all
  # other panels are unaffected.
  use_md_labels <- !is.null(group_labels)

  # ---- assemble the long data.frame ---------------------------------------
  if (is.null(df)) {
    stopifnot(!is.null(seurat))
    so <- seurat
    keep <- so$group %in% groups
    if (!is.null(cell_type) && is.null(facet_celltype)) {
      keep <- keep & so$cell_type %in% cell_type
    }
    so <- subset(so, cells = colnames(so)[keep])
    expr_mat <- Seurat::GetAssayData(so, layer = "data")
    avail <- intersect(genes, rownames(expr_mat))
    missing <- setdiff(genes, avail)
    if (length(missing)) warning("ff_violin_compare: genes not in matrix: ",
                                 paste(missing, collapse = ", "))
    pieces <- lapply(avail, function(g) {
      d <- data.frame(group = as.character(so$group),
                      expr  = as.numeric(expr_mat[g, ]),
                      gene  = g, stringsAsFactors = FALSE)
      if (!is.null(facet_celltype)) d[[facet_celltype]] <- so@meta.data[[facet_celltype]]
      d
    })
    df <- do.call(rbind, pieces)
    genes <- avail
  } else {
    df <- df[df$group %in% groups, , drop = FALSE]
    if (!("gene" %in% colnames(df)))
      stop("ff_violin_compare: pre-built df must have a 'gene' column")
    genes <- intersect(genes, unique(df$gene))
  }

  df$group <- factor(df$group, levels = groups)

  # ---- subtitle handling (Nature-minimal: default NULL = no subtitle) --------
  # When subtitle is explicitly supplied and cnf_note is TRUE, append CN-Female
  # caveat. When NULL (default), no subtitle block is rendered.
  if (!is.null(subtitle) && cnf_note && "CN-Female" %in% groups) {
    subtitle <- paste0(subtitle,
      " CN-Female n=2, cell-level only (no pseudobulk stat).")
  }

  # ---- per-gene panel builder ---------------------------------------------
  # x labels carry per-group %-positive so a near-flat violin reads as
  # "near zero", not "not measured".
  pct_labels <- function(sub) {
    pct <- tapply(sub$expr, sub$group, function(x) round(100 * mean(x > 0)))
    vapply(levels(sub$group), function(g) {
      p <- ifelse(is.na(pct[[g]]), 0L, pct[[g]])
      if (use_md_labels) {
        disp <- if (g %in% names(group_labels)) group_labels[[g]] else g
        sprintf("%s<br>%d%%+", disp, p)   # ggtext markdown (superscript + %+)
      } else {
        sprintf("%s\n%d%%+", g, p)
      }
    }, character(1))
  }

  if (is.null(facet_celltype)) {
    # ---- mode 1: one violin+box per gene, patchwork grid ------------------
    panels <- lapply(genes, function(g) {
      sub <- df[df$gene == g, , drop = FALSE]
      sub$group <- factor(sub$group, levels = groups)
      x_labs <- pct_labels(sub)
      ttl_col <- if (!is.null(title_colors) && g %in% names(title_colors))
                   title_colors[[g]] else "black"
      p <- ggplot(sub, aes(x = group, y = expr, fill = group)) +
        geom_violin(scale = "width", trim = TRUE, alpha = 0.9, linewidth = 0.25) +
        geom_boxplot(width = 0.13, outlier.shape = NA, alpha = 1, linewidth = 0.25) +
        scale_fill_manual(values = group_cols) +
        scale_x_discrete(labels = x_labs)
      if (sig_brackets) {
        # Cap the y-axis at a robust upper bound (99th pct) so a few high
        # outlier cells don't stretch the axis and leave low-expression panels
        # mostly empty. Brackets are placed compactly just above the data and
        # the rare tail beyond the cap is clipped (density still computed on the
        # full data). pairwise compare-means stars (p.signif).
        d_hi <- as.numeric(stats::quantile(sub$expr, 0.99, na.rm = TRUE))
        if (!is.finite(d_hi) || d_hi <= 0) d_hi <- max(sub$expr, na.rm = TRUE)
        if (!is.finite(d_hi) || d_hi <= 0) d_hi <- 1
        nc  <- length(sig_comparisons)
        by  <- d_hi * (1.03 + 0.13 * (seq_len(nc) - 1))   # stacked brackets
        top <- max(by) + d_hi * 0.10
        p <- p +
          ggpubr::stat_compare_means(
            comparisons   = sig_comparisons,
            method        = sig_method,
            label         = "p.signif",
            label.y       = by,           # explicit, compact stacking
            size          = 2.0,          # star glyph size (small for Nature panel)
            bracket.size  = 0.2,
            tip.length    = 0.01,
            vjust         = 0.4
          ) +
          scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
          coord_cartesian(ylim = c(0, top))
      }
      # brackets replace the subtitle stat; otherwise show log2FC/padj.
      p +
        labs(title = g,
             subtitle = if (sig_brackets) NULL else .ff_stat_label(deg_deseq2, g),
             y = y_lab, x = NULL) +
        theme_nature_contract(base_size = base_size) +
        theme(
          legend.position = "none",
          axis.text.x  = if (use_md_labels)
                           ggtext::element_markdown(size = base_size - 2, colour = "black", lineheight = 1.0)
                         else
                           element_text(size = base_size - 2, colour = "black", lineheight = 0.9),
          axis.text.y  = element_text(size = base_size - 2, colour = "black"),
          axis.title.y = element_text(size = base_size - 1.5),
          plot.title   = element_text(size = base_size, face = "italic", colour = ttl_col),
          plot.subtitle = element_text(size = base_size - 2, colour = "grey35"),
          plot.margin  = margin(3, 4, 2, 3)
        )
    })
    combined <- patchwork::wrap_plots(panels, ncol = ncol)
    # Only add plot_annotation when at least one of title/subtitle is non-NULL
    # (Nature-minimal: no floating text above the panel grid by default).
    if (!is.null(title) || !is.null(subtitle)) {
      combined <- combined +
        patchwork::plot_annotation(
          title = title, subtitle = subtitle,
          theme = theme(
            plot.title    = element_text(size = base_size + 0.5, face = "bold", family = "Arial"),
            plot.subtitle = element_text(size = base_size - 1, colour = "grey35",
                                         family = "Arial", lineheight = 1.05)
          )
        )
    }
    return(combined)
  }

  # ---- mode 2: genes x cell_type grid (single faceted ggplot) -------------
  df$gene <- factor(df$gene, levels = genes)
  ct_levels <- levels(factor(df[[facet_celltype]]))
  df[[facet_celltype]] <- factor(df[[facet_celltype]], levels = ct_levels)

  # %-positive annotation: one label per group x gene x cell-type, placed
  # just under the panel. Compact: rendered as the x-axis tick text instead of
  # geom_text to avoid clutter in a tall grid.
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.9, linewidth = 0.22) +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 1, linewidth = 0.22) +
    scale_fill_manual(values = group_cols) +
    facet_grid(gene ~ .data[[facet_celltype]], scales = "free_y") +
    labs(title = title, subtitle = subtitle, x = NULL, y = y_lab) +
    theme_nature_contract(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.text.x   = element_text(angle = 40, hjust = 1, size = base_size - 2, colour = "black"),
      axis.text.y   = element_text(size = base_size - 2.5, colour = "black"),
      strip.text.x  = element_text(size = base_size - 1, face = "bold"),
      strip.text.y  = element_text(size = base_size - 2, face = "italic"),
      # Nature-minimal: suppress title/subtitle space when NULL
      plot.title    = if (is.null(title)) element_blank()
                      else element_text(size = base_size + 0.5, face = "bold", family = "Arial"),
      plot.subtitle = if (is.null(subtitle)) element_blank()
                      else element_text(size = base_size - 1.5, colour = "grey35", family = "Arial",
                                        lineheight = 1.05),
      plot.margin   = margin(3, 3, 2, 3)
    )
  p
}

# ---------------------------------------------------------------------------
# ff_umap — UMAP plotting component
# ---------------------------------------------------------------------------

#' Canonical UMAP panel: rasterised points (600 dpi) + vector axes/legend.
#'
#' Point size auto-scales to cell count (clamped 0.1-0.4) so it reads crisply at
#' panel size across PNG/PDF/TIFF; alpha ~0.6. Uses coord_equal and a minimal
#' theme_nature_contract. Legend keys are sized and labels given room so they are
#' fully visible. When color_by is "group", COLORS$group supplies the palette.
#'
#' @param obj Seurat object with a UMAP reduction.
#' @param color_by meta column to colour by (e.g. "group", "cell_type",
#'   "assigned_state").
#' @param palette named colour vector. When NULL and color_by == "group",
#'   COLORS$group is used; when color_by maps to CD8 states, COLORS$cd8_state is
#'   used; otherwise a Set2 palette is generated.
#' @param title optional panel title (method-level / neutral).
#' @param pt optional point size override; NULL -> auto (clamped 0.1-0.4).
#' @param reduction name of the reduction to plot (default "umap").
#' @param legend_title legend title (default = color_by).
#' @param label logical; draw group labels on the embedding (default FALSE).
#' @param label_repel logical; when label=TRUE, use ggrepel::geom_text_repel
#'   with halo instead of plain geom_text (default TRUE). When labelling,
#'   the colour legend is hidden because labels replace it.
#' @param alpha point opacity (default 0.6). Set to 1 for full-opacity,
#'   high-contrast overview UMAPs (e.g. F1A) where washed-out pastel clouds
#'   read as low quality; keep <1 for dense panels where density shading helps.
#' @param keep_legend logical; when label=TRUE, also draw the colour legend
#'   instead of letting labels replace it (default FALSE). Used by F1A where the
#'   panel carries both on-plot labels and a colour key.
#' @return a ggplot ready for save_pub_r().
ff_umap <- function(obj, color_by, palette = NULL, title = NULL, pt = NULL,
                    reduction = "umap", legend_title = NULL,
                    label = FALSE, label_repel = TRUE, alpha = 0.6,
                    keep_legend = FALSE) {
  if (!reduction %in% names(obj@reductions)) {
    alt <- intersect(c("umap", "wnn.umap", "ref.umap"), names(obj@reductions))
    if (length(alt)) reduction <- alt[1]
  }
  emb <- Seurat::Embeddings(obj, reduction = reduction)
  meta_col <- as.character(obj@meta.data[[color_by]])

  df <- data.frame(
    UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
    value  = meta_col, stringsAsFactors = FALSE
  )

  n <- nrow(df)
  if (is.null(pt)) {
    # auto: more cells -> smaller points, clamp 0.1-0.4
    pt <- max(0.1, min(0.4, 8000 / n))
  }

  # ---- resolve palette -----------------------------------------------------
  is_numeric_val <- suppressWarnings(!any(is.na(as.numeric(meta_col)))) &&
    !color_by %in% c("group", "cell_type", "assigned_state")
  if (is_numeric_val) {
    df$value <- as.numeric(meta_col)
  } else {
    lv <- sort(unique(df$value))
    df$value <- factor(df$value, levels = lv)
    if (is.null(palette)) {
      if (color_by == "group") {
        palette <- COLORS$group[intersect(names(COLORS$group), lv)]
      } else if (all(lv %in% names(COLORS$cd8_state))) {
        palette <- COLORS$cd8_state[lv]
      } else {
        palette <- setNames(
          RColorBrewer::brewer.pal(max(3, length(lv)), "Set2")[seq_along(lv)], lv)
      }
    }
  }

  if (is.null(legend_title)) legend_title <- color_by

  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = value)) +
    ggrastr::rasterise(
      geom_point(size = pt, alpha = alpha, stroke = 0), dpi = 600) +
    coord_equal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = legend_title) +
    theme_nature_contract() +
    theme(
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),
      axis.text   = element_blank(),
      axis.title  = element_text(size = 5.5),
      legend.key.size = unit(3, "mm"),
      legend.text  = element_text(size = 5, margin = margin(l = 1)),
      legend.title = element_text(size = 5.5),
      # Nature-minimal: suppress title space when title=NULL
      plot.title   = if (is.null(title)) element_blank()
                     else element_text(size = 7, face = "bold")
    )

  if (is_numeric_val) {
    p <- p + scale_colour_viridis_c(option = "viridis")
  } else {
    if (label) {
      # By default labels replace the legend; keep_legend=TRUE shows both
      # (matches the F1A overview panel, which carries on-plot labels AND a key).
      p <- p + scale_colour_manual(values = palette, drop = FALSE)
      if (keep_legend) {
        p <- p + guides(colour = guide_legend(
          override.aes = list(size = 2, alpha = 1), ncol = 1))
      } else {
        p <- p + guides(colour = "none")
      }
      cen <- aggregate(cbind(UMAP_1, UMAP_2) ~ value, data = df, FUN = median)
      if (label_repel) {
        p <- p + ggrepel::geom_text_repel(
          data = cen, aes(label = value),
          colour = "black", size = 2, fontface = "bold",
          segment.size = 0.2, segment.color = "grey50",
          min.segment.length = 0.2,
          box.padding = 0.35, point.padding = 0.1,
          max.overlaps = Inf,
          bg.color = "white", bg.r = 0.08,
          show.legend = FALSE
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = cen, aes(label = value), colour = "black", size = 2,
          fontface = "bold", show.legend = FALSE)
      }
    } else {
      p <- p + scale_colour_manual(values = palette, drop = FALSE) +
        guides(colour = guide_legend(
          override.aes = list(size = 2, alpha = 1), ncol = 1))
    }
  }
  p
}

# ---------------------------------------------------------------------------
# Fig 1d — leukemic-blast directional pathway butterfly (GO-BP ORA)
# ---------------------------------------------------------------------------
# Bars show -log10(adjusted P): down-in-LOY terms extend left, and
# up-in-LOY terms extend right. The input method is specified by the caller.
#
# With no term list, selection uses the smallest adjusted P values. The
# explicit lists below instead display a manually selected subset of terms,
# intersected with the ORA results and optionally backfilled. They are not
# a comprehensive summary of all significant enrichment results.
#
# Up-in-LOY display: cell-cycle and proliferation terms. Other significant
# immune-activating and antigen-receptor signalling terms are excluded
# from this display list; consult the full ORA tables for all results.
FF_F1D_UP_TERMS <- c(
  "chromosome segregation",
  "mitotic cell cycle phase transition",
  "regulation of mitotic cell cycle phase transition",
  "cell cycle G1/S phase transition",
  "mitotic sister chromatid cohesion",
  "regulation of mitotic cell cycle",
  "regulation of chromosome segregation",
  "regulation of chromosome organization",
  "spindle assembly",
  "regulation of mitotic metaphase/anaphase transition",
  "regulation of ubiquitin-dependent protein catabolic process",
  "protein polyubiquitination",
  "hematopoietic progenitor cell differentiation"
)

# Down-in-LOY display: antigen-presentation and immune-effector terms.
# Cytoplasmic-translation and ribosome-biogenesis terms are excluded
# from this display list; consult the full ORA tables for all results.
FF_F1D_DOWN_TERMS <- c(
  "antigen processing and presentation of peptide antigen",
  "MHC protein complex assembly",
  "positive regulation of cell killing",
  "positive regulation of immune effector process",
  "positive regulation of cytokine production",
  "positive regulation of T cell mediated cytotoxicity",
  "positive regulation of lymphocyte activation",
  "positive regulation of T cell activation",
  "leukocyte mediated cytotoxicity",
  "cell killing",
  "regulation of innate immune response",
  "type II interferon production",
  "positive regulation of defense response"
)

# ---------------------------------------------------------------------------
# Fixed display term sets for comparison across enrichment tables.
# With keep_missing = TRUE, all requested terms remain in the panel;
# unavailable or non-significant terms are visually distinguished.
# ---------------------------------------------------------------------------
FF_F1E_ORIG_UP_TERMS <- c(
  "mitotic cell cycle phase transition",
  "chromosome segregation",
  "regulation of mitotic cell cycle",
  "cell cycle checkpoint signaling",
  "proteasome-mediated ubiquitin-dependent protein catabolic process",
  "cell cycle G2/M phase transition",
  "regulation of mitotic sister chromatid separation",
  "mRNA processing",
  "double-strand break repair",
  "nuclear division",
  "centromere complex assembly",
  "regulation of DNA repair",
  "mitotic spindle assembly checkpoint signaling",
  "mitotic G2/M transition checkpoint",
  "stem cell division"
)
FF_F1E_ORIG_DOWN_TERMS <- c(
  "positive regulation of inflammatory response",
  "T cell mediated immunity",
  "leukocyte mediated immunity",
  "MHC protein complex assembly",
  "positive regulation of defense response",
  "tumor necrosis factor superfamily cytokine production",
  "phagocytosis",
  "regulation of immune effector process",
  "antigen processing and presentation",
  "immune response-regulating signaling pathway",
  "ATP biosynthetic process",
  "respiratory electron transport chain",
  "electron transport chain",
  "aerobic respiration",
  "oxidative phosphorylation"
)

#' Normalise a GO term for case-insensitive matching
.ff_norm_term <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

#' Capitalise first letter for display, preserving leading acronyms (mRNA, ATP,
#' MHC, DNA, TNF, T cell ...).
.ff_cap1 <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^[a-z][A-Z]", x) | grepl("^[A-Z]", x),
         x, paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
}

#' Pick the figure terms for one direction: curated whitelist (kept significant)
#' first, then back-fill to `n` with the most-significant remaining terms.
#'
#' @return list(df = tibble(term, nlp, padj, source), dropped = chr, added = chr)
.ff_select_terms <- function(ora_tbl, whitelist = NULL, n = 15, backfill = TRUE,
                             keep_missing = FALSE) {
  if (is.null(ora_tbl) || nrow(ora_tbl) == 0) {
    return(list(df = NULL,
                dropped = if (is.null(whitelist)) character(0) else whitelist,
                added = character(0)))
  }
  tbl <- ora_tbl
  tbl$.key <- .ff_norm_term(tbl$Description)
  tbl$nlp  <- -log10(pmax(tbl$p.adjust, 1e-300))

  # keep_missing -> "same pathways, refreshed data": plot EVERY whitelist term in
  # its given order, looking up its new p.adjust (NA/0 if it has no enrichment in
  # the new data), and flag significance. No dropping, no backfill, no re-sorting.
  if (keep_missing && !is.null(whitelist)) {
    wl_key <- .ff_norm_term(whitelist)
    idx    <- match(wl_key, tbl$.key)
    padj   <- ifelse(is.na(idx), NA_real_, tbl$p.adjust[idx])
    nlp    <- ifelse(is.na(padj), 0, -log10(pmax(padj, 1e-300)))
    disp   <- ifelse(is.na(idx), .ff_cap1(whitelist), .ff_cap1(tbl$Description[idx]))
    sig    <- !is.na(padj) & padj < 0.05
    src    <- ifelse(is.na(idx), "not enriched", ifelse(sig, "match", "ns"))
    df <- tibble::tibble(term = disp, nlp = nlp, padj = padj,
                         sig = sig, source = src)
    return(list(df = df,
                dropped = whitelist[is.na(idx)],          # term absent from ORA
                added = character(0)))
  }

  # whitelist = NULL -> pure top-N by significance (transparent, no curation)
  if (is.null(whitelist)) {
    keep <- utils::head(tbl[order(-tbl$nlp), , drop = FALSE], n)
    df <- tibble::tibble(term = .ff_cap1(keep$Description),
                         nlp = keep$nlp, padj = keep$p.adjust, source = "top-N")
    return(list(df = df, dropped = character(0), added = character(0)))
  }

  wl_key   <- .ff_norm_term(whitelist)

  kept_idx <- match(wl_key, tbl$.key)          # whitelist order, NA = not sig now
  dropped  <- whitelist[is.na(kept_idx)]
  keep     <- tbl[stats::na.omit(kept_idx), , drop = FALSE]
  keep$source <- "curated"

  added <- character(0)
  if (backfill && nrow(keep) < n) {
    pool <- tbl[!(tbl$.key %in% keep$.key), , drop = FALSE]
    pool <- pool[order(-pool$nlp), , drop = FALSE]
    take <- utils::head(pool, n - nrow(keep))
    if (nrow(take) > 0) { take$source <- "backfill"; added <- take$Description }
    keep <- dplyr::bind_rows(keep, take)
  }
  keep <- utils::head(keep[order(-keep$nlp), , drop = FALSE], n)
  df <- tibble::tibble(term = .ff_cap1(keep$Description),
                       nlp = keep$nlp, padj = keep$p.adjust, source = keep$source)
  list(df = df, dropped = dropped, added = added)
}

#' Build one side (panel) of the butterfly.
#'
#' @param df tibble(term, nlp) ordered any way (re-ordered here: most-sig at top)
#' @param side "down" (blue, reversed x, labels left) or "up" (red, labels right)
#' @param xmax shared axis maximum
#' @param base_size theme size
.ff_butterfly_side <- function(df, side, xmax, wrap = 32, base_size = 6.5,
                               cap = Inf, title = NULL, sig_line = NULL,
                               family = "Arial") {
  base_family <- if (is.null(family)) "" else family
  df <- df[order(df$nlp), , drop = FALSE]          # smallest first -> bottom
  df$term <- factor(df$term, levels = df$term)
  df$lab  <- stringr::str_wrap(as.character(df$term), width = wrap)
  df$lab  <- factor(df$lab, levels = df$lab)
  df$nlp_plot <- pmin(df$nlp, cap)                 # cap bar length (n-inflation)
  df$capped   <- df$nlp > cap
  if (!("sig" %in% names(df))) df$sig <- TRUE       # default: all significant
  fill <- if (side == "up") COLORS$diverging[["high"]] else COLORS$diverging[["low"]]
  if (is.null(title)) title <- if (side == "up") "Up in LOY" else "Up in CN-Male"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = nlp_plot, y = lab, alpha = sig)) +
    ggplot2::geom_col(fill = fill, width = 0.74) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.32),
                                guide = "none") +
    theme_nature_contract(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = base_size - 0.5,
                                           colour = "black", lineheight = 0.85),
      plot.title   = ggplot2::element_text(size = base_size, face = "bold",
                                           colour = fill, hjust = 0.5)
    )

  # dashed significance threshold (e.g. sig_line = 0.05 -> x = -log10(0.05))
  if (!is.null(sig_line)) {
    p <- p + ggplot2::geom_vline(xintercept = -log10(sig_line),
                                 linetype = "dashed", linewidth = 0.3,
                                 colour = "grey50")
  }

  # annotate true value on capped bars (so the cap never hides magnitude)
  if (any(df$capped)) {
    cd <- df[df$capped, , drop = FALSE]
    cd$txt <- sprintf("%.0f", cd$nlp)
    p <- p + ggplot2::geom_text(
      data = cd, ggplot2::aes(x = cap, y = lab, label = txt),
      hjust = if (side == "down") -0.15 else 1.15, colour = "white",
      fontface = "bold", size = base_size * 0.32, inherit.aes = FALSE)
  }

  if (side == "down") {
    p + ggplot2::scale_x_reverse(limits = c(xmax, 0),
                                 labels = function(x) {
                                   lab <- scales::label_number()(x)
                                   lab[abs(x) < .Machine$double.eps^0.5] <- ""
                                   lab
                                 },
                                 expand = ggplot2::expansion(mult = c(0.02, 0))) +
      ggplot2::scale_y_discrete(position = "left") +
      ggplot2::labs(title = title, x = NULL)
  } else {
    p + ggplot2::scale_x_continuous(limits = c(0, xmax),
                                    expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::scale_y_discrete(position = "right") +
      ggplot2::labs(title = title, x = NULL)
  }
}

#' Leukemic-blast directional pathway butterfly (Fig 1d)
#'
#' @param ora output of run_blast_ora()
#' @param out_basepath path WITHOUT extension; save_pub_r adds .pdf/.tiff/.png
#' @param n_per_side terms per wing (default 15)
#' @param up_terms,down_terms curated GO:BP whitelists
#' @param width_mm,height_mm panel size (2-col default)
#' @return list(paths, plot, selection = list(up, down)) — selection tables carry
#'   the fresh p.adjust + curated/backfill provenance for the figure legend.
ff_pathway_butterfly_leukemia <- function(
  ora,
  out_basepath = file.path(loy_paths("single_cell")$results, "figures",
                           "fig1", "F1D_pathway_butterfly_leukemia"),
  n_per_side = 15,
  up_terms = NULL, down_terms = NULL,
  width_mm = 183, height_mm = 95, base_size = 6.5,
  nlp_cap = Inf, backfill = TRUE, balance = TRUE, keep_missing = FALSE,
  sig_line = NULL, family = "Arial",
  up_title = "Up in LOY", down_title = "Up in CN-Male"
) {
  up_sel   <- .ff_select_terms(ora$up_in_loy, up_terms,  n = n_per_side,
                               backfill = backfill, keep_missing = keep_missing)
  down_sel <- .ff_select_terms(ora$up_in_cn,  down_terms, n = n_per_side,
                               backfill = backfill, keep_missing = keep_missing)
  if (is.null(up_sel$df) || is.null(down_sel$df))
    stop("ff_pathway_butterfly_leukemia: ORA returned no terms for one direction.")

  # balance wings to equal row count so the butterfly aligns vertically
  if (balance) {
    k <- min(nrow(up_sel$df), nrow(down_sel$df))
    up_sel$df   <- utils::head(up_sel$df[order(-up_sel$df$nlp), , drop = FALSE], k)
    down_sel$df <- utils::head(down_sel$df[order(-down_sel$df$nlp), , drop = FALSE], k)
  }

  xmax <- max(pmin(c(up_sel$df$nlp, down_sel$df$nlp), nlp_cap), na.rm = TRUE)
  xmax <- ceiling(xmax * 1.02)

  p_down <- .ff_butterfly_side(down_sel$df, "down", xmax, base_size = base_size,
                               cap = nlp_cap, title = down_title, sig_line = sig_line,
                               family = family)
  p_up   <- .ff_butterfly_side(up_sel$df,   "up",   xmax, base_size = base_size,
                               cap = nlp_cap, title = up_title, sig_line = sig_line,
                               family = family)

  xlab <- patchwork::wrap_elements(grid::textGrob(
    expression(-log[10] ~ "adjusted " * italic(P)),
    gp = grid::gpar(fontsize = base_size,
                    fontfamily = if (is.null(family)) "" else family)))

  body  <- patchwork::wrap_plots(p_down, p_up, nrow = 1, widths = c(1, 1))
  final <- patchwork::wrap_plots(body, xlab, ncol = 1, heights = c(1, 0.045))

  dir.create(dirname(out_basepath), recursive = TRUE, showWarnings = FALSE)
  paths <- save_pub_r(final, out_basepath, width_mm = width_mm,
                      height_mm = height_mm, family = family)

  # source-data CSV (DESIGN sec 0.1)
  sd <- dplyr::bind_rows(
    transform(down_sel$df, direction = "Up in CN-Male (down in LOY)"),
    transform(up_sel$df,   direction = "Up in LOY")
  )
  utils::write.csv(sd, paste0(out_basepath, "_sourcedata.csv"), row.names = FALSE)

  list(paths = paths, plot = final,
       selection = list(up = up_sel, down = down_sel), xmax = xmax)
}

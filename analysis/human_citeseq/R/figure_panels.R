# Manuscript panels: Y-linked expression, blast analyses and CD8 states.

library(ggplot2)
.FF_PANELS_DIR <- file.path(loy_paths("single_cell")$results, "figures")

.ff_prep_out <- function(out_basepath) {
  dir.create(dirname(out_basepath), recursive = TRUE, showWarnings = FALSE)
  out_basepath
}

.ff_render_paired_yx_bubble <- function(pairs, dosage_df, ds_blast, out_basepath) {
  GROUPS3 <- c("LOY", "CN-Male", "CN-Female")
  ctd <- dosage_df[dosage_df$cell_type == "Leukemia", ]

  # Order gene pairs by the largest Y-gene detection fraction in either male
  # group. Differential-expression significance does not determine ordering.
  ydetect <- vapply(pairs, function(p) {
    g <- p[1]
    if (!(g %in% colnames(ctd))) return(-Inf)
    max(mean(ctd[[g]][ctd$group == "LOY"]     > 0),
        mean(ctd[[g]][ctd$group == "CN-Male"] > 0))
  }, numeric(1))
  pairs       <- pairs[order(-ydetect)]
  pair_labels <- vapply(pairs, function(p) paste0(p[1], " / ", p[2]), character(1))
  # Per-facet gene labels: each dot is labelled by its OWN gene (not the combined
  # "Y / X" pair), so the Y-linked and X-linked facets are unambiguous. Rows stay
  # aligned by homolog pair via a shared rank: reversing within each homolog set
  # puts pair-rank-1 at the TOP of BOTH facets (free_y discrete axis = first level
  # at bottom).
  genes_y     <- vapply(pairs, `[`, character(1), 1)   # Y genes, detectability order
  genes_x     <- vapply(pairs, `[`, character(1), 2)   # X homologs, same order
  gene_levels <- c(rev(genes_y), rev(genes_x))

  # per (gene, group) summary: % expressing + mean log-norm expr.
  rec <- list()
  for (i in seq_along(pairs)) {
    p <- pairs[[i]]
    for (j in 1:2) {
      g  <- p[j]
      hl <- c("Y-linked", "X-linked")[j]
      if (!(g %in% colnames(ctd))) { warning("F2cd bubble: gene absent in dosage_df: ", g); next }
      for (grp in GROUPS3) {
        e <- ctd[[g]][ctd$group == grp]
        rec[[length(rec) + 1]] <- data.frame(
          pair = pair_labels[i], homolog = hl, gene = g, group = grp,
          pct = 100 * mean(e > 0, na.rm = TRUE), mean_expr = mean(e, na.rm = TRUE),
          stringsAsFactors = FALSE)
      }
    }
  }
  bub <- do.call(rbind, rec)

  # z-score mean expr WITHIN each gene so low-absolute genes still show contrast.
  bub$z <- stats::ave(bub$mean_expr, bub$gene, FUN = function(x) {
    s <- stats::sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s
  })
  bub$group   <- factor(bub$group, levels = GROUPS3)
  bub$gene    <- factor(bub$gene, levels = gene_levels)
  bub$homolog <- factor(bub$homolog, levels = c("Y-linked", "X-linked"))

  p <- ggplot(bub, aes(x = group, y = gene)) +
    geom_point(aes(size = pct, fill = z), shape = 21, stroke = 0.2, colour = "grey30") +
    facet_wrap(~ homolog, scales = "free_y", nrow = 1) +
    scale_size_area(name = "% expressing", max_size = 5.5, limits = c(0, 100)) +
    scale_fill_gradient2(name = "mean expr\n(z, per gene)",
                         low = COLORS$diverging[["low"]], mid = "grey95",
                         high = COLORS$diverging[["high"]], midpoint = 0) +
    scale_x_discrete(labels = function(b) parse(text = group_label(b))) +
    labs(x = NULL, y = NULL) +
    theme_nature_contract(base_size = 8) +
    theme(
      panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
      axis.text.y      = element_text(size = 7.5, face = "italic"),
      axis.text.x      = element_text(size = 7.5),
      strip.text       = element_text(size = 8, face = "bold"),
      legend.title     = element_text(size = 7.5),
      legend.text      = element_text(size = 7),
      legend.position  = "right",
      legend.box       = "vertical"
    )
  # Export the paired-gene panel at 89 x 120 mm.
  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 89, height_mm = 120)
}

fig_final_F1A_umap_celltype <- function(
    seurat_filtered,
    out_basepath = file.path(.FF_PANELS_DIR, "fig1/F1A_overview_umap_celltype")) {
  so <- seurat_filtered
  ct_col <- if ("cell_type" %in% colnames(so@meta.data)) "cell_type" else "Broad_cell_identity"
  so@meta.data[[ct_col]][is.na(so@meta.data[[ct_col]])] <- "Other"
  red <- intersect(c("umap", "wnn.umap", "ref.umap"), names(so@reductions))[1]

  # Display annotated cell types in the overview.
  drop_lvls <- c("Unidentified", "Other")
  keep <- !(as.character(so@meta.data[[ct_col]]) %in% drop_lvls)
  so <- subset(so, cells = colnames(so)[keep])
  so@meta.data[[ct_col]] <- droplevels(factor(so@meta.data[[ct_col]]))
  Seurat::Idents(so) <- so@meta.data[[ct_col]]

  # Project palette (COLORS$cell_type), named by identity; grey for any extra.
  lv   <- levels(Seurat::Idents(so))
  cols <- COLORS$cell_type[lv]
  cols[is.na(cols)] <- "#BDBDBD"
  names(cols) <- lv

  # Coloured cloud + legend + simple corner-arrow axes (NO ticks / numbers).
  emb <- Seurat::Embeddings(so, reduction = red)
  xr <- range(emb[, 1]); yr <- range(emb[, 2])
  dx <- diff(xr); dy <- diff(yr); sp <- max(dx, dy)
  alen <- 0.14 * sp
  ox <- xr[1] - 0.02 * dx
  oy <- yr[1] - 0.02 * dy
  gap <- 0.045 * sp
  arr <- grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")

  p_ct <- Seurat::DimPlot(so, reduction = red, label = FALSE, cols = cols,
                          pt.size = 0.8, raster = TRUE, raster.dpi = c(600, 600)) +
    Seurat::NoAxes() +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    annotate("segment", x = ox, xend = ox + alen, y = oy, yend = oy,
             arrow = arr, linewidth = 0.4, colour = "black") +
    annotate("segment", x = ox, xend = ox, y = oy, yend = oy + alen,
             arrow = arr, linewidth = 0.4, colour = "black") +
    annotate("text", x = ox + alen / 2, y = oy - gap, label = "UMAP 1",
             size = 2.6, hjust = 0.5, vjust = 1) +
    annotate("text", x = ox - gap, y = oy + alen / 2, label = "UMAP 2",
             size = 2.6, angle = 90, hjust = 0.5, vjust = 1) +
    coord_equal(clip = "off") +
    theme(plot.margin = margin(3, 3, 11, 11, "mm"),
          legend.position = "right",
          legend.text = element_text(size = 8))

  save_pub_r(p_ct, .ff_prep_out(out_basepath), width_mm = 140, height_mm = 120)
}

fig_final_F1A_umap_group <- function(
    seurat_filtered,
    out_basepath = file.path(.FF_PANELS_DIR, "fig1/F1A_overview_umap_group")) {
  so <- seurat_filtered
  group_cols <- COLORS$group

  p_grp <- Seurat::DimPlot(so, group.by = "group", cols = group_cols,
                           pt.size = 0.2, raster = TRUE, raster.dpi = c(600, 600)) +
    labs(title = NULL, x = "UMAP 1", y = "UMAP 2") +
    scale_color_manual(values = group_cols) +
    theme_nature_contract() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_text(size = 5.5)) +
    guides(color = guide_legend(override.aes = list(size = 2.5), ncol = 1,
                                title = "Group"))

  save_pub_r(p_grp, .ff_prep_out(out_basepath), width_mm = 70, height_mm = 60)
}

fig_final_F1B_blast_volcano <- function(
    degs_by_celltype,
    out_basepath = file.path(.FF_PANELS_DIR, "fig1/F1B_blast_volcano")) {
  deg_tbl <- degs_by_celltype[["blast"]][["deseq2"]]

  force_label <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-G",
                   "B2M", "TAPBP", "PSMB8", "PSMB9", "TAP1", "TAP2", "CIITA")
  fc_cut <- 1; padj_cut <- 0.05; top_label <- 15

  df <- deg_tbl[!is.na(deg_tbl$log2FC) & !is.na(deg_tbl$padj), ]
  df$neglog10p_raw <- -log10(pmax(df$padj, 1e-300))
  p99   <- stats::quantile(df$neglog10p_raw, 0.99, na.rm = TRUE)
  y_cap <- min(p99, 50)
  df$neglog10p <- pmin(df$neglog10p_raw, y_cap)
  df$capped    <- df$neglog10p_raw > y_cap

  df$signif <- dplyr::case_when(
    df$padj < padj_cut & df$log2FC >  fc_cut ~ "UP in LOY",
    df$padj < padj_cut & df$log2FC < -fc_cut ~ "DOWN in LOY",
    TRUE ~ "ns"
  )
  sig_df <- df[df$signif != "ns", ]
  df_top <- dplyr::slice_head(
    dplyr::arrange(sig_df, dplyr::desc(abs(sig_df$log2FC) * sig_df$neglog10p_raw)),
    n = top_label)
  df_force <- df[df$gene %in% force_label, ]
  df_label <- dplyr::distinct(dplyr::bind_rows(df_top, df_force), gene, .keep_all = TRUE)

  col_map <- c("UP in LOY" = "#E41A1C", "DOWN in LOY" = "#377EB8", "ns" = "grey70")

  p <- ggplot(df, aes(x = log2FC, y = neglog10p, color = signif)) +
    ggrastr::rasterise(geom_point(alpha = 0.5, size = 0.25), dpi = 600) +
    {
      capped_df <- df[df$capped & df$signif != "ns", ]
      if (nrow(capped_df) > 0)
        geom_point(data = capped_df, shape = 17, size = 1.2, alpha = 0.8)
      else NULL
    } +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed",
               linewidth = 0.35, color = "grey40") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed",
               linewidth = 0.35, color = "grey40") +
    annotate("text", x = Inf, y = y_cap, hjust = 1.05, vjust = -0.3,
             label = paste0("y capped at ", round(y_cap, 0), "; ▲ = off-scale"),
             size = 1.6, color = "grey40") +
    ggrepel::geom_text_repel(
      data = df_label, aes(label = gene), size = 1.8,
      max.overlaps = 60, force = 2, force_pull = 0.5, color = "black",
      min.segment.length = 0, segment.size = 0.25,
      box.padding = 0.35, point.padding = 0.1, seed = 42) +
    scale_color_manual(values = col_map, name = NULL) +
    coord_cartesian(ylim = c(0, y_cap * 1.06)) +
    labs(x = "log₂ fold change (LOY / CN-Male)", y = "-log₁₀ adjusted p") +
    theme_nature_contract() +
    theme(legend.position = "top")

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 89, height_mm = 80)
}

fig_final_F1C_genepair_leukemia <- function(
    gene_dosage_results, degs_by_celltype,
    out_basepath = file.path(.FF_PANELS_DIR, "fig1/F1C_genepair_leukemia")) {
  ds_blast <- degs_by_celltype[["blast"]][["deseq2"]]

  gene_pairs_8 <- list(
    c("UTY",    "KDM6A"),  c("USP9Y",  "USP9X"),
    c("ZFY",    "ZFX"),    c("EIF1AY", "EIF1AX"),
    c("RPS4Y1", "RPS4X"),  c("DDX3Y",  "DDX3X"),
    c("TMSB4Y", "TMSB4X"), c("KDM5D",  "KDM5C")
  )
  gp_genes <- unlist(gene_pairs_8)   # Y, X, Y, X, ... order

  y_genes <- vapply(gene_pairs_8, `[`, character(1), 1)
  x_genes <- vapply(gene_pairs_8, `[`, character(1), 2)
  gene_title_cols <- c(
    setNames(rep("#002060", length(y_genes)), y_genes),   # Y-linked navy
    setNames(rep("#C2561E", length(x_genes)), x_genes)    # X-linked rust
  )

  dosage_df <- gene_dosage_results$dosage_df
  ct_leuk   <- dosage_df[dosage_df$cell_type == "Leukemia", ]
  f1c_long  <- do.call(rbind, lapply(gp_genes, function(g) {
    if (!(g %in% colnames(ct_leuk))) return(NULL)
    data.frame(group = ct_leuk$group, expr = ct_leuk[[g]], gene = g,
               stringsAsFactors = FALSE)
  }))

  p_f1c <- ff_violin_compare(
    df         = f1c_long,
    genes      = gp_genes,
    deg_deseq2 = ds_blast,
    groups     = c("LOY", "CN-Male", "CN-Female"),
    ncol       = 4,
    title      = NULL,
    group_labels = c("LOY"       = "M<sup>LOY</sup>",
                     "CN-Male"   = "M<sup>CN</sup>",
                     "CN-Female" = "F<sup>CN</sup>"),
    sig_brackets = TRUE,
    title_colors = gene_title_cols
  )

  save_pub_r(p_f1c, .ff_prep_out(out_basepath), width_mm = 183, height_mm = 160)
}

fig_final_F1D_pathway_butterfly <- function(
    blast_ora_f1d,
    out_basepath = file.path(.FF_PANELS_DIR, "fig1/F1D_pathway_butterfly_leukemia")) {
  res <- ff_pathway_butterfly_leukemia(
    blast_ora_f1d,
    out_basepath = .ff_prep_out(out_basepath),
    n_per_side = 13,
    up_terms = FF_F1D_UP_TERMS, down_terms = FF_F1D_DOWN_TERMS,
    backfill = FALSE, balance = TRUE, nlp_cap = Inf,
    family = NULL,                       # default PDF device
    up_title = "Up in LOY", down_title = "Up in CN-Male",
    width_mm = 183, height_mm = 105, base_size = 8
  )
  res$paths
}

fig_final_F1E_pathway_origterms <- function(
    blast_ora_f1e,
    out_basepath = file.path(.FF_PANELS_DIR, "fig1/F1E_pathway_butterfly_origterms")) {
  res <- ff_pathway_butterfly_leukemia(
    blast_ora_f1e,
    out_basepath = .ff_prep_out(out_basepath),
    up_terms = FF_F1E_ORIG_UP_TERMS, down_terms = FF_F1E_ORIG_DOWN_TERMS,
    keep_missing = TRUE, balance = TRUE, sig_line = 0.05, nlp_cap = Inf,
    family = NULL,                       # default PDF device
    up_title = "Up in LOY", down_title = "Up in CN-Male",
    width_mm = 183, height_mm = 95
  )
  res$paths
}

fig_final_F2cd_paired_YX <- function(
    gene_dosage_results, degs_by_celltype,
    out_basepath = file.path(.FF_PANELS_DIR, "fig2/F2cd_paired_YX")) {
  .ff_render_paired_yx_bubble(
    pairs = list(
      c("RPS4Y1", "RPS4X"), c("DDX3Y", "DDX3X"), c("UTY", "KDM6A"),
      c("USP9Y", "USP9X"),  c("ZFY", "ZFX"),     c("EIF1AY", "EIF1AX"),
      c("KDM5D", "KDM5C"),  c("TMSB4Y", "TMSB4X")
    ),
    dosage_df = gene_dosage_results$dosage_df,
    ds_blast  = degs_by_celltype[["blast"]][["deseq2"]],
    out_basepath = out_basepath
  )
}

fig_final_F2CD_dosage_compensation <- function(
    gene_dosage_results,
    out_basepath = file.path(.FF_PANELS_DIR, "fig2/F2CD_dosage_compensation")) {
  res <- gene_dosage_results

  # group_means matrix: rows = groups, cols = genes (Leukemia blasts)
  gm <- as.data.frame(res$group_means)  # rownames = LOY / CN-Male / CN-Female

  # Gene pairs (same 7-pair subset)
  ALL_PAIRS <- list(
    c("UTY", "KDM6A"), c("KDM5D", "KDM5C"), c("ZFY", "ZFX"),
    c("USP9Y", "USP9X"), c("DDX3Y", "DDX3X"), c("EIF1AY", "EIF1AX"),
    c("RPS4Y1", "RPS4X")
  )

  # Compute ratio of group means: LOY mean(X-homolog) / CN-Male mean(X-homolog)
  comp_rows <- list()
  for (pair in ALL_PAIRS) {
    y_gene <- pair[1]; x_gene <- pair[2]
    pair_label <- paste0(y_gene, "/", x_gene)

    cnm_x <- gm["CN-Male", x_gene]
    loy_x <- gm["LOY",     x_gene]
    x_fold <- if (!is.na(cnm_x) && cnm_x > 0) loy_x / cnm_x else NA_real_

    category <- if (is.na(x_fold)) "Unknown"
                else if (x_fold >= 2.5)  "Overcompensation"
                else if (x_fold >= 1.3)  "Partial compensation"
                else if (x_fold >= 0.95) "No compensation"
                else                     "Dosage loss"

    comp_rows[[length(comp_rows) + 1]] <- data.frame(
      pair     = pair_label,
      x_fold   = x_fold,
      category = category,
      stringsAsFactors = FALSE
    )
  }

  cdf <- dplyr::bind_rows(comp_rows)
  cdf$pair <- factor(cdf$pair, levels = cdf$pair[order(cdf$x_fold)])
  cdf$category <- factor(cdf$category,
                          levels = c("Dosage loss", "No compensation",
                                     "Partial compensation", "Overcompensation"))

  cat_cols <- c(
    "Dosage loss"          = "#E41A1C",
    "No compensation"      = "#95A5A6",
    "Partial compensation" = "#E67E22",
    "Overcompensation"     = "#8E44AD"
  )

  p <- ggplot(cdf, aes(x = pair, y = x_fold, color = category)) +
    geom_hline(yintercept = 1.0, linetype = "dashed",
               color = "grey50", linewidth = 0.35) +
    geom_segment(aes(xend = pair, y = 1.0, yend = x_fold),
                 linewidth = 0.6) +
    geom_point(size = 1.8) +
    geom_text(aes(label = paste0(round(x_fold, 1), "x")),
              hjust = -0.4, size = 2.0, fontface = "bold") +
    scale_color_manual(values = cat_cols, name = NULL) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(
      y = "X-homolog expression ratio (LOY mean / CN-Male mean)",
      x = NULL
    ) +
    theme_nature_contract() +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      axis.text.y      = element_text(size = 6, face = "bold")
    )

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 80, height_mm = 65)
}

fig_final_F3a_ysig_umap <- function(
    seurat_filtered, lineage_annotation,
    out_basepath = file.path(.FF_PANELS_DIR, "fig3/F3a_ysig_umap"),
    annotation_mode = c("all", "blast_only"),
    style = c("original", "teal", "dashed", "shaded", "minimal", "navy_bold",
              "violet_dotted", "no_arrows", "grey_shaded", "shade_only", "all_regions"),
    layout = c("horizontal", "vertical"), return_plot = FALSE) {
  annotation_mode <- match.arg(annotation_mode)
  style <- match.arg(style)
  layout <- match.arg(layout)
  if (style == "minimal") annotation_mode <- "blast_only"
  so <- seurat_filtered

  # Y-signature score from the nine listed Y-linked genes.
  Y_GENES <- c("UTY", "USP9Y", "ZFY", "EIF1AY", "DDX3Y",
               "RPS4Y1", "KDM5D", "TMSB4Y", "PRKY")
  y_present <- intersect(Y_GENES, rownames(so))
  if (length(setdiff(Y_GENES, y_present)))
    warning("Y genes absent: ", paste(setdiff(Y_GENES, y_present), collapse = ", "))
  set.seed(42)
  so <- AddModuleScore(so, features = list(y_present), name = "Ysig", seed = 42)
  so$Ysig_score <- so$Ysig1

  red <- intersect(c("umap", "wnn.umap", "ref.umap"), names(so@reductions))[1]

  # ---- assemble plotting frame; high-score cells drawn on top ----------------
  # Facet order = positive control (CN-Male) -> depleted (LOY). CN-Female omitted.
  # Facet strips use plotmath display labels M^CN / M^LOY (COLORS/GROUP_LABELS).
  GROUPS <- if (layout == "vertical") c("LOY", "CN-Male") else c("CN-Male", "LOY")
  LAB_LEVELS <- group_label(GROUPS)        # c("M^{plain(CN)}", "M^{plain(LOY)}")
  emb <- Seurat::Embeddings(so, reduction = red)
  lineage <- setNames(lineage_annotation$labels$lineage,
                      lineage_annotation$labels$barcode)[colnames(so)]
  df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                   Ysig   = so$Ysig_score,
                   group  = factor(as.character(so$group), levels = GROUPS),
                   cell_type = unname(lineage),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$group), ]
  df$group_lab <- factor(group_label(df$group), levels = LAB_LEVELS)

  # Rough, data-derived cell-type regions. The 92% KDE boundary suppresses
  # scattered annotation outliers while retaining disconnected populations.
  CT_LEVELS <- c("Leukemia", "CD8 T", "CD4 T", "B cell", "NK")
  rough_outline <- function(cell_type) {
    d <- df[df$cell_type == cell_type, c("UMAP_1", "UMAP_2")]
    if (style == "original") {
      kde <- MASS::kde2d(d$UMAP_1, d$UMAP_2, n = 160)
    } else {
      # Extend the density grid to close contours at the observed extrema.
      pad_range <- function(x) range(x) + c(-1, 1) * diff(range(x)) * 0.15
      kde <- MASS::kde2d(d$UMAP_1, d$UMAP_2, n = 160,
                         lims = c(pad_range(d$UMAP_1), pad_range(d$UMAP_2)))
    }
    z <- sort(as.vector(kde$z), decreasing = TRUE)
    level <- z[which(cumsum(z) / sum(z) >= 0.92)[1]]
    pieces <- contourLines(kde$x, kde$y, kde$z, levels = level)
    areas <- vapply(pieces, function(p) abs(sum(
      p$x * c(p$y[-1], p$y[1]) - c(p$x[-1], p$x[1]) * p$y
    ) / 2), numeric(1))
    nested <- vapply(seq_along(pieces), function(i) any(vapply(
      which(areas > areas[i]),
      function(j) sp::point.in.polygon(
        pieces[[i]]$x[1], pieces[[i]]$y[1], pieces[[j]]$x, pieces[[j]]$y
      ) > 0,
      logical(1)
    )), logical(1))
    pieces <- pieces[!nested]
    do.call(rbind, lapply(seq_along(pieces), function(i) data.frame(
      UMAP_1 = pieces[[i]]$x, UMAP_2 = pieces[[i]]$y,
      cell_type = cell_type, piece = i
    )))
  }
  outline_base <- rough_outline("Leukemia")
  stopifnot(identical(unique(outline_base$cell_type), "Leukemia"))
  outline_df <- do.call(rbind, lapply(LAB_LEVELS, function(group_lab) {
    transform(outline_base,
              group_lab = factor(group_lab, levels = LAB_LEVELS))
  }))

  label_df <- data.frame(
    cell_type = CT_LEVELS,
    label = c("Leukemic blasts", "CD8 T", "CD4 T", "B cells", "NK"),
    UMAP_1 = c(-2.43, 11.37, 10.75, 5.09, 12.10),
    UMAP_2 = c(0.28, 0.80, -2.47, -8.03, 2.17),
    label_x = c(-7.2, 8.6, 8.4, 2.7, 9.2),
    label_y = c(6.4, 1.0, -4.2, -8.9, 4.5)
  )
  if (style != "original") {
    # Place the blast label directly above its region, without a leader arrow.
    label_df$label_x[1] <- -5.5
    label_df$label_y[1] <- 7.6
  }
  annotated_types <- if (annotation_mode == "all") CT_LEVELS else "Leukemia"
  label_df <- label_df[label_df$cell_type %in% annotated_types, ]
  label_df <- do.call(rbind, lapply(LAB_LEVELS, function(group_lab) {
    transform(label_df,
              group_lab = factor(group_lab, levels = LAB_LEVELS))
  }))

  border <- switch(style, original = COLORS$cell_type[["Leukemia"]],
                   teal = "#007D83", dashed = "#303030",
                   shaded = "#007D83", minimal = "#303030",
                   navy_bold = "#163E72", violet_dotted = "#64428A",
                   no_arrows = "#007D83", grey_shaded = "#404040",
                   shade_only = "#007D83", all_regions = "#007D83")
  line_type <- switch(style, dashed = "33", violet_dotted = "11", "solid")
  line_width <- switch(style, shaded = 0.35, grey_shaded = 0.3,
                       navy_bold = 0.85, violet_dotted = 0.7, 0.55)
  # Boundary styles share the same data-derived coordinates.
  # A fill sits UNDER the point layer; expression colors are never tinted.
  shade_layers <- if (style %in% c("shaded", "grey_shaded", "shade_only")) list(
    geom_polygon(data = outline_df,
                 aes(UMAP_1, UMAP_2, group = interaction(group_lab, piece)),
                 inherit.aes = FALSE,
                 fill = if (style == "grey_shaded") "#E8E8E8" else "#DCEFED",
                 colour = NA)
  ) else list()
  outline_layers <- list(
    geom_path(data = outline_df,
              aes(UMAP_1, UMAP_2,
                  group = interaction(group_lab, piece)),
              inherit.aes = FALSE, colour = "white", linewidth = line_width + 0.45,
              linetype = line_type),
    geom_path(data = outline_df,
              aes(UMAP_1, UMAP_2,
                  group = interaction(group_lab, piece)),
              inherit.aes = FALSE,
              colour = border, linewidth = line_width, linetype = line_type)
  )
  if (style == "shade_only") outline_layers <- list()
  if (style == "all_regions") {
    # Show the pre-existing immune lineage regions as an alternative to arrows.
    immune_layers <- lapply(setdiff(CT_LEVELS, "Leukemia"), function(ct) {
      region <- rough_outline(ct)
      region <- do.call(rbind, lapply(LAB_LEVELS, function(group_lab)
        transform(region, group_lab = factor(group_lab, levels = LAB_LEVELS))))
      geom_path(data = region,
                aes(UMAP_1, UMAP_2, group = interaction(group_lab, piece)),
                inherit.aes = FALSE, colour = COLORS$cell_type[[ct]],
                linewidth = 0.3, linetype = "33")
    })
    outline_layers <- c(outline_layers, immune_layers)
  }
  leader_types <- if (style == "original") annotated_types else setdiff(annotated_types, "Leukemia")
  if (style %in% c("no_arrows", "all_regions")) leader_types <- character()
  leader_layers <- lapply(leader_types, function(cell_type) {
    d <- label_df[label_df$cell_type == cell_type, ]
    geom_segment(data = d,
                 aes(x = label_x, y = label_y,
                     xend = UMAP_1, yend = UMAP_2),
                 inherit.aes = FALSE,
                 colour = if (style == "grey_shaded") "#555555" else COLORS$cell_type[[cell_type]],
                 linewidth = 0.4,
                 arrow = if (style == "grey_shaded") NULL else
                   grid::arrow(length = grid::unit(0.8, "mm"), type = "closed"))
  })
  df <- df[order(df$Ysig), ]

  # colour clamp (1-99%) so the depletion reads with contrast
  lo <- as.numeric(quantile(df$Ysig, 0.01))
  hi <- as.numeric(quantile(df$Ysig, 0.99))

  # ---- corner arrows: drawn in the LEFTMOST facet only (group == GROUPS[1]) ---
  xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)
  dx <- diff(xr); dy <- diff(yr); sp <- max(dx, dy)
  alen <- 0.16 * sp
  ox <- xr[1] - 0.02 * dx
  oy <- yr[1] - 0.02 * dy
  gap <- 0.05 * sp
  arr <- grid::arrow(length = grid::unit(1.3, "mm"), type = "closed")
  g1  <- factor(group_label("CN-Male"), levels = LAB_LEVELS)
  seg <- data.frame(group_lab = g1, x = c(ox, ox), y = c(oy, oy),
                    xend = c(ox + alen, ox), yend = c(oy, oy + alen))
  lab <- data.frame(group_lab = g1,
                    x = c(ox + alen / 2, ox - gap),
                    y = c(oy - gap, oy + alen / 2),
                    label = c("UMAP 1", "UMAP 2"), angle = c(0, 90))

  p_f3a <- ggplot(df, aes(UMAP_1, UMAP_2, colour = Ysig)) +
    shade_layers +
    ggrastr::rasterise(geom_point(size = 0.4, stroke = 0), dpi = 600) +
    outline_layers +
    leader_layers +
    geom_label(data = label_df,
               aes(x = label_x, y = label_y, label = label),
               inherit.aes = FALSE, family = "Arial", fontface = "bold",
               colour = "black", fill = "white", linewidth = 0,
               size = if (style == "original") 1.8 else 2.1,
               label.padding = unit(0.08, "lines")) +
    facet_wrap(~ group_lab, nrow = if (layout == "vertical") 2 else 1,
               labeller = label_parsed) +
    scale_colour_gradient(low = "lightgrey", high = "darkred",
                          limits = c(lo, hi), oob = scales::squish,
                          name = "Y signature") +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, arrow = arr, linewidth = 0.4, colour = "black") +
    geom_text(data = lab, aes(x = x, y = y, label = label, angle = angle),
              inherit.aes = FALSE, size = 2.3, hjust = 0.5, vjust = 1) +
    coord_equal(clip = "off") +
    theme_nature_contract() +
    theme(
      axis.line  = element_blank(), axis.ticks = element_blank(),
      axis.text  = element_blank(), axis.title = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      legend.position   = "right",
      legend.title      = element_text(size = 8),
      legend.text       = element_text(size = 7),
      legend.key.height = unit(5, "mm"), legend.key.width = unit(3, "mm"),
      panel.spacing     = unit(2, "mm"),
      plot.margin       = margin(3, 3, 9, 9, "mm")
    )

  if (layout == "vertical") {
    p_f3a <- p_f3a +
      guides(colour = guide_colourbar(direction = "horizontal",
             title.position = "top", title.hjust = 0.5,
             barwidth = unit(25, "mm"), barheight = unit(2.5, "mm"))) +
      theme(legend.position = "bottom", panel.spacing = unit(5, "mm"),
            plot.margin = margin(2, 2, 3, 4, "mm"))
  }
  if (return_plot) return(p_f3a)
  save_pub_r(p_f3a, .ff_prep_out(out_basepath),
             width_mm = if (layout == "vertical") 70 else 130,
             height_mm = if (layout == "vertical") 150 else 75)
}

fig_final_F3b_ygene_dotplot <- function(
    seurat_filtered, lineage_annotation,
    out_basepath = file.path(.FF_PANELS_DIR, "fig3/F3b_ygene_dotplot")) {
  so <- seurat_filtered
  # NK-separated lineage annotation gives NK cells their own
  # column rather than lumped into CD8 T. Fallback to cell_type for any barcode
  # missing from the lineage map.
  lin <- setNames(lineage_annotation$labels$lineage, lineage_annotation$labels$barcode)[colnames(so)]
  lin[is.na(lin)] <- as.character(so$cell_type)[is.na(lin)]
  # unname(): lin carries NA names for cells absent from the lineage map (they
  # fail Seurat's named-vector assignment check); values stay positionally aligned.
  so$cell_type <- factor(unname(lin), levels = c("Leukemia", "CD8 T", "CD4 T", "B cell", "NK"))

  F3B_GENES <- c("RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "USP9Y", "UTY")
  CT_KEEP <- c("Leukemia", "CD8 T", "CD4 T", "B cell", "NK")
  cells_keep <- so$group %in% c("LOY", "CN-Male") & !is.na(so$cell_type)
  so <- subset(so, cells = colnames(so)[cells_keep])
  so[["Cell type"]] <- factor(as.character(so$cell_type), levels = CT_KEEP)
  so$Karyotype <- factor(ifelse(so$group == "LOY", "LOY", "CN"),
                         levels = c("LOY", "CN"))

  p_f3b <- scplotter::FeatureStatPlot(
    object = so,
    features = F3B_GENES,
    plot_type = "dot",
    ident = "Karyotype",
    columns_split_by = "Cell type",
    assay = "RNA",
    column_name_annotation = FALSE,
    column_annotation = "Karyotype",
    column_annotation_type = list(Karyotype = "simple"),
    column_annotation_palcolor = list(
      Karyotype = c("LOY" = "#FF1D25", "CN" = "#34358E")
    ),
    column_annotation_agg = list(Karyotype = dplyr::first),
    column_annotation_params = list(
      Karyotype = list(height = grid::unit(2.5, "mm"), border = TRUE,
                       show_legend = TRUE)
    ),
    row_name_annotation = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    add_reticle = TRUE,
    flip = FALSE,
    column_names_rot = 0
  )
  # scplotter builds split and label annotations before returning the
  # ComplexHeatmap object. Keep the cell-type split titles, remove its redundant
  # colour strip, and update both stored label surfaces to true plotmath text.
  h <- p_f3b@ht_list[[1]]
  h@top_annotation <- h@top_annotation[, "Karyotype"]
  h@column_title <- CT_KEEP
  h@column_title_param$gp <- grid::gpar(fontfamily = "Arial", fontsize = 11)
  group_labels <- rep(
    expression(M^{plain(LOY)}, M^{plain(CN)}), length(CT_KEEP)
  )
  h@column_names_param$labels <- group_labels
  assign("value", group_labels, envir = h@column_names_param$anno@var_env)
  p_f3b@ht_list[[1]] <- h
  save_pub_r(function() ComplexHeatmap::draw(p_f3b), .ff_prep_out(out_basepath),
             width_mm = 165, height_mm = 95)
}

fig_final_F3c_ygene_per_patient <- function(
    seurat_filtered,
    out_basepath = file.path(.FF_PANELS_DIR, "fig3/F3c_ygene_per_patient")) {
  so <- seurat_filtered

  # Canonical 9 Y-linked genes (same as F3a; recomputed inline, same seed).
  Y_GENES <- c("UTY", "USP9Y", "ZFY", "EIF1AY", "DDX3Y",
               "RPS4Y1", "KDM5D", "TMSB4Y", "PRKY")
  GROUPS3 <- c("LOY", "CN-Male", "CN-Female")

  y_present <- intersect(Y_GENES, rownames(so))
  if (length(setdiff(Y_GENES, y_present)))
    warning("Y genes absent from object: ", paste(setdiff(Y_GENES, y_present), collapse = ", "))
  set.seed(42)
  so <- AddModuleScore(so, features = list(y_present), name = "Ysig", seed = 42)
  so$Ysig_score <- so$Ysig1

  donor <- get_donor_id(so)
  pp_df <- data.frame(
    donor = donor,
    group = factor(as.character(so$group), levels = GROUPS3),
    cell_type = as.character(so$cell_type),
    Ysig  = so$Ysig_score,
    stringsAsFactors = FALSE
  )
  # LOY S1-S6 + CN-Male contrast; keep the 4 lineages present in both male groups
  # (NK absent in LOY/CN-Male — all NK are CN-Female; Unidentified dropped).
  PP_CT <- c("Leukemia", "CD8 T", "CD4 T", "B cell")
  pp_df <- pp_df[pp_df$group %in% c("LOY", "CN-Male") & pp_df$cell_type %in% PP_CT, ]
  pp_df$cell_type <- factor(pp_df$cell_type, levels = PP_CT)

  # Order donors: LOY S1..S6 first, then CN-Male AML* (group-coloured strips)
  loy_donors <- sort(unique(pp_df$donor[pp_df$group == "LOY"]))
  cnm_donors <- sort(unique(pp_df$donor[pp_df$group == "CN-Male"]))
  pp_df$donor <- factor(pp_df$donor, levels = c(loy_donors, cnm_donors))

  # donor strip fill by group (LOY red / CN-Male blue) via a group label in facet
  pp_df$donor_lab <- paste0(pp_df$donor,
                            ifelse(pp_df$group == "LOY", " (LOY)", " (CN-M)"))
  lab_levels <- c(paste0(loy_donors, " (LOY)"), paste0(cnm_donors, " (CN-M)"))
  pp_df$donor_lab <- factor(pp_df$donor_lab, levels = lab_levels)

  p_f3c <- ggplot(pp_df, aes(x = cell_type, y = Ysig, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, alpha = 0.9) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.18, alpha = 1) +
    facet_wrap(~ donor_lab, ncol = 4) +
    scale_fill_manual(values = COLORS$group[c("LOY", "CN-Male")], guide = "none") +
    labs(x = NULL, y = "Y-signature score") +
    theme_nature_contract(base_size = 7.5) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1, size = 6.2, colour = "black"),
      axis.text.y = element_text(size = 6.2, colour = "black"),
      axis.title.y = element_text(size = 7),
      strip.text  = element_text(size = 6.8, face = "bold"),
      plot.title  = element_blank()
    )

  save_pub_r(p_f3c, .ff_prep_out(out_basepath), width_mm = 183, height_mm = 120)
}

fig_final_F3d_ysig_persample_violin <- function(
    seurat_filtered, lineage_annotation,
    out_basepath = file.path(.FF_PANELS_DIR, "fig3/F3d_ysig_persample_violin")) {
  so <- seurat_filtered
  la <- lineage_annotation

  # Original Y-signature gene set (unchanged from the CITEseq script).
  y_signature_genes <- c(
    "SRY", "TSPYL2", "RPS4Y1", "ZFY", "PRKY", "USP9Y", "DDX3Y",
    "UTY", "TMSB4Y", "KDM5D", "EIF1AY"
  )
  CT_LEVELS <- c("Leukemia", "CD8 T", "CD4 T", "B cell", "NK")  # Leukemia, CD8, CD4, B, NK
  LOY_SAMP  <- paste0("S", 1:6)

  # Y-signature score calculated with AddModuleScore.
  y_present <- intersect(y_signature_genes, rownames(so))
  so <- AddModuleScore(so, features = list(Y_signature = y_present),
                       name = "Y_Signature", assay = "RNA")

  # clean lineage (NK separated) + donor id
  lin <- setNames(la$labels$lineage, la$labels$barcode)
  don <- get_donor_id(so); names(don) <- colnames(so)

  plot_df <- data.frame(
    Y_Signature1 = so$Y_Signature1,
    cell_type    = lin[colnames(so)],
    sample_id    = don[colnames(so)],
    group        = as.character(so$group),
    stringsAsFactors = FALSE
  )
  # change 2: LOY S1-S6 only; change 3: 5 cell types incl NK
  plot_df <- plot_df[plot_df$group == "LOY" &
                     plot_df$sample_id %in% LOY_SAMP &
                     plot_df$cell_type %in% CT_LEVELS, ]
  plot_df$cell_type <- factor(plot_df$cell_type, levels = CT_LEVELS)
  plot_df$sample_id <- factor(plot_df$sample_id, levels = LOY_SAMP)

  # Publication typography and export dimensions.
  # Distribution-only panel. The former four cell-level Wilcoxon brackets per
  # donor were visually redundant (all ****) and n-inflated; the donor-resolved
  # pattern is shown directly without a forest of 24 brackets.
  plot2_by_sample <- ggplot(plot_df, aes(x = cell_type, y = Y_Signature1)) +
    geom_violin(aes(fill = cell_type), scale = "width", trim = TRUE,
                linewidth = 0.22, show.legend = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", linewidth = 0.22) +
    scale_fill_manual(values = COLORS$cell_type[CT_LEVELS]) +
    facet_wrap(~ sample_id, nrow = 1) +
    labs(
      x = NULL,
      y = "Y-signature score"
    ) +
    theme_nature_contract(base_size = 7.5) +
    theme(
      strip.text  = element_text(size = 7, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6.2, colour = "black"),
      axis.text.y = element_text(size = 6.2, colour = "black"),
      axis.title.y = element_text(size = 7)
    )

  save_pub_r(plot2_by_sample, .ff_prep_out(out_basepath), width_mm = 183, height_mm = 80)
}

fig_final_F3e_ygene_dotplot_by_group <- function(
    seurat_filtered, lineage_annotation,
    out_basepath = file.path(.FF_PANELS_DIR, "fig3/F3e_ygene_dotplot_by_group")) {
  so <- seurat_filtered
  la <- lineage_annotation   # clean lineage (NK separated)

  # Display the six listed genes; low-expression ZFY and undetected PRKY
  # are excluded from this panel.
  y_genes_in_data <- c("RPS4Y1", "DDX3Y", "EIF1AY", "UTY", "KDM5D", "USP9Y")

  # 2-group comparison (LOY vs CN-Male); LOY first.
  so <- subset(so, cells = colnames(so)[so$group %in% c("LOY", "CN-Male")])
  so$group <- factor(as.character(so$group), levels = c("LOY", "CN-Male"))

  # Clean lineage so NK is split out of coarse "CD8 T"; safety fallback to cell_type.
  lin <- setNames(la$labels$lineage, la$labels$barcode)[colnames(so)]
  lin[is.na(lin)] <- as.character(so$cell_type)[is.na(lin)]
  ct_levels <- c("Leukemia", "CD8 T", "CD4 T", "B cell", "NK")
  so$cell_type <- factor(lin, levels = ct_levels)

  # per (gene, cell_type, group): % expressing + mean log-norm expr.
  expr   <- as.matrix(Seurat::GetAssayData(so, assay = "RNA", layer = "data")[y_genes_in_data, , drop = FALSE])
  grp_ct <- interaction(as.character(so$group), as.character(so$cell_type),
                        drop = TRUE, sep = "@@")
  rec <- list()
  for (g in y_genes_in_data) {
    ev  <- expr[g, ]
    pct <- tapply(ev, grp_ct, function(x) 100 * mean(x > 0))
    avg <- tapply(ev, grp_ct, mean)
    for (k in names(pct)) {
      parts <- strsplit(k, "@@", fixed = TRUE)[[1]]
      rec[[length(rec) + 1]] <- data.frame(
        gene = g, group = parts[1], cell_type = parts[2],
        pct = as.numeric(pct[[k]]), avg = as.numeric(avg[[k]]),
        stringsAsFactors = FALSE)
    }
  }
  dp <- do.call(rbind, rec)

  # z-score mean expr WITHIN each gene (fill = relative level, not absolute).
  dp$z <- stats::ave(dp$avg, dp$gene, FUN = function(x) {
    s <- stats::sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s
  })
  dp$gene      <- factor(dp$gene, levels = rev(y_genes_in_data))   # RPS4Y1 on top
  dp$cell_type <- factor(dp$cell_type, levels = ct_levels)
  dp$group     <- factor(dp$group, levels = c("LOY", "CN-Male"))

  p <- ggplot(dp, aes(x = group, y = gene)) +
    geom_point(aes(size = pct, fill = z), shape = 21, stroke = 0.2, colour = "grey30") +
    facet_grid(. ~ cell_type) +
    scale_size_area(name = "% expressing", max_size = 6, limits = c(0, 100)) +
    scale_fill_gradient2(name = "mean expr\n(z, per gene)",
                         low = COLORS$diverging[["low"]], mid = "grey95",
                         high = COLORS$diverging[["high"]], midpoint = 0) +
    scale_x_discrete(labels = function(b) parse(text = group_label(b))) +
    labs(x = NULL, y = NULL) +
    theme_nature_contract() +
    theme(
      panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
      axis.text.y      = element_text(face = "italic"),
      panel.spacing    = unit(1.2, "mm"),
      legend.position  = "right"
    )

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 165, height_mm = 62)
}

fig_final_F3supp_y_by_celltype <- function(
    gene_dosage_results,
    out_basepath = file.path(.FF_PANELS_DIR, "fig3/F3supp_y_by_celltype")) {
  # Display the six listed genes; low-expression ZFY, PRKY and TMSB4Y
  # are excluded from this panel.
  Y_GENES     <- c("UTY", "USP9Y", "EIF1AY", "DDX3Y", "RPS4Y1", "KDM5D")
  dosage_df   <- gene_dosage_results$dosage_df
  dosage_ct   <- dosage_df[dosage_df$cell_type %in% c("Leukemia", "CD8 T"), ]
  dosage_ct$group     <- factor(dosage_ct$group, levels = c("LOY", "CN-Male", "CN-Female"))
  dosage_ct$cell_type <- factor(dosage_ct$cell_type, levels = c("Leukemia", "CD8 T"))

  avail_y <- intersect(Y_GENES, colnames(dosage_ct))

  long_supp <- do.call(rbind, lapply(avail_y, function(g) {
    data.frame(group = dosage_ct$group, cell_type = dosage_ct$cell_type,
               gene = g, expr = dosage_ct[[g]])
  }))
  long_supp$gene <- factor(long_supp$gene, levels = avail_y)

  p_f3supp <- ggplot(long_supp, aes(x = group, y = expr, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.22) +
    geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9, linewidth = 0.22) +
    scale_fill_manual(values = COLORS$group) +
    scale_x_discrete(labels = function(b) parse(text = group_label(b))) +
    facet_grid(gene ~ cell_type, scales = "free_y") +
    labs(x = NULL, y = "log-norm expr") +
    theme_nature_contract() +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(angle = 40, hjust = 1, size = 4.5),
      strip.text.x = element_text(size = 5.5, face = "bold"),
      strip.text.y = element_text(size = 4.5, face = "italic"),
      plot.margin  = margin(3, 3, 2, 3)
    )

  save_pub_r(p_f3supp, .ff_prep_out(out_basepath), width_mm = 100, height_mm = 135)
}

fig_final_F4B_blast_gsea_hallmark <- function(
    deg_pathway_enrichment,
    out_basepath = file.path(.FF_PANELS_DIR, "fig4/F4B_blast_gsea_hallmark")) {
  tbl <- deg_pathway_enrichment$combined

  # Filter: Leukemia blast, deseq2, hallmark — top 20 by padj
  df <- tbl[tbl$cell_type == "Leukemia (blast)" &
             tbl$method == "deseq2" &
             tbl$collection == "hallmark" &
             !is.na(tbl$padj), ]
  df <- dplyr::slice_head(dplyr::arrange(df, padj), n = 20)

  stopifnot(nrow(df) > 0)

  df$pathway <- gsub("^(HALLMARK_|GOBP_)", "", df$pathway)
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))
  df$neglog10padj <- -log10(pmax(df$padj, 1e-300))

  p <- ggplot(df, aes(x = NES, y = pathway, size = size, color = neglog10padj)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.35) +
    scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "-log10 padj") +
    scale_size_continuous(range = c(1, 4), name = "set size") +
    labs(
      x = expression("NES (" * M^{plain(LOY)} * " vs " * M^{plain(CN)} * ")"),
      y = NULL
    ) +
    theme_nature_contract(base_size = 7.5) +
    theme(axis.text.y = element_text(size = 7))

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 110, height_mm = 90)
}

fig_final_F4D_blast_mhc_ap <- function(
    seurat_filtered, degs_by_celltype,
    out_basepath = file.path(.FF_PANELS_DIR, "fig4/F4D_blast_mhc_ap")) {
  so           <- seurat_filtered
  blast_deseq2 <- degs_by_celltype[["blast"]]$deseq2

  AP_GENES <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "CIITA",
                "HLA-DPA1", "HLA-DQA2", "HLA-DQB2", "HLA-DMA", "HLA-DRA", "CD74")

  p_f4d <- ff_violin_compare(
    seurat     = so,
    genes      = AP_GENES,
    cell_type  = "Leukemia",
    deg_deseq2 = blast_deseq2,
    groups     = c("LOY", "CN-Male"),
    ncol       = 6,
    base_size  = 8,
    group_labels = c("LOY" = "M<sup>LOY</sup>",
                     "CN-Male" = "M<sup>CN</sup>")
  )

  save_pub_r(p_f4d, .ff_prep_out(out_basepath), width_mm = 183, height_mm = 70)
}

fig_final_F4E_nlrc5_ciita <- function(
    seurat_filtered, degs_by_celltype,
    out_basepath = file.path(.FF_PANELS_DIR, "fig4/F4E_nlrc5_ciita")) {
  so           <- seurat_filtered
  blast_deseq2 <- degs_by_celltype[["blast"]]$deseq2

  p_f4e <- ff_violin_compare(
    seurat     = so,
    genes      = c("NLRC5", "CIITA"),
    cell_type  = "Leukemia",
    deg_deseq2 = blast_deseq2,
    groups     = c("LOY", "CN-Male"),
    ncol       = 2,
    base_size  = 8,
    group_labels = c("LOY" = "M<sup>LOY</sup>",
                     "CN-Male" = "M<sup>CN</sup>")
  )

  save_pub_r(p_f4e, .ff_prep_out(out_basepath), width_mm = 70, height_mm = 55)
}

fig_final_F4F_blast_apm_score <- function(
    blast_cue_scores,
    out_basepath = file.path(.FF_PANELS_DIR, "fig4/F4F_blast_apm_score")) {
  scores_df <- blast_cue_scores$scores
  stopifnot(all(c("donor", "group", "blast_apm_score") %in% colnames(scores_df)))

  scores_plot <- scores_df
  scores_plot$group <- factor(scores_plot$group, levels = c("LOY", "CN-Male", "CN-Female"))
  scores_plot <- scores_plot[!is.na(scores_plot$blast_apm_score), ]

  loy_vals <- scores_plot$blast_apm_score[scores_plot$group == "LOY"]
  cnm_vals <- scores_plot$blast_apm_score[scores_plot$group == "CN-Male"]
  wt <- wilcox.test(loy_vals, cnm_vals, exact = FALSE)
  wt_p <- wt$p.value
  p_label <- if (wt_p < 0.001) formatC(wt_p, format = "e", digits = 1) else
               formatC(wt_p, format = "f", digits = 3)

  group_cols_f <- COLORS$group[c("LOY", "CN-Male", "CN-Female")]

  y_min   <- min(scores_plot$blast_apm_score, na.rm = TRUE)
  y_max   <- max(scores_plot$blast_apm_score, na.rm = TRUE)
  y_range <- y_max - y_min
  ann_y   <- y_max + 0.08 * y_range
  y_top   <- y_max + 0.30 * y_range

  p_F <- ggplot(scores_plot,
                aes(x = group, y = blast_apm_score, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.75, linewidth = 0.3,
                colour = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 1, linewidth = 0.3,
                 fill = "white", colour = "grey30") +
    ggbeeswarm::geom_beeswarm(size = 1.2, cex = 2.5, alpha = 0.9, shape = 21,
                              stroke = 0.4, colour = "grey20") +
    annotate("segment",
             x = 1, xend = 2,
             y = ann_y, yend = ann_y,
             linewidth = 0.35, colour = "black") +
    annotate("text",
             x = 1.5, y = ann_y + 0.04 * y_range,
             label = paste0("p = ", p_label),
             size = 1.8, hjust = 0.5, vjust = 0, colour = "black", family = "Arial") +
    scale_fill_manual(values = group_cols_f, na.value = "grey70") +
    coord_cartesian(ylim = c(NA, y_top)) +
    labs(
      x = NULL,
      y = "Blast APM module score (donor-level)"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(size = 5.5, colour = "black"),
      axis.text.y  = element_text(size = 5,   colour = "black"),
      axis.title.y = element_text(size = 5.5),
      plot.margin  = margin(4, 5, 3, 4)
    )

  save_pub_r(p_F, .ff_prep_out(out_basepath), width_mm = 55, height_mm = 55)
}

fig_final_F4opt_scenic_blast_tf <- function(
    scenic_blast_diff_file,
    out_basepath = file.path(.FF_PANELS_DIR, "fig4/F4opt_scenic_blast_tf")) {
  requireNamespace("data.table", quietly = TRUE)

  diff_file <- scenic_blast_diff_file
  n_top <- 15
  cell_type <- "Leukemia (blast)"

  d <- data.table::fread(diff_file)
  sig <- d[!is.na(padj) & padj < 0.05]
  sig <- sig[grepl("_act$", regulon)]

  top_up <- sig[log2FC > 0][order(-log2FC)][seq_len(min(n_top, .N))]
  top_dn <- sig[log2FC < 0][order(log2FC)][seq_len(min(n_top, .N))]
  top <- rbind(top_up, top_dn)

  top[, tf_name := gsub("_act$", "", regulon)]
  top[, direction := ifelse(log2FC > 0, "UP in LOY", "DOWN in LOY")]
  top$tf_name <- factor(top$tf_name, levels = top$tf_name[order(top$log2FC)])

  total_sig <- nrow(sig)

  p <- ggplot(top, aes(x = log2FC, y = tf_name, fill = direction)) +
    geom_col(width = 0.65) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.35) +
    scale_fill_manual(
      values = c("UP in LOY" = COLORS$group[["LOY"]],
                 "DOWN in LOY" = COLORS$group[["CN-Male"]]),
      name = NULL
    ) +
    labs(
      x = "log2FC (AUCell regulon activity)",
      y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      legend.position = "top",
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = 6)
    )

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 80, height_mm = 100)
}

fig_final_F5C1_state_proportions <- function(
    cd8_states,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5C1_state_proportions")) {
  STATE_LVL  <- c("Naive", "Effector", "Dysfunctional")
  STATE_COLS <- COLORS$cd8_state

  props <- cd8_states$proportions
  props <- props[props$group %in% c("LOY", "CN-Male"), , drop = FALSE]
  props$group          <- factor(props$group, levels = c("CN-Male", "LOY"))
  props$assigned_state <- factor(props$assigned_state, levels = STATE_LVL)

  p_f5c1 <- ggplot(props,
            aes(x = group, y = prop * 100, fill = assigned_state)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prop * 100)),
              position = position_dodge(width = 0.75),
              vjust = -0.4, size = 2.3) +
    scale_fill_manual(values = STATE_COLS, name = "State") +
    scale_x_discrete(labels = function(b) parse(text = group_label(b))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "% of CD8 T cells", title = NULL) +
    theme_nature_contract(base_size = 8) +
    theme(legend.position = "right")

  save_pub_r(p_f5c1, .ff_prep_out(out_basepath),
             width_mm = 70, height_mm = 55)
}

fig_final_F5E1_bach2_violin_celltype <- function(
    seurat_filtered, lineage_annotation,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5E1_bach2_violin_celltype")) {
  so <- seurat_filtered

  obj_b2 <- subset(so, group %in% c("LOY", "CN-Male"))
  meta_b2 <- obj_b2@meta.data
  # Compare clean CD8 T cells with B cells. Clean lineage labels exclude
  # marker-defined NK cells from the CD8 compartment.
  lin <- setNames(lineage_annotation$labels$lineage,
                  lineage_annotation$labels$barcode)[colnames(obj_b2)]
  lin[is.na(lin)] <- as.character(meta_b2$cell_type)[is.na(lin)]
  meta_b2$BACH2 <- as.numeric(LayerData(obj_b2[["RNA"]], layer = "data")["BACH2", ])
  meta_b2$celltype <- factor(unname(lin), levels = c("CD8 T", "B cell"))
  meta_b2 <- meta_b2[!is.na(meta_b2$celltype), ]
  meta_b2$celltype <- droplevels(meta_b2$celltype)
  meta_b2$group <- factor(meta_b2$group, levels = c("LOY", "CN-Male"))

  B2_COLS <- c("LOY" = "#E41A1C", "CN-Male" = "#377EB8")

  p_f5e1 <- ggplot(meta_b2, aes(x = celltype, y = BACH2, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.25,
                position = position_dodge(0.8)) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6, linewidth = 0.25,
                 position = position_dodge(0.8)) +
    scale_fill_manual(values = B2_COLS,
                      labels = function(b) parse(text = group_label(b))) +
    labs(x = NULL, y = expression(italic(BACH2) ~ "(log-norm)"), fill = "Group") +
    theme_nature_contract(base_size = 8) +
    theme(
      axis.text.x     = element_text(angle = 0, hjust = 0.5, size = 7.5),
      legend.position = "top",
      legend.key.size = unit(2.5, "mm")
    )

  save_pub_r(p_f5e1, .ff_prep_out(out_basepath), width_mm = 75, height_mm = 65)
}

fig_final_F5E2_bach2_violin_state <- function(
    cd8_states, degs_by_cd8_state,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5E2_bach2_violin_state")) {
  # BACH2 expression in naive/stem-like CD8 cells, LOY versus male controls.
  cd8_b <- cd8_states$cd8_obj
  keep  <- cd8_b$group %in% c("LOY", "CN-Male") &
           as.character(cd8_b$assigned_state) == "Naive"
  cd8_b <- subset(cd8_b, cells = colnames(cd8_b)[keep])

  df <- data.frame(
    group = factor(as.character(cd8_b$group), levels = c("LOY", "CN-Male")),
    expr  = as.numeric(LayerData(cd8_b[["RNA"]], layer = "data")["BACH2", ]),
    stringsAsFactors = FALSE
  )

  # Naive-state pseudobulk DESeq2 stat for BACH2
  lab <- NULL
  d <- degs_by_cd8_state$Naive$deseq2
  if (!is.null(d) && "BACH2" %in% d$gene) {
    r <- d[d$gene == "BACH2", ]
    lab <- sprintf("log2FC %+.2f, padj %.0e", r$log2FC[1], r$padj[1])
  }

  p_f5e2 <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = expr, fill = group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.25) +
    ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.95, linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = COLORS$group[c("LOY", "CN-Male")], guide = "none") +
    ggplot2::scale_x_discrete(labels = function(b) parse(text = group_label(b))) +
    {
      if (!is.null(lab))
        ggplot2::annotate("text", x = 1.5, y = max(df$expr, na.rm = TRUE) * 1.03,
                          label = lab, size = 2.0, family = "Arial",
                          fontface = "plain", colour = "black")
      else NULL
    } +
    ggplot2::labs(x = NULL, y = "BACH2 (log-norm)") +
    theme_nature_contract(base_size = 8, base_family = "Arial") +
    ggplot2::theme(legend.position = "none")

  save_pub_r(p_f5e2, .ff_prep_out(out_basepath), width_mm = 52, height_mm = 62,
             family = "Arial")
}

fig_final_F5E3_bach2_regulon <- function(
    scenic_cd8_by_state_imported, cd8_states,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5E3_bach2_regulon")) {
  # BACH2 activating-regulon AUCell score in naive/stem-like CD8 cells.
  # Match precomputed AUCell rows to the clean CD8 metadata before plotting.
  naive_res <- scenic_cd8_by_state_imported$Naive
  auc_nm <- grep("auc_per_cell", names(naive_res), value = TRUE)
  if (length(auc_nm) != 1) stop("Expected one Naive AUCell per-cell table.")
  auc <- naive_res[[auc_nm]]
  if (!("Cell" %in% colnames(auc))) stop("Naive AUCell table lacks Cell column.")

  bcols <- grep("^BACH2", colnames(auc), value = TRUE)
  if (length(bcols) == 0) stop("No BACH2 regulon column in SCENIC CSV.")
  bcol  <- if ("BACH2_act" %in% bcols) "BACH2_act" else bcols[1]   # the (+) regulon

  cd8_meta <- cd8_states$cd8_obj@meta.data
  cells <- intersect(as.character(auc$Cell), rownames(cd8_meta))
  cells <- cells[
    cd8_meta[cells, "group"] %in% c("LOY", "CN-Male") &
    as.character(cd8_meta[cells, "assigned_state"]) == "Naive"
  ]
  if (!length(cells)) stop("No clean Naive CD8 cells overlap the AUCell table.")
  auc_idx <- match(cells, as.character(auc$Cell))

  df <- data.frame(
    group = factor(as.character(cd8_meta[cells, "group"]), levels = c("LOY", "CN-Male")),
    auc   = as.numeric(auc[auc_idx, bcol]),
    stringsAsFactors = FALSE
  )

  p_f5e3 <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = auc, fill = group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.25) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = COLORS$group[c("LOY", "CN-Male")], guide = "none") +
    ggplot2::scale_x_discrete(labels = function(b) parse(text = group_label(b))) +
    ggplot2::labs(x = NULL, y = "BACH2 regulon (AUCell)") +
    theme_nature_contract(base_size = 8, base_family = "Arial") +
    ggplot2::theme(legend.position = "none")

  save_pub_r(p_f5e3, .ff_prep_out(out_basepath), width_mm = 52, height_mm = 62,
             family = "Arial")
}

fig_final_F5b_naive_volcano <- function(
    degs_by_cd8_state_clean,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5b_naive_volcano")) {
  v <- degs_by_cd8_state_clean$Naive$deseq2
  v <- v[!is.na(v$padj) & !is.na(v$log2FC), ]
  v$nlp <- -log10(pmax(v$padj, 1e-300))
  fc_thr <- 1; p_thr <- 0.05
  v$dir <- ifelse(v$padj < p_thr & v$log2FC >=  fc_thr, "Up in LOY",
           ifelse(v$padj < p_thr & v$log2FC <= -fc_thr, "Down in LOY", "NS"))
  v$dir <- factor(v$dir, levels = c("Up in LOY", "Down in LOY", "NS"))
  cols <- c("Up in LOY"   = unname(COLORS$group["LOY"]),
            "Down in LOY" = unname(COLORS$group["CN-Male"]), "NS" = "grey80")
  # cap y so a few astronomically-significant genes don't flatten the rest
  ycap <- as.numeric(stats::quantile(v$nlp[v$dir != "NS"], 0.99, na.rm = TRUE))
  if (!is.finite(ycap) || ycap <= 0) ycap <- max(v$nlp, na.rm = TRUE)
  v$nlp_d <- pmin(v$nlp, ycap)

  # Label biologically interesting CD8 genes, not the MT/ribosomal housekeeping
  # genes that dominate the top-by-significance list. Curated set = naive/stem &
  # memory TFs, effector/cytotoxic, exhaustion, signaling/TGFb; plus the top
  # non-housekeeping hits each direction to fill in.
  GOI <- c("BACH2","LEF1","FOXO1","TCF7","CCR7","SELL","IL7R","KLF2","KLF3","MYB",
           "ID3","SATB1","TXNIP","STAT1","STAT4","TGFBR2","SMAD7","IFITM1","KLF12",
           "TOX","TNFRSF9","BATF","GZMA","GZMK","GZMB","GZMH","PRF1","NKG7","GNLY",
           "IFNG","TBX21","EOMES","KLRG1","CX3CR1","RUNX3","ZNF683","ITGAE")
  hk  <- "^(MT-|RP[SL]|MRP[SL]|EEF|EIF[0-9]|FAU|MALAT1|XIST|HB[AB])"
  sig <- v[v$dir != "NS", ]
  # Keep curated labels only when they clear a stronger visual threshold. This
  # preserves the key biology while avoiding the low-y label pile-up around the
  # significance line at the final slide-panel size.
  curated <- sig[sig$gene %in% GOI & sig$nlp >= 2.5, ]
  nonhk   <- sig[!grepl(hk, sig$gene), ]
  up <- utils::head(nonhk[nonhk$log2FC > 0, ][order(-nonhk$nlp[nonhk$log2FC > 0]), ], 6)
  dn <- utils::head(nonhk[nonhk$log2FC < 0, ][order(-nonhk$nlp[nonhk$log2FC < 0]), ], 6)
  lab <- unique(rbind(curated, up, dn))

  p <- ggplot(v, aes(log2FC, nlp_d, colour = dir)) +
    geom_point(data = v[v$dir == "NS", ], size = 0.8, stroke = 0, alpha = 0.75) +
    geom_point(data = v[v$dir != "NS", ], size = 1.2, stroke = 0, alpha = 0.95) +
    geom_point(data = v[v$gene == "BACH2", ], aes(x = log2FC, y = nlp_d),
               inherit.aes = FALSE, shape = 21, size = 2.0, stroke = 0.35,
               fill = unname(cols["Up in LOY"]), colour = "black") +
    geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed",
               linewidth = 0.2, colour = "grey60") +
    geom_hline(yintercept = -log10(p_thr), linetype = "dashed",
               linewidth = 0.2, colour = "grey60") +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label = gene,
          size = ifelse(gene == "BACH2", 3.2, 2.5),
          fontface = ifelse(gene == "BACH2", "bold.italic", "italic")),
                             max.overlaps = Inf, segment.size = 0.15,
                             family = "Arial", colour = "black", min.segment.length = 0,
                             box.padding = 0.45, point.padding = 0.18,
                             force = 2, seed = 42, show.legend = FALSE) +
    scale_size_identity() +
    scale_colour_manual(values = cols, name = NULL,
                        breaks = c("Up in LOY", "Down in LOY")) +
    labs(x = expression("log2FC  (" * M^{plain(LOY)} * " vs " * M^{plain(CN)} * ")"),
         y = expression(-log[10] ~ padj)) +
    theme_nature_contract(base_size = 9.5, base_family = "Arial") +
    theme(legend.position = "bottom",
          plot.margin = margin(2, 3, 1, 2, unit = "mm"))

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 115, height_mm = 82,
             family = "Arial")
}

fig_final_F5g_naive_pathway <- function(
    deg_pathway_enrichment_by_cd8_state_clean,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5g_naive_pathway")) {
  pw <- deg_pathway_enrichment_by_cd8_state_clean$by_state$Naive
  pw <- pw[pw$collection == "hallmark" & !is.na(pw$padj) & pw$padj < 0.05, ]
  if (!nrow(pw)) stop("F5g: no significant Hallmark pathways for Naive state")
  pw <- pw[order(-abs(pw$NES)), ]
  pw <- pw[!duplicated(pw$pathway), ]          # one row per pathway
  pw <- utils::head(pw, 15)
  pw$lab <- gsub("_", " ", sub("^HALLMARK_", "", pw$pathway))
  pw$lab <- factor(pw$lab, levels = pw$lab[order(pw$NES)])
  pw$dir <- factor(ifelse(pw$NES > 0, "Up in LOY", "Down in LOY"),
                   levels = c("Up in LOY", "Down in LOY"))
  cols <- c("Up in LOY"   = unname(COLORS$group["LOY"]),
            "Down in LOY" = unname(COLORS$group["CN-Male"]))

  p <- ggplot(pw, aes(NES, lab, fill = dir)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = expression("NES  (" * M^{plain(LOY)} * " vs " * M^{plain(CN)} * ")"), y = NULL) +
    theme_nature_contract(base_size = 8, base_family = "") +
    theme(axis.text.y = element_text(size = 7, colour = "black"),
          legend.position = "bottom")

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 115, height_mm = 85, family = NULL)
}

fig_final_F5g2_naive_pathway_gobp <- function(
    deg_pathway_enrichment_by_cd8_state_clean,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5g2_naive_pathway_gobp")) {
  pw <- deg_pathway_enrichment_by_cd8_state_clean$by_state$Naive
  pw <- pw[pw$collection == "gobp" & !is.na(pw$padj) & pw$padj < 0.05, ]
  keep <- "IMMUN|T_CELL|LYMPHOCYTE|LEUKOCYTE|CYTOKINE|INTERFERON|INTERLEUKIN|CHEMOKINE|TUMOR_NECROSIS|CYTOTOX|CELL_KILLING|KILLING_OF|DIFFERENTIATION|MIGRATION|CHEMOTAXIS|ACTIVATION|TGF|SMAD|WNT|STEM_CELL|MEMORY|EXHAUST|ADHESION|APOPTO|TRANSLATION|OXIDATIVE_PHOSPHOR|RESPIRATION"
  drop <- "SKELET|NEURON|NEURAL|CRANIAL|EMBRYONIC|EPITHELI|MUSCLE|CARDIAC|HEART|RENAL|KIDNEY|EYE|AXON|SYNAP|MORPHOGENESIS|MESENCHYM|OSSIF|BONE|CARTILAGE|LIMB|HAIR|PIGMENT"
  pw <- pw[grepl(keep, pw$pathway) & !grepl(drop, pw$pathway), ]
  if (!nrow(pw)) stop("F5g2: no curated significant GOBP pathways for Naive")
  pw <- pw[!duplicated(pw$pathway), ]
  # take the strongest up and down (by NES), up to 8 each, so both directions show
  up <- utils::head(pw[order(-pw$NES), ], 8)
  dn <- utils::head(pw[order(pw$NES), ], 8)
  pw <- unique(rbind(up, dn))
  pw$lab <- gsub("_", " ", sub("^GOBP_", "", pw$pathway))
  # Full publication labels: wrap long GO terms rather than truncating them.
  pw$lab <- stringr::str_wrap(pw$lab, width = 38)
  pw$lab <- factor(pw$lab, levels = pw$lab[order(pw$NES)])
  pw$dir <- factor(ifelse(pw$NES > 0, "Up in LOY", "Down in LOY"),
                   levels = c("Up in LOY", "Down in LOY"))
  cols <- c("Up in LOY"   = unname(COLORS$group["LOY"]),
            "Down in LOY" = unname(COLORS$group["CN-Male"]))

  p <- ggplot(pw, aes(NES, lab, fill = dir)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = expression("NES  (" * M^{plain(LOY)} * " vs " * M^{plain(CN)} * ")"), y = NULL) +
    theme_nature_contract(base_size = 8, base_family = "") +
    theme(axis.text.y = element_text(size = 6.2, colour = "black", lineheight = 0.86),
          legend.position = "bottom")

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 140, height_mm = 105, family = NULL)
}

fig_final_F5d_naive_scenic <- function(
    scenic_cd8_by_state_imported,
    out_basepath = file.path(.FF_PANELS_DIR, "fig5/F5d_naive_scenic")) {
  d <- scenic_cd8_by_state_imported$Naive[["CD8_T_Naive_LOY_vs_CN_diff_regulons.csv"]]
  required_columns <- c("regulon", "padj", "log2FC")
  if (is.null(d) || !all(required_columns %in% names(d))) {
    stop("Supply Naive/CD8_T_Naive_LOY_vs_CN_diff_regulons.csv with ",
         "regulon, padj and log2FC columns in the CD8 regulon directory.")
  }
  d <- d[grepl("_act$", d$regulon) & !is.na(d$padj) & d$padj < 0.05, ]
  d$tf <- sub("_act$", "", d$regulon)
  # drop housekeeping / uninformative TF families
  d <- d[!grepl("^(RP[SL]|MRP[SL]|MT-|EEF|EIF[0-9]|FAU|ZNF[0-9]|ZSCAN|ZFP|ZBTB[0-9])", d$tf), ]
  up <- utils::head(d[order(-d$log2FC), ], 10)
  dn <- utils::head(d[order(d$log2FC), ], 10)
  reg <- unique(rbind(up, dn))
  reg$tf  <- factor(reg$tf, levels = reg$tf[order(reg$log2FC)])
  reg$dir <- factor(ifelse(reg$log2FC > 0, "Up in LOY", "Up in CN-Male"),
                    levels = c("Up in LOY", "Up in CN-Male"))
  cols <- c("Up in LOY"     = unname(COLORS$group["LOY"]),
            "Up in CN-Male" = unname(COLORS$group["CN-Male"]))

  p <- ggplot(reg, aes(log2FC, tf, fill = dir)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = expression("regulon activity log2FC  (" * M^{plain(LOY)} * " vs " * M^{plain(CN)} * ")"), y = NULL) +
    theme_nature_contract(base_size = 8, base_family = "") +
    theme(axis.text.y = element_text(size = 7, colour = "black"),
          legend.position = "bottom")

  save_pub_r(p, .ff_prep_out(out_basepath), width_mm = 100, height_mm = 95, family = NULL)
}

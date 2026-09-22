# Marker-based NK annotation within the lymphoid compartment.

annotate_nk_lymphoid <- function(seurat_filtered) {
  set.seed(1)
  NK_RESTRICTED <- c("KLRF1", "NCR1", "NCR3", "KLRD1")
  T_MARKERS     <- c("CD3D", "CD3E", "CD3G")
  CD8_MARKERS   <- c("CD8A", "CD8B")
  LYMPHOID_CT   <- c("B cell", "CD4 T", "CD8 T", "NK", "Unidentified")

  cat("=== MODULE NK: re-cluster lymphoid compartment (LOY vs CN-Male) ===\n")
  obj <- seurat_filtered
  Seurat::DefaultAssay(obj) <- "RNA"

  # Lymphoid (non-blast) compartment, ALL THREE GROUPS. CN-Female is included
  # ONLY to ANCHOR the NK cluster boundary during clustering: it carries the
  # bulk of fine-validated NK, and dropping it BEFORE clustering collapses NK
  # purity (validation precision ~98% -> ~65%, the cluster absorbs NKR+ CD8).
  # CN-Female is removed immediately after the NK call (below) and enters NO
  # subsequent male-group comparison, proportion, DEG or pathway analysis.
  lym <- subset(obj, subset = cell_type %in% LYMPHOID_CT)
  cat(sprintf("Lymphoid (all 3 groups, non-blast; CN-Female = anchor only): %d cells\n",
              ncol(lym)))
  print(table(lym$cell_type, lym$group))

  # Deterministic re-cluster recipe.
  lym <- Seurat::FindVariableFeatures(lym, nfeatures = 2000, verbose = FALSE)
  lym <- Seurat::ScaleData(lym, verbose = FALSE)
  lym <- Seurat::RunPCA(lym, npcs = 30, verbose = FALSE)
  lym <- harmony::RunHarmony(lym, group.by.vars = "orig.ident.x", verbose = FALSE)
  lym <- Seurat::FindNeighbors(lym, reduction = "harmony", dims = 1:30, verbose = FALSE)
  lym <- Seurat::FindClusters(lym, resolution = 0.5, verbose = FALSE)
  lym <- Seurat::RunUMAP(lym, reduction = "harmony", dims = 1:30, verbose = FALSE)
  cat(sprintf("\nClusters at res=0.5: %d\n", length(unique(lym$seurat_clusters))))

  # Per-cluster mean expression of NK-restricted / T / CD8 panels.
  expr <- Seurat::GetAssayData(lym, layer = "data")
  mean_in <- function(genes, cells) {
    genes <- intersect(genes, rownames(expr))
    if (!length(genes)) return(NA_real_)
    round(mean(Matrix::colMeans(expr[genes, cells, drop = FALSE])), 3)
  }
  cl <- sort(unique(as.integer(as.character(lym$seurat_clusters))))
  score_tbl <- data.frame()
  for (c in cl) {
    cells <- colnames(lym)[lym$seurat_clusters == as.character(c)]
    score_tbl <- rbind(score_tbl, data.frame(
      cluster = c, n = length(cells),
      NK_restricted = mean_in(NK_RESTRICTED, cells),
      Tcell = mean_in(T_MARKERS, cells),
      CD8 = mean_in(CD8_MARKERS, cells)
    ))
  }
  cat("\n--- per-cluster mean expression (NK_restricted / Tcell / CD8) ---\n")
  print(score_tbl, row.names = FALSE)

  # Programmatic NK call.
  score_tbl$is_NK <- score_tbl$NK_restricted > score_tbl$Tcell &
    score_tbl$NK_restricted > 0.2 &
    score_tbl$CD8 < score_tbl$NK_restricted
  nk_clusters <- score_tbl$cluster[score_tbl$is_NK]
  cat(sprintf("\nNK-called cluster(s) (NK_restricted>Tcell & >0.2 & >CD8): %s\n",
              paste(nk_clusters, collapse = ", ")))

  lym$nk_call <- ifelse(lym$seurat_clusters %in% as.character(nk_clusters),
                        "NK", "non-NK")
  cat("\n=== NK membership x group (all 3; CN-Female anchored clustering) ===\n")
  print(table(lym$nk_call, lym$group))

  # Drop the CN-Female anchor now. The UMAP/cluster/nk_call labels are stored
  # per-cell, so the NK cluster stays well-anchored, but CN-Female enters no
  # comparison hereafter (everything below is LOY vs CN-Male).
  lym <- subset(lym, subset = group %in% c("LOY", "CN-Male"))
  lym$group <- droplevels(factor(lym$group))
  cat(sprintf("After dropping CN-Female anchor: %d cells (LOY+CN-Male)\n", ncol(lym)))

  # Canonical unified donor id (R/plot_style.R get_donor_id): donor_id -> sample_id
  # -> orig.ident with an NA guard. LOY donor_id is NA so it resolves to
  # sample_id S1..S6; CN-Male resolves via donor_id.
  lym$donor <- get_donor_id(lym)

  # Per-group NK x group table (LOY vs CN-Male).
  per_group <- table(lym$nk_call, lym$group)
  cat("\n=== NK cluster membership x group ===\n")
  print(per_group)
  cat("\n=== NK as % of lymphoid, per group ===\n")
  for (g in colnames(per_group)) {
    cat(sprintf("  %-9s NK=%4d / %5d  (%.1f%%)\n",
                g, per_group["NK", g], sum(per_group[, g]),
                100 * per_group["NK", g] / sum(per_group[, g])))
  }

  # Per-donor data.frame.
  donors <- unique(lym$donor)
  per_donor <- do.call(rbind, lapply(donors, function(d) {
    idx <- lym$donor == d
    grp <- unique(lym$group[idx])[1]
    n_lymph <- sum(idx)
    n_nk <- sum(idx & lym$nk_call == "NK")
    data.frame(donor = d, group = grp, n_nk = n_nk, n_lymph = n_lymph,
               pct = round(100 * n_nk / n_lymph, 1), stringsAsFactors = FALSE)
  }))
  # Order: LOY donors (S1..S6) then CN-Male donors.
  per_donor <- per_donor[order(per_donor$group != "LOY", per_donor$donor), ]
  rownames(per_donor) <- NULL
  cat("\n=== NK per donor (raw + % of donor lymphoid) ===\n")
  print(per_donor, row.names = FALSE)

  # Abundance summary — donor-level medians (robust to the AML0114 outlier).
  med_loy <- median(per_donor$pct[per_donor$group == "LOY"])
  med_cn  <- median(per_donor$pct[per_donor$group == "CN-Male"])
  cat(sprintf(
    "\nDonor-level median NK fraction: LOY %.1f%% vs CN-Male %.1f%%.\n",
    med_loy, med_cn))
  cat("CAVEAT: CN-Male NK is single-donor-dominated (AML0114) and LOY is fully\n")
  cat("  pooled (library confound). Abundance is NOT robustly testable in either\n")
  cat("  direction. Report as 'NK PRESENT IN BOTH GROUPS; abundance not reliably\n")
  cat("  comparable here' — NOT depletion and NOT equality.\n")

  # Validation vs CN-Male fine Cell_type_identity NK (only ground truth available).
  cm <- lym$group == "CN-Male"
  fine_nk <- !is.na(lym$Cell_type_identity) & grepl("NK", lym$Cell_type_identity)
  my_nk <- lym$nk_call == "NK"
  tp <- sum(cm & fine_nk & my_nk)
  fn <- sum(cm & fine_nk & !my_nk)
  fp <- sum(cm & !fine_nk & my_nk)
  n_fine_nk <- sum(cm & fine_nk)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  cat("\n=== Validation vs CN-Male fine Cell_type_identity NK ===\n")
  cat(sprintf("  fine-annotated NK (CN-Male): %d\n", n_fine_nk))
  cat(sprintf("  sensitivity (recovered): %d/%d = %.0f%%\n",
              tp, tp + fn, 100 * sensitivity))
  cat(sprintf("  precision: %d/%d = %.0f%%\n",
              tp, tp + fp, 100 * precision))

  list(
    obj = lym,
    nk_cluster_id = nk_clusters,
    per_group = per_group,
    per_donor = per_donor,
    validation = list(sensitivity = sensitivity, precision = precision,
                      n_fine_nk = n_fine_nk),
    markers_used = list(NK_restricted = NK_RESTRICTED, T = T_MARKERS,
                        CD8 = CD8_MARKERS),
    score_table = score_tbl,
    generated_at = Sys.time()
  )
}

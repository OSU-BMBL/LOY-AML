# NK-aware lineage annotation and clean CD8 state classification.

annotate_lineage <- function(seurat_filtered, nk_annotated) {
  nkbc <- colnames(nk_annotated$obj)[nk_annotated$obj$nk_call == "NK"]
  obj <- seurat_filtered
  keep <- obj$group %in% c("LOY", "CN-Male")
  bc <- colnames(obj)[keep]
  coarse <- as.character(obj$cell_type[keep])
  lineage <- ifelse(bc %in% nkbc, "NK", coarse)  # only NK cells move
  labels <- data.frame(barcode = bc, group = as.character(obj$group[keep]),
                       coarse = coarse, lineage = lineage, stringsAsFactors = FALSE)

  cat("=== Stage 1: lineage = coarse cell_type, NK-cluster cells -> NK ===\n")
  cat("Reassigned (coarse -> NK):\n")
  print(table(coarse_was = coarse[bc %in% nkbc], group = labels$group[bc %in% nkbc]))
  cat("\n=== lineage x group ===\n")
  print(table(labels$lineage, labels$group))

  # Validate the NK removal against CN-Male fine annotation (ground truth).
  fine <- ifelse(is.na(obj$Cell_type_identity[keep]), "<NA>",
                 as.character(obj$Cell_type_identity[keep]))
  is_cm <- labels$group == "CN-Male"
  fine_nk <- grepl("NK", fine)
  cd8_now <- labels$lineage == "CD8 T"   # clean CD8 after NK removal
  cat("\n=== Validation (CN-Male): clean 'CD8 T' vs fine annotation ===\n")
  cat(sprintf("  our NK removed from CD8 (CN-Male): %d | fine NK in coarse CD8 T: %d\n",
      sum(is_cm & labels$coarse == "CD8 T" & labels$lineage == "NK"),
      sum(is_cm & labels$coarse == "CD8 T" & fine_nk)))
  cm_cd8 <- is_cm & cd8_now & labels$coarse == "CD8 T"
  fine_cm_cd8 <- ifelse(grepl("NK", fine[cm_cd8]), "NK",
                 ifelse(grepl("^CD8", fine[cm_cd8]) & !grepl("NK-like", fine[cm_cd8]), "CD8 T", "other"))
  cat(sprintf("  clean CD8 (CN-Male) purity vs fine: %.0f%% CD8, %.0f%% residual NK, %.0f%% other\n",
      100 * mean(fine_cm_cd8 == "CD8 T"), 100 * mean(fine_cm_cd8 == "NK"),
      100 * mean(fine_cm_cd8 == "other")))

  list(labels = labels,
       validation_note = "NK removed uniformly via re-cluster nk_call; validated vs CN-Male fine")
}

classify_cd8_3state_clean <- function(seurat_filtered, lineage_annotation) {
  state_order <- c("Naive", "Effector", "Dysfunctional")
  markers <- cd8_3state_markers()
  clean_bc <- lineage_annotation$labels$barcode[
    lineage_annotation$labels$lineage == "CD8 T"]
  cat(sprintf("=== Stage 2: CD8 3-state on CLEAN CD8 (%d cells) ===\n", length(clean_bc)))

  cd8 <- subset(seurat_filtered, cells = clean_bc)
  cd8 <- Seurat::RunUMAP(cd8, reduction = "harmony", dims = 1:20, verbose = FALSE)
  cd8 <- Seurat::FindNeighbors(cd8, reduction = "harmony", dims = 1:20, verbose = FALSE)
  cd8 <- Seurat::FindClusters(cd8, resolution = 0.5, verbose = FALSE)

  avail <- lapply(markers, function(g) intersect(g, rownames(cd8)))
  set.seed(42)
  cd8 <- Seurat::AddModuleScore(cd8, features = avail, name = "state3_")
  score_mat <- as.matrix(cd8@meta.data[, paste0("state3_", 1:3)])
  colnames(score_mat) <- names(markers)
  cd8$assigned_state <- factor(names(markers)[apply(score_mat, 1, which.max)],
                               levels = state_order)
  cd8$new_Naive <- score_mat[, "Naive"]
  cd8$new_Effector <- score_mat[, "Effector"]
  cd8$new_Dysfunctional <- score_mat[, "Dysfunctional"]
  s_sorted <- apply(score_mat, 1, function(x) sort(x, decreasing = TRUE))
  cd8$classification_margin <- s_sorted[1, ] - s_sorted[2, ]

  props <- cd8@meta.data %>%
    dplyr::count(group, assigned_state) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(total = sum(n), prop = n / total,
                  pct = paste0(round(prop * 100, 1), "%")) %>%
    dplyr::ungroup()
  cat("\n=== CLEAN CD8 state proportions ===\n")
  print(props %>% dplyr::select(group, assigned_state, n, pct))

  list(cd8_obj = cd8, proportions = props, markers = markers,
       n_per_group = table(cd8$group),
       n_cells = ncol(cd8), source = "clean_lineage")
}

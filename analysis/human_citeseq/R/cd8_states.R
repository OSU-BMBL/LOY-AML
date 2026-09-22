# R/cd8_states.R
# CD8 T cell 3-state classification: Naive / Effector / Dysfunctional

#' CD8 3-state marker lists used by the classifier and validation targets
#'
#' @return named list of marker genes for Naive / Effector / Dysfunctional
cd8_3state_markers <- function() {
  list(
    Naive = c("IL7R", "CCR7", "SELL", "FOXO1", "KLF2", "KLF3",
              "LEF1", "TCF7", "ACTN1", "FOXP1"),
    Effector = c("GZMA", "GZMB", "GZMH", "GZMK", "GNLY", "PRF1",
                 "IFNG", "TNF", "NKG7", "FGFBP2", "CX3CR1", "KLRD1",
                 "KLRG1", "FCGR3A", "HOPX", "ZEB2", "TBX21", "PLEK",
                 "EOMES", "RUNX3"),
    Dysfunctional = c(
      "PDCD1", "LAYN", "HAVCR2", "LAG3", "CD244", "CTLA4",
      "LILRB1", "TIGIT", "TOX", "VSIR", "BTLA", "ENTPD1",
      "CD160", "LAIR1",
      "ANXA1", "LMNA", "GPR183", "ITGB1", "KLRB1",
      "S1PR1", "CD27", "BCL2"
    )
  )
}

#' Classify CD8 T cells into 3 states using AddModuleScore
#'
#' Steps:
#' 1. Subset CD8 T from full filtered object (keeps Harmony reduction)
#' 2. RunUMAP on Harmony dims 1:20
#' 3. FindNeighbors + FindClusters (res=0.5)
#' 4. AddModuleScore with Naive/Effector/Dysfunctional marker gene sets
#' 5. Winner-take-all classification, compute margins
#' 6. Compute proportions and chi-square test
#'
#' @param obj Seurat object (full filtered, must have Harmony reduction)
#' @return List with cd8 Seurat object, proportions, chi_sq, margin stats
classify_cd8_3state <- function(obj) {
  state_order <- c("Naive", "Effector", "Dysfunctional")

  # 3-state markers (Linghua-based, Exhaustion + Memory merged into Dysfunctional)
  markers <- cd8_3state_markers()

  # Step 1: Subset CD8 T (keeps Harmony embeddings from full object)
  cat("=== Subsetting CD8 T cells ===\n")
  cd8 <- subset(obj, subset = cell_type == "CD8 T")
  cat("CD8 T cells:", ncol(cd8), "\n")
  rm(obj); gc()

  # Step 2: Re-run UMAP on full-object Harmony dims
  cat("=== Running UMAP on Harmony dims 1:20 ===\n")
  cd8 <- Seurat::RunUMAP(cd8, reduction = "harmony", dims = 1:20, verbose = FALSE)

  # Step 3: Clustering
  cat("=== Clustering ===\n")
  cd8 <- Seurat::FindNeighbors(cd8, reduction = "harmony", dims = 1:20, verbose = FALSE)
  cd8 <- Seurat::FindClusters(cd8, resolution = 0.5, verbose = FALSE)

  # Step 4: AddModuleScore
  avail <- lapply(markers, function(g) intersect(g, rownames(cd8)))
  cat("Gene sets: Naive", length(avail$Naive), "/ Effector", length(avail$Effector),
      "/ Dysfunctional", length(avail$Dysfunctional), "\n")

  cat("=== Running AddModuleScore ===\n")
  set.seed(42)
  cd8 <- Seurat::AddModuleScore(cd8, features = avail, name = "state3_")
  score_cols <- paste0("state3_", 1:3)
  score_mat <- as.matrix(cd8@meta.data[, score_cols])
  colnames(score_mat) <- names(markers)

  # Step 5: Winner-take-all assignment
  cd8$assigned_state <- factor(
    names(markers)[apply(score_mat, 1, which.max)],
    levels = state_order
  )

  # Store named scores
  cd8$new_Naive <- score_mat[, "Naive"]
  cd8$new_Effector <- score_mat[, "Effector"]
  cd8$new_Dysfunctional <- score_mat[, "Dysfunctional"]

  # Classification margin
  s_max <- apply(score_mat, 1, max)
  s_2nd <- apply(score_mat, 1, function(x) sort(x, decreasing = TRUE)[2])
  cd8$classification_margin <- s_max - s_2nd

  # Step 6: Proportions and chi-square
  cat("\n=== State proportions ===\n")
  props <- cd8@meta.data %>%
    dplyr::count(group, assigned_state) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(total = sum(n), prop = n / total,
                  pct = paste0(round(prop * 100, 1), "%")) %>%
    dplyr::ungroup()
  print(props %>% dplyr::select(group, assigned_state, n, pct))

  # Chi-square
  tbl <- props %>%
    dplyr::select(group, assigned_state, n) %>%
    tidyr::pivot_wider(names_from = assigned_state, values_from = n, values_fill = 0)
  mat <- as.matrix(tbl[, -1])
  rownames(mat) <- tbl$group
  chi_sq <- chisq.test(mat)
  cat(sprintf("Chi-square: X2=%.2f, df=%d, p=%.2e\n",
              chi_sq$statistic, chi_sq$parameter, chi_sq$p.value))

  # Margin stats
  overall_margin <- median(cd8$classification_margin)
  low_conf <- sum(cd8$classification_margin < 0.1) / ncol(cd8)
  cat(sprintf("Overall margin median: %.3f, low-confidence: %.1f%%\n",
              overall_margin, 100 * low_conf))

  n_per_group <- table(cd8$group)

  list(
    cd8_obj = cd8,
    proportions = props,
    chi_sq = chi_sq,
    margin_median = overall_margin,
    low_confidence_pct = low_conf,
    markers = markers,
    n_per_group = n_per_group
  )
}

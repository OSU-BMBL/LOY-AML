# Blast cytokine expression and antigen-presentation module scores.

check_blast_cytokines <- function(degs_by_celltype, obj) {
  CYTOKINES <- c("IL12A", "IL12B", "IL18", "IL15", "IL15RA",
                 "IL6", "IL10", "CXCL9", "CXCL10", "IFNB1")
  ROLE <- c(
    IL12A = "Richer: miR → BACH2 degradation",
    IL12B = "Richer: miR → BACH2 degradation",
    IL18 = "Richer: miR → BACH2 degradation",
    IL15 = "Richer: miR → BACH2 degradation",
    IL15RA = "IL-15 trans-presentation receptor",
    IL6 = "Immunosuppressive control",
    IL10 = "Immunosuppressive control",
    CXCL9 = "IFNg-induced chemokine",
    CXCL10 = "IFNg-induced chemokine",
    IFNB1 = "Type I IFN"
  )

  blast_deseq <- degs_by_celltype$blast$deseq2
  blast_wilcox <- degs_by_celltype$blast$wilcox

  blast <- subset(obj, cell_type == "Leukemia")
  expr <- Seurat::GetAssayData(blast, layer = "data")

  results <- do.call(rbind, lapply(CYTOKINES, function(g) {
    deseq_row <- blast_deseq[blast_deseq$gene == g, ]
    wilcox_row <- blast_wilcox[blast_wilcox$gene == g, ]

    in_assay <- g %in% rownames(expr)
    det_all <- if (in_assay) mean(expr[g, ] > 0) * 100 else 0

    loy_mask <- blast$group == "LOY"
    det_loy <- if (in_assay) mean(expr[g, loy_mask] > 0) * 100 else 0
    det_cn <- if (in_assay) mean(expr[g, !loy_mask] > 0) * 100 else 0

    data.frame(
      gene = g,
      role = ROLE[g],
      in_assay = in_assay,
      det_all_pct = round(det_all, 1),
      det_loy_pct = round(det_loy, 1),
      det_cn_pct = round(det_cn, 1),
      deseq2_log2FC = if (nrow(deseq_row) > 0) round(deseq_row$log2FC, 3) else NA,
      deseq2_padj = if (nrow(deseq_row) > 0) deseq_row$padj else NA,
      deseq2_baseMean = if (nrow(deseq_row) > 0) round(deseq_row$baseMean, 1) else NA,
      floor_status = if (det_all < 5) "Below floor (<5%)" else "Expressed",
      stringsAsFactors = FALSE
    )
  }))

  rm(blast, expr); gc(verbose = FALSE)

  cat("=== Blast Inflammatory Cytokines (Richer 2016 Model) ===\n")
  print(results[, c("gene", "role", "det_loy_pct", "det_cn_pct",
                     "deseq2_log2FC", "deseq2_padj", "floor_status")])

  richer_genes <- c("IL12A", "IL12B", "IL18", "IL15")
  richer_rows <- results[results$gene %in% richer_genes, ]
  any_loy_specific_down <- any(!is.na(richer_rows$deseq2_padj) &
                                richer_rows$deseq2_padj < 0.05 &
                                richer_rows$deseq2_log2FC < 0)

  cat(sprintf("\nRicher miRNA model support: %s\n",
    if (any_loy_specific_down) "SUPPORTED — key cytokines DOWN in LOY"
    else "NOT SUPPORTED — IL-12/IL-18/IL-15 not LOY-specifically reduced"
  ))

  results
}

compute_blast_cue_score <- function(obj, cd8_states, blast_cytokines) {
  APM_GENES <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
                 "TAPBP", "PSMB8", "PSMB9")
  INFLAM_CANDIDATES <- c("IL12A", "IL12B", "IL15", "IL18",
                         "CXCL9", "CXCL10", "IFNB1")

  expressed <- blast_cytokines$gene[blast_cytokines$det_all_pct >= 5 &
                                      blast_cytokines$gene %in% INFLAM_CANDIDATES]
  cat(sprintf("Inflammatory cue genes passing 5%% floor: %s\n",
              paste(expressed, collapse = ", ")))

  blast <- subset(obj, cell_type == "Leukemia")
  blast_expr <- Seurat::GetAssayData(blast, layer = "data")

  apm_avail <- intersect(APM_GENES, rownames(blast_expr))
  inflam_avail <- intersect(expressed, rownames(blast_expr))

  set.seed(42)
  blast <- Seurat::AddModuleScore(blast, features = list(apm_avail), name = "APM_")
  if (length(inflam_avail) >= 2) {
    blast <- Seurat::AddModuleScore(blast, features = list(inflam_avail), name = "Inflam_")
  }

  combined_genes <- unique(c(apm_avail, inflam_avail))
  if (length(combined_genes) >= 3) {
    blast <- Seurat::AddModuleScore(blast, features = list(combined_genes), name = "BlastCue_")
  }

  blast_donors <- get_donor_id(blast)
  donor_group <- tapply(as.character(blast$group), blast_donors,
                        function(x) names(which.max(table(x))))

  scores <- data.frame(
    donor = unique(blast_donors),
    stringsAsFactors = FALSE
  )
  scores$group <- as.character(donor_group[scores$donor])
  scores$blast_apm_score <- tapply(blast$APM_1, blast_donors, mean)[scores$donor]
  if ("Inflam_1" %in% colnames(blast@meta.data)) {
    scores$blast_inflam_score <- tapply(blast$Inflam_1, blast_donors, mean)[scores$donor]
  } else {
    scores$blast_inflam_score <- NA_real_
  }
  if ("BlastCue_1" %in% colnames(blast@meta.data)) {
    scores$blast_cue_score <- tapply(blast$BlastCue_1, blast_donors, mean)[scores$donor]
  } else {
    scores$blast_cue_score <- NA_real_
  }
  scores$blast_ncount <- tapply(blast$nCount_RNA, blast_donors, mean)[scores$donor]

  rm(blast, blast_expr); gc(verbose = FALSE)

  cd8 <- cd8_states$cd8_obj
  cd8_donors <- get_donor_id(cd8)
  cd8_expr <- as.matrix(Seurat::GetAssayData(cd8, layer = "data"))

  scores$cd8_bach2 <- tapply(cd8_expr["BACH2", ], cd8_donors, mean)[scores$donor]
  CYTO <- intersect(c("GZMB", "GZMH", "GNLY", "NKG7", "PRF1"), rownames(cd8))
  set.seed(42)
  cd8 <- Seurat::AddModuleScore(cd8, features = list(CYTO), name = "Cyto_")
  scores$cd8_cytotoxic <- tapply(cd8$Cyto_1, cd8_donors, mean)[scores$donor]
  scores$cd8_ncount <- tapply(cd8$nCount_RNA, cd8_donors, mean)[scores$donor]

  cors <- list()
  test_pairs <- list(
    c("blast_cue_score", "cd8_bach2"),
    c("blast_cue_score", "cd8_cytotoxic"),
    c("blast_apm_score", "cd8_bach2"),
    c("blast_inflam_score", "cd8_bach2")
  )
  for (pair in test_pairs) {
    x <- scores[[pair[1]]]; y <- scores[[pair[2]]]
    valid <- !is.na(x) & !is.na(y) & length(x) >= 5
    if (sum(valid) >= 5) {
      ct <- suppressWarnings(cor.test(x[valid], y[valid], method = "spearman"))
      r_x <- rank(x[valid]); r_y <- rank(y[valid])
      r_bn <- rank(scores$blast_ncount[valid]); r_cn <- rank(scores$cd8_ncount[valid])
      res_x <- residuals(lm(r_x ~ r_bn + r_cn))
      res_y <- residuals(lm(r_y ~ r_bn + r_cn))
      pt <- cor.test(res_x, res_y, method = "pearson")

      cors[[length(cors) + 1]] <- data.frame(
        from = pair[1], to = pair[2],
        raw_rho = unname(ct$estimate), raw_p = ct$p.value,
        partial_r = unname(pt$estimate), partial_p = pt$p.value,
        n = sum(valid), stringsAsFactors = FALSE
      )
      cat(sprintf("%s vs %s: raw rho=%.3f (p=%.3f), partial r=%.3f (p=%.3f)\n",
                  pair[1], pair[2], ct$estimate, ct$p.value,
                  pt$estimate, pt$p.value))
    }
  }

  list(
    scores = scores,
    correlations = if (length(cors) > 0) do.call(rbind, cors) else NULL,
    apm_genes = apm_avail,
    inflam_genes = inflam_avail,
    combined_genes = combined_genes
  )
}

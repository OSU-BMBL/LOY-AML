# Differential expression: LOY versus male cytogenetically normal controls.
# Cell-level Wilcoxon tests and donor-level pseudobulk DESeq2 are returned
# separately for each cell type or CD8 state. Female controls are excluded
# from these contrasts.
#' Aggregate per-cell counts into a (donor × cell_type) pseudobulk matrix
#'
#' @param obj Seurat object
#' @param celltype Cell type label in obj$cell_type
#' @param min_cells minimum cells per donor to keep that donor
#' @return list(counts matrix genes x donors, colData data.frame donor, group)
make_pseudobulk <- function(obj, celltype, min_cells = 10) {
  sub <- subset(obj, subset = cell_type == celltype & group %in% c("LOY", "CN-Male"))
  donor <- get_donor_id(sub)
  sub$donor <- donor

  keep_donors <- names(which(table(donor) >= min_cells))
  sub <- subset(sub, subset = donor %in% keep_donors)

  counts <- Seurat::GetAssayData(sub, assay = "RNA", layer = "counts")
  donor_v <- factor(sub$donor)
  # indicator matrix cells x donors; pb = genes x donors via counts %*% indicator
  indic <- Matrix::sparse.model.matrix(~ 0 + donor_v)
  colnames(indic) <- levels(donor_v)
  pb <- as.matrix(counts %*% indic)

  donor_to_group <- tapply(as.character(sub$group), as.character(sub$donor),
                            function(x) x[1])
  coldata <- data.frame(
    donor = colnames(pb),
    group = factor(donor_to_group[colnames(pb)],
                   levels = c("CN-Male", "LOY")),
    n_cells = as.integer(table(donor_v)[colnames(pb)]),
    row.names = colnames(pb),
    stringsAsFactors = FALSE
  )
  list(counts = pb, coldata = coldata)
}

#' Run pseudobulk DESeq2 LOY vs CN-Male on one cell type
#'
#' @return tibble (gene, baseMean, log2FC, stat, pvalue, padj, method, cell_type,
#'         n_donors_LOY, n_donors_CN)
run_deseq2_pseudobulk <- function(obj, celltype, min_cells = 10,
                                   min_counts = 10, min_samples = 3) {
  cat(sprintf("  [DESeq2] %s pseudobulk ...\n", celltype))
  pb <- make_pseudobulk(obj, celltype, min_cells = min_cells)
  cd <- pb$coldata
  cnt <- round(pb$counts)
  storage.mode(cnt) <- "integer"

  n_loy <- sum(cd$group == "LOY")
  n_cn <- sum(cd$group == "CN-Male")
  if (n_loy < 2 || n_cn < 2) {
    warning("Too few donors for DESeq2 ", celltype, ": LOY=", n_loy, " CN=", n_cn)
    return(NULL)
  }
  keep <- rowSums(cnt >= min_counts) >= min_samples
  cnt <- cnt[keep, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt, colData = cd,
                                        design = ~ group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("group", "LOY", "CN-Male"),
                         independentFiltering = TRUE)

  out <- tibble::tibble(
    gene = rownames(res),
    baseMean = res$baseMean,
    log2FC = res$log2FoldChange,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj,
    method = "DESeq2_pseudobulk",
    cell_type = celltype,
    n_donors_LOY = n_loy,
    n_donors_CN = n_cn
  )
  out <- out[!is.na(out$log2FC), ]
  out <- dplyr::arrange(out, pvalue)
  cat(sprintf("    n_genes_tested=%d  n_padj<0.05=%d  |log2FC|>1_padj<0.05=%d\n",
              nrow(out), sum(out$padj < 0.05, na.rm = TRUE),
              sum(out$padj < 0.05 & abs(out$log2FC) > 1, na.rm = TRUE)))
  out
}

#' Run Wilcoxon (cell-level) LOY vs CN-Male on one cell type via presto
#'
#' @return tibble (gene, log2FC, auc, pvalue, padj, pct_LOY, pct_CN,
#'         method, cell_type, n_cells_LOY, n_cells_CN)
run_wilcox_cell <- function(obj, celltype) {
  cat(sprintf("  [Wilcoxon] %s cell-level ...\n", celltype))
  sub <- subset(obj, subset = cell_type == celltype & group %in% c("LOY", "CN-Male"))
  counts <- Seurat::GetAssayData(sub, assay = "RNA", layer = "data")  # log-normalized
  grp <- as.character(sub$group)

  res <- presto::wilcoxauc(counts, grp)
  res_loy <- res[res$group == "LOY", ]

  n_loy <- sum(grp == "LOY")
  n_cn <- sum(grp == "CN-Male")

  out <- tibble::tibble(
    gene = res_loy$feature,
    log2FC = res_loy$logFC,   # presto logFC is already on log scale (diff of log-norm means)
    auc = res_loy$auc,
    pvalue = res_loy$pval,
    padj = res_loy$padj,
    pct_LOY = res_loy$pct_in / 100,
    pct_CN = res_loy$pct_out / 100,
    method = "Wilcoxon_cell",
    cell_type = celltype,
    n_cells_LOY = n_loy,
    n_cells_CN = n_cn
  )
  out <- dplyr::arrange(out, pvalue)
  cat(sprintf("    n_genes=%d  n_padj<0.05=%d  |log2FC|>0.25_padj<0.05=%d\n",
              nrow(out), sum(out$padj < 0.05, na.rm = TRUE),
              sum(out$padj < 0.05 & abs(out$log2FC) > 0.25, na.rm = TRUE)))
  out
}

#' Top-level DEG pipeline: Wilcoxon + DESeq2 for blast, CD8 T, CD4 T
#'
#' @param obj Seurat object (full filtered)
#' @return list with elements blast/cd8/cd4, each containing $wilcox + $deseq2 tibbles
run_degs_by_celltype <- function(obj) {
  cat("=== Module I: DEG analysis LOY vs CN-Male ===\n")
  cat("Cell types: Leukemia (blasts), CD8 T, CD4 T\n")
  cat("Methods: Wilcoxon (cell-level, presto) + DESeq2 (pseudobulk per donor)\n\n")

  celltype_map <- c(blast = "Leukemia", cd8 = "CD8 T", cd4 = "CD4 T")
  out <- list()
  for (nm in names(celltype_map)) {
    ct <- celltype_map[[nm]]
    cat(sprintf("\n--- %s (%s) ---\n", nm, ct))
    out[[nm]] <- list(
      wilcox = run_wilcox_cell(obj, ct),
      deseq2 = run_deseq2_pseudobulk(obj, ct)
    )
  }
  out$meta <- list(
    contrast = "LOY vs CN-Male",
    cn_female_note = "CN-Female n=2 excluded from primary DEG (used only for specific gene checks)",
    generated_at = Sys.time()
  )
  cat("\n=== Module I complete ===\n")
  out
}

#' Flatten degs list to one long tibble for downstream tests/tables
#'
#' @param degs_list output of run_degs_by_celltype
#' @return tibble with all DEGs, unified schema
flatten_degs <- function(degs_list) {
  dfs <- list()
  cell_keys <- setdiff(names(degs_list), "meta")
  for (nm in cell_keys) {
    for (m in c("wilcox", "deseq2")) {
      d <- degs_list[[nm]][[m]]
      if (is.null(d) || nrow(d) == 0) next
      d$cell_group <- nm
      dfs[[paste(nm, m, sep = "_")]] <- d
    }
  }
  dplyr::bind_rows(dfs)
}

# ══════════════════════════════════════════════════════════════════════════════
# State-specific CD8 contrasts compare LOY and male control cells within
# each state, separating these comparisons from pooled CD8 composition.
# ---------------------------------------------------------------------------
#' Run DESeq2 pseudobulk + Wilcoxon on a pre-subset Seurat object
#'
#' Helper used by run_degs_by_cd8_state. Skips the cell_type filter (caller
#' has already subsetted).
#'
#' @param sub Seurat object already subset to the desired cell population
#' @param label string for logs (e.g. "CD8 Effector")
#' @param min_cells minimum cells per donor to retain donor
#' @return list(deseq2, wilcox, n_donors_LOY, n_donors_CN, n_cells_total)
run_deg_pair_on_subset <- function(sub, label, min_cells = 5) {
  cat(sprintf("  [%s] n_cells_total=%d\n", label, ncol(sub)))
  donor <- get_donor_id(sub)
  sub$donor <- donor
  keep_donors <- names(which(table(donor) >= min_cells))
  if (length(keep_donors) < 4) {
    warning(label, ": too few donors passing min_cells=", min_cells,
            " (kept ", length(keep_donors), ")")
    return(list(deseq2 = NULL, wilcox = NULL,
                n_donors_LOY = 0, n_donors_CN = 0,
                n_cells_total = ncol(sub)))
  }
  sub <- subset(sub, subset = donor %in% keep_donors)

  # Pseudobulk
  counts <- Seurat::GetAssayData(sub, assay = "RNA", layer = "counts")
  donor_v <- factor(sub$donor)
  indic <- Matrix::sparse.model.matrix(~ 0 + donor_v)
  colnames(indic) <- levels(donor_v)
  pb <- as.matrix(counts %*% indic)

  donor_to_group <- tapply(as.character(sub$group), as.character(sub$donor),
                            function(x) x[1])
  cd <- data.frame(
    donor = colnames(pb),
    group = factor(donor_to_group[colnames(pb)], levels = c("CN-Male", "LOY")),
    n_cells = as.integer(table(donor_v)[colnames(pb)]),
    row.names = colnames(pb),
    stringsAsFactors = FALSE
  )
  n_loy <- sum(cd$group == "LOY")
  n_cn <- sum(cd$group == "CN-Male")
  cat(sprintf("    donors post-filter: LOY=%d  CN-Male=%d\n", n_loy, n_cn))
  if (n_loy < 2 || n_cn < 2) {
    warning(label, ": DESeq2 needs >=2 donors per group (LOY=",
            n_loy, " CN=", n_cn, ")")
    return(list(deseq2 = NULL, wilcox = NULL,
                n_donors_LOY = n_loy, n_donors_CN = n_cn,
                n_cells_total = ncol(sub)))
  }

  cnt <- round(pb); storage.mode(cnt) <- "integer"
  keep <- rowSums(cnt >= 10) >= 3
  cnt <- cnt[keep, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt, colData = cd,
                                        design = ~ group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("group", "LOY", "CN-Male"),
                         independentFiltering = TRUE)
  deseq2_tbl <- tibble::tibble(
    gene = rownames(res), baseMean = res$baseMean,
    log2FC = res$log2FoldChange, stat = res$stat,
    pvalue = res$pvalue, padj = res$padj,
    method = "DESeq2_pseudobulk", subset_label = label,
    n_donors_LOY = n_loy, n_donors_CN = n_cn
  )
  deseq2_tbl <- deseq2_tbl[!is.na(deseq2_tbl$log2FC), ]
  deseq2_tbl <- dplyr::arrange(deseq2_tbl, pvalue)
  n_sig <- sum(deseq2_tbl$padj < 0.05, na.rm = TRUE)
  n_sig_fc <- sum(deseq2_tbl$padj < 0.05 & abs(deseq2_tbl$log2FC) > 1,
                  na.rm = TRUE)
  cat(sprintf("    DESeq2: n_tested=%d  padj<0.05=%d  |FC|>1 & padj<0.05=%d\n",
              nrow(deseq2_tbl), n_sig, n_sig_fc))

  # Wilcoxon on cell-level normalized data
  counts_log <- Seurat::GetAssayData(sub, assay = "RNA", layer = "data")
  grp <- as.character(sub$group)
  res_w <- presto::wilcoxauc(counts_log, grp)
  res_loy <- res_w[res_w$group == "LOY", ]
  wilcox_tbl <- tibble::tibble(
    gene = res_loy$feature,
    log2FC = res_loy$logFC,
    auc = res_loy$auc, pvalue = res_loy$pval, padj = res_loy$padj,
    pct_LOY = res_loy$pct_in / 100, pct_CN = res_loy$pct_out / 100,
    method = "Wilcoxon_cell", subset_label = label,
    n_cells_LOY = sum(grp == "LOY"),
    n_cells_CN = sum(grp == "CN-Male")
  )
  wilcox_tbl <- dplyr::arrange(wilcox_tbl, pvalue)

  list(deseq2 = deseq2_tbl, wilcox = wilcox_tbl,
       n_donors_LOY = n_loy, n_donors_CN = n_cn,
       n_cells_total = ncol(sub))
}

#' Subset-matched DEG: LOY vs CN-Male within each CD8 T 3-state
#'
#' @param cd8_states output of classify_cd8_3state (list containing $cd8_obj)
#' @param min_cells minimum cells per donor per state
#' @return list with elements Naive/Effector/Dysfunctional each containing
#'         $deseq2 and $wilcox tibbles; plus $meta and $key_genes summary.
run_degs_by_cd8_state <- function(cd8_states, min_cells = 5) {
  cat("=== Module I-2: Subset-matched CD8 T DEG (LOY vs CN-Male) ===\n")
  cat("Rationale: pooled CD8 T DEG conflates composition shift\n")
  cat("           (Naive 15.7%->44.8%, Effector 65.2%->29.4%)\n")
  cat("           with per-cell reprogramming. Running state-by-state.\n\n")

  cd8 <- cd8_states$cd8_obj
  cd8 <- subset(cd8, subset = group %in% c("LOY", "CN-Male"))

  states <- c("Naive", "Effector", "Dysfunctional")
  out <- list()
  for (st in states) {
    cat(sprintf("\n--- CD8 %s state ---\n", st))
    sub <- subset(cd8, subset = assigned_state == st)
    out[[st]] <- run_deg_pair_on_subset(sub,
                                        label = paste("CD8", st),
                                        min_cells = min_cells)
  }

  # Key-gene summary across all 3 states + legacy pooled CD8 (for comparison)
  key_genes <- c(
    # TGFb axis
    "TGFBR1", "TGFBR2", "SMAD3", "SMAD4", "SMAD7", "TGFB1",
    # Effector machinery
    "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "IFNG", "TNF", "NKG7",
    # Activation / TCR / costim
    "CD8A", "CD8B", "CD69", "CD28", "CTLA4", "PDCD1", "LAG3", "TIGIT", "TOX",
    # Naive / memory / TRM TFs
    "LEF1", "TCF7", "FOXO1", "IL7R", "CCR7", "SELL", "RUNX3",
    # Effector TFs
    "TBX21", "EOMES",
    # Inflammatory
    "NFKB1", "NFKBIA", "JUN", "FOS",
    # Other
    "PTGER4", "BCL2", "KLRG1"
  )

  summary_rows <- list()
  for (st in states) {
    d <- out[[st]]$deseq2
    if (is.null(d)) next
    k <- d[d$gene %in% key_genes, c("gene", "log2FC", "padj"), drop = FALSE]
    if (nrow(k) == 0) next
    k$cd8_state <- st
    summary_rows[[st]] <- k
  }
  summary_long <- dplyr::bind_rows(summary_rows)

  summary_wide_log2FC <- if (nrow(summary_long) > 0) {
    summary_long %>%
      dplyr::select(gene, cd8_state, log2FC) %>%
      tidyr::pivot_wider(names_from = cd8_state, values_from = log2FC,
                         names_prefix = "log2FC_")
  } else NULL

  summary_wide_padj <- if (nrow(summary_long) > 0) {
    summary_long %>%
      dplyr::select(gene, cd8_state, padj) %>%
      tidyr::pivot_wider(names_from = cd8_state, values_from = padj,
                         names_prefix = "padj_")
  } else NULL

  cat("\n=== Key-gene summary across states (log2FC) ===\n")
  if (!is.null(summary_wide_log2FC)) {
    print(as.data.frame(summary_wide_log2FC), row.names = FALSE)
  }

  out$meta <- list(
    contrast = "LOY vs CN-Male",
    method_primary = "DESeq2 pseudobulk (donor-level), subset-matched by CD8 state",
    rationale = paste("Pooled CD8 T differential expression can conflate",
                      "composition shift with within-state expression changes"),
    states = states,
    min_cells_per_donor = min_cells,
    generated_at = Sys.time()
  )
  out$key_genes <- list(
    long = summary_long,
    wide_log2FC = summary_wide_log2FC,
    wide_padj = summary_wide_padj
  )
  cat("\n=== Module I-2 complete ===\n")
  out
}

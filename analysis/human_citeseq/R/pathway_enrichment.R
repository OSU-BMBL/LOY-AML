# Hallmark and GO biological-process enrichment for LOY-versus-control
# differential expression. Rankings from DESeq2 and cell-level Wilcoxon
# results are analysed separately.
#' Fetch MSigDB collection as list(pathway_name -> character vector of genes)
#'
#' @param category "H" (Hallmark) or "C5" with subcategory "GO:BP"
#' @param subcategory Optional subcategory
#' @return named list of gene symbol vectors
msigdb_list <- function(category, subcategory = NULL) {
  df <- msigdbr::msigdbr(species = "Homo sapiens", collection = category,
                         subcollection = subcategory)
  split(df$gene_symbol, df$gs_name)
}

#' Build signed ranking vector for fgsea from DESeq2 results
#'
#' Uses wald stat directly (ties handled by jitter for stability).
#'
#' @param deg_tbl tibble from run_deseq2_pseudobulk
#' @return named numeric vector (names = gene symbol)
rank_by_deseq2 <- function(deg_tbl) {
  x <- deg_tbl[!is.na(deg_tbl$stat), c("gene", "stat")]
  v <- x$stat + runif(nrow(x), -1e-6, 1e-6)   # jitter for ties
  names(v) <- x$gene
  sort(v, decreasing = TRUE)
}

#' Build signed ranking vector from Wilcoxon presto output
#'
#' Signed by logFC, magnitude by -log10(pvalue) — standard ranking for fgsea.
rank_by_wilcox <- function(deg_tbl) {
  x <- deg_tbl[!is.na(deg_tbl$pvalue) & !is.na(deg_tbl$log2FC), ]
  v <- sign(x$log2FC) * (-log10(pmax(x$pvalue, 1e-300)))
  v <- v + runif(length(v), -1e-6, 1e-6)
  names(v) <- x$gene
  sort(v, decreasing = TRUE)
}

#' Run fgsea on one ranking against one gene set list
#'
#' @param ranks named numeric vector (decreasing)
#' @param pathways list of pathways
#' @param label string to tag the output
#' @return tibble with pathway, pval, padj, NES, size, leadingEdge, collection, method, cell_type
fgsea_one <- function(ranks, pathways, collection_label,
                      method_label, cell_type_label, minsize = 15, maxsize = 500) {
  res <- fgsea::fgsea(pathways = pathways, stats = ranks,
                       minSize = minsize, maxSize = maxsize,
                       eps = 0)
  res <- tibble::as_tibble(res)
  if (nrow(res) == 0) return(res)
  res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ","),
                              character(1))
  res$collection <- collection_label
  res$method <- method_label
  res$cell_type <- cell_type_label
  dplyr::arrange(res, padj)
}

#' Top-level pathway enrichment pipeline
#'
#' @param degs_list output of run_degs_by_celltype
#' @param collections character vec of MSigDB collection tags
#' @return list with $hallmark, $gobp each a tibble; $combined = all rows rbind
run_deg_pathway_enrichment <- function(degs_list,
                                        collections = c("hallmark", "gobp")) {
  cat("=== Module J: DEG pathway enrichment (fgsea) ===\n")
  cat("Loading MSigDB collections ...\n")

  pathway_sets <- list()
  if ("hallmark" %in% collections) {
    pathway_sets$hallmark <- msigdb_list("H")
    cat(sprintf("  Hallmark: %d sets\n", length(pathway_sets$hallmark)))
  }
  if ("gobp" %in% collections) {
    pathway_sets$gobp <- msigdb_list("C5", "GO:BP")
    cat(sprintf("  GO BP:    %d sets\n", length(pathway_sets$gobp)))
  }

  cell_keys <- setdiff(names(degs_list), "meta")
  cell_to_label <- c(blast = "Leukemia (blast)", cd8 = "CD8 T", cd4 = "CD4 T")

  rows <- list()
  for (nm in cell_keys) {
    for (mth in c("deseq2", "wilcox")) {
      dtbl <- degs_list[[nm]][[mth]]
      if (is.null(dtbl) || nrow(dtbl) == 0) next
      ranks <- if (mth == "deseq2") rank_by_deseq2(dtbl) else rank_by_wilcox(dtbl)
      for (coll in names(pathway_sets)) {
        key <- paste(nm, mth, coll, sep = "_")
        cat(sprintf("\n[%s] ranks=%d  pathways=%d\n",
                    key, length(ranks), length(pathway_sets[[coll]])))
        rows[[key]] <- fgsea_one(ranks, pathway_sets[[coll]],
                                  collection_label = coll,
                                  method_label = mth,
                                  cell_type_label = cell_to_label[[nm]])
      }
    }
  }
  combined <- dplyr::bind_rows(rows)
  cat(sprintf("\nTotal pathway rows: %d  padj<0.05: %d\n",
              nrow(combined), sum(combined$padj < 0.05, na.rm = TRUE)))

  out <- list(
    combined = combined,
    by_cell = split(combined, combined$cell_type),
    meta = list(
      generated_at = Sys.time(),
      contrast = "LOY vs CN-Male",
      collections = collections
    )
  )
  cat("=== Module J complete ===\n")
  out
}

# ══════════════════════════════════════════════════════════════════════════════
# Enrichment of state-specific CD8 differential expression.
#
# Companion to Module I-2 (per-state DEG). Runs fgsea on each (state, method,
# collection) DEG ranking so pathway calls are composition-matched. Mirrors
# Module A-2 (per-state PROGENy) on the DEG-ranked-list side.
# ══════════════════════════════════════════════════════════════════════════════

#' Pathway enrichment on per-CD8-state DEG rankings
#'
#' @param degs_by_cd8_state output of run_degs_by_cd8_state (list with
#'   Naive/Effector/Dysfunctional elements each containing $deseq2 + $wilcox)
#' @param collections character vec of MSigDB collection tags
#' @return list with $combined (all rows), $by_state (split by state), $meta
run_deg_pathway_enrichment_by_cd8_state <- function(
  degs_by_cd8_state,
  collections = c("hallmark", "gobp")
) {
  cat("=== Module J-2: Per-state DEG pathway enrichment (fgsea) ===\n")
  cat("Input: run_degs_by_cd8_state output (Naive / Effector / Dysfunctional)\n")

  cat("\nLoading MSigDB collections ...\n")
  pathway_sets <- list()
  if ("hallmark" %in% collections) {
    pathway_sets$hallmark <- msigdb_list("H")
    cat(sprintf("  Hallmark: %d sets\n", length(pathway_sets$hallmark)))
  }
  if ("gobp" %in% collections) {
    pathway_sets$gobp <- msigdb_list("C5", "GO:BP")
    cat(sprintf("  GO BP:    %d sets\n", length(pathway_sets$gobp)))
  }

  states <- c("Naive", "Effector", "Dysfunctional")
  rows <- list()
  for (st in states) {
    st_entry <- degs_by_cd8_state[[st]]
    if (is.null(st_entry)) {
      cat(sprintf("\n[%s] missing — skipping\n", st))
      next
    }
    for (mth in c("deseq2", "wilcox")) {
      dtbl <- st_entry[[mth]]
      if (is.null(dtbl) || nrow(dtbl) == 0) next
      ranks <- if (mth == "deseq2") rank_by_deseq2(dtbl) else rank_by_wilcox(dtbl)
      for (coll in names(pathway_sets)) {
        key <- paste(st, mth, coll, sep = "_")
        cat(sprintf("\n[%s] ranks=%d  pathways=%d\n",
                    key, length(ranks), length(pathway_sets[[coll]])))
        res <- fgsea_one(ranks, pathway_sets[[coll]],
                          collection_label = coll,
                          method_label = mth,
                          cell_type_label = paste0("CD8 ", st))
        if (is.null(res) || nrow(res) == 0) next
        res$state <- st
        rows[[key]] <- res
      }
    }
  }
  combined <- dplyr::bind_rows(rows)
  cat(sprintf("\nTotal per-state pathway rows: %d  padj<0.05: %d\n",
              nrow(combined), sum(combined$padj < 0.05, na.rm = TRUE)))

  out <- list(
    combined = combined,
    by_state = if (nrow(combined) > 0) split(combined, combined$state) else list(),
    meta = list(
      generated_at = Sys.time(),
      contrast = "LOY vs CN-Male (within state)",
      states = states,
      collections = collections,
      rationale = paste(
        "Pooled fgsea on CD8 T DEG can conflate composition shift",
        "with per-cell reprogramming. Per-state fgsea on subset-matched DEG",
        "rankings provides state-specific pathway comparisons."
      )
    )
  )
  cat("=== Module J-2 complete ===\n")
  out
}

# ══════════════════════════════════════════════════════════════════════════════
# GO biological-process over-representation among up- and down-regulated
# blast genes. The caller supplies the differential-expression table.
# The figure workflow uses the cell-level Wilcoxon result; term-display
# subsets are defined explicitly in R/plot_helpers.R.
# ---------------------------------------------------------------------------
#' GO-BP ORA on the blast LOY-vs-CN-Male DEG
#'
#' Uses enrichGO(ont="BP") with BH multiple-testing adjustment
#' on significant up-in-LOY and down-in-LOY genes separately, with the tested
#' genes as the universe. Returned terms depend on the supplied P-value
#' and q-value cutoffs; the figure workflow includes an unfiltered branch.
#'
#' @param blast_deg data.frame with columns gene, log2FC, padj (the
#'   cell-level Wilcoxon or donor-level DESeq2 result).
#' @param padj_cut adjusted-p cutoff for calling a gene DE (default 0.05).
#' @param lfc_cut absolute log2FC cutoff for the gene list (default 0).
#' @param pvalueCutoff,qvalueCutoff enrichGO significance cutoffs.
#' @return list(up_in_loy, up_in_cn, meta). Each table is a tibble of the full
#'   enrichGO result (Description, p.adjust, Count, GeneRatio, geneID, ...).
run_blast_ora <- function(blast_deg,
                          padj_cut = 0.05, lfc_cut = 0,
                          pvalueCutoff = 0.05, qvalueCutoff = 0.10,
                          simplify_cutoff = NULL,
                          minGSSize = 10, maxGSSize = 500) {
  stopifnot(all(c("gene", "log2FC", "padj") %in% colnames(blast_deg)))
  cat("=== Module J-3: GO-BP ORA on blast DEG (Fig 1d) ===\n")

  universe <- unique(as.character(blast_deg$gene))
  ok  <- !is.na(blast_deg$padj) & !is.na(blast_deg$log2FC) & blast_deg$padj < padj_cut
  up   <- unique(as.character(blast_deg$gene[ok & blast_deg$log2FC >  lfc_cut]))
  down <- unique(as.character(blast_deg$gene[ok & blast_deg$log2FC < -lfc_cut]))
  cat(sprintf("  universe=%d  up-in-LOY=%d  down-in-LOY(=up-in-CN)=%d  (padj<%.3g, |log2FC|>%.3g)\n",
              length(universe), length(up), length(down), padj_cut, lfc_cut))

  go <- function(g, lab) {
    if (length(g) < 10) {
      cat(sprintf("  [%s] <10 genes — skipping\n", lab)); return(NULL)
    }
    r <- clusterProfiler::enrichGO(
      gene = g, universe = universe,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL",
      ont = "BP", pAdjustMethod = "BH",
      pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize, maxGSSize = maxGSSize
    )
    if (is.null(r) || nrow(as.data.frame(r)) == 0) {
      cat(sprintf("  [%s] enrichGO returned 0 terms\n", lab)); return(NULL)
    }
    n_raw <- nrow(as.data.frame(r))
    if (!is.null(simplify_cutoff)) {
      r <- tryCatch(
        clusterProfiler::simplify(r, cutoff = simplify_cutoff,
                                  by = "p.adjust", select_fun = min),
        error = function(e) {
          cat(sprintf("  [%s] simplify failed (%s) — keeping raw\n",
                      lab, conditionMessage(e))); r })
    }
    out <- tibble::as_tibble(as.data.frame(r))
    cat(sprintf("  [%s] %d significant GO:BP terms%s\n", lab, nrow(out),
                if (!is.null(simplify_cutoff))
                  sprintf(" (simplified from %d, cutoff %.2f)", n_raw, simplify_cutoff)
                else ""))
    out
  }

  res <- list(
    up_in_loy = go(up,   "Up in LOY"),
    up_in_cn  = go(down, "Up in CN-Male / down in LOY"),
    meta = list(
      generated_at = Sys.time(),
      contrast = "Leukemia blast: LOY vs CN-Male",
      n_up = length(up), n_down = length(down), universe = length(universe),
      padj_cut = padj_cut, lfc_cut = lfc_cut, simplify_cutoff = simplify_cutoff,
      enrichGO = list(ont = "BP", pAdjustMethod = "BH",
                      pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
    )
  )
  cat("=== Module J-3 complete ===\n")
  res
}

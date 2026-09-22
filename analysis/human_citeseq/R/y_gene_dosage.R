# R/y_gene_dosage.R
# Module: Y/X gene pair total dosage analysis
# Compares total paired expression across LOY, CN-Male, CN-Female
#
# Biological rationale:
#   CN-Male: 1X + 1Y  -> total dosage = X_gene + Y_gene
#   CN-Female: 2X (XCI-escape genes expressed from both X) -> total dosage = X_gene (represents 2X)
#   LOY: 1X, Y lost  -> total dosage = X_gene (represents 1X only)
#   Hypothesis: CN-Male(X+Y) ≈ CN-Female(2X) > LOY(1X)

#' Run gene dosage analysis for X-Y homolog pairs
#'
#' Computes total paired expression per cell across three groups (LOY, CN-Male,
#' CN-Female) and tests whether dosage compensation holds:
#' CN-Male ≈ CN-Female > LOY.
#'
#' @param obj Seurat object (full filtered, must contain $group and $cell_type)
#' @return List with dosage_df, stat_results, all_stats, individual_gene_stats
run_gene_dosage_analysis <- function(obj) {
  # All X-Y gene pairs (Y-expressed genes + X homologs)
  GENE_PAIRS <- list(
    c("UTY", "KDM6A"),       # H3K27 demethylase
    c("KDM5D", "KDM5C"),     # H3K4me3 demethylase (MHC-I mechanism)
    c("ZFY", "ZFX"),         # Zinc finger TF
    c("DDX3Y", "DDX3X"),     # RNA helicase
    c("EIF1AY", "EIF1AX"),   # Translation initiation
    c("USP9Y", "USP9X"),     # Deubiquitinase
    c("RPS4Y1", "RPS4X"),    # Ribosomal protein
    c("TMSB4Y", "TMSB4X"),   # Thymosin beta-4
    c("PRKY", "PRKX")        # Protein kinase
  )
  CELL_TYPES <- c("Leukemia", "CD8 T")
  GROUPS     <- c("LOY", "CN-Male", "CN-Female")

  # ── 1. Extract expression matrix ─────────────────────────────────────────
  all_genes       <- unique(unlist(GENE_PAIRS))
  available_genes <- all_genes[all_genes %in% rownames(obj)]
  missing_genes   <- setdiff(all_genes, available_genes)
  if (length(missing_genes) > 0) {
    warning("Missing genes: ", paste(missing_genes, collapse = ", "))
  }

  expr_mat <- Seurat::GetAssayData(obj, layer = "data")[available_genes, , drop = FALSE]
  dosage_df <- as.data.frame(t(as.matrix(expr_mat)))
  dosage_df$cell_barcode <- rownames(dosage_df)
  dosage_df$group        <- obj$group
  dosage_df$cell_type    <- obj$cell_type

  # Keep only the three groups of interest
  dosage_df <- dosage_df[dosage_df$group %in% GROUPS, ]

  # ── 2. Compute total_paired_expr per gene pair ────────────────────────────
  #   CN-Male:   Y_gene + X_gene
  #   CN-Female: X_gene  (2X active copies; single-cell value represents 2X output)
  #   LOY:       X_gene  (1X only)
  for (pair in GENE_PAIRS) {
    y_gene  <- pair[1]
    x_gene  <- pair[2]
    col_nm  <- paste0("dosage_", y_gene, "_", x_gene)

    if (!(y_gene %in% available_genes) || !(x_gene %in% available_genes)) next

    dosage_df[[col_nm]] <- ifelse(
      dosage_df$group == "CN-Male",
      dosage_df[[y_gene]] + dosage_df[[x_gene]],
      dosage_df[[x_gene]]   # CN-Female and LOY: X only
    )
  }

  # ── 3. Statistical testing ────────────────────────────────────────────────
  stat_results        <- list()
  individual_gene_stats <- list()

  for (ct in CELL_TYPES) {
    ct_data <- dosage_df[dosage_df$cell_type == ct, ]
    if (nrow(ct_data) == 0) next

    present_groups <- intersect(GROUPS, ct_data$group)
    if (length(present_groups) < 2) next

    results_rows <- list()

    for (pair in GENE_PAIRS) {
      y_gene  <- pair[1]
      x_gene  <- pair[2]
      col_nm  <- paste0("dosage_", y_gene, "_", x_gene)
      pair_lbl <- paste0(y_gene, "/", x_gene)

      if (!(col_nm %in% colnames(ct_data))) next

      # Subsets per group
      loy_vals <- ct_data[[col_nm]][ct_data$group == "LOY"]
      cnm_vals <- ct_data[[col_nm]][ct_data$group == "CN-Male"]
      cnf_vals <- ct_data[[col_nm]][ct_data$group == "CN-Female"]

      # a. Kruskal-Wallis 3-group test (only if all three groups present)
      kw_p <- NA_real_
      if (all(c("LOY", "CN-Male", "CN-Female") %in% ct_data$group)) {
        kw_data  <- ct_data[ct_data$group %in% c("LOY", "CN-Male", "CN-Female"),
                             c(col_nm, "group")]
        kw_test  <- kruskal.test(as.formula(paste0("`", col_nm, "` ~ group")),
                                 data = kw_data)
        kw_p     <- kw_test$p.value
      }

      # b. Pairwise Wilcoxon tests; collect p-values for BH correction
      comparisons <- list(
        list(g1 = "CN-Male",   g2 = "CN-Female", v1 = cnm_vals, v2 = cnf_vals),
        list(g1 = "CN-Male",   g2 = "LOY",       v1 = cnm_vals, v2 = loy_vals),
        list(g1 = "CN-Female", g2 = "LOY",       v1 = cnf_vals, v2 = loy_vals)
      )

      pair_rows <- list()
      for (comp in comparisons) {
        if (length(comp$v1) < 2 || length(comp$v2) < 2) next
        wt <- wilcox.test(comp$v1, comp$v2, alternative = "two.sided")

        med1 <- median(comp$v1)
        med2 <- median(comp$v2)
        fc   <- if (med2 != 0) med1 / med2 else NA_real_

        pair_rows[[length(pair_rows) + 1]] <- data.frame(
          cell_type       = ct,
          gene_pair       = pair_lbl,
          comparison      = paste0(comp$g1, " vs ", comp$g2),
          median_group1   = med1,
          median_group2   = med2,
          fold_change     = fc,
          W_statistic     = wt$statistic,
          p_value         = wt$p.value,
          kw_p_value      = kw_p,
          stringsAsFactors = FALSE
        )
      }

      if (length(pair_rows) > 0) {
        pr_df  <- dplyr::bind_rows(pair_rows)
        pr_df$p_adj <- p.adjust(pr_df$p_value, method = "BH")
        results_rows <- c(results_rows, list(pr_df))
      }
    }

    if (length(results_rows) > 0) {
      stat_results[[ct]] <- dplyr::bind_rows(results_rows)
    }

    # ── 4. Individual gene expression stats per group ─────────────────────
    ind_rows <- list()
    for (gene in available_genes) {
      for (grp in present_groups) {
        vals <- ct_data[[gene]][ct_data$group == grp]
        ind_rows[[length(ind_rows) + 1]] <- data.frame(
          cell_type = ct,
          gene      = gene,
          group     = grp,
          median    = median(vals),
          mean      = mean(vals),
          pct_expr  = mean(vals > 0) * 100,
          n_cells   = length(vals),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(ind_rows) > 0) {
      individual_gene_stats[[ct]] <- dplyr::bind_rows(ind_rows)
    }
  }

  all_stats <- dplyr::bind_rows(stat_results)

  # ── 5. group_means matrix (Leukemia blasts, mean expression per group) ───
  blast_df <- dosage_df[dosage_df$cell_type == "Leukemia", ]
  gm_rows <- list()
  for (grp in GROUPS) {
    grp_data <- blast_df[blast_df$group == grp, ]
    if (nrow(grp_data) == 0) next
    means <- vapply(available_genes, function(g) mean(grp_data[[g]]), numeric(1))
    gm_rows[[grp]] <- means
  }
  group_means <- do.call(rbind, gm_rows)

  # ── 6. pair_dosage_comparison (CN-Male Y+X vs CN-Female 2X, Leukemia) ───
  pdc_rows <- list()
  for (pair in GENE_PAIRS) {
    y_gene <- pair[1]
    x_gene <- pair[2]
    if (!(y_gene %in% available_genes) || !(x_gene %in% available_genes)) next
    cnm_blast <- blast_df[blast_df$group == "CN-Male", ]
    cnf_blast <- blast_df[blast_df$group == "CN-Female", ]
    cn_male_total  <- mean(cnm_blast[[y_gene]]) + mean(cnm_blast[[x_gene]])
    cn_female_total <- mean(cnf_blast[[x_gene]])  # 2X represented by single-cell X value
    pdc_rows[[length(pdc_rows) + 1]] <- data.frame(
      pair = paste0(y_gene, "/", x_gene),
      cn_male_total = cn_male_total,
      cn_female_total = cn_female_total,
      stringsAsFactors = FALSE
    )
  }
  pair_dosage_comparison <- dplyr::bind_rows(pdc_rows)

  # ── 7. Print summary ──────────────────────────────────────────────────────
  cat("Gene dosage analysis complete.\n")
  cat(sprintf("  Groups present: %s\n",
              paste(sort(unique(dosage_df$group)), collapse = ", ")))
  cat(sprintf("  Cell types analysed: %s\n", paste(names(stat_results), collapse = ", ")))

  for (ct in names(stat_results)) {
    res <- stat_results[[ct]]
    cat(sprintf("\n  [%s]\n", ct))
    for (pair in unique(res$gene_pair)) {
      sub <- res[res$gene_pair == pair, ]
      cat(sprintf("    %s:\n", pair))
      for (i in seq_len(nrow(sub))) {
        sig <- if (!is.na(sub$p_adj[i]) && sub$p_adj[i] < 0.05) "*" else " "
        cat(sprintf("      %s%s  median1=%.3f  median2=%.3f  FC=%.2f  p_adj=%.3g\n",
                    sig, sub$comparison[i],
                    sub$median_group1[i], sub$median_group2[i],
                    ifelse(is.na(sub$fold_change[i]), NA, sub$fold_change[i]),
                    sub$p_adj[i]))
      }
    }
  }

  list(
    dosage_df             = dosage_df,
    stat_results          = stat_results,
    all_stats             = all_stats,
    individual_gene_stats = dplyr::bind_rows(individual_gene_stats),
    group_means           = group_means,
    pair_dosage_comparison = pair_dosage_comparison
  )
}

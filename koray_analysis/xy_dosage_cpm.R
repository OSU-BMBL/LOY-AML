# Purpose: compare X-Y homolog dosage from processed expression data.
# Input: sample-by-gene CPM and sample metadata; see xy_dosage_cpm() arguments.
# Output: homolog summaries and the existing dosage plots.
# Author: Koray Dogan Kaya.
# Sourcing defines the helper only; direct Rscript execution runs the original
# patient-pseudobulk example with three optional file/output path arguments.

xy_dosage_cpm <- function(
    cpm,
    gender,
    pairs,
    group = NULL,
    condition = NULL,
    reference_group = NULL,
    female_label = "Female",
    facet_ncol = 4,
    homolog_colors = c(X = "lightpink", Y = "#301934"),
    y_breaks = c(0, 0.5, 1),
    pair_label_size = 14,
    group_label_size = 12,
    group_label_angle = 0,
    condition_label_size = 12,
    out_prefix = NULL,
    width = 14,
    height = NULL,
    dpi = 300
) {
  
  cpm <- as.matrix(cpm)
  storage.mode(cpm) <- "double"
  
  if (is.null(rownames(cpm)) || is.null(colnames(cpm)))
    stop("cpm must have sample row names and gene column names.")
  
  if (any(!is.finite(cpm)) || any(cpm < 0))
    stop("cpm must contain finite, non-negative values.")
  
  n <- nrow(cpm)
  
  align_vector <- function(x, default, label) {
    if (is.null(x))
      return(rep(default, n))
    
    if (!is.null(names(x))) {
      x <- unname(x[rownames(cpm)])
      if (anyNA(x))
        stop(label, " is missing values for one or more matrix samples.")
    } else if (length(x) != n) {
      stop(label, " must have one value per matrix row.")
    }
    
    as.character(x)
  }
  
  gender <- align_vector(gender, NA_character_, "gender")
  group <- align_vector(group, "All", "group")
  condition <- align_vector(condition, "All", "condition")
  
  gender <- dplyr::case_when(
    tolower(trimws(gender)) %in% c("m", "male") ~ "Male",
    tolower(trimws(gender)) %in% c("f", "female") ~ "Female",
    TRUE ~ NA_character_
  )
  
  if (anyNA(gender))
    stop("gender must contain only Male/Female or M/F.")
  
  # Accept either c(Y_gene = "X_gene") or a geneY/geneX data frame.
  if (is.atomic(pairs) && !is.null(names(pairs))) {
    pair_df <- tibble::tibble(
      geneY = names(pairs),
      geneX = unname(as.character(pairs))
    )
  } else {
    pair_df <- as.data.frame(pairs)
    
    if (!all(c("geneY", "geneX") %in% names(pair_df)))
      stop("pairs must be a named Y=X vector or contain geneY and geneX columns.")
    
    pair_df <- pair_df |>
      dplyr::transmute(
        geneY = as.character(geneY),
        geneX = as.character(geneX)
      )
  }
  
  pair_df <- pair_df |>
    dplyr::distinct(geneY, geneX) |>
    dplyr::mutate(Pair = paste0(geneY, "/", geneX))
  
  genes <- unique(c(pair_df$geneY, pair_df$geneX))
  missing_genes <- setdiff(genes, colnames(cpm))
  
  if (length(missing_genes))
    stop(
      "Genes absent from cpm: ",
      paste(missing_genes, collapse = ", ")
    )
  
  male_groups <- unique(group[gender == "Male"])
  
  if (!length(male_groups))
    stop("At least one male sample is required for normalization.")
  
  if (is.null(reference_group))
    reference_group <- male_groups[1]
  
  if (!reference_group %in% male_groups)
    stop("reference_group is absent from the male samples.")
  
  if (female_label %in% male_groups)
    stop("female_label duplicates an existing male group name.")
  
  plot_group <- ifelse(
    gender == "Female",
    female_label,
    group
  )
  
  group_levels <- c(
    male_groups,
    if (any(gender == "Female")) female_label
  )
  
  condition_levels <- unique(condition)
  
  sample_meta <- tibble::tibble(
    sample = rownames(cpm),
    gender,
    group,
    condition,
    plot_group
  )
  
  pair_map <- dplyr::bind_rows(
    pair_df |>
      dplyr::transmute(
        gene = geneX,
        Pair,
        Homolog = "X"
      ),
    pair_df |>
      dplyr::transmute(
        gene = geneY,
        Pair,
        Homolog = "Y"
      )
  )
  
  sample_values <- as.data.frame(
    cpm[, genes, drop = FALSE],
    check.names = FALSE
  ) |>
    tibble::rownames_to_column("sample") |>
    dplyr::left_join(sample_meta, by = "sample") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(genes),
      names_to = "gene",
      values_to = "CPM"
    ) |>
    dplyr::left_join(pair_map, by = "gene")
  
  summary_df <- sample_values |>
    dplyr::group_by(
      condition,
      plot_group,
      Pair,
      Homolog
    ) |>
    dplyr::summarise(
      n_samples = sum(is.finite(CPM)),
      mean_CPM = mean(CPM, na.rm = TRUE),
      sd_CPM = stats::sd(CPM, na.rm = TRUE),
      se_CPM = sd_CPM / sqrt(n_samples),
      .groups = "drop"
    )
  
  # Pair-specific normalization denominator within each condition panel.
  reference <- summary_df |>
    dplyr::filter(plot_group == reference_group) |>
    dplyr::select(
      condition,
      Pair,
      Homolog,
      mean_CPM
    ) |>
    tidyr::pivot_wider(
      names_from = Homolog,
      values_from = mean_CPM,
      names_prefix = "ref_"
    ) |>
    dplyr::mutate(
      reference_total = ref_X + ref_Y
    )
  
  expected_references <-
    length(condition_levels) * nrow(pair_df)
  
  if (
    nrow(reference) != expected_references ||
    any(!is.finite(reference$reference_total)) ||
    any(reference$reference_total <= 0)
  ) {
    stop(
      "Every condition must contain male samples from reference_group = ",
      reference_group,
      " with a positive X+Y total."
    )
  }
  
  plot_data <- summary_df |>
    dplyr::left_join(
      reference |>
        dplyr::select(
          condition,
          Pair,
          reference_total
        ),
      by = c("condition", "Pair")
    ) |>
    dplyr::mutate(
      Relative_expression = mean_CPM / reference_total,
      condition = factor(condition, levels = condition_levels),
      plot_group = factor(plot_group, levels = group_levels),
      Pair = factor(Pair, levels = pair_df$Pair),
      Homolog = factor(Homolog, levels = c("Y", "X"))
    )
  
  metrics <- summary_df |>
    dplyr::select(
      condition,
      plot_group,
      Pair,
      Homolog,
      mean_CPM
    ) |>
    tidyr::pivot_wider(
      names_from = Homolog,
      values_from = mean_CPM,
      values_fill = 0
    ) |>
    dplyr::left_join(
      reference,
      by = c("condition", "Pair")
    ) |>
    dplyr::mutate(
      Total = X + Y,
      Y_fraction = dplyr::if_else(
        Total > 0,
        Y / Total,
        NA_real_
      ),
      Relative_X = X / reference_total,
      Relative_Y = Y / reference_total,
      Retained_dosage = Total / reference_total,
      X_log2FC = log2((X + 1e-6) / (ref_X + 1e-6)),
      Pair_log2FC = log2(
        (Total + 1e-6) /
          (reference_total + 1e-6)
      ),
      Compensation_efficiency = dplyr::case_when(
        plot_group == reference_group ~ NA_real_,
        plot_group == female_label ~ NA_real_,
        ref_Y <= Y ~ NA_real_,
        TRUE ~ (X - ref_X) / (ref_Y - Y)
      )
    ) |>
    dplyr::arrange(
      condition,
      factor(plot_group, levels = group_levels),
      factor(Pair, levels = pair_df$Pair)
    )
  # -----------------------------------------------------------------------
  # Plot: all homolog pairs in one horizontal row per condition
  # -----------------------------------------------------------------------
  
  group_axis_labels <- function(x) {
    as.expression(lapply(x, function(z) {
      if (z %in% c("M-CN", "CN-Male", "CN_Male")) {
        quote(M^CN)
      } else if (z %in% c("LOY", "M-LOY", "M_LOY")) {
        quote(M^LOY)
      } else {
        bquote(plain(.(z)))
      }
    }))
  }
  
  y_axis_label <- if (reference_group %in%
                      c("M-CN", "CN-Male", "CN_Male")) {
    bquote("Expression relative to " * M^CN * " X+Y total")
  } else if (reference_group %in%
             c("LOY", "M-LOY", "M_LOY")) {
    bquote("Expression relative to " * M^LOY * " X+Y total")
  } else {
    paste0(
      "Expression relative to ",
      reference_group,
      " X+Y total"
    )
  }
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = plot_group,
      y = Relative_expression,
      fill = Homolog
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::scale_fill_manual(
      values = homolog_colors,
      breaks = c("Y", "X"),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      labels = group_axis_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = \(x) as.character(x),
      minor_breaks = NULL,
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_axis_label,
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      
      # X/Y pair names above each facet
      strip.text.x = ggplot2::element_text(
        face = "bold.italic",
        size = pair_label_size
      ),
      
      # Condition labels when condition has multiple levels
      strip.text.y = ggplot2::element_text(
        face = "bold",
        size = condition_label_size
      ),
      
      strip.background = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(0.7, "lines"),
      
      # M^CN, M^LOY and Female labels below the bars
      axis.text.x = ggplot2::element_text(
        size = group_label_size,
        angle = group_label_angle,
        hjust = if (group_label_angle == 0) 0.5 else 1,
        vjust = if (group_label_angle == 0) 0.5 else 1
      )
    )
  
  if (length(condition_levels) > 1) {
    p <- p +
      ggplot2::facet_grid(
        rows = ggplot2::vars(condition),
        cols = ggplot2::vars(Pair),
        scales = "free_x",
        space = "free_x"
      )
  } else {
    p <- p +
      ggplot2::facet_grid(
        cols = ggplot2::vars(Pair),
        scales = "free_x",
        space = "free_x"
      )
  }
  
  
  # -----------------------------------------------------------------------
  # Save wide plot
  # -----------------------------------------------------------------------
  
  if (!is.null(out_prefix)) {
    dir.create(
      dirname(out_prefix),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    plot_width <- max(
      width,
      2 * nrow(pair_df)
    )
    
    if (is.null(height)) {
      height <- 3.2 * length(condition_levels) + 2
    }
    
    ggplot2::ggsave(
      paste0(out_prefix, ".pdf"),
      plot = p,
      width = plot_width,
      height = height,
      units = "in",
      limitsize = FALSE,
      bg = "white"
    )
    
    ggplot2::ggsave(
      paste0(out_prefix, ".png"),
      plot = p,
      width = plot_width,
      height = height,
      units = "in",
      dpi = dpi,
      limitsize = FALSE,
      bg = "white"
    )
    list(
      sample_values = sample_values,
      summary = summary_df,
      metrics = metrics,
      plot_data = plot_data,
      plot = p,
      reference_group = reference_group
    )
  }

  }
# Standalone example. HSC.Rmd sources the function without running this block.
if (sys.nframe() == 0L) {
  paths <- commandArgs(trailingOnly = TRUE)
  if (length(paths) == 0L) {
    paths <- c(
      "data/patient_xy/CITEseq_blasts_pseudobulk_counts.csv",
      "data/patient_xy/CITEseq_sample_metadata.csv",
      "results/patient_xy/xypair_dosage_TMM_CPM"
    )
  }
  if (length(paths) != 3L) {
    stop("Usage: Rscript xy_dosage_cpm.R COUNTS.csv METADATA.csv OUTPUT_PREFIX")
  }
  ipf <- paths[[1L]]
  input <- read.csv(ipf,header = T, check.names = F, stringsAsFactors = F)
  #gcn <- names(input)[1]
  genes <- trimws(as.character(input[[1]]))
  sample_data <- input[-1]

  counts <- as.matrix(sample_data)
  rownames(counts) <- genes
  storage.mode(counts) <- "double"

  dge <- edgeR::DGEList(counts = counts)
  edgeR::calcNormFactors(dge,method = "TMM")
  tmm_cpm <- edgeR::cpm(dge, normalized.lib.sizes = T, log = F)

  # Metadata in one line
  mfn <- paths[[2L]]
  meta <- readr::read_csv(mfn, show_col_types = F)

  all_pairs <- c(
    EIF1AY = "EIF1AX",
    USP9Y  = "USP9X",
    TMSB4Y = "TMSB4X",
    PRKY   = "PRKX",
    ZFY    = "ZFX",
    RPS4Y1 = "RPS4X",
    DDX3Y  = "DDX3X",
    KDM5D  = "KDM5C",
    UTY    = "KDM6A"
  )

  # Align metadata and matrix samples, placing the reference group first.
  meta2 <- meta |>
    dplyr::filter(sample_id %in% colnames(tmm_cpm)) |>
    dplyr::arrange(
      factor(group, levels = c("CN-Male", "LOY")),
      sample_id
    )

  sample_cpm <- t(tmm_cpm[, meta2$sample_id, drop = FALSE])
  pairs <- all_pairs[
    names(all_pairs) %in% colnames(sample_cpm) &
      unname(all_pairs) %in% colnames(sample_cpm)
  ]
  xy_result <- xy_dosage_cpm(
    cpm = sample_cpm,
    gender = stats::setNames(meta2$sex, meta2$sample_id),
    group = stats::setNames(meta2$group, meta2$sample_id),
    condition = NULL,
    pairs = pairs,
    reference_group = "CN-Male",
    out_prefix = paths[[3L]],
    width = 14,
    height = 8
  )



  writeLines(capture.output(sessionInfo()),
             file.path(dirname(paths[[3L]]), "sessionInfo.txt"))
}

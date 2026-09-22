#!/usr/bin/env Rscript
# Render the Y-signature UMAP from saved expression and lineage objects.
# Run from the repository root; optional arguments override the input paths.

source("config/paths.R")
paths <- loy_paths("single_cell")
args <- getOption("loy.script_args", commandArgs(trailingOnly = TRUE))
if (length(args) > 3L) {
  stop("Usage: Rscript scripts/run.R y-signature-umap ",
       "[cohort.qs] [lineage.qs] [output_prefix]")
}
cohort_file <- if (length(args) >= 1L) args[[1L]] else
  file.path(paths$data, "figure_inputs", "y_signature_cohort.qs")
lineage_file <- if (length(args) >= 2L) args[[2L]] else
  file.path(paths$data, "figure_inputs", "y_signature_lineage.qs")
output_prefix <- if (length(args) >= 3L) args[[3L]] else
  file.path(paths$results, "figures", "fig3", "F3a_ysig_umap")
for (input_file in c(cohort_file, lineage_file)) {
  if (!file.exists(input_file)) stop("Input object not found: ", input_file)
}

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(ggrastr)
})
source("analysis/human_citeseq/R/plot_style.R")
source("analysis/human_citeseq/R/figure_panels.R")

# Preserve the saved expression values for this direct rendering workflow.
seurat_filtered <- qs::qread(cohort_file)
lineage_object <- qs::qread(lineage_file)
if (!"lineage" %in% colnames(lineage_object[[]])) {
  stop("The lineage object must contain a lineage metadata column.")
}
lineage_annotation <- list(labels = data.frame(
  barcode = colnames(lineage_object),
  lineage = as.character(lineage_object$lineage)
))
rm(lineage_object)
invisible(gc())
files <- fig_final_F3a_ysig_umap(
  seurat_filtered, lineage_annotation, out_basepath = output_prefix
)
cat("Files written:\n", paste(files, collapse = "\n"), "\n")

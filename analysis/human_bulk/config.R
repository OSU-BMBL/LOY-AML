# Shared inputs and outputs for the human bulk RNA-seq workflows.
source("config/paths.R")
bulk_paths <- loy_paths("bulk")

bulk_output <- function(...) {
  path <- file.path(bulk_paths$results, ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

bulk_files <- list(
  cohort_counts = file.path(bulk_paths$data, "counts", "study.csv"),
  cohort_metadata = file.path(bulk_paths$data, "metadata", "study.csv"),
  alliance_counts = file.path(bulk_paths$data, "counts", "alliance.tsv"),
  alliance_metadata = file.path(bulk_paths$data, "metadata", "alliance.xlsx"),
  montri_counts = file.path(bulk_paths$data, "counts", "montri.csv"),
  montri_metadata = file.path(bulk_paths$data, "metadata", "montri.csv"),
  counts = file.path(bulk_paths$results, "prepared", "counts.csv"),
  metadata = file.path(bulk_paths$results, "prepared", "samples.csv"),
  fractions = file.path(bulk_paths$data, "deconvolution", "cell_fractions.csv"),
  gene_sets = file.path(bulk_paths$data, "gene_sets", "msigdb_human.qs"),
  signature = file.path(bulk_paths$root, "analysis", "human_bulk", "resources", "loy_signature.txt"),
  zfy_reference = file.path(bulk_paths$data, "references", "zfy_knockdown.xlsx"),
  volcano = file.path(bulk_paths$data, "figure_inputs", "volcano"),
  deg_overlap = file.path(bulk_paths$data, "figure_inputs", "deg_overlap"),
  zfy_deg = file.path(bulk_paths$data, "figure_inputs", "zfy_overlap", "loy_vs_cn_male.csv"),
  hallmark = file.path(bulk_paths$data, "figure_inputs", "pathways", "hallmark.csv"),
  go_pathways = file.path(bulk_paths$data, "figure_inputs", "pathways", "go_selected.xlsx"),
  go_development = file.path(bulk_paths$data, "figure_inputs", "pathways", "go_development.xlsx"),
  go_immune = file.path(bulk_paths$data, "figure_inputs", "pathways", "go_immune.xlsx"),
  zfy_reactome = file.path(bulk_paths$data, "figure_inputs", "pathways", "zfy_overlap_reactome.tsv")
)

bulk_dirs <- list(
  prepared = bulk_output("prepared"),
  differential_expression = bulk_output("tables", "differential_expression"),
  enrichment = bulk_output("tables", "enrichment")
)

bulk_contrast_name <- function(groups) {
  contrasts <- c(
    "LOY-Male_vs_CN_AML-Male" = "loy_vs_cn_male",
    "LOY-Male_vs_CN_AML-Female" = "loy_vs_cn_female",
    "CN_AML-Male_vs_CN_AML-Female" = "cn_male_vs_cn_female"
  )
  unname(contrasts[paste(groups, collapse = "_vs_")])
}

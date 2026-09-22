# R/input.R
# Load and prepare filtered Seurat object for pipeline

#' Load the cohort-filtered Seurat object
#'
#' @param path Path to the filtered .qs file
#' @return Seurat object with RNA assay set and normalized
load_filtered_seurat <- function(
    path = file.path(loy_paths("single_cell")$results, "objects", "cohort.qs")) {
  if (!file.exists(path)) {
    stop("Filtered Seurat object not found: ", path,
         "\nCreate this file with the preprocessing workflow or supply its path.")
  }
  cat("Loading Seurat object from", path, "...\n")
  obj <- qs::qread(path)
  DefaultAssay(obj) <- "RNA"
  # Join Seurat v5 RNA layers before accessing counts and normalizing.
  if (inherits(obj[["RNA"]], "Assay5")) {
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    cat("Joined RNA assay layers (Seurat v5)\n")
  }
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  cat("Loaded:", ncol(obj), "cells,", nrow(obj), "genes\n")
  cat("Groups:", paste(names(table(obj$group)), collapse = ", "), "\n")
  cat("Cell types:", paste(names(table(obj$cell_type)), collapse = ", "), "\n")
  expected_groups <- c("LOY", "CN-Male", "CN-Female")
  missing <- setdiff(expected_groups, unique(obj$group))
  if (length(missing) > 0) {
    warning("Missing expected groups: ", paste(missing, collapse = ", "))
  }
  cat("Group x Cell type:\n")
  print(table(obj$group, obj$cell_type))
  obj
}

#' Subset Seurat object to a specific cell type
#'
#' @param obj Seurat object
#' @param celltype Character string matching cell_type metadata column
#' @return Subsetted Seurat object
subset_celltype <- function(obj, celltype) {
  sub <- subset(obj, subset = cell_type == celltype)
  cat("Subset", celltype, ":", ncol(sub), "cells\n")
  sub
}

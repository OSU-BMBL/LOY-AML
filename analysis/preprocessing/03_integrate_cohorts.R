#!/usr/bin/env Rscript

source("config/paths.R")
source("analysis/preprocessing/config.R")
preprocessing <- preprocessing_paths()

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(qs)
  library(dplyr)
})

# -------------------------------
# Parameters
# -------------------------------
LOY_FILE <- preprocessing$loy
OUTPUT_FILE <- preprocessing$integrated
NPCS_TO_CALCULATE <- 50
NPCS_TO_USE <- 30

# Load reference objects exported by 02_preprocess_reference.Rmd.
loy_file <- LOY_FILE
sample_manifest <- read.csv(preprocessing$reference_samples, stringsAsFactors = FALSE)
if (!"sample_id" %in% names(sample_manifest)) stop("Reference sample manifest requires sample_id.")
cn_files <- file.path(
  preprocessing$reference_objects, paste0(sample_manifest$sample_id, ".qs")
)
if (!file.exists(loy_file)) stop("Missing LOY object: ", loy_file)
if (length(cn_files) == 0) stop("No reference objects in: ", preprocessing$reference_objects)
if (!all(file.exists(cn_files))) stop("Run reference preparation for every sample in the manifest.")

# Load LOY group
loy_obj <- qread(loy_file)
loy_obj@meta.data$group <- "LOY"

# Load CN-Male group
cn_objs <- lapply(cn_files, function(f) {
  obj <- qread(f)
  obj@meta.data$group <- "CN-Male"
  return(obj)
})

# Merge all objects
all_objs <- c(list(loy_obj), cn_objs)
cat("Merging", length(all_objs), "Seurat objects...\n")

integrated_seurat_obj <- merge(all_objs[[1]], y = all_objs[-1])

# -------------------------------
# RNA processing & Harmony integration
# -------------------------------
cat("Running RNA integration pipeline...\n")

DefaultAssay(integrated_seurat_obj) <- "RNA"

integrated_seurat_obj <- NormalizeData(integrated_seurat_obj, verbose = FALSE)
integrated_seurat_obj <- FindVariableFeatures(
  integrated_seurat_obj, selection.method = "vst", nfeatures = 2000,
  verbose = FALSE
)
integrated_seurat_obj <- ScaleData(integrated_seurat_obj, verbose = FALSE)
integrated_seurat_obj <- RunPCA(integrated_seurat_obj, npcs = NPCS_TO_CALCULATE, verbose = FALSE)

cat("Running Harmony to correct for batch effects...\n")
integrated_seurat_obj <- RunHarmony(
  object = integrated_seurat_obj,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  verbose = FALSE
)

cat("--- Harmony integration complete ---\n")

# Downstream analysis
integrated_seurat_obj <- FindNeighbors(integrated_seurat_obj, reduction = "harmony", dims = 1:NPCS_TO_USE)
integrated_seurat_obj <- RunUMAP(integrated_seurat_obj, reduction = "harmony", dims = 1:NPCS_TO_USE)
integrated_seurat_obj <- FindClusters(integrated_seurat_obj, graph.name = "RNA_snn")

# -------------------------------
# ADT processing (if present)
# -------------------------------
if ("ADT" %in% Assays(integrated_seurat_obj)) {
  cat("Processing ADT assay...\n")
  integrated_seurat_obj <- NormalizeData(integrated_seurat_obj, assay = "ADT", normalization.method = "CLR")
  integrated_seurat_obj <- ScaleData(integrated_seurat_obj, assay = "ADT")
}


# Attach supplied annotations by cell barcode.
annotations <- read.csv(preprocessing$cell_annotations, stringsAsFactors = FALSE)
if (!"NAME" %in% names(annotations)) stop("Cell annotations require a NAME column.")
if (anyNA(annotations$NAME) || any(!nzchar(annotations$NAME)) ||
    anyDuplicated(annotations$NAME)) {
  stop("Cell annotation NAME values must be nonempty and unique.")
}

# Keep both origin fields for downstream sample assignment and visualization.
integrated_seurat_obj$orig.ident.x <- as.character(integrated_seurat_obj$orig.ident)
annotation_columns <- setdiff(
  names(annotations), c("NAME", names(integrated_seurat_obj@meta.data))
)
matched_rows <- match(colnames(integrated_seurat_obj), annotations$NAME)
for (column in annotation_columns) {
  integrated_seurat_obj[[column]] <- annotations[[column]][matched_rows]
}

integrated_seurat_obj <- FindClusters(
  integrated_seurat_obj, graph.name = "RNA_snn", resolution = 0.6
)
DefaultAssay(integrated_seurat_obj) <- "RNA"
Idents(integrated_seurat_obj) <- "seurat_clusters"
qsave(integrated_seurat_obj, OUTPUT_FILE)

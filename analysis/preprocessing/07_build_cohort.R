#!/usr/bin/env Rscript

source("config/paths.R")
source("analysis/preprocessing/config.R")
preprocessing <- preprocessing_paths()
# Construct the strict adult CN/LOY cohort.
# Rebuild Seurat object with 3 strict groups: CN-Male (46,XY adult), CN-Female (46,XX adult), LOY
#
# CN-Male: 7 verified adult 46,XY donors (AML0160, AML3133, AML0114, AML0310, AML0361, AML2451, AML2910)
# CN-Female: 2 adult 46,XX donors (AML2123, AML4897) — processed from raw 10x data
# LOY: all existing LOY cells from the master object
#
# Removes: AML3762 (trisomy 8), AML4340+AML1371 (partial LOY), AML0024 (complex karyotype),
#          pediatric samples, tiny pseudobulk samples

library(Seurat)
library(qs)
library(harmony)
library(dplyr)

# --- Configuration ---
master_path   <- preprocessing$shared_genes
metadata_path <- preprocessing$reference_metadata
cluster_path  <- preprocessing$reference_annotations
raw_data_dir  <- preprocessing$reference_counts
out_path      <- preprocessing$cohort
out_dir       <- dirname(out_path)

NPCS_TO_CALCULATE <- 50
NPCS_TO_USE       <- 30

# Strict CN-Male adult donors to keep
cn_male_donors <- c("AML0160", "AML3133", "AML0114", "AML0310",
                     "AML0361", "AML2451", "AML2910")

# CN-Female donors and their raw 10x pool directories
cn_female_donors <- c("AML2123", "AML4897")
cn_female_pools  <- c("2020-07-29-AML2123-c1",
                       "2020-07-29-AML2123-c2",
                       "2020-10-06-AML4897")

# QC thresholds (matching original preprocessing)
NCOUNT_RNA_MIN <- 1000
NCOUNT_RNA_MAX <- 30000
PERCENT_MT_MAX <- 15

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: print summary ---
print_summary <- function(obj, label) {
  cat("\n=== ", label, " ===\n", sep = "")
  cat("Total cells:", ncol(obj), "\n")

  cat("\n-- Cells per group --\n")
  print(table(obj$group))

  if ("cell_type" %in% colnames(obj@meta.data)) {
    cat("\n-- Cells per cell_type --\n")
    print(table(obj$cell_type))

    cat("\n-- Cross-tab: group x cell_type --\n")
    print(table(obj$group, obj$cell_type))
  }

  ident_col <- if ("orig.ident.x" %in% colnames(obj@meta.data)) "orig.ident.x" else "orig.ident"
  cat("\n-- Unique", ident_col, "per group --\n")
  for (g in sort(unique(obj$group))) {
    vals <- unique(obj@meta.data[[ident_col]][obj$group == g])
    cat("  ", g, ":", length(vals), "samples —", paste(vals, collapse = ", "), "\n")
  }
  cat("\n")
}

# ============================================================
# PART 1: Load master object and filter to LOY + strict CN-Male
# ============================================================
cat("========================================\n")
cat("PART 1: Load and filter master object\n")
cat("========================================\n")

cat("Loading master Seurat object from", master_path, "...\n")
obj <- qread(master_path)
cat("Loaded:", ncol(obj), "cells,", nrow(obj), "genes\n")

# Determine orig.ident column
ident_col <- if ("orig.ident.x" %in% colnames(obj@meta.data)) "orig.ident.x" else "orig.ident"
cat("Using ident column:", ident_col, "\n")

# Load cell-to-donor metadata (skip TYPE row at row 2)
cat("Loading combined metadata...\n")
meta_raw <- read.csv(metadata_path, stringsAsFactors = FALSE)
required_metadata <- c("NAME", "biosample_id", "donor_id", "sex", "disease__ontology_label")
if (!all(required_metadata %in% names(meta_raw))) {
  stop("Reference metadata must contain: ", paste(required_metadata, collapse = ", "))
}
# Remove TYPE row (biosample_id == "group")
meta <- meta_raw[meta_raw$biosample_id != "group", ]
cat("Metadata: ", nrow(meta), " cells with donor mapping\n")

# Map each cell in the master object to a donor_id
cat("Mapping cells to donor_id...\n")
cell_names <- colnames(obj)
donor_map <- setNames(meta$donor_id, meta$NAME)
obj@meta.data$donor_id <- unname(donor_map[cell_names])

# LOY cells won't be in GSE185381 metadata — they have NA donor_id
n_loy   <- sum(obj$group == "LOY", na.rm = TRUE)
n_mapped <- sum(!is.na(obj$donor_id))
n_na     <- sum(is.na(obj$donor_id))
cat("  LOY cells:", n_loy, "\n")
cat("  Cells with donor mapping:", n_mapped, "\n")
cat("  Cells without donor mapping (expect ~LOY count):", n_na, "\n")

# Show donor_id distribution for CN-Male cells
cn_cells <- obj$group == "CN-Male"
cat("\n-- CN-Male donor_id distribution (before filtering) --\n")
print(table(obj$donor_id[cn_cells], useNA = "ifany"))

# Filter: keep LOY cells + cells from strict CN-Male donors
keep_loy    <- obj$group == "LOY"
keep_cn     <- !is.na(obj$donor_id) & obj$donor_id %in% cn_male_donors
keep_cells  <- keep_loy | keep_cn

cat("\n-- Filtering decisions --\n")
cat("  Keeping LOY cells:", sum(keep_loy), "\n")
cat("  Keeping strict CN-Male cells:", sum(keep_cn), "\n")
cat("  Removing:", sum(!keep_cells), "cells\n")

# Show what's being removed
removed_donors <- unique(obj$donor_id[!keep_cells & !is.na(obj$donor_id)])
cat("  Removed donors:", paste(removed_donors, collapse = ", "), "\n")

obj_filtered <- subset(obj, cells = colnames(obj)[keep_cells])

print_summary(obj_filtered, "FILTERED (LOY + strict CN-Male)")

# Free master object
rm(obj)
gc()

# ============================================================
# PART 2: Process CN-Female from raw 10x data
# ============================================================
cat("========================================\n")
cat("PART 2: Process CN-Female from raw 10x\n")
cat("========================================\n")

# Load clustering metadata for cell_type mapping
cat("Loading clustering metadata...\n")
cluster_meta_raw <- read.csv(cluster_path, stringsAsFactors = FALSE)
if (!all(c("NAME", "Broad_cell_identity") %in% names(cluster_meta_raw))) {
  stop("Reference annotations require NAME and Broad_cell_identity.")
}
cluster_meta <- cluster_meta_raw[cluster_meta_raw$Broad_cell_identity != "group", ]
cat("Clustering metadata:", nrow(cluster_meta), "cells\n")

# Build Broad_cell_identity -> cell_type mapping
# Match the labels used in the existing pipeline
broad_to_celltype <- function(broad) {
  malignant_types <- c("HSC", "MPP", "GMP", "Granulocyte",
                        "CD14+ monocyte", "HLA-II+ monocyte", "CD16+ monocyte",
                        "CD11c+", "DC precursor", "Ery", "MEP")
  dplyr::case_when(
    broad %in% malignant_types ~ "Leukemia",
    broad == "CD8+ T"          ~ "CD8 T",
    broad == "CD4+ T"          ~ "CD4 T",
    broad == "NK"              ~ "NK",
    broad == "B"               ~ "B cell",
    broad == "gd T"            ~ "Unidentified",
    broad == "MAIT"            ~ "Unidentified",
    broad == "Plasmablast"     ~ "Unidentified",
    broad == "Plasma cell"     ~ "Unidentified",
    broad == "Pre-B"           ~ "Unidentified",
    broad == "Pro-B"           ~ "Unidentified",
    broad == "LymP"            ~ "Unidentified",
    broad == "pDC"             ~ "Unidentified",
    broad == "cDC1"            ~ "Unidentified",
    broad == "cDC2"            ~ "Unidentified",
    broad == "Megakaryocyte"   ~ "Unidentified",
    broad == "Perivascular cell" ~ "Unidentified",
    TRUE                       ~ "Unidentified"
  )
}

# Build lookup: NAME -> Broad_cell_identity
cluster_lookup <- setNames(cluster_meta$Broad_cell_identity, cluster_meta$NAME)

# Process each female pool
female_objs <- list()

for (pool_name in cn_female_pools) {
  cat("\n--- Processing pool:", pool_name, "---\n")

  pool_path <- file.path(raw_data_dir, pool_name)
  if (!dir.exists(pool_path)) {
    stop("Pool directory not found: ", pool_path)
  }

  # Read 10x data (returns list with Gene Expression + Antibody Capture)
  input_data <- Read10X(pool_path)

  # Create Seurat object from Gene Expression only
  if (is.list(input_data)) {
    counts <- input_data$`Gene Expression`
  } else {
    counts <- input_data
  }
  sobj <- CreateSeuratObject(counts = counts, project = pool_name)
  cat("  Raw barcodes:", ncol(sobj), "\n")

  # Rename barcodes to pool:barcode format (strip -1 suffix)
  new_names <- paste0(pool_name, ":", gsub("-1$", "", colnames(sobj)))
  sobj <- RenameCells(sobj, new.names = new_names)

  # Filter to CN-Female cells using combined_metadata
  female_cells_in_meta <- meta$NAME[meta$donor_id %in% cn_female_donors &
                                     meta$sex == "female" &
                                     meta$disease__ontology_label == "acute myeloid leukemia"]
  keep_female <- intersect(colnames(sobj), female_cells_in_meta)
  cat("  Cells matching CN-Female metadata:", length(keep_female), "\n")

  if (length(keep_female) == 0) {
    cat("  WARNING: No female cells found in this pool, skipping\n")
    next
  }

  sobj <- subset(sobj, cells = keep_female)
  cat("  After metadata filter:", ncol(sobj), "\n")

  # QC: percent mitochondrial
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

  # QC filter
  n_before <- ncol(sobj)
  sobj <- subset(sobj,
                 subset = nCount_RNA > NCOUNT_RNA_MIN &
                          nCount_RNA < NCOUNT_RNA_MAX &
                          percent.mt < PERCENT_MT_MAX)
  cat("  After QC filter:", ncol(sobj), "(removed", n_before - ncol(sobj), ")\n")

  # Map cell types from clustering metadata
  matched_broad <- cluster_lookup[colnames(sobj)]
  sobj@meta.data$Broad_cell_identity <- unname(matched_broad)
  sobj@meta.data$cell_type <- broad_to_celltype(unname(matched_broad))

  n_unmapped <- sum(is.na(matched_broad))
  if (n_unmapped > 0) {
    cat("  WARNING:", n_unmapped, "cells not found in clustering metadata -> Unidentified\n")
    sobj$cell_type[is.na(sobj$cell_type)] <- "Unidentified"
  }

  cat("  Cell type breakdown:\n")
  print(table(sobj$cell_type))

  # Set group and donor metadata
  sobj$group <- "CN-Female"
  sobj@meta.data$donor_id <- unname(donor_map[colnames(sobj)])

  female_objs[[pool_name]] <- sobj
  cat("  Done:", ncol(sobj), "cells\n")
}

# Merge female pools
cat("\nMerging", length(female_objs), "CN-Female pools...\n")
if (length(female_objs) == 1) {
  obj_female <- female_objs[[1]]
} else {
  obj_female <- merge(female_objs[[1]], y = female_objs[-1])
}
cat("Total CN-Female cells:", ncol(obj_female), "\n")

cat("\n-- CN-Female cell_type summary --\n")
print(table(obj_female$cell_type))
cat("\n-- CN-Female donor_id summary --\n")
print(table(obj_female$donor_id))

# Free individual pools
rm(female_objs)
gc()

# ============================================================
# PART 3: Merge and re-integrate
# ============================================================
cat("\n========================================\n")
cat("PART 3: Merge and re-integrate\n")
cat("========================================\n")

cat("Merging filtered male+LOY (", ncol(obj_filtered), " cells) with female (",
    ncol(obj_female), " cells)...\n", sep = "")
obj_merged <- merge(obj_filtered, y = obj_female)
cat("Merged object:", ncol(obj_merged), "cells,", nrow(obj_merged), "genes\n")

# Free intermediate objects
rm(obj_filtered, obj_female)
gc()

# Standard Seurat + Harmony pipeline
cat("NormalizeData...\n")
DefaultAssay(obj_merged) <- "RNA"
obj_merged <- NormalizeData(obj_merged, verbose = FALSE)

cat("FindVariableFeatures...\n")
obj_merged <- FindVariableFeatures(obj_merged, selection.method = "vst",
                                    nfeatures = 2000, verbose = FALSE)

cat("ScaleData...\n")
obj_merged <- ScaleData(obj_merged, verbose = FALSE)

cat("RunPCA...\n")
obj_merged <- RunPCA(obj_merged, npcs = NPCS_TO_CALCULATE, verbose = FALSE)

cat("RunHarmony (batch correction across pools)...\n")
# Use orig.ident for LOY cells, project name for CN-Female
# The harmony group.by.vars should be the pool/sample identifier
# For the filtered object, orig.ident.x is the pool name
# For the female object, orig.ident is the pool name
# After merge, check which column is available
harmony_col <- if ("orig.ident.x" %in% colnames(obj_merged@meta.data)) "orig.ident.x" else "orig.ident"
# For CN-Female cells, orig.ident.x may be NA — fill with orig.ident
if (harmony_col == "orig.ident.x") {
  na_mask <- is.na(obj_merged@meta.data[["orig.ident.x"]])
  if (any(na_mask)) {
    obj_merged@meta.data[["orig.ident.x"]][na_mask] <- as.character(
      obj_merged@meta.data[["orig.ident"]][na_mask]
    )
    cat("  Filled", sum(na_mask), "NA orig.ident.x values from orig.ident\n")
  }
}
cat("  Harmony group.by.vars:", harmony_col, "\n")
cat("  Unique batches:", length(unique(obj_merged@meta.data[[harmony_col]])), "\n")

obj_merged <- RunHarmony(
  object        = obj_merged,
  group.by.vars = harmony_col,
  reduction.use = "pca",
  verbose       = FALSE
)
cat("Harmony complete.\n")

cat("FindNeighbors...\n")
obj_merged <- FindNeighbors(obj_merged, reduction = "harmony", dims = 1:NPCS_TO_USE)

cat("RunUMAP...\n")
obj_merged <- RunUMAP(obj_merged, reduction = "harmony", dims = 1:NPCS_TO_USE)

cat("FindClusters...\n")
obj_merged <- FindClusters(obj_merged, graph.name = "RNA_snn")

print_summary(obj_merged, "INTEGRATED 3-GROUP OBJECT")

# Save
cat("Saving to", out_path, "...\n")
qsave(obj_merged, out_path)
cat("Saved successfully.\n")

# ============================================================
# PART 4: Validation
# ============================================================
cat("\n========================================\n")
cat("PART 4: Validation\n")
cat("========================================\n")

# Total cells per group
cat("\n-- Total cells per group --\n")
print(table(obj_merged$group))

# Group x cell_type cross-tab
if ("cell_type" %in% colnames(obj_merged@meta.data)) {
  cat("\n-- Group x cell_type cross-tab --\n")
  print(table(obj_merged$group, obj_merged$cell_type))
}

# Verify zero AML3762 cells
ident_col_check <- if ("orig.ident.x" %in% colnames(obj_merged@meta.data)) "orig.ident.x" else "orig.ident"
aml3762_count <- sum(grepl("AML3762", obj_merged@meta.data[[ident_col_check]]), na.rm = TRUE)
cat("\nAML3762 cells (should be 0):", aml3762_count, "\n")
if (aml3762_count > 0) warning("AML3762 cells found! Trisomy 8 samples not fully removed.")

# Verify zero AML3762 via donor_id
aml3762_donor <- sum(obj_merged$donor_id == "AML3762", na.rm = TRUE)
cat("AML3762 by donor_id (should be 0):", aml3762_donor, "\n")

# Verify zero pediatric cells
if ("ap_aml_age" %in% colnames(obj_merged@meta.data)) {
  n_ped <- sum(obj_merged$ap_aml_age == "pediatric_AML", na.rm = TRUE)
  cat("Pediatric cells (should be 0):", n_ped, "\n")
  if (n_ped > 0) warning("Pediatric cells found!")
} else {
  cat("No ap_aml_age column — pediatric check N/A (filtered by donor list)\n")
}

# Unique orig.ident per group
cat("\n-- Unique orig.ident per group --\n")
for (g in sort(unique(obj_merged$group))) {
  vals <- unique(obj_merged@meta.data[[ident_col_check]][obj_merged$group == g])
  cat("  ", g, ":", length(vals), "samples\n")
  for (v in sort(vals)) cat("    ", v, "\n")
}

# CN-Female cell count
n_female <- sum(obj_merged$group == "CN-Female")
cat("\nCN-Female total cells:", n_female, "\n")

# Donor breakdown
cat("\n-- Donor breakdown --\n")
print(table(obj_merged$group, obj_merged$donor_id, useNA = "ifany"))

cat("\n====== ALL DONE ======\n")

# Free memory
rm(obj_merged)
gc()

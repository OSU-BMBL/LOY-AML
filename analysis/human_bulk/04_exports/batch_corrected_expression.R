#!/usr/bin/env Rscript

source("analysis/human_bulk/config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
})

counts_path <- bulk_files$counts
metadata_path <- bulk_files$metadata
output_dir <- bulk_output("exports", "batch_corrected_expression")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- as.matrix(read.csv(
  counts_path,
  row.names = 1,
  check.names = FALSE
))
metadata <- read.csv(
  metadata_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(metadata$sample_id),
  all(metadata$sample_id %in% colnames(counts))
)

counts <- counts[, metadata$sample_id, drop = FALSE]
counts_input <- round(pmax(counts, 0))
storage.mode(counts_input) <- "integer"

metadata$Group <- ifelse(
  metadata$Cyto.Group == "LOY",
  "MLOY",
  ifelse(metadata$Sex == "Male", "MCN", "FCN")
)
metadata$Group <- factor(metadata$Group, levels = c("MCN", "MLOY", "FCN"))
metadata$Batch <- factor(metadata$Batch)

stopifnot(
  identical(colnames(counts_input), metadata$sample_id),
  identical(
    as.integer(table(metadata$Group)[c("MLOY", "MCN", "FCN")]),
    c(23L, 218L, 203L)
  )
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_input,
  colData = metadata,
  design = ~ Batch + Group
)
vsd <- vst(dds, blind = FALSE)

group_design <- model.matrix(~ Group, data = as.data.frame(colData(vsd)))
batch_corrected <- removeBatchEffect(
  assay(vsd),
  batch = colData(vsd)$Batch,
  design = group_design
)

stopifnot(
  identical(dim(batch_corrected), dim(counts_input)),
  identical(rownames(batch_corrected), rownames(counts_input)),
  identical(colnames(batch_corrected), metadata$sample_id),
  all(is.finite(batch_corrected))
)

metadata_out <- data.frame(
  sample_id = metadata$sample_id,
  Group = as.character(metadata$Group),
  Cyto.Group = metadata$Cyto.Group,
  Sex = metadata$Sex,
  Batch = as.character(metadata$Batch),
  stringsAsFactors = FALSE
)

counts_out <- data.frame(
  gene = rownames(counts_input),
  counts_input,
  check.names = FALSE
)
corrected_out <- data.frame(
  gene = rownames(batch_corrected),
  round(batch_corrected, 6),
  check.names = FALSE
)

group_columns <- split(metadata$sample_id, metadata$Group)
lookup <- data.frame(gene = rownames(batch_corrected), check.names = FALSE)
for (group_name in c("MLOY", "MCN", "FCN")) {
  group_matrix <- batch_corrected[, group_columns[[group_name]], drop = FALSE]
  lookup[[paste0(group_name, "_mean_batch_corrected_VST")]] <- rowMeans(group_matrix)
  lookup[[paste0(group_name, "_median_batch_corrected_VST")]] <- apply(
    group_matrix,
    1,
    median
  )
}
lookup$MLOY_minus_MCN <- (
  lookup$MLOY_mean_batch_corrected_VST -
    lookup$MCN_mean_batch_corrected_VST
)
lookup$MLOY_minus_FCN <- (
  lookup$MLOY_mean_batch_corrected_VST -
    lookup$FCN_mean_batch_corrected_VST
)
lookup$MCN_minus_FCN <- (
  lookup$MCN_mean_batch_corrected_VST -
    lookup$FCN_mean_batch_corrected_VST
)
lookup[, -1] <- round(lookup[, -1], 6)

write.csv(
  counts_out,
  file.path(output_dir, "input_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  corrected_out,
  file.path(output_dir, "batch_corrected_vst.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  metadata_out,
  file.path(output_dir, "sample_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  lookup,
  file.path(output_dir, "gene_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO.txt"))

cat("Output directory:", output_dir, "\n")
cat("Genes:", nrow(counts_input), "\n")
cat("Samples:", ncol(counts_input), "\n")
print(table(metadata$Group, metadata$Batch))

source("analysis/human_bulk/config.R")

counts_path <- bulk_files$counts
meta_path <- bulk_files$metadata
ciber_path <- bulk_files$fractions
out_dir <- bulk_output("exports", "hsc_weighted_expression")

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)
ciber <- read.csv(ciber_path)
names(ciber)[1] <- "sample_id"

hsc_columns <- c("HSC", "HSC.like", "Prog.like")
stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(meta$sample_id),
  !anyDuplicated(ciber$sample_id),
  all(hsc_columns %in% names(ciber)),
  setequal(colnames(counts), meta$sample_id),
  setequal(colnames(counts), ciber$sample_id)
)

sample_ids <- colnames(counts)
meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]
ciber <- ciber[match(sample_ids, ciber$sample_id), , drop = FALSE]
stopifnot(identical(sample_ids, meta$sample_id), identical(sample_ids, ciber$sample_id))

counts[counts < 0] <- 0
library_sizes <- colSums(counts)
if (any(library_sizes <= 0)) stop("At least one sample has a non-positive library size.")

# Library-size CPM, without TMM normalization.
linear_cpm <- sweep(as.matrix(counts), 2, library_sizes / 1e6, "/")
hsc_total_fraction <- rowSums(ciber[, hsc_columns, drop = FALSE])
hsc_weighted_cpm <- sweep(linear_cpm, 2, hsc_total_fraction, "*")
hsc_weighted_log2cpm <- log2(hsc_weighted_cpm + 1)

group <- ifelse(
  meta$Cyto.Group == "LOY" & meta$Sex == "Male", "LOY",
  ifelse(
    meta$Cyto.Group == "CN_AML" & meta$Sex == "Male", "CN-Male",
    ifelse(meta$Cyto.Group == "CN_AML" & meta$Sex == "Female", "CN-Female", "Other")
  )
)

matrix_out <- cbind(
  gene = rownames(hsc_weighted_log2cpm),
  as.data.frame(hsc_weighted_log2cpm, check.names = FALSE)
)
matrix_nonlog_out <- cbind(
  gene = rownames(hsc_weighted_cpm),
  as.data.frame(hsc_weighted_cpm, check.names = FALSE)
)
metadata_out <- data.frame(
  sample_id = sample_ids,
  group = group,
  Cyto.Group = meta$Cyto.Group,
  Sex = meta$Sex,
  Batch = meta$Batch,
  HSC = ciber$HSC,
  HSC_like = ciber$HSC.like,
  Prog_like = ciber$Prog.like,
  HSC_total_fraction = hsc_total_fraction,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(digits = 15)
write.csv(
  matrix_out,
  file.path(out_dir, "hsc_weighted_log2cpm.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  matrix_nonlog_out,
  file.path(out_dir, "hsc_weighted_cpm.csv"),
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  metadata_out,
  file.path(out_dir, "sample_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf(
  "Exported %d genes across %d samples to %s\n",
  nrow(hsc_weighted_log2cpm), ncol(hsc_weighted_log2cpm), out_dir
))

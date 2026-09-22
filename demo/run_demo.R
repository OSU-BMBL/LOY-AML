#!/usr/bin/env Rscript
# Synthetic example only. Run from the repository root.
required <- c("Seurat", "Matrix", "ggplot2", "DESeq2", "tibble", "dplyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
if (!file.exists("demo/data/counts.csv")) stop("Run from the repository root.")

suppressPackageStartupMessages(library(Seurat))
source("analysis/human_citeseq/R/plot_style.R")
source("analysis/human_citeseq/R/differential_expression.R")

counts <- as.matrix(read.csv("demo/data/counts.csv", row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer"
metadata <- read.csv("demo/data/cell_metadata.csv", row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(counts), rownames(metadata)))
obj <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE),
                                  meta.data = metadata, project = "SYNTHETIC_DEMO")
results <- run_deseq2_pseudobulk(obj, celltype = "CD8 T", min_cells = 10)
if (is.null(results)) stop("No result returned.")
out <- "demo/output"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(out, "synthetic_pseudobulk_DE.csv"), row.names = FALSE)

plot_data <- results[is.finite(results$log2FC) & !is.na(results$padj), ]
p <- ggplot2::ggplot(plot_data, ggplot2::aes(log2FC, -log10(pmax(padj, 1e-300)))) +
  ggplot2::geom_point(size = 1, alpha = 0.6, color = "#377EB8") +
  ggplot2::theme_classic() +
  ggplot2::labs(title = "Synthetic example — no patient data", x = "log2 fold change",
                y = "−log10 adjusted P value")
ggplot2::ggsave(file.path(out, "synthetic_volcano.pdf"), p, width = 5, height = 4)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("Synthetic example outputs written to ", out)

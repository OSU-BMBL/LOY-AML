# Installation

## Software environments

| Workflow | Software | Version record |
|---|---|---|
| Primary bulk RNA-seq processing | nf-core/rnaseq 3.14.0, Nextflow 23.10.1, STAR 2.6.1d, Salmon 1.10.1 | [Pipeline versions](../environments/bulk_primary_software_versions.yml) |
| Bulk differential expression and enrichment | R 4.4.1, DESeq2 1.44.0, fgsea 1.30.0, msigdbr 7.5.1 | [R session](../environments/bulk_hallmark_saved_session.txt) |
| Single-cell annotation | R 4.4.1, Seurat 5.1.0, SeuratObject 5.0.2 | [R session](../environments/citeseq_annotation_saved_session.txt) |
| Single-cell analysis export | R 4.4.1, Seurat 5.4.0, SeuratObject 5.3.0 | [R session](../environments/citeseq_export_session_20260803.txt) |
| Y-signature UMAP rendering | R 4.4.3, Seurat 5.5.1, SeuratObject 5.4.0, ggplot2 4.0.3 | [Package versions](../environments/figure2a_R_packages_20260903.tsv) |

Use the environment corresponding to the workflow. The environment records include the operating system and additional package versions. The [dependency inventory](../environments/source_dependencies.tsv) lists packages referenced by the source code.

## Hardware and installation time

The downstream workflows use standard CPUs and require no specialized hardware. Allow 16 GB RAM for a typical desktop setup; 32 GB or more is recommended for full single-cell datasets.

On a desktop with 4 CPU cores, 16 GB RAM and an SSD, estimated installation time is **20–60 minutes for the demo environment** and **1–2 hours for the full workflow dependencies**. These are planning estimates; compilation and network speed affect the time required.

## Synthetic demo

The demo requires R, Seurat, Matrix, DESeq2, ggplot2, dplyr and tibble. For an R 4.4 environment:

```r
install.packages(c("BiocManager", "Seurat", "ggplot2", "dplyr", "tibble"))
BiocManager::install(version = "3.19")
BiocManager::install("DESeq2")
```

For specific package versions, use the package records above and the corresponding CRAN or Bioconductor archives. See the official [Bioconductor installation guide](https://bioconductor.org/install/).

From the repository root:

```sh
Rscript demo/run_demo.R
```

## Analysis dependencies

Notebook workflows use `rmarkdown`, `knitr` and Pandoc. The command runner uses `targets` and `tidyselect` for single-cell workflows.

- **Bulk RNA-seq:** DESeq2, edgeR, limma, GSVA, fgsea, msigdbr, qs, tidyverse, EnhancedVolcano, readxl, openxlsx, ggpubr, enrichR and biomaRt, together with the plotting packages imported by each script.
- **Single-cell RNA-seq / CITE-seq:** Seurat, SeuratObject, harmony, targets, qs, Matrix, DESeq2, presto, fgsea, msigdbr, dplyr, tidyr and tibble.
- **Single-cell figures:** ggplot2, ggrepel, ggrastr, ragg, patchwork, ggbeeswarm, data.table, stringr, clusterProfiler, org.Hs.eg.db, ComplexHeatmap, scplotter, scales, RColorBrewer, ggpubr, ggtext, sp and MASS.
- **CD8 trajectory:** monocle3, SeuratWrappers and igraph.

Bulk enrichment uses an MSigDB table with `gs_cat` and `gs_subcat` columns. The single-cell enrichment functions use the newer `collection` and `subcollection` interface to msigdbr. Keep the corresponding package environments separate.

PDF and raster figure exports use Cairo, ragg and Arial. Source-package installation may require a compiler and platform libraries. Downstream R workflows run on CPU.

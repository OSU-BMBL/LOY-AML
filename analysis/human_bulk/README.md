# Human bulk RNA-seq

Differential expression, pathway enrichment, cell-fraction analysis and Y-chromosome expression in the integrated AML cohort.

## Run a workflow

Run commands from the repository root:

```sh
Rscript scripts/run.R bulk-prepare
Rscript scripts/run.R bulk-de
Rscript scripts/run.R bulk-hallmark
Rscript scripts/run.R plot-volcano
```

`bulk-prepare` combines cohort count matrices and sample metadata. The analysis notebooks read its outputs from `results/bulk/prepared/`. Figure workflows read the count matrices, cell fractions or curated result tables specified in [config.R](config.R).

## Organisation

- **01_preprocessing** — cohort harmonisation and sample summaries.
- **02_analysis** — DESeq2 contrasts, GO/Reactome/Hallmark enrichment and ZFY comparisons.
- **03_figures** — volcano plots, enrichment plots, expression heatmaps and cell-fraction figures.
- **04_exports** — batch-corrected VST and HSC-fraction-weighted expression matrices.
- **resources** — the seven-gene Y-chromosome signature.

## Inputs and outputs

Input paths are defined in [config.R](config.R). Place cohort data, reference annotations and curated plotting tables under `data/bulk/`; see [input formats](../../docs/inputs.md#bulk-rna-seq) for the directory layout and required columns. The signature resource is read directly from this repository.

Results are written to `results/bulk/`: prepared matrices in `prepared/`, statistical tables in `tables/`, plots in `figures/` and expression exports in `exports/`. Each figure workflow has its own output folder.

Set `LOY_DATA_DIR` and `LOY_RESULTS_DIR` to use other data and output locations. Both settings refer to parent directories containing a `bulk/` subdirectory.

## Figure and export commands

| Output | Command |
| --- | --- |
| Volcano plots | `plot-volcano` |
| Hallmark bar and dot plots | `plot-hallmark-bars`, `plot-hallmark-dots` |
| GO pathway plots | `plot-go` |
| DEG overlap diagrams | `plot-deg-overlap` |
| Cell-fraction plots and correlations | `plot-cell-fractions`, `plot-deconvolution` |
| Y-gene expression and signature heatmaps | `plot-y-genes`, `plot-y-signature` |
| Developmental, immune and ZFY-overlap pathways | `plot-zfy` |
| ZFY differential-expression comparison | `bulk-zfy` |
| Expression matrices | `export-vst`, `export-hsc` |

Use `Rscript scripts/run.R <command>` for each entry. Package requirements are listed in [software requirements](../../docs/installation.md).

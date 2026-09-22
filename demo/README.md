# Synthetic pseudobulk example

A small synthetic dataset illustrates donor-level count aggregation and differential-expression analysis with the single-cell workflow's `run_deseq2_pseudobulk()` function.

## Dataset

| File | Content |
|---|---|
| `data/counts.csv` | Count matrix with 400 synthetic genes and 120 cells. |
| `data/cell_metadata.csv` | Six donors, three per group, with 20 cells per donor. |

Gene and cell identifiers begin with `SIMGENE` and `DEMO_`. The dataset was generated using seed `20260922`, gamma-distributed counts, donor scaling and simulated group effects.

## Run

Install the [demo dependencies](../docs/installation.md), then run from the repository root:

```sh
Rscript demo/run_demo.R
```

**Estimated runtime:** 1–5 minutes on a desktop with 4 CPU cores, 16 GB RAM and an SSD, excluding installation. This is an estimate, not a measured benchmark.

## Expected outputs

The script writes to `demo/output/`:

| File | Content |
|---|---|
| `synthetic_pseudobulk_DE.csv` | Differential-expression results for genes retained after filtering. |
| `synthetic_volcano.pdf` | Gene fold changes plotted against adjusted P values. |
| `sessionInfo.txt` | R and package versions used for execution. |

The results table contains `gene`, `baseMean`, `log2FC`, `stat`, `pvalue`, `padj`, `method`, `cell_type`, `n_donors_LOY` and `n_donors_CN`. Each retained row should report three donors per group and method `DESeq2_pseudobulk`.

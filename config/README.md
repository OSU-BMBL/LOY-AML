# Configuration

## Shared paths

All commands use the same project layout:

```text
data/
  bulk/
  single_cell/
results/
  bulk/
  single_cell/
```

[paths.R](paths.R) resolves these directories. To use data or result directories outside the repository, set `LOY_DATA_DIR` and `LOY_RESULTS_DIR` before running a command. Each directory contains separate `bulk/` and `single_cell/` subdirectories.

For example, in R:

```r
Sys.setenv(LOY_DATA_DIR = "/path/to/data", LOY_RESULTS_DIR = "/path/to/results")
```

The command runner sets the repository root automatically. When sourcing a script interactively, start R in the repository root.

## Workflow inputs

- [Bulk configuration](../analysis/human_bulk/config.R) defines count matrices, metadata, gene sets and curated figure tables.
- [Preprocessing configuration](../analysis/preprocessing/config.R) defines sequencing inputs, reference metadata and intermediate Seurat objects.
- Single-cell analysis reads the cohort object from `results/single_cell/objects/cohort.qs` and regulon tables from `data/single_cell/regulons/`.

Change individual paths in the relevant configuration file when using differently named inputs. File schemas are described in [Data inputs](../docs/inputs.md).

## Commands

[workflows.csv](workflows.csv) maps command names to scripts and notebooks. [single_cell_targets.csv](single_cell_targets.csv) provides short names for single-cell targets. List commands with:

```sh
Rscript scripts/run.R --list
```

Examples:

```sh
Rscript scripts/run.R bulk-prepare
Rscript scripts/run.R bulk-de
Rscript scripts/run.R plot-volcano
Rscript scripts/run.R single-cell cd8-de
Rscript scripts/run.R single-cell bach2-expression bach2-regulon
```

Notebook reports are written to each workflow's `reports/` directory under `results/`. Single-cell target objects are stored in `results/single_cell/_targets/`.

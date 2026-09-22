# LOY–AML

Human bulk RNA-seq, single-cell RNA-seq, CITE-seq and figure-generation code accompanying *Non-concomitant Y chromosome loss shapes the genetic and immune landscape of acute myeloid leukemia and predicts clinical outcome*.

[Bulk RNA-seq](analysis/human_bulk/) · [Single-cell analysis](analysis/human_citeseq/) · [Preprocessing](analysis/preprocessing/) · [Data inputs](docs/inputs.md)

## 1. System requirements

- **R:** 4.4.1 for the recorded analysis environments; 4.4.3 for the Y-signature renderer.
- **Software:** workflow-specific dependencies and versions are listed in the [installation guide](docs/installation.md), [package inventory](environments/source_dependencies.tsv) and [environment records](environments/).
- **Operating system:** recorded analysis sessions used macOS 15.6 on Apple Silicon; additional session records include macOS 26.5.1. These records describe the original analysis environments.
- **Hardware:** a standard CPU desktop; no GPU or other specialized hardware is required by the downstream workflows.

Timing estimates below assume a desktop with **4 CPU cores, 16 GB RAM and an SSD**. They are estimates, not measured benchmarks.

## 2. Installation

```sh
git clone https://github.com/OSU-BMBL/LOY-AML.git
cd LOY-AML
```

For the demo, install the following in R 4.4:

```r
install.packages(c("BiocManager", "Seurat", "ggplot2", "dplyr", "tibble"))
BiocManager::install(version = "3.19")
BiocManager::install("DESeq2")
```

For study workflows, install the additional packages in the [installation guide](docs/installation.md). **Estimated installation time:** 20–60 minutes for the demo environment; 1–2 hours for the full workflow dependencies. Source compilation and download speed affect installation time.

## 3. Demo

The repository includes synthetic counts for **400 genes and 120 cells from six donors** in `demo/data/`. Run from the repository root:

```sh
Rscript demo/run_demo.R
```

Expected files in `demo/output/`:

- `synthetic_pseudobulk_DE.csv` — gene-level fold changes, test statistics, P values and adjusted P values.
- `synthetic_volcano.pdf` — a volcano plot of the synthetic results.
- `sessionInfo.txt` — R and package versions.

**Estimated runtime:** 1–5 minutes, excluding installation. See the [demo guide](demo/README.md) for input and output schemas.

## 4. Use with study data

Arrange count matrices, metadata and reference inputs according to [Data inputs](docs/inputs.md). Set file locations in the [workflow configuration](config/README.md), then select a command:

```sh
Rscript scripts/run.R --list
Rscript scripts/run.R bulk-de
Rscript scripts/run.R single-cell bach2-expression
```

Inputs default to `data/bulk/` and `data/single_cell/`; outputs are written to the corresponding directories under `results/`. The workflow READMEs describe preprocessing order, analysis commands and figure inputs. Study data are distributed separately from the code.

## Code availability and license

Source code and the synthetic demo are available in this repository and as a [single ZIP archive](https://github.com/OSU-BMBL/LOY-AML/archive/refs/heads/main.zip). Project code is released under the [MIT License](LICENSE). Third-party dependencies retain their respective licenses.

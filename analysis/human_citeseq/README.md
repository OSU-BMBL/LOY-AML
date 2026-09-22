# Single-cell RNA-seq and CITE-seq

Human AML transcriptome analysis with cell-lineage annotation, CD8 T-cell state
classification, differential expression, pathway enrichment and figure export.
The workflow uses Seurat for cell-level analysis, DESeq2 for donor pseudobulks
and `targets` to manage analysis dependencies.

## Run an analysis

Install the [software dependencies](../../docs/installation.md), prepare the
[study inputs](../../docs/inputs.md) and run commands from the repository root.
The [preprocessing workflow](../preprocessing/README.md) produces the cohort
object consumed by the analysis.

```sh
# CD8-state differential expression and its upstream dependencies
Rscript scripts/run.R single-cell cd8-de

# BACH2 expression panel
Rscript scripts/run.R single-cell bach2-expression

# Complete analysis and figure workflow
Rscript scripts/run.R single-cell
```

Each named target builds only its required dependencies. The complete workflow
also uses the precomputed regulon tables described below.

## Input data

The analysis reads `results/single_cell/objects/cohort.qs`, a Seurat object with
an RNA counts layer, UMAP and Harmony reductions, and these metadata fields:

- `group`: `LOY`, `CN-Male` or `CN-Female`.
- `cell_type` and `Cell_type_identity`: broad and detailed cell annotations.
- `orig.ident` and `orig.ident.x`: sample and integration-batch identifiers.
- `donor_id` and `sample_id`: donor identifiers; either may be missing where the
  other identifies the donor.

Regulon figures read CSV files under `data/single_cell/regulons/`:

```text
regulons/
├── cd8_states/
│   ├── Naive/
│   ├── Effector/
│   └── Dysfunctional/
└── blasts/
    └── Leukemia_LOY_vs_CN_diff_regulons.csv
```

The BACH2 panel uses one naive-state CSV with `auc_per_cell` in its filename,
a `Cell` column and a `BACH2_act` or other BACH2 regulon column. The naive-state
regulon summary uses `Naive/CD8_T_Naive_LOY_vs_CN_diff_regulons.csv`. Differential
regulon tables contain `regulon`, `log2FC` and `padj` columns.

## Standalone figures

The Y-signature renderer reads saved expression and lineage objects directly.
Its default inputs are `data/single_cell/figure_inputs/y_signature_cohort.qs`
and `y_signature_lineage.qs`. The first contains the normalized RNA data and
UMAP embedding; the second supplies cell barcodes and a `lineage` metadata
column. Optional arguments override both input paths and the output prefix.

```sh
Rscript scripts/run.R y-signature-umap
Rscript scripts/run.R y-signature-umap cohort.qs lineage.qs results/y_signature
```

The CD8 trajectory uses the `cd8_states` target, performs RPCA integration and
learns a Monocle3 trajectory rooted in naive cells. It additionally requires
`monocle3`, `SeuratWrappers` and `igraph`.

```sh
Rscript scripts/run.R single-cell trajectory-input
Rscript scripts/run.R cd8-trajectory
```

## Source layout

| Component | Source |
| --- | --- |
| Input preparation and annotations | `R/input.R`, `R/nk_annotation.R`, `R/lineage_annotation.R`, `R/cd8_states.R` |
| Differential expression and enrichment | `R/differential_expression.R`, `R/pathway_enrichment.R` |
| Gene dosage and antigen presentation | `R/y_gene_dosage.R`, `R/blast_apm.R` |
| Regulon table import | `R/scenic_import.R` |
| Figure panels, styles and exports | `R/figure_panels.R`, `R/plot_helpers.R`, `R/plot_style.R` |
| Workflow dependencies | `_targets.R` |

Figures are written to `results/single_cell/figures/` as PDF, TIFF and PNG;
`targets` stores intermediate objects in `results/single_cell/_targets/`.
Configure the data and results roots with `LOY_DATA_DIR` and `LOY_RESULTS_DIR`
as described in the [input guide](../../docs/inputs.md).

# Data inputs

Place study inputs under `data/bulk/` and `data/single_cell/`. Generated tables and objects are written under `results/`. Paths can be changed in the [configuration files](../config/README.md).

## Bulk RNA-seq

```text
data/bulk/
├── counts/            study.csv, alliance.tsv, montri.csv
├── metadata/          study.csv, alliance.xlsx, montri.csv
├── deconvolution/     cell_fractions.csv
├── gene_sets/         msigdb_human.qs
├── references/        zfy_knockdown.xlsx
└── figure_inputs/
    ├── volcano/       Contrast-specific DE tables
    ├── deg_overlap/   Contrast-specific DE tables
    ├── zfy_overlap/   loy_vs_cn_male.csv
    └── pathways/     Curated enrichment tables
```

### Counts and sample metadata

`bulk-prepare` combines the three cohort inputs. Count matrices contain gene identifiers and sample columns. Study and integrated metadata use `sample_id`, `Cyto.Group`, `Race`, `Age`, `Sex` and `Batch`. The Alliance workbook supplies the same fields except `Batch`, which is assigned by the script. Montri metadata identifies samples by `Lab_ID`.

The Alliance count table retains its source layout: column 2 contains gene identifiers, columns 1 and 3 contain annotations, and remaining columns contain sample counts. The study matrix uses its first column for gene identifiers; the Montri matrix uses gene row names. The preparation script harmonizes cohort-specific group and sample labels.

Prepared inputs are written to:

```text
results/bulk/prepared/
├── counts.csv
└── samples.csv
```

The count matrix has gene row names and sample columns. `samples.csv` contains matching sample identifiers. Group values are `LOY` and `CN_AML`; sex values are `Male` and `Female`. Already prepared datasets can be placed at these paths, or configured through `bulk_files$counts` and `bulk_files$metadata` in the [bulk configuration](../analysis/human_bulk/config.R).

### Deconvolution and gene sets

- **cell_fractions.csv:** CIBERSORTx sample identifiers in `Mixture`, compartment fractions and QC columns. Relevant compartments include `HSC`, `HSC-like`, `Prog-like`, `T` and `CTL`; QC columns include `P-value`, `Correlation` and `RMSE`.
- **msigdb_human.qs:** an MSigDB table serialized with `qs`, containing `gs_cat`, `gs_subcat`, `gs_name` and `gene_symbol`.
- **zfy_knockdown.xlsx:** the published ZFY comparison table, sheet `S11A`, with three header rows before columns `Gene`, `padj` and `log2FoldChange`.

The seven-gene Y-chromosome signature is supplied in the repository's `analysis/human_bulk/resources/` directory.

### Figure tables

Volcano and DEG-overlap inputs use three filenames:

```text
loy_vs_cn_male.csv
loy_vs_cn_female.csv
cn_male_vs_cn_female.csv
```

These tables contain gene identifiers in the first column and DESeq2 columns such as `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and `padj`. Place the tables selected for each figure in its corresponding input directory. The ZFY comparison uses `zfy_overlap/loy_vs_cn_male.csv`.

The pathway directory contains:

| File | Columns |
|---|---|
| `hallmark.csv` | `pathway`, `pval`, `padj`, `NES`, `size`, `leadingEdge`; permutation tables also contain `nMoreExtreme` |
| `go_selected.xlsx` | `pathway`, `NES`, `padj`, `size`, `leadingEdge` |
| `go_development.xlsx`, `go_immune.xlsx` | `pathway`, `NES`, `padj`, `size`, `leadingEdge` |
| `zfy_overlap_reactome.tsv` | Enrichr table containing `Term`, `Adjusted P-value`, `Odds Ratio`, `Combined Score` and `Genes` |

Figure inputs preserve the selected genes and pathways for each display. Analysis tables are written to `results/bulk/tables/`; individual file locations can be assigned in the bulk configuration.

## Single-cell RNA-seq and CITE-seq

### Sequencing inputs and annotations

```text
data/single_cell/
├── cellranger/          Pooled LOY count matrices
├── sample_barcodes/     Demultiplexed sample barcodes
├── metadata/           cell_annotations.csv
├── reference/          Reference counts, metadata and sample manifest
├── regulons/           Precomputed SCENIC tables
└── figure_inputs/      Saved objects for direct figure rendering
```

The [preprocessing guide](../analysis/preprocessing/README.md#inputs) defines the Cell Ranger layout, reference sample manifest and cell-annotation columns. Preprocessing creates distinct intermediate objects and writes the downstream cohort to `results/single_cell/objects/cohort.qs`.

### Analysis object

The Seurat cohort object contains raw and normalized RNA layers, UMAP and Harmony reductions, and cell metadata including `group`, `cell_type`, `orig.ident`, `orig.ident.x`, `sample_id`, `donor_id` and `Cell_type_identity`. Group values are `LOY`, `CN-Male` and `CN-Female`.

### Regulon tables and direct figure inputs

SCENIC state-specific tables are placed under `regulons/cd8_states/{Naive,Effector,Dysfunctional}/`. The BACH2 display reads an `auc_per_cell` CSV with cell identifiers in `Cell` and BACH2 activity in `BACH2_act` or a matching regulon column. The blast differential-regulon table is `regulons/blasts/Leukemia_LOY_vs_CN_diff_regulons.csv`.

The direct Y-signature renderer reads `figure_inputs/y_signature_cohort.qs` and `figure_inputs/y_signature_lineage.qs`. The second object supplies a `lineage` metadata column. Both paths can also be provided as command arguments. See the [single-cell guide](../analysis/human_citeseq/README.md) for commands.

## Synthetic example

The repository includes `demo/data/counts.csv` and `demo/data/cell_metadata.csv`, with matching synthetic cell identifiers. See the [demo guide](../demo/README.md).

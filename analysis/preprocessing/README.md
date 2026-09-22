# Single-cell preprocessing

Prepare pooled LOY CITE-seq samples and the GSE185381 reference cohort for downstream single-cell analysis.

## Workflow

Run the numbered stages from the repository root:

```bash
Rscript scripts/run.R preprocess-loy
Rscript scripts/run.R preprocess-reference
Rscript scripts/run.R integrate-cohorts
Rscript scripts/run.R assign-samples
Rscript scripts/run.R annotate-cells
Rscript scripts/run.R filter-genes
Rscript scripts/run.R build-cohort
```

1. **LOY preprocessing:** RNA quality control, Harmony integration and ADT normalization.
2. **Reference preparation:** extract male AML cells into per-sample Seurat objects.
3. **Cohort integration:** combine objects, integrate RNA profiles and attach cell annotations.
4. **Sample assignment:** map LOY barcodes to demultiplexed samples within each pool.
5. **Cell annotation:** inspect markers, apply cluster labels and summarize composition.
6. **Gene selection:** retain shared expressed genes and Y-chromosome genes.
7. **Cohort construction:** select adult CN donors, add female CN reference cells and integrate the three groups.

## Inputs

Paths are defined in [config.R](config.R), using the shared [path configuration](../../config/paths.R).

```text
data/single_cell/
├── cellranger/             Pool1_multi_results, Pool2_multi_results, Pool3_multi_results
├── sample_barcodes/        S1–S6 demultiplexed barcode folders
├── metadata/
│   └── cell_annotations.csv
└── reference/
    ├── counts/            One 10x count-matrix folder per reference sample
    ├── sample_manifest.csv
    ├── metadata.csv
    └── cell_annotations.csv
```

- Each Cell Ranger pool contains `sample_filtered_feature_bc_matrix.h5` with RNA and antibody counts.
- Each demultiplexed sample contains `sample_filtered_feature_bc_matrix/barcodes.tsv.gz`.
- The integration annotation table uses `NAME` as the unique cell barcode; additional columns supply cell annotations. Existing Seurat metadata fields are retained.
- Reference `metadata.csv` contains `NAME`, `biosample_id`, `donor_id`, `sex` and `disease__ontology_label`.
- Reference `sample_manifest.csv` has one `sample_id` per selected sample, matching its folder name under `counts/`. Sample IDs must be unique and present among male AML cells in the reference metadata.
- Reference `cell_annotations.csv` contains `NAME` and `Broad_cell_identity`. Reference cell names use `sample:barcode`, without the 10x `-1` suffix.

## Outputs

Intermediate objects are written to `results/single_cell/preprocessing/`: `loy.qs`, per-sample reference objects, `integrated.qs`, `sample_assigned.qs`, `annotated.qs` and `shared_genes.qs`.

The final object is `results/single_cell/objects/cohort.qs`, used by the [single-cell analysis workflow](../human_citeseq/). All stages retain their inputs and write a separate output.

Gene selection queries Ensembl through `biomaRt`. Software dependencies are listed in the [installation guide](../../docs/installation.md).

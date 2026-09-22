# KG-1a, mouse CITE-seq and X-Y dosage analyses

## Contents

| File | Analysis |
| --- | --- |
| `bulkxy_story.Rmd` | KG-1a bulk RNA-seq, modules, pathways and X-Y dosage |
| `LOY_chemo_CITEseq.Rmd` | Mouse CITE-seq QC, LOY classification and pseudobulk analysis |
| `HSC.Rmd` | X-Y dosage from the supplied patient expression spreadsheet |
| `xy_dosage_cpm.R` | Shared dosage function and a separately executable patient-pseudobulk example |
| `build_tx2gene2_all_primary.py` | Ensembl 116 transcript-to-gene mapping from the kallisto-index FASTA |

## Setup

Keep these five files together. Use the R environment in which the analyses
were verified, with the packages loaded by each script and `rmarkdown` installed.
The R syntax requires R 4.1 or newer. The transcript-mapping builder uses Python
3 and either a frozen primary-map TSV or a MySQL/MariaDB client for its Ensembl
query; PyYAML is needed only when reading its optional YAML configuration.
This component does not install or upgrade packages.

## Input and output locations

Edit file locations under `params:` at the top of each Rmd. Absolute paths may
point to data already on your computer; relative paths are resolved from the
folder containing these scripts when run as instructed below. Processed study
data are supplied separately and are not included here.

| Workflow | Paths to set |
| --- | --- |
| Bulk | `kallisto_out`, `tx2gene_file`, `alignment_rate`, `pairs_root`, `root_path` |
| Mouse | `p1_cellranger_dir`, `p2_cellranger_dir`, `feature_reference_file`, `reference_dir`, `analysis_root` |
| HSC | `xlsx_file`, `xy_helper`, `output_prefix` |

`root_path` is the bulk output root. `pairs_root` is the original input directory
containing the pair files. Keep its original subdirectory structure: module
analysis checks `results/figures/pairs.csv`, `figures/pairs.csv`, then `pairs.csv`;
the dosage analysis reads `pairs.csv`. This preserves the original lookup order.

Use the complete original alignment summary: bulk sample selection retains the
original index-based order. Preserve kallisto sample directory names, the mouse
Cell Ranger `outs/per_sample_outs/S*/sample_filtered_feature_bc_matrix.h5`
structure, the pool-level filtered H5 files, and the reference `genes/genes.gtf`
or `genes/genes.gtf.gz` file. The actual tx2gene mapping is tab-delimited even
when its filename ends in `.csv`; select the correct existing file, not a
renamed or converted substitute.

Choose NEW output directories for validation runs. Existing exports may be
overwritten if their output directory is reused.

## Run

Open a terminal in this folder. Run only the desired analysis, one at a time:

```bash
Rscript --vanilla -e 'rmarkdown::render("bulkxy_story.Rmd", output_format="html_document", knit_root_dir=getwd())'
Rscript --vanilla -e 'rmarkdown::render("LOY_chemo_CITEseq.Rmd", output_format="html_document", knit_root_dir=getwd())'
Rscript --vanilla -e 'rmarkdown::render("HSC.Rmd", output_format="html_document", knit_root_dir=getwd())'
```

In RStudio, open the chosen Rmd and use Knit with Knit Directory set to Document
Directory. Run Console chunks only after setting the working directory to this
folder.

To run the patient-pseudobulk example separately:

```bash
Rscript --vanilla xy_dosage_cpm.R data/patient_xy/CITEseq_blasts_pseudobulk_counts.csv data/patient_xy/CITEseq_sample_metadata.csv results/patient_xy/xypair_dosage_TMM_CPM
```

The helper's example does not run when it is sourced by `HSC.Rmd`.

Rebuilding the transcript mapping is optional when the existing mapping table
is supplied. Use the SAME FASTA and frozen mapping as the original build:

```bash
python3 build_tx2gene2_all_primary.py --fasta data/reference/transcriptome.fa.gz --primary-map data/reference/ensembl116_all_alt_to_primary.tsv --out results/tx2gene/tx2gene2.tsv --qc-out results/tx2gene/tx2gene2.mapping_qc.tsv
```

Those are example path locations; set them to the actual input files. Without
`--primary-map`, the builder queries the release-locked Ensembl database.
MSigDB and KEGG steps may also access their external resources. This is not a
fully offline environment bundle.

## Outputs and versions

Each Rmd produces its original plots and tables plus an HTML report. Bulk exports
include `DE_LOY_vs_WT.xlsx` and `enrichment_LOY_vs_WT.xlsx`; mouse exports include
`LOY_chemo_review_data.xlsx` and saved R objects. HSC and the standalone X-Y
example save their dosage plots. Bulk, HSC and the standalone example also save
`sessionInfo.txt`; mouse records `Session_info` in its review workbook.

The analytical calculations and original figure settings are retained. Compare
the reported numerical outputs with the verified reference results after changing
paths. Upstream signature weighting is not reconstructed by the HSC spreadsheet
reader; it uses the supplied processed values.


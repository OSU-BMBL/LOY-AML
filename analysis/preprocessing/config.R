# File interfaces for human CITE-seq preprocessing.
# Paths are resolved from the repository root by config/paths.R.
preprocessing_paths <- function() {
  paths <- loy_paths("single_cell")
  intermediate <- file.path(paths$results, "preprocessing")
  reference_objects <- file.path(intermediate, "reference_samples")
  objects <- file.path(paths$results, "objects")
  for (directory in c(intermediate, reference_objects, objects)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }

  list(
    root = paths$root,
    cellranger = file.path(paths$data, "cellranger"),
    sample_barcodes = file.path(paths$data, "sample_barcodes"),
    reference_counts = file.path(paths$data, "reference", "counts"),
    reference_samples = file.path(paths$data, "reference", "sample_manifest.csv"),
    reference_metadata = file.path(paths$data, "reference", "metadata.csv"),
    reference_annotations = file.path(paths$data, "reference", "cell_annotations.csv"),
    cell_annotations = file.path(paths$data, "metadata", "cell_annotations.csv"),
    reference_objects = reference_objects,
    loy = file.path(intermediate, "loy.qs"),
    integrated = file.path(intermediate, "integrated.qs"),
    sample_assigned = file.path(intermediate, "sample_assigned.qs"),
    annotated = file.path(intermediate, "annotated.qs"),
    shared_genes = file.path(intermediate, "shared_genes.qs"),
    cohort = file.path(objects, "cohort.qs")
  )
}

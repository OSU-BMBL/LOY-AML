#!/usr/bin/env Rscript
# Command-line entry point for analysis notebooks, scripts and targets.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Use Rscript scripts/run.R <command>.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- dirname(dirname(script_file))
Sys.setenv(LOY_REPO_ROOT = root)
setwd(root)
source("config/paths.R")

commands <- read.csv("config/workflows.csv", stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[[1L]] %in% c("--help", "-h", "--list")) {
  cat("Usage: Rscript scripts/run.R <command> [arguments]\n\n")
  for (i in seq_len(nrow(commands))) {
    cat(sprintf("  %-23s %s\n", commands$command[[i]], commands$description[[i]]))
  }
  cat("\nSingle-cell targets: single-cell [name ...]\n")
  cell_targets <- read.csv("config/single_cell_targets.csv", stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cell_targets))) {
    cat(sprintf("  %-23s %s\n", cell_targets$name[[i]], cell_targets$description[[i]]))
  }
  cat("Y-signature renderer: y-signature-umap [cohort.qs lineage.qs output_prefix]\n")
  cat("Trajectory renderer: cd8-trajectory [targets_store output_directory]\n")
  quit(status = 0)
}

index <- match(args[[1L]], commands$command)
if (is.na(index)) stop("Unknown command: ", args[[1L]], ". Use --list.", call. = FALSE)
entry <- commands[index, , drop = FALSE]
extra <- args[-1L]
if (length(extra) && !entry$command %in% c("single-cell", "y-signature-umap", "cd8-trajectory")) {
  stop("This command takes no additional arguments. Set input paths in its configuration file.", call. = FALSE)
}
options(loy.script_args = extra)
paths <- loy_paths(entry$workflow)

if (entry$type == "rmd") {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install rmarkdown to render notebooks.", call. = FALSE)
  report_dir <- file.path(paths$results, "reports")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  rmarkdown::render(
    entry$script, knit_root_dir = root, output_dir = report_dir,
    envir = new.env(parent = globalenv())
  )
} else if (entry$type == "targets") {
  if (!requireNamespace("targets", quietly = TRUE)) stop("Install targets to run this workflow.", call. = FALSE)
  store <- file.path(paths$results, "_targets")
  if (length(extra)) {
    cell_targets <- read.csv("config/single_cell_targets.csv", stringsAsFactors = FALSE)
    matches <- match(extra, cell_targets$name)
    extra[!is.na(matches)] <- cell_targets$target[matches[!is.na(matches)]]
    targets::tar_make(script = entry$script, store = store, names = tidyselect::all_of(extra))
  } else {
    targets::tar_make(script = entry$script, store = store)
  }
} else {
  sys.source(entry$script, envir = globalenv(), chdir = FALSE)
}

# Shared project locations. Paths can be redirected with environment variables.

loy_root <- function() {
  root <- normalizePath(
    Sys.getenv("LOY_REPO_ROOT", unset = getwd()),
    winslash = "/", mustWork = TRUE
  )
  if (!file.exists(file.path(root, "config", "paths.R"))) {
    stop("Run from the repository root or set LOY_REPO_ROOT.", call. = FALSE)
  }
  root
}

loy_absolute_path <- function(path, root = loy_root()) {
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\|~)", path)
  if (!absolute) path <- file.path(root, path)
  normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
}

loy_paths <- function(workflow = c("bulk", "single_cell")) {
  workflow <- match.arg(workflow)
  root <- loy_root()
  data_root <- loy_absolute_path(Sys.getenv("LOY_DATA_DIR", unset = "data"), root)
  result_root <- loy_absolute_path(Sys.getenv("LOY_RESULTS_DIR", unset = "results"), root)
  list(
    root = root,
    data = file.path(data_root, workflow),
    results = file.path(result_root, workflow)
  )
}

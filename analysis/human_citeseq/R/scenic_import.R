# Import precomputed state-specific SCENIC result tables.

import_scenic_cd8_by_state_results <- function(
    dir = file.path(loy_paths("single_cell")$data, "regulons", "cd8_states")) {
  states <- c("Naive", "Effector", "Dysfunctional")
  result <- list(Naive = NULL, Effector = NULL, Dysfunctional = NULL,
                 status = "missing_input", dir = dir)
  if (!dir.exists(dir)) {
    warning("SCENIC input directory not found: ", dir,
            ". Place state-specific CSV files in Naive, Effector and ",
            "Dysfunctional subdirectories.")
    return(result)
  }
  for (state in states) {
    csv_files <- list.files(file.path(dir, state), pattern = "\\.csv$",
                            full.names = TRUE)
    if (length(csv_files)) {
      result[[state]] <- lapply(csv_files, function(path) {
        tryCatch(
          read.csv(path),
          error = function(error) {
            warning("Cannot read SCENIC table ", path, ": ", conditionMessage(error))
            NULL
          }
        )
      })
      names(result[[state]]) <- basename(csv_files)
      result[[state]] <- Filter(Negate(is.null), result[[state]])
      if (!length(result[[state]])) result[[state]] <- NULL
    }
  }
  if (any(!vapply(result[states], is.null, logical(1)))) {
    result$status <- "populated"
  } else {
    warning("No readable state-specific SCENIC CSV files found in ", dir)
  }
  result
}

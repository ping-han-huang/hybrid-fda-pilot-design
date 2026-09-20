# Shared file I/O and FPCA workflow helpers.

ensure_parent_directory <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

load_rdata_environment <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env
}

load_rdata_object <- function(path, object_name) {
  env <- load_rdata_environment(path)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Object '", object_name, "' not found in ", path)
  }
  get(object_name, envir = env, inherits = FALSE)
}

estimate_fpca_collection <- function(pilot_data, subject_sizes, grid_size = 25L) {
  fpca_estimates <- lapply(seq_along(pilot_data), function(dataset_id) {
    estimates_by_size <- lapply(subject_sizes, function(subject_size) {
      size_name <- paste0("N=", subject_size)
      sparse_replicates <- pilot_data[[dataset_id]][[size_name]]
      if (is.null(sparse_replicates)) {
        stop("Sparse pilot data not found for ", size_name,
             " in dataset ", dataset_id)
      }

      lapply(sparse_replicates, function(sparse_data) {
        fdapace::FPCA(
          sparse_data$Ly,
          sparse_data$Lt,
          optns = list(
            nRegGrid = grid_size,
            methodBwCov = "GCV",
            methodBwMu = "GCV",
            FVEthreshold = 1,
            kernel = "gauss",
            dataType = "Sparse"
          )
        )
      })
    })
    names(estimates_by_size) <- paste0("N=", subject_sizes)
    estimates_by_size
  })
  names(fpca_estimates) <- paste0("dataset=", seq_along(pilot_data))
  fpca_estimates
}

# Evaluate the FPCA estimates under exhaustive and PSS search.
source("paper-analysis/functions/evaluation_functions.R")

set.seed(24857)
subject_sizes <- seq(50, 240, by = 10)
grid_size <- 25L
num_observations <- 5L
num_components <- 3L
grid <- seq(1, grid_size, length.out = grid_size)

fpca_estimates <- load_rdata_object(
  "paper-analysis/intermediate/case-study-medfly/fpca_BIBD.Rdata",
  "fpca_estimates"
)
functional_data <- load_rdata_environment(
  "paper-analysis/data/case-study-medfly/MedrealfunctionalDataSets.Rdata"
)
subject_index_offset <- length(functional_data$dataset) - length(subject_sizes)
if (subject_index_offset < 0) {
  stop("Medfly functional-data object does not contain all requested subject sizes")
}

true_curve_getter <- function(dataset_id, subject_size) {
  subject_index <- match(subject_size, subject_sizes) + subject_index_offset
  functional_data$dataset[[subject_index]][[dataset_id]][["Xtrue.X"]]
}

evaluation_results <- evaluate_fpca_collection(
  fpca_estimates = fpca_estimates,
  subject_sizes = subject_sizes,
  true_curve_getter = true_curve_getter,
  grid = grid,
  num_observations = num_observations,
  num_components = num_components,
  true_eigenfunctions = functional_data$Xtrue.ef,
  true_eigenvalues = functional_data$Xtrue.ev,
  true_sigma = functional_data$Xtrue.sigma,
  true_optimal_value = functional_data$FxTrue
)

output_file <- "paper-analysis/intermediate/case-study-medfly/evaluation_BIBD.Rdata"
ensure_parent_directory(output_file)
save(evaluation_results, file = output_file)

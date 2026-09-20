# Evaluate the FPCA estimates under exhaustive and PSS search.
source("paper-analysis/functions/evaluation_functions.R")

set.seed(24857)
subject_sizes <- seq(50, 240, by = 10)
grid_size <- 25L
num_observations <- 5L
num_components <- 5L
grid <- seq(0, 1, length.out = grid_size)

fpca_estimates <- load_rdata_object(
  "paper-analysis/intermediate/simulation-fourier/fpca_BIBD.Rdata",
  "fpca_estimates"
)
functional_data <- load_rdata_environment(
  "paper-analysis/data/simulation-fourier/FourierfunctionalDataSets.Rdata"
)

true_curve_getter <- function(dataset_id, subject_size) {
  true_curves <- functional_data$dataset[[dataset_id]][["Xtrue.X"]]
  true_curves[seq_len(subject_size), , drop = FALSE]
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

output_file <- "paper-analysis/intermediate/simulation-fourier/evaluation_BIBD.Rdata"
ensure_parent_directory(output_file)
save(evaluation_results, file = output_file)

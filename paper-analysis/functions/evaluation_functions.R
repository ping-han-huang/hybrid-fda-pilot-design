# Shared evaluation helpers.

source("paper-analysis/functions/workflow_functions.R")
source("paper-analysis/functions/pss_functions.R")

evaluate_relative_rmse <- function(fpca, true_curves, num_subjects) {
  estimated_curves <- fpca$xiEst %*% t(fpca$phi) +
    matrix(fpca$mu, nrow = num_subjects, ncol = length(fpca$mu), byrow = TRUE)
  relative_mse <- vapply(seq_len(num_subjects), function(i) {
    sum((estimated_curves[i, ] - true_curves[i, ])^2) / sum(true_curves[i, ]^2)
  }, numeric(1))
  sqrt(sum(relative_mse) / num_subjects)
}

find_grid_positions <- function(points, grid) {
  match(points, grid)
}

evaluate_single_design <- function(indices, grid, eigenfunctions, eigenvalues,
                                   sigma_squared) {
  num_candidates <- length(grid)
  if (nrow(eigenfunctions) != num_candidates) {
    stop("length(grid) != nrow(eigenfunctions)")
  }

  eigenvalue_matrix <- diag(eigenvalues, nrow = length(eigenvalues))
  h_matrix <- eigenvalue_matrix %*% t(eigenfunctions)
  gamma_matrix <- eigenfunctions %*% eigenvalue_matrix %*% t(eigenfunctions) +
    diag(sigma_squared, num_candidates)

  sum(diag(
    h_matrix[, indices, drop = FALSE] %*%
      solve(
        gamma_matrix[indices, indices, drop = FALSE],
        t(h_matrix[, indices, drop = FALSE])
      )
  ))
}

exhaustive_search <- function(num_points, grid, eigenfunctions, eigenvalues,
                              sigma_squared) {
  candidate_designs <- combn(length(grid), num_points)
  criterion_values <- apply(
    candidate_designs,
    2,
    evaluate_single_design,
    grid = grid,
    eigenfunctions = eigenfunctions,
    eigenvalues = eigenvalues,
    sigma_squared = sigma_squared
  )
  designs <- apply(candidate_designs, 2, function(indices) grid[indices])

  list(design = designs, criVal = criterion_values)
}

optimal_exhaustive <- function(exhaustive_result) {
  best_index <- which.max(exhaustive_result$criVal)
  list(
    design = exhaustive_result$design[, best_index],
    criVal = exhaustive_result$criVal[best_index]
  )
}

evaluate_are <- function(fpca, true_eigenfunctions, true_eigenvalues,
                         true_sigma, true_optimal_value, grid,
                         num_points, num_components) {
  estimated_pss <- pss_x(
    num_points = num_points,
    candidate_set = grid,
    eigenfunctions = fpca$phi,
    eigenvalues = fpca$lambda,
    sigma_squared = fpca$sigma2
  )

  positions <- find_grid_positions(estimated_pss$design, grid)
  selected_eigenfunctions <- true_eigenfunctions[
    positions, seq_len(num_components), drop = FALSE
  ]
  eigenvalue_matrix <- diag(true_eigenvalues, nrow = length(true_eigenvalues))
  h_matrix <- eigenvalue_matrix %*% t(selected_eigenfunctions)
  gamma_matrix <- selected_eigenfunctions %*% eigenvalue_matrix %*%
    t(selected_eigenfunctions) + diag(true_sigma^2, num_points)
  true_value_at_estimated_design <- sum(
    diag(h_matrix %*% solve(gamma_matrix, t(h_matrix)))
  )

  (true_optimal_value - true_value_at_estimated_design) / true_optimal_value
}

evaluate_fpca_collection <- function(fpca_estimates, subject_sizes,
                                     true_curve_getter, grid, num_observations,
                                     num_components, true_eigenfunctions,
                                     true_eigenvalues, true_sigma,
                                     true_optimal_value) {
  lapply(seq_along(fpca_estimates), function(dataset_id) {
    lapply(subject_sizes, function(subject_size) {
      size_name <- paste0("N=", subject_size)
      fpca_replicates <- fpca_estimates[[dataset_id]][[size_name]]
      if (is.null(fpca_replicates)) {
        stop(
          "FPCA estimates not found for ", size_name,
          " in dataset ", dataset_id
        )
      }
      true_curves <- true_curve_getter(dataset_id, subject_size)

      lapply(fpca_replicates, function(fpca) {
        exhaustive <- exhaustive_search(
          num_points = num_observations,
          grid = grid,
          eigenfunctions = fpca$phi,
          eigenvalues = fpca$lambda,
          sigma_squared = fpca$sigma2
        )
        are <- evaluate_are(
          fpca = fpca,
          true_eigenfunctions = true_eigenfunctions,
          true_eigenvalues = true_eigenvalues,
          true_sigma = true_sigma,
          true_optimal_value = true_optimal_value,
          grid = grid,
          num_points = num_observations,
          num_components = num_components
        )
        relative_rmse <- evaluate_relative_rmse(
          fpca,
          true_curves,
          length(fpca$inputData$Ly)
        )

        list(
          exhaustive = exhaustive,
          relative_rmse = relative_rmse,
          are = are
        )
      })
    })
  })
}

# Among designs attaining the requested estimated-efficiency threshold,
# summarize their true criterion values by either the median or the minimum.
evaluate_design_at_true_efficiency <- function(estimated_exhaustive,
                                               true_exhaustive,
                                               efficiency_level,
                                               efficiency_case = c(
                                                 "median", "worst"
                                               )) {
  efficiency_case <- match.arg(efficiency_case)
  estimated_optimum <- max(estimated_exhaustive$criVal)
  estimated_efficiency <- estimated_exhaustive$criVal / estimated_optimum
  candidate_indices <- which(estimated_efficiency >= efficiency_level)
  true_candidate_values <- true_exhaustive$criVal[candidate_indices]

  switch(
    efficiency_case,
    median = median(true_candidate_values),
    worst = min(true_candidate_values)
  )
}

collect_evaluation_results <- function(results, group_label, subject_sizes,
                                       true_exhaustive, efficiency_levels,
                                       efficiency_cases = c(
                                         "median", "worst"
                                       )) {
  efficiency_cases <- match.arg(
    efficiency_cases,
    choices = c("median", "worst"),
    several.ok = TRUE
  )

  efficiency_rows <- list()
  metric_rows <- list()
  efficiency_index <- 1L
  metric_index <- 1L
  true_optimum <- max(true_exhaustive$criVal)

  for (dataset_id in seq_along(results)) {
    for (subject_index in seq_along(subject_sizes)) {
      current_subject_size <- subject_sizes[subject_index]
      simulations <- results[[dataset_id]][[subject_index]]

      for (simulation_id in seq_along(simulations)) {
        current_result <- simulations[[simulation_id]]

        metric_rows[[metric_index]] <- data.frame(
          dataset = dataset_id,
          group = group_label,
          subject_size = current_subject_size,
          simulation = simulation_id,
          relative_rmse = current_result$relative_rmse,
          are = current_result$are
        )
        metric_index <- metric_index + 1L

        for (efficiency_level in efficiency_levels) {
          for (efficiency_case in efficiency_cases) {
            true_value <- evaluate_design_at_true_efficiency(
              estimated_exhaustive = current_result$exhaustive,
              true_exhaustive = true_exhaustive,
              efficiency_level = efficiency_level,
              efficiency_case = efficiency_case
            )

            efficiency_rows[[efficiency_index]] <- data.frame(
              dataset = dataset_id,
              group = group_label,
              subject_size = current_subject_size,
              simulation = simulation_id,
              efficiency_level = efficiency_level,
              efficiency_case = efficiency_case,
              true_efficiency = true_value / true_optimum
            )
            efficiency_index <- efficiency_index + 1L
          }
        }
      }
    }
  }

  list(
    efficiency = do.call(rbind, efficiency_rows),
    metrics = do.call(rbind, metric_rows)
  )
}

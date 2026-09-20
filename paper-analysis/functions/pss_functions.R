# PSS search used in the simulation and case-study evaluations.

pss_x <- function(num_points, candidate_set, eigenfunctions, eigenvalues,
                  sigma_squared, subset_size = 3L, max_iterations = 3000L,
                  stopping_length = 30L) {
  num_candidates <- length(candidate_set)
  if (nrow(eigenfunctions) != num_candidates) {
    stop("length(candidate_set) != nrow(eigenfunctions)")
  }

  eigenvalue_matrix <- diag(eigenvalues, nrow = length(eigenvalues))
  h_matrix <- eigenvalue_matrix %*% t(eigenfunctions)
  gamma_matrix <- eigenfunctions %*% eigenvalue_matrix %*% t(eigenfunctions) +
    diag(sigma_squared, num_candidates)

  weights <- (eigenfunctions^2) %*% eigenvalues
  optimal_indices <- sort(
    sample.int(num_candidates, num_points, prob = weights)
  )
  efficiency_curve <- rep(NA_real_, max_iterations)

  iteration <- 0L
  repeat {
    iteration <- iteration + 1L
    remaining <- setdiff(seq_len(num_candidates), optimal_indices)
    sampled_remaining <- sort(
      sample.int(length(remaining), subset_size, prob = weights[remaining])
    )
    candidate_indices <- sort(c(remaining[sampled_remaining], optimal_indices))
    subsets <- combn(num_points + subset_size, num_points)

    criterion_values <- vapply(seq_len(ncol(subsets)), function(k) {
      indices <- candidate_indices[subsets[, k]]
      sum(diag(h_matrix[, indices, drop = FALSE] %*%
                 solve(gamma_matrix[indices, indices, drop = FALSE],
                       t(h_matrix[, indices, drop = FALSE]))))
    }, numeric(1))

    best <- which.max(criterion_values)
    optimal_indices <- candidate_indices[subsets[, best]]
    efficiency_curve[iteration] <- criterion_values[best]

    converged <- iteration >= stopping_length &&
      length(unique(
        efficiency_curve[(iteration - stopping_length + 1L):iteration]
      )) == 1L
    if (converged || iteration == max_iterations) break
  }

  list(
    design = candidate_set[optimal_indices],
    efficiency_curve = efficiency_curve[seq_len(iteration)]
  )
}

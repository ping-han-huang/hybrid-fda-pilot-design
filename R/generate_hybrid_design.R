# Example: generate the proposed Hybrid design for a new study.
#
# Run from the repository root:
#   Rscript R/generate_hybrid_design.R
#
# Edit only the design parameters below for your application.

source("R/hybrid_design.R")

# Design parameters ------------------------------------------------------------
N <- 30          # number of subjects
gsize <- 25       # number of candidate observation locations
num_obs <- 5      # observations collected per subject
W <- 0.5          # proportion assigned to the adjacent-pair component
seed <- 3429      # optional seed for reproducibility

# Generate the design ----------------------------------------------------------
design <- generate_hybrid_design(
  N = N,
  gsize = gsize,
  num_obs = num_obs,
  W = W,
  seed = seed
)

# A successful design is an N x gsize incidence matrix.
if (is.matrix(design)) {
  print(design)
  cat("\nRow sums (observations per subject):\n")
  print(table(rowSums(design)))
} else {
  print(design)
}

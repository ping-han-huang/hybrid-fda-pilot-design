# Estimate FPCA from the supplied sparse pilot datasets.
source("paper-analysis/functions/workflow_functions.R")

set.seed(9853)
subject_sizes <- seq(50, 240, by = 10)
grid_size <- 25L

pilot_data <- load_rdata_object("paper-analysis/data/case-study-medfly/pilotMedrealRandom.Rdata", "pilot_randomDesign")
fpca_estimates <- estimate_fpca_collection(
  pilot_data = pilot_data,
  subject_sizes = subject_sizes,
  grid_size = grid_size
)

output_file <- "paper-analysis/intermediate/case-study-medfly/fpca_random.Rdata"
ensure_parent_directory(output_file)
save(fpca_estimates, file = output_file)

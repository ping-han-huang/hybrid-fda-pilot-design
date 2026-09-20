# Generate the supplementary cluster-size RRMSE and ARE plots.

source("paper-analysis/functions/plotting_functions.R")

subject_sizes <- seq(50, 240, by = 10)

group_levels <- c("q = 2", "q = 3", "q = 4", "q = 5")
# Keep the original supplementary cluster-size color palette.
group_colors <- c(
  "q = 2" = "#0072B2",
  "q = 3" = "#E69F00",
  "q = 4" = "#009E73",
  "q = 5" = "#CC79A7"
)

metric_data <- dplyr::bind_rows(
  collect_metric_file(
    "paper-analysis/intermediate/supplementary-cluster-size/evaluation_cluster2.Rdata",
    "q = 2",
    subject_sizes
  ),
  collect_metric_file(
    "paper-analysis/intermediate/supplementary-cluster-size/evaluation_cluster3.Rdata",
    "q = 3",
    subject_sizes
  ),
  collect_metric_file(
    "paper-analysis/intermediate/supplementary-cluster-size/evaluation_cluster4.Rdata",
    "q = 4",
    subject_sizes
  ),
  collect_metric_file(
    "paper-analysis/intermediate/supplementary-cluster-size/evaluation_cluster5.Rdata",
    "q = 5",
    subject_sizes
  )
)

metric_summary <- summarize_metric_results(
  metric_data,
  include_composite = FALSE
)

rm(metric_data)
gc()

save_paper_boxplot(
  data = metric_summary,
  value_column = "rrmse_normalized",
  y_label = "RRMSE",
  group_levels = group_levels,
  group_colors = group_colors,
  subject_sizes = subject_sizes,
  output_file = "paper-analysis/results/supplementary-cluster-size/cluster_size_rrmse.eps",
  y_breaks = seq(0, 1, by = 0.1),
  legend_position = c(0.88, 0.8)
)

save_paper_boxplot(
  data = metric_summary,
  value_column = "are_normalized",
  y_label = "ARE",
  group_levels = group_levels,
  group_colors = group_colors,
  subject_sizes = subject_sizes,
  output_file = "paper-analysis/results/supplementary-cluster-size/cluster_size_are.eps",
  y_breaks = seq(0, 1, by = 0.1),
  legend_position = c(0.88, 0.8)
)

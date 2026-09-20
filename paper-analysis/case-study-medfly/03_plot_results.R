# Generate the paper-format plots for the Medfly case study.

source("paper-analysis/functions/plotting_functions.R")

subject_sizes <- seq(50, 240, by = 10)
efficiency_levels <- c(0.99, 0.97, 0.95)
efficiency_cases <- c("median", "worst")

group_levels <- c("BIBD", "Random", "Hybrid")
group_colors <- c(
  "BIBD" = "#88CCEE",
  "Random" = "#CC6677",
  "Hybrid" = "#DDCC77"
)

true_exhaustive <- load_rdata_object(
  "paper-analysis/data/case-study-medfly/Medreal_XtrueExh.Rdata",
  "Xtrue.Exh"
)

bibd_results <- collect_evaluation_file(
  "paper-analysis/intermediate/case-study-medfly/evaluation_BIBD.Rdata",
  "BIBD",
  subject_sizes,
  true_exhaustive,
  efficiency_levels,
  efficiency_cases
)
random_results <- collect_evaluation_file(
  "paper-analysis/intermediate/case-study-medfly/evaluation_random.Rdata",
  "Random",
  subject_sizes,
  true_exhaustive,
  efficiency_levels,
  efficiency_cases
)
hybrid_results <- collect_evaluation_file(
  "paper-analysis/intermediate/case-study-medfly/evaluation_hybrid.Rdata",
  "Hybrid",
  subject_sizes,
  true_exhaustive,
  efficiency_levels,
  efficiency_cases
)

efficiency_data <- dplyr::bind_rows(
  bibd_results$efficiency,
  random_results$efficiency,
  hybrid_results$efficiency
)
metric_data <- dplyr::bind_rows(
  bibd_results$metrics,
  random_results$metrics,
  hybrid_results$metrics
)

rm(bibd_results, random_results, hybrid_results, true_exhaustive)
gc()

efficiency_summary <- summarize_efficiency_results(efficiency_data)
metric_summary <- summarize_metric_results_median_scaled(metric_data)

rm(efficiency_data, metric_data)
gc()

save_efficiency_plots(
  efficiency_summary = efficiency_summary,
  efficiency_levels = efficiency_levels,
  efficiency_cases = efficiency_cases,
  group_levels = group_levels,
  group_colors = group_colors,
  subject_sizes = subject_sizes,
  output_directory = "paper-analysis/results/case-study-medfly",
  file_prefix = "medfly"
)

save_paper_boxplot(
  data = metric_summary,
  value_column = "rrmse_median_scaled",
  y_label = "RRMSE",
  group_levels = group_levels,
  group_colors = group_colors,
  subject_sizes = subject_sizes,
  output_file = "paper-analysis/results/case-study-medfly/medfly_rrmse.eps",
  y_breaks = pretty(metric_summary$rrmse_median_scaled),
  legend_position = c(0.88, 0.8)
)

save_paper_boxplot(
  data = metric_summary,
  value_column = "are_median_scaled",
  y_label = "ARE",
  group_levels = group_levels,
  group_colors = group_colors,
  subject_sizes = subject_sizes,
  output_file = "paper-analysis/results/case-study-medfly/medfly_are.eps",
  y_breaks = pretty(metric_summary$are_median_scaled),
  legend_position = c(0.88, 0.8)
)

save_paper_boxplot(
  data = metric_summary,
  value_column = "median_scaled_criterion",
  y_label = "Composite Criterion",
  group_levels = group_levels,
  group_colors = group_colors,
  subject_sizes = subject_sizes,
  output_file = "paper-analysis/results/case-study-medfly/medfly_composite_criterion.eps",
  y_breaks = pretty(metric_summary$median_scaled_criterion),
  legend_position = c(0.88, 0.8)
)

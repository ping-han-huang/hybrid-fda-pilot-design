# Shared data-reduction and plotting helpers for paper figures.

source("paper-analysis/functions/evaluation_functions.R")

collect_evaluation_file <- function(path, group_label, subject_sizes,
                                    true_exhaustive, efficiency_levels,
                                    efficiency_cases) {
  evaluation_results <- load_rdata_object(path, "evaluation_results")
  collected_results <- collect_evaluation_results(
    results = evaluation_results,
    group_label = group_label,
    subject_sizes = subject_sizes,
    true_exhaustive = true_exhaustive,
    efficiency_levels = efficiency_levels,
    efficiency_cases = efficiency_cases
  )
  rm(evaluation_results)
  gc()
  collected_results
}

collect_metric_file <- function(path, group_label, subject_sizes) {
  evaluation_results <- load_rdata_object(path, "evaluation_results")
  metric_rows <- list()
  metric_index <- 1L

  for (dataset_id in seq_along(evaluation_results)) {
    for (subject_index in seq_along(subject_sizes)) {
      current_subject_size <- subject_sizes[subject_index]
      simulations <- evaluation_results[[dataset_id]][[subject_index]]

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
      }
    }
  }

  rm(evaluation_results)
  gc()
  do.call(rbind, metric_rows)
}

normalize_metric <- function(values) {
  value_range <- range(values, na.rm = TRUE)
  if (diff(value_range) == 0) {
    return(rep(0, length(values)))
  }
  (values - value_range[1]) / diff(value_range)
}

# We use median scaling when computing ARE and RRMSE for the Mediterranean fruit fly case study.
summarize_metric_results_median_scaled <- function(metric_data,
                                                    w_are = 0.5,
                                                    w_rrmse = 0.5) {
  if (!"dataset" %in% names(metric_data)) {
    stop("metric_data must contain a dataset column for within-dataset scaling.")
  }
  if (!isTRUE(all.equal(w_are + w_rrmse, 1))) {
    stop("w_are and w_rrmse must sum to 1.")
  }

  metric_data |>
    dplyr::group_by(dataset) |>
    dplyr::mutate(
      are_median_dataset = median(are, na.rm = TRUE),
      rrmse_median_dataset = median(relative_rmse, na.rm = TRUE),
      are_median_scaled = are / are_median_dataset,
      rrmse_median_scaled = relative_rmse / rrmse_median_dataset,
      median_scaled_criterion =
        w_are * are_median_scaled + w_rrmse * rrmse_median_scaled
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(group, subject_size, simulation) |>
    dplyr::summarize(
      are_median_scaled = mean(are_median_scaled, na.rm = TRUE),
      rrmse_median_scaled = mean(rrmse_median_scaled, na.rm = TRUE),
      median_scaled_criterion = mean(median_scaled_criterion, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_efficiency_results <- function(efficiency_data) {
  efficiency_data |>
    dplyr::group_by(
      group,
      subject_size,
      simulation,
      efficiency_level,
      efficiency_case
    ) |>
    dplyr::summarize(
      efficiency = mean(true_efficiency, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_metric_results <- function(metric_data, include_composite = TRUE) {
  metric_summary <- metric_data |>
    dplyr::group_by(group, subject_size, simulation) |>
    dplyr::summarize(
      are_mean = mean(are, na.rm = TRUE),
      rrmse_mean = mean(relative_rmse, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      are_normalized = normalize_metric(are_mean),
      rrmse_normalized = normalize_metric(rrmse_mean)
    )

  if (include_composite) {
    metric_summary <- metric_summary |>
      dplyr::mutate(
        composite_criterion =
          0.5 * are_normalized + 0.5 * rrmse_normalized
      )
  }

  metric_summary
}

save_paper_boxplot <- function(data, value_column, y_label, group_levels,
                               group_colors, subject_sizes, output_file,
                               y_breaks, legend_position) {
  
  plot_data <- data
  
  plot_data$plot_value <- plot_data[[value_column]]
  
  plot_data$group_factor <- factor(
    plot_data$group,
    levels = group_levels
  )
  
  ensure_parent_directory(output_file)
  
  grDevices::cairo_ps(
    filename = output_file,
    width = 6.57,
    height = 3.93,
    onefile = FALSE,
    fallback_resolution = 300
  )
  
  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = factor(subject_size, levels = subject_sizes),
      y = plot_value,
      fill = group_factor
    )
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::xlab("Number of Subjects") +
    ggplot2::ylab(y_label) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::scale_fill_manual(
      name = "Designs",
      values = group_colors
    ) +
    
    ggplot2::theme_gray() +
    
    ggplot2::theme(
      
      text = ggplot2::element_text(
        family = "Times New Roman"
      ),
      
      axis.title.x = ggplot2::element_text(
        size = 14,
        family = "Times New Roman",
        margin = grid::unit(c(3, 0, 0, 0), "mm")
      ),
      
      axis.title.y = ggplot2::element_text(
        size = 14,
        family = "Times New Roman",
        margin = grid::unit(c(0, 3, 0, 0), "mm")
      ),
      
      axis.text = ggplot2::element_text(
        size = 12,
        family = "Times New Roman"
      ),
      
      legend.position = legend_position,
      
      legend.title = ggplot2::element_text(
        size = 14,
        family = "Times New Roman"
      ),
      
      legend.text = ggplot2::element_text(
        size = 14,
        family = "Times New Roman"
      ),
      
      legend.margin = ggplot2::margin(3, 10, 5, 5)
    )
  
  print(plot_object)
  
  dev.off()
  
  invisible(output_file)
}

save_efficiency_plots <- function(efficiency_summary, efficiency_levels,
                                  efficiency_cases, group_levels,
                                  group_colors, subject_sizes,
                                  output_directory, file_prefix) {
  for (efficiency_case in efficiency_cases) {
    y_label <- switch(
      efficiency_case,
      median = "Median-Case Efficiency",
      worst = "Worst-Case Efficiency"
    )

    for (efficiency_level in efficiency_levels) {
      level_percent <- round(100 * efficiency_level)
      plot_data <- efficiency_summary[
        efficiency_summary$efficiency_case == efficiency_case &
          abs(efficiency_summary$efficiency_level - efficiency_level) < 1e-10,
        ,
        drop = FALSE
      ]
      output_file <- file.path(
        output_directory,
        sprintf(
          "%s_%s_case_efficiency_%d.eps",
          file_prefix,
          efficiency_case,
          level_percent
        )
      )

      save_paper_boxplot(
        data = plot_data,
        value_column = "efficiency",
        y_label = y_label,
        group_levels = group_levels,
        group_colors = group_colors,
        subject_sizes = subject_sizes,
        output_file = output_file,
        y_breaks = seq(0.45, 1, by = 0.05),
        legend_position = c(0.89, 0.179)
      )
    }
  }
}

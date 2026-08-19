# ============================================================
# 18_compare_diversity_methods_final.R
#
# Purpose:
# Compare the Primary M3 framework with pooled and audio-only
# probability sensitivities and five deterministic workflows:
# naive retention, hard top-1, and global thresholds 0.50, 0.70,
# and 0.90.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()

calibrated_file <- file.path(
  paths$processed_root,
  "calibration",
  "calibrated_candidates.rds"
)
mc_file <- file.path(
  paths$results_root,
  "05_monte_carlo",
  "models",
  "monte_carlo_results.rds"
)

if (!file.exists(calibrated_file)) {
  stop("Calibrated candidate RDS not found: ", calibrated_file)
}
if (!file.exists(mc_file)) {
  stop("Monte Carlo result object not found: ", mc_file)
}

calibrated <- readRDS(calibrated_file) %>%
  dplyr::mutate(
    predicted_species = as.character(.data$predicted_species),
    year_month = as.character(.data$year_month),
    event_id = as.character(.data$event_id),
    id = as.character(.data$id)
  )
mc_results <- readRDS(mc_file)

month_levels <- mc_results$month_levels
species_levels <- mc_results$species_levels

method_metadata <- tibble::tribble(
  ~method, ~method_label, ~method_class, ~result_type,
  "Primary_M3", "Primary M3: species-aware", "Reference", "Monte Carlo",
  "Primary_M2", "Primary M2: pooled", "Probability sensitivity", "Monte Carlo",
  "Audio_only", "Audio-only sensitivity", "Probability sensitivity", "Monte Carlo",
  "Naive", "Naive: all candidates", "Deterministic", "Point estimate",
  "Hard_top1", "Hard top-1", "Deterministic", "Point estimate",
  "Threshold_0.50", "Threshold >= 0.50", "Deterministic", "Point estimate",
  "Threshold_0.70", "Threshold >= 0.70", "Deterministic", "Point estimate",
  "Threshold_0.90", "Threshold >= 0.90", "Deterministic", "Point estimate"
)

method_data <- list(
  Naive = calibrated,
  Hard_top1 = calibrated %>%
    dplyr::group_by(.data$event_id) %>%
    dplyr::arrange(
      dplyr::desc(.data$birdnet_score),
      .data$predicted_species,
      .data$id,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = 1L) %>%
    dplyr::ungroup()
)

for (threshold in config$thresholds) {
  threshold_name <- paste0("Threshold_", sprintf("%.2f", threshold))
  method_data[[threshold_name]] <- calibrated %>%
    dplyr::filter(.data$birdnet_score >= threshold)
}

deterministic_method_ids <- purrr::map(
  method_data,
  ~as.character(.x$id)
)

calculate_deterministic_method <- function(data, method_name) {
  count_matrix <- build_species_month_count_matrix(
    data,
    month_levels = month_levels,
    species_levels = species_levels
  )

  monthly <- monthly_diversity_from_matrix(
    count_matrix,
    month_levels,
    config$q_values
  ) %>%
    dplyr::mutate(
      method = .env$method_name,
      estimate = .data$diversity,
      lower_95 = .data$diversity,
      upper_95 = .data$diversity
    )

  partition <- purrr::map_dfr(
    c("equal_month", "abundance_weighted"),
    function(weighting) {
      purrr::map_dfr(
        config$q_values,
        function(q_value) {
          partition_hill_diversity(count_matrix, q_value, weighting) %>%
            tidyr::pivot_longer(
              cols = c(alpha, gamma, beta),
              names_to = "component",
              values_to = "estimate"
            ) %>%
            dplyr::mutate(
              method = .env$method_name,
              weighting = .env$weighting,
              q = .env$q_value,
              lower_95 = .data$estimate,
              upper_95 = .data$estimate,
              component_label = component_label(.data$component),
              q_label = q_label(.data$q)
            )
        }
      )
    }
  )

  quality <- tibble::tibble(
    method = method_name,
    retained_candidate_records = nrow(data),
    retained_events = dplyr::n_distinct(data$event_id),
    retained_species = dplyr::n_distinct(data$predicted_species),
    retained_months = dplyr::n_distinct(data$year_month)
  )

  list(monthly = monthly, partition = partition, quality = quality)
}

method_results <- purrr::imap(
  method_data,
  ~calculate_deterministic_method(.x, .y)
)

deterministic_monthly <- dplyr::bind_rows(
  purrr::map(method_results, "monthly")
)
deterministic_partition <- dplyr::bind_rows(
  purrr::map(method_results, "partition")
)
method_quality <- dplyr::bind_rows(
  purrr::map(method_results, "quality")
)

probabilistic_monthly <- mc_results$monthly_summary %>%
  dplyr::transmute(
    method = .data$scenario,
    year_month,
    q,
    q_label,
    estimate = .data$median,
    lower_95,
    upper_95,
    accepted_event_species = NA_real_
  )

probabilistic_partition <- mc_results$partition_summary %>%
  dplyr::transmute(
    method = .data$scenario,
    weighting,
    q,
    q_label,
    component,
    component_label,
    estimate = .data$median,
    lower_95,
    upper_95,
    n_months_used = NA_real_
  )

all_monthly <- dplyr::bind_rows(
  probabilistic_monthly,
  deterministic_monthly %>%
    dplyr::mutate(q_label = q_label(.data$q))
) %>%
  dplyr::left_join(method_metadata, by = "method") %>%
  dplyr::arrange(
    factor(.data$method, levels = method_metadata$method),
    .data$year_month,
    .data$q
  )

all_partition <- dplyr::bind_rows(
  probabilistic_partition,
  deterministic_partition
) %>%
  dplyr::left_join(method_metadata, by = "method") %>%
  dplyr::arrange(
    factor(.data$method, levels = method_metadata$method),
    .data$weighting,
    .data$q,
    factor(.data$component, levels = c("alpha", "gamma", "beta"))
  )

reference_monthly <- all_monthly %>%
  dplyr::filter(.data$method == "Primary_M3") %>%
  dplyr::select(
    year_month,
    q,
    reference_estimate = estimate,
    reference_lower_95 = lower_95,
    reference_upper_95 = upper_95
  )

monthly_difference <- all_monthly %>%
  dplyr::left_join(reference_monthly, by = c("year_month", "q")) %>%
  dplyr::mutate(
    absolute_difference = .data$estimate - .data$reference_estimate,
    percentage_difference = 100 * (
      .data$estimate / .data$reference_estimate - 1
    )
  )

reference_partition <- all_partition %>%
  dplyr::filter(.data$method == "Primary_M3") %>%
  dplyr::select(
    weighting,
    q,
    component,
    reference_estimate = estimate,
    reference_lower_95 = lower_95,
    reference_upper_95 = upper_95
  )

partition_difference <- all_partition %>%
  dplyr::left_join(
    reference_partition,
    by = c("weighting", "q", "component")
  ) %>%
  dplyr::mutate(
    absolute_difference = .data$estimate - .data$reference_estimate,
    percentage_difference = 100 * (
      .data$estimate / .data$reference_estimate - 1
    )
  )

monthly_bias_summary <- monthly_difference %>%
  dplyr::filter(.data$method != "Primary_M3") %>%
  dplyr::group_by(
    .data$method,
    .data$method_label,
    .data$method_class,
    .data$q,
    .data$q_label
  ) %>%
  dplyr::summarise(
    median_percentage_difference = stats::median(
      .data$percentage_difference,
      na.rm = TRUE
    ),
    lower_quartile_percentage_difference = stats::quantile(
      .data$percentage_difference,
      0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    upper_quartile_percentage_difference = stats::quantile(
      .data$percentage_difference,
      0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )

model_dir <- ensure_dir(
  silwood_result_dir(paths, "06_method_comparison", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "06_method_comparison", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "06_method_comparison", "figures")
)

readr::write_csv(
  all_monthly,
  file.path(table_dir, "method_comparison_monthly.csv"),
  na = ""
)
readr::write_csv(
  all_partition,
  file.path(table_dir, "method_comparison_partition.csv"),
  na = ""
)
readr::write_csv(
  monthly_difference,
  file.path(table_dir, "method_difference_from_m3_monthly.csv"),
  na = ""
)
readr::write_csv(
  partition_difference,
  file.path(table_dir, "method_difference_from_m3_partition.csv"),
  na = ""
)
readr::write_csv(
  monthly_bias_summary,
  file.path(table_dir, "monthly_bias_summary.csv"),
  na = ""
)
readr::write_csv(
  method_quality,
  file.path(table_dir, "deterministic_method_quality.csv")
)

saveRDS(
  list(
    method_metadata = method_metadata,
    monthly = all_monthly,
    partition = all_partition,
    monthly_difference = monthly_difference,
    partition_difference = partition_difference,
    method_quality = method_quality,
    deterministic_method_ids = deterministic_method_ids
  ),
  file.path(model_dir, "method_comparison_results.rds"),
  compress = "xz"
)

manuscript_table_dir <- ensure_dir(
  silwood_result_dir(paths, "06_method_comparison", "tables")
)
method_comparison_publication <- partition_difference %>%
  dplyr::filter(
    .data$weighting == "equal_month",
    .data$method != "Primary_M3"
  ) %>%
  dplyr::transmute(
    Workflow = .data$method_label,
    Workflow_class = .data$method_class,
    Diversity_order = .data$q_label,
    Component = .data$component_label,
    Estimate = round(.data$estimate, 3),
    Primary_M3_reference = round(.data$reference_estimate, 3),
    Absolute_difference = round(.data$absolute_difference, 3),
    Percentage_difference = round(.data$percentage_difference, 2)
  )
readr::write_csv(
  method_comparison_publication,
  file.path(manuscript_table_dir, "Table_method_comparison_from_M3.csv"),
  na = ""
)

heatmap_data <- partition_difference %>%
  dplyr::filter(
    .data$weighting == "equal_month",
    .data$method != "Primary_M3"
  ) %>%
  dplyr::mutate(
    method_label = factor(
      .data$method_label,
      levels = rev(method_metadata$method_label[-1])
    ),
    component_label = factor(
      .data$component_label,
      levels = c("Alpha", "Gamma", "Beta")
    ),
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    ),
    display_label = sprintf("%+.1f%%", .data$percentage_difference),
    label_colour = dplyr::if_else(
      abs(.data$percentage_difference) >= 35,
      "white",
      "black"
    )
  )

heatmap_limit <- max(abs(heatmap_data$percentage_difference), na.rm = TRUE)

bias_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = .data$component_label,
    y = .data$method_label,
    fill = .data$percentage_difference
  )
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = .data$display_label,
      colour = .data$label_colour
    ),
    size = 3.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~q_label, nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#2C5AA0",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-heatmap_limit, heatmap_limit),
    oob = scales::squish,
    labels = scales::label_number(suffix = "%", accuracy = 1)
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::labs(
    title = "Percentage difference from Primary M3",
    subtitle = "Equal-month alpha, gamma and multiplicative beta diversity",
    x = "Partition component",
    y = NULL,
    fill = "Difference\nfrom M3",
    caption = "Blue values are lower than Primary M3; red values are higher."
  ) +
  theme_silwood(10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  )

save_figure(
  bias_heatmap,
  file.path(figure_dir, "figure_07_method_bias_heatmap"),
  width = 13,
  height = 6.8
)

forest_data <- all_partition %>%
  dplyr::filter(.data$weighting == "equal_month") %>%
  dplyr::mutate(
    method_label = factor(
      .data$method_label,
      levels = rev(method_metadata$method_label)
    ),
    component_label = factor(
      .data$component_label,
      levels = c("Alpha", "Gamma", "Beta")
    ),
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    )
  )

partition_forest <- ggplot2::ggplot(
  forest_data,
  ggplot2::aes(
    x = .data$estimate,
    y = .data$method_label,
    colour = .data$method_label
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$lower_95, xmax = .data$upper_95),
    orientation = "y",
    width = 0,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(size = 2.7) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$component_label),
    cols = ggplot2::vars(.data$q_label),
    scales = "free_x"
  ) +
  ggplot2::scale_colour_manual(values = palette$method) +
  ggplot2::labs(
    title = "Diversity estimates across identification workflows",
    subtitle = "Primary M3 intervals reflect full Monte Carlo uncertainty; deterministic methods are point estimates",
    x = "Effective diversity",
    y = NULL,
    colour = "Workflow"
  ) +
  theme_silwood(8.5) +
  ggplot2::theme(
    legend.position = "none",
    strip.text.y = ggplot2::element_text(angle = 0)
  )

save_figure(
  partition_forest,
  file.path(figure_dir, "figure_08_method_comparison_forest"),
  width = 14,
  height = 10
)

monthly_plot_data <- all_monthly %>%
  dplyr::filter(
    .data$method %in% c(
      "Primary_M3",
      "Naive",
      "Hard_top1",
      "Threshold_0.50",
      "Threshold_0.70",
      "Threshold_0.90"
    )
  ) %>%
  dplyr::mutate(
    month_date = as.Date(paste0(.data$year_month, "-01")),
    method_label = factor(
      .data$method_label,
      levels = method_metadata$method_label
    ),
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    )
  )

monthly_methods_figure <- ggplot2::ggplot(
  monthly_plot_data,
  ggplot2::aes(
    x = .data$month_date,
    y = .data$estimate,
    colour = .data$method_label,
    group = .data$method_label
  )
) +
  ggplot2::geom_ribbon(
    data = monthly_plot_data %>%
      dplyr::filter(.data$method == "Primary_M3"),
    ggplot2::aes(
      ymin = .data$lower_95,
      ymax = .data$upper_95,
      fill = .data$method_label
    ),
    inherit.aes = TRUE,
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::facet_wrap(~q_label, scales = "free_y", ncol = 1) +
  ggplot2::scale_colour_manual(values = palette$method) +
  ggplot2::scale_fill_manual(values = palette$method) +
  ggplot2::scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
  ggplot2::labs(
    title = "Monthly diversity under deterministic workflows",
    subtitle = "Primary M3 is shown with its full Monte Carlo interval",
    x = NULL,
    y = "Effective species number",
    colour = "Workflow",
    fill = "Workflow"
  ) +
  theme_silwood(9)

save_figure(
  monthly_methods_figure,
  file.path(figure_dir, "figure_S03_monthly_method_comparison"),
  width = 12.5,
  height = 10.5
)

message("Method comparison completed.")
message("Tables: ", table_dir)
message("Figures: ", figure_dir)

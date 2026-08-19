# ============================================================
# 20_low_coverage_month_sensitivity.R
#
# Purpose:
# Recalculate Primary M3 and deterministic workflow partitions
# after excluding low represented-coverage months.
#
# The script reuses the exact Primary M3 species-by-month count
# arrays generated in script 17. It does not rerun Monte Carlo
# classification draws, so all-month and higher-coverage results
# are paired at the simulation level.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()

mc_file <- file.path(
  paths$results_root,
  "05_monte_carlo",
  "models",
  "monte_carlo_results.rds"
)
method_file <- file.path(
  paths$results_root,
  "06_method_comparison",
  "models",
  "method_comparison_results.rds"
)
coverage_file <- file.path(
  paths$results_root,
  "07_sampling_coverage",
  "models",
  "sampling_coverage_results.rds"
)
calibrated_file <- file.path(
  paths$processed_root,
  "calibration",
  "calibrated_candidates.rds"
)

for (required_file in c(mc_file, method_file, coverage_file, calibrated_file)) {
  if (!file.exists(required_file)) {
    stop("Required result object not found: ", required_file)
  }
}

mc_results <- readRDS(mc_file)
method_results <- readRDS(method_file)
coverage_results <- readRDS(coverage_file)
calibrated <- readRDS(calibrated_file) %>%
  dplyr::mutate(
    id = as.character(.data$id),
    year_month = as.character(.data$year_month),
    predicted_species = as.character(.data$predicted_species)
  )

coverage <- coverage_results$coverage
higher_months <- coverage %>%
  dplyr::filter(!.data$low_coverage) %>%
  dplyr::pull(year_month)
low_months <- coverage %>%
  dplyr::filter(.data$low_coverage) %>%
  dplyr::pull(year_month)

month_indices <- match(higher_months, mc_results$month_levels)
if (any(is.na(month_indices))) {
  stop("At least one higher-coverage month is absent from the Monte Carlo array.")
}

primary_count_array <- mc_results$primary_count_array
n_simulations <- dim(primary_count_array)[1]

high_partition_draws <- vector(
  "list",
  n_simulations * length(config$q_values)
)
output_index <- 1L

for (simulation_id in seq_len(n_simulations)) {
  high_count_matrix <- primary_count_array[
    simulation_id,
    month_indices,
    ,
    drop = FALSE
  ]
  high_count_matrix <- matrix(
    high_count_matrix,
    nrow = length(month_indices),
    ncol = length(mc_results$species_levels),
    dimnames = list(higher_months, mc_results$species_levels)
  )

  for (q_value in config$q_values) {
    high_partition_draws[[output_index]] <- partition_hill_diversity(
      high_count_matrix,
      q_value,
      "equal_month"
    ) %>%
      tidyr::pivot_longer(
        cols = c(alpha, gamma, beta),
        names_to = "component",
        values_to = "diversity"
      ) %>%
      dplyr::mutate(
        simulation_id = .env$simulation_id,
        q = .env$q_value,
        q_label = q_label(.env$q_value),
        component_label = component_label(.data$component),
        month_set = "Higher represented coverage"
      )
    output_index <- output_index + 1L
  }
}

high_partition_draws <- dplyr::bind_rows(high_partition_draws)
high_m3_summary <- high_partition_draws %>%
  dplyr::group_by(
    .data$q,
    .data$q_label,
    .data$component,
    .data$component_label
  ) %>%
  dplyr::group_modify(~summarise_distribution(.x$diversity)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    method = "Primary_M3",
    method_label = "Primary M3: species-aware",
    weighting = "equal_month",
    estimate = .data$median
  )

high_deterministic_results <- purrr::imap_dfr(
  method_results$deterministic_method_ids,
  function(method_ids, method_name) {
    high_data <- calibrated %>%
      dplyr::filter(
        .data$id %in% method_ids,
        .data$year_month %in% higher_months
      )

    count_matrix <- build_species_month_count_matrix(
      high_data,
      month_levels = higher_months,
      species_levels = mc_results$species_levels
    )

    purrr::map_dfr(
      config$q_values,
      function(q_value) {
        partition_hill_diversity(
          count_matrix,
          q_value,
          "equal_month"
        ) %>%
          tidyr::pivot_longer(
            cols = c(alpha, gamma, beta),
            names_to = "component",
            values_to = "estimate"
          ) %>%
          dplyr::mutate(
            method = .env$method_name,
            weighting = "equal_month",
            q = .env$q_value,
            q_label = q_label(.env$q_value),
            component_label = component_label(.data$component),
            lower_95 = .data$estimate,
            upper_95 = .data$estimate
          )
      }
    )
  }
)

method_metadata <- method_results$method_metadata
high_partition <- dplyr::bind_rows(
  high_m3_summary %>%
    dplyr::select(
      method,
      method_label,
      weighting,
      q,
      q_label,
      component,
      component_label,
      estimate,
      lower_95,
      upper_95
    ),
  high_deterministic_results %>%
    dplyr::left_join(
      method_metadata %>% dplyr::select(method, method_label),
      by = "method"
    ) %>%
    dplyr::select(
      method,
      method_label,
      weighting,
      q,
      q_label,
      component,
      component_label,
      estimate,
      lower_95,
      upper_95
    )
)

high_reference <- high_partition %>%
  dplyr::filter(.data$method == "Primary_M3") %>%
  dplyr::select(
    q,
    component,
    high_reference = estimate,
    high_reference_lower = lower_95,
    high_reference_upper = upper_95
  )

high_difference <- high_partition %>%
  dplyr::left_join(high_reference, by = c("q", "component")) %>%
  dplyr::mutate(
    absolute_difference = .data$estimate - .data$high_reference,
    percentage_difference = 100 * (
      .data$estimate / .data$high_reference - 1
    ),
    month_set = "Higher represented coverage"
  )

all_difference <- method_results$partition_difference %>%
  dplyr::filter(
    .data$weighting == "equal_month",
    .data$method %in% names(method_results$deterministic_method_ids)
  ) %>%
  dplyr::select(
    method,
    method_label,
    q,
    q_label,
    component,
    component_label,
    all_estimate = estimate,
    all_reference = reference_estimate,
    all_absolute_difference = absolute_difference,
    all_percentage_difference = percentage_difference
  )

robustness <- all_difference %>%
  dplyr::left_join(
    high_difference %>%
      dplyr::filter(.data$method != "Primary_M3") %>%
      dplyr::select(
        method,
        q,
        component,
        high_estimate = estimate,
        high_reference,
        high_absolute_difference = absolute_difference,
        high_percentage_difference = percentage_difference
      ),
    by = c("method", "q", "component")
  ) %>%
  dplyr::mutate(
    all_direction = sign(.data$all_percentage_difference),
    high_direction = sign(.data$high_percentage_difference),
    direction_stable = .data$all_direction == .data$high_direction,
    absolute_bias_change_points = abs(
      .data$high_percentage_difference - .data$all_percentage_difference
    )
  )

robustness_summary <- robustness %>%
  dplyr::group_by(.data$method, .data$method_label) %>%
  dplyr::summarise(
    comparisons = dplyr::n(),
    direction_stable_n = sum(.data$direction_stable, na.rm = TRUE),
    direction_stable_proportion = mean(.data$direction_stable, na.rm = TRUE),
    mean_absolute_bias_change_points = mean(
      .data$absolute_bias_change_points,
      na.rm = TRUE
    ),
    maximum_absolute_bias_change_points = max(
      .data$absolute_bias_change_points,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

all_primary <- mc_results$partition_summary %>%
  dplyr::filter(
    .data$scenario == "Primary_M3",
    .data$weighting == "equal_month"
  ) %>%
  dplyr::transmute(
    q,
    q_label,
    component,
    component_label,
    all_median = median,
    all_lower_95 = lower_95,
    all_upper_95 = upper_95
  )

primary_all_vs_high <- all_primary %>%
  dplyr::left_join(
    high_m3_summary %>%
      dplyr::select(
        q,
        component,
        high_median = median,
        high_lower_95 = lower_95,
        high_upper_95 = upper_95
      ),
    by = c("q", "component")
  ) %>%
  dplyr::mutate(
    absolute_change = .data$high_median - .data$all_median,
    percentage_change = 100 * (
      .data$high_median / .data$all_median - 1
    )
  )

quality_checks <- tibble::tibble(
  check = c(
    "All represented months",
    "Low-coverage months",
    "Higher-coverage months",
    "Primary simulations reused",
    "High-coverage Primary partition rows",
    "Deterministic method comparisons"
  ),
  value = c(
    length(mc_results$month_levels),
    length(low_months),
    length(higher_months),
    n_simulations,
    nrow(high_m3_summary),
    nrow(robustness)
  )
)

assert_expected(length(mc_results$month_levels), 23, "All represented months")
assert_expected(length(low_months), 9, "Low-coverage months")
assert_expected(length(higher_months), 14, "Higher-coverage months")

model_dir <- ensure_dir(
  silwood_result_dir(paths, "08_coverage_robustness", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "08_coverage_robustness", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "08_coverage_robustness", "figures")
)

readr::write_csv(
  high_partition,
  file.path(table_dir, "higher_coverage_partition.csv"),
  na = ""
)
readr::write_csv(
  high_difference,
  file.path(table_dir, "higher_coverage_difference_from_m3.csv"),
  na = ""
)
readr::write_csv(
  robustness,
  file.path(table_dir, "method_bias_robustness.csv"),
  na = ""
)
readr::write_csv(
  robustness_summary,
  file.path(table_dir, "method_bias_robustness_summary.csv"),
  na = ""
)
readr::write_csv(
  primary_all_vs_high,
  file.path(table_dir, "primary_m3_all_vs_higher_coverage.csv"),
  na = ""
)
readr::write_csv(
  quality_checks,
  file.path(table_dir, "coverage_robustness_quality_checks.csv")
)

saveRDS(
  list(
    high_partition_draws = high_partition_draws,
    high_partition = high_partition,
    high_difference = high_difference,
    robustness = robustness,
    robustness_summary = robustness_summary,
    primary_all_vs_high = primary_all_vs_high,
    low_months = low_months,
    higher_months = higher_months
  ),
  file.path(model_dir, "coverage_robustness_results.rds"),
  compress = "xz"
)

manuscript_table_dir <- ensure_dir(
  silwood_result_dir(paths, "08_coverage_robustness", "tables")
)
robustness_publication <- robustness_summary %>%
  dplyr::transmute(
    Workflow = .data$method_label,
    Comparisons = .data$comparisons,
    Stable_direction_n = .data$direction_stable_n,
    Stable_direction_proportion = round(
      .data$direction_stable_proportion,
      3
    ),
    Mean_absolute_bias_change_points = round(
      .data$mean_absolute_bias_change_points,
      2
    ),
    Maximum_absolute_bias_change_points = round(
      .data$maximum_absolute_bias_change_points,
      2
    )
  )
readr::write_csv(
  robustness_publication,
  file.path(manuscript_table_dir, "Table_low_coverage_robustness.csv"),
  na = ""
)

robustness_plot_data <- robustness %>%
  dplyr::select(
    method,
    method_label,
    q,
    q_label,
    component,
    component_label,
    all_percentage_difference,
    high_percentage_difference
  ) %>%
  tidyr::pivot_longer(
    cols = c(all_percentage_difference, high_percentage_difference),
    names_to = "month_set",
    values_to = "percentage_difference"
  ) %>%
  dplyr::mutate(
    month_set = dplyr::recode(
      .data$month_set,
      all_percentage_difference = "All represented months",
      high_percentage_difference = "Higher-coverage months"
    ),
    method_label = factor(
      .data$method_label,
      levels = rev(unique(method_results$method_metadata$method_label[
        method_results$method_metadata$method %in%
          names(method_results$deterministic_method_ids)
      ]))
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

robustness_figure <- ggplot2::ggplot(
  robustness_plot_data,
  ggplot2::aes(
    x = .data$percentage_difference,
    y = .data$method_label
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_line(
    ggplot2::aes(
      group = interaction(
        .data$method_label,
        .data$q_label,
        .data$component_label
      )
    ),
    colour = "grey65",
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = .data$month_set),
    size = 2.6
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$component_label),
    cols = ggplot2::vars(.data$q_label),
    scales = "free_x"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "All represented months" = "#0072B2",
      "Higher-coverage months" = "#D55E00"
    )
  ) +
  ggplot2::labs(
    title = "Robustness of deterministic-method bias to low-coverage months",
    subtitle = "Points show percentage differences from the corresponding Primary M3 reference",
    x = "Difference from Primary M3 (%)",
    y = NULL,
    colour = "Month set"
  ) +
  theme_silwood(8.5) +
  ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))

save_figure(
  robustness_figure,
  file.path(figure_dir, "figure_11_method_bias_coverage_robustness"),
  width = 14,
  height = 9.5
)

primary_plot_data <- primary_all_vs_high %>%
  tidyr::pivot_longer(
    cols = c(all_median, high_median),
    names_to = "month_set",
    values_to = "median"
  ) %>%
  dplyr::mutate(
    lower_95 = dplyr::if_else(
      .data$month_set == "all_median",
      .data$all_lower_95,
      .data$high_lower_95
    ),
    upper_95 = dplyr::if_else(
      .data$month_set == "all_median",
      .data$all_upper_95,
      .data$high_upper_95
    ),
    month_set = dplyr::recode(
      .data$month_set,
      all_median = "All represented months",
      high_median = "Higher-coverage months"
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

primary_figure <- ggplot2::ggplot(
  primary_plot_data,
  ggplot2::aes(
    x = .data$median,
    y = .data$month_set,
    colour = .data$month_set
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$lower_95, xmax = .data$upper_95),
    orientation = "y",
    width = 0,
    linewidth = 0.7
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$component_label),
    cols = ggplot2::vars(.data$q_label),
    scales = "free_x"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "All represented months" = "#0072B2",
      "Higher-coverage months" = "#D55E00"
    )
  ) +
  ggplot2::labs(
    title = "Primary M3 partitions before and after excluding low-coverage months",
    x = "Effective diversity",
    y = NULL,
    colour = "Month set"
  ) +
  theme_silwood(8.5) +
  ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))

save_figure(
  primary_figure,
  file.path(figure_dir, "figure_S05_primary_m3_all_vs_higher_coverage"),
  width = 13,
  height = 8.5
)

message("Low-coverage robustness analysis completed.")
message("Higher-coverage months: ", length(higher_months))

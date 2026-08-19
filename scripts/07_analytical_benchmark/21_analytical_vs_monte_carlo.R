# ============================================================
# 21_analytical_vs_monte_carlo.R
#
# Output 2: analytical benchmark for Primary M3.
#
# The script calculates:
#   - exact expected candidate-record counts;
#   - expected richness under conditional independence;
#   - plug-in q = 1 and q = 2 Hill diversity;
#   - analytical / plug-in temporal partitions;
#   - fixed-probability Monte Carlo;
#   - comparisons with full Monte Carlo from script 17.
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
full_mc_file <- file.path(
  paths$results_root,
  "05_monte_carlo",
  "models",
  "monte_carlo_results.rds"
)

if (!file.exists(calibrated_file)) {
  stop("Calibrated candidate RDS not found: ", calibrated_file)
}
if (!file.exists(full_mc_file)) {
  stop("Full Monte Carlo result object not found: ", full_mc_file)
}

calibrated <- readRDS(calibrated_file) %>%
  dplyr::mutate(
    year_month = as.character(.data$year_month),
    predicted_species = as.character(.data$predicted_species),
    p_primary = clip_probability(.data$p_primary)
  ) %>%
  dplyr::arrange(.data$year_month, .data$predicted_species, .data$id)
full_mc <- readRDS(full_mc_file)

month_levels <- full_mc$month_levels
species_levels <- full_mc$species_levels
aggregation <- build_candidate_aggregation(calibrated)

probability_at_least_one <- function(probability) {
  probability <- clip_probability(probability)
  1 - exp(sum(log1p(-probability)))
}

species_month_expectations <- calibrated %>%
  dplyr::group_by(.data$year_month, .data$predicted_species) %>%
  dplyr::summarise(
    candidate_n = dplyr::n(),
    expected_correct_candidates = sum(.data$p_primary),
    probability_at_least_one_correct = probability_at_least_one(
      .data$p_primary
    ),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    year_month = month_levels,
    predicted_species = species_levels,
    fill = list(
      candidate_n = 0L,
      expected_correct_candidates = 0,
      probability_at_least_one_correct = 0
    )
  ) %>%
  dplyr::mutate(
    year_month = factor(.data$year_month, levels = month_levels),
    predicted_species = factor(
      .data$predicted_species,
      levels = species_levels
    )
  ) %>%
  dplyr::arrange(.data$year_month, .data$predicted_species)

expected_count_matrix <- matrix(
  species_month_expectations$expected_correct_candidates,
  nrow = length(month_levels),
  ncol = length(species_levels),
  byrow = TRUE,
  dimnames = list(month_levels, species_levels)
)

presence_probability_matrix <- matrix(
  species_month_expectations$probability_at_least_one_correct,
  nrow = length(month_levels),
  ncol = length(species_levels),
  byrow = TRUE,
  dimnames = list(month_levels, species_levels)
)

expected_monthly_richness <- rowSums(presence_probability_matrix)
expected_monthly_total <- rowSums(expected_count_matrix)

monthly_plugin <- purrr::map_dfr(
  seq_along(month_levels),
  function(month_index) {
    counts <- expected_count_matrix[month_index, ]
    tibble::tibble(
      year_month = month_levels[[month_index]],
      q = config$q_values,
      analytical_estimate = c(
        expected_monthly_richness[[month_index]],
        hill_number_from_counts(counts, 1),
        hill_number_from_counts(counts, 2)
      ),
      analytical_type = c(
        "Expected richness",
        "Plug-in Hill number",
        "Plug-in Hill number"
      )
    )
  }
) %>%
  dplyr::mutate(q_label = q_label(.data$q))

species_gamma_probability <- calibrated %>%
  dplyr::group_by(.data$predicted_species) %>%
  dplyr::summarise(
    probability_at_least_one_correct = probability_at_least_one(
      .data$p_primary
    ),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    predicted_species = species_levels,
    fill = list(probability_at_least_one_correct = 0)
  )

expected_gamma_richness <- sum(
  species_gamma_probability$probability_at_least_one_correct
)

calculate_analytical_partition <- function(weighting) {
  weighting_value <- as.character(weighting)
  if (weighting_value == "equal_month") {
    month_weights <- rep(1 / length(month_levels), length(month_levels))
  } else if (weighting_value == "abundance_weighted") {
    month_weights <- expected_monthly_total / sum(expected_monthly_total)
  } else {
    stop("Unknown analytical weighting: ", weighting_value)
  }

  q0_alpha <- sum(month_weights * expected_monthly_richness)
  q0_gamma <- expected_gamma_richness
  q0_beta <- q0_gamma / q0_alpha

  q0_result <- tibble::tibble(
    weighting = weighting_value,
    q = 0,
    component = c("alpha", "gamma", "beta"),
    analytical_estimate = c(q0_alpha, q0_gamma, q0_beta),
    analytical_type = c(
      "Expected-richness alpha component",
      "Expected-richness gamma component",
      "Ratio of expected-richness components"
    )
  )

  q12_result <- purrr::map_dfr(
    c(1, 2),
    function(q_value) {
      partition_hill_diversity(
        expected_count_matrix,
        q_value,
        weighting_value
      ) %>%
        tidyr::pivot_longer(
          cols = c(alpha, gamma, beta),
          names_to = "component",
          values_to = "analytical_estimate"
        ) %>%
        dplyr::mutate(
          weighting = .env$weighting_value,
          q = .env$q_value,
          analytical_type = "Plug-in Hill partition"
        )
    }
  )

  dplyr::bind_rows(q0_result, q12_result) %>%
    dplyr::mutate(
      q_label = q_label(.data$q),
      component_label = component_label(.data$component)
    )
}

analytical_partition <- dplyr::bind_rows(
  calculate_analytical_partition("equal_month"),
  calculate_analytical_partition("abundance_weighted")
)

set.seed(config$seed + 100000L)
fixed_monthly_draws <- vector(
  "list",
  config$n_fixed_probability_simulations
)
fixed_partition_draws <- vector(
  "list",
  config$n_fixed_probability_simulations
)

message(
  "Running fixed-probability Monte Carlo: ",
  config$n_fixed_probability_simulations,
  " simulations."
)

for (simulation_id in seq_len(config$n_fixed_probability_simulations)) {
  if (
    simulation_id == 1L ||
      simulation_id %% max(
        1L,
        floor(config$n_fixed_probability_simulations / 10L)
      ) == 0L
  ) {
    message(
      "  fixed-probability simulation ",
      simulation_id,
      "/",
      config$n_fixed_probability_simulations
    )
  }

  candidate_state <- stats::rbinom(
    nrow(calibrated),
    size = 1L,
    prob = calibrated$p_primary
  )
  count_matrix <- aggregate_candidate_states(candidate_state, aggregation)

  fixed_monthly_draws[[simulation_id]] <- monthly_diversity_from_matrix(
    count_matrix,
    month_levels,
    config$q_values
  ) %>%
    dplyr::mutate(simulation_id = .env$simulation_id)

  fixed_partition_draws[[simulation_id]] <- purrr::map_dfr(
    c("equal_month", "abundance_weighted"),
    function(weighting) {
      purrr::map_dfr(
        config$q_values,
        function(q_value) {
          partition_hill_diversity(
            count_matrix,
            q_value,
            weighting
          ) %>%
            tidyr::pivot_longer(
              cols = c(alpha, gamma, beta),
              names_to = "component",
              values_to = "diversity"
            ) %>%
            dplyr::mutate(
              simulation_id = .env$simulation_id,
              weighting = .env$weighting,
              q = .env$q_value,
              q_label = q_label(.env$q_value),
              component_label = component_label(.data$component)
            )
        }
      )
    }
  )
}

fixed_monthly_draws <- dplyr::bind_rows(fixed_monthly_draws)
fixed_partition_draws <- dplyr::bind_rows(fixed_partition_draws)

fixed_monthly_summary <- fixed_monthly_draws %>%
  dplyr::group_by(.data$year_month, .data$q) %>%
  dplyr::group_modify(~summarise_distribution(.x$diversity)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q_label = q_label(.data$q))

fixed_partition_summary <- fixed_partition_draws %>%
  dplyr::group_by(
    .data$weighting,
    .data$q,
    .data$q_label,
    .data$component,
    .data$component_label
  ) %>%
  dplyr::group_modify(~summarise_distribution(.x$diversity)) %>%
  dplyr::ungroup()

full_monthly_summary <- full_mc$monthly_summary %>%
  dplyr::filter(.data$scenario == "Primary_M3") %>%
  dplyr::select(
    year_month,
    q,
    full_mean = mean,
    full_median = median,
    full_lower_95 = lower_95,
    full_upper_95 = upper_95
  )

full_partition_summary <- full_mc$partition_summary %>%
  dplyr::filter(.data$scenario == "Primary_M3") %>%
  dplyr::select(
    weighting,
    q,
    component,
    full_mean = mean,
    full_median = median,
    full_lower_95 = lower_95,
    full_upper_95 = upper_95
  )

monthly_comparison <- monthly_plugin %>%
  dplyr::left_join(
    fixed_monthly_summary %>%
      dplyr::select(
        year_month,
        q,
        fixed_mean = mean,
        fixed_median = median,
        fixed_lower_95 = lower_95,
        fixed_upper_95 = upper_95
      ),
    by = c("year_month", "q")
  ) %>%
  dplyr::left_join(full_monthly_summary, by = c("year_month", "q")) %>%
  dplyr::mutate(
    analytical_minus_fixed_mean = .data$analytical_estimate - .data$fixed_mean,
    analytical_percent_error = 100 * (
      .data$analytical_estimate / .data$fixed_mean - 1
    ),
    analytical_inside_fixed_95 =
      .data$analytical_estimate >= .data$fixed_lower_95 &
      .data$analytical_estimate <= .data$fixed_upper_95,
    fixed_minus_full_mean = .data$fixed_mean - .data$full_mean,
    fixed_percent_difference_from_full = 100 * (
      .data$fixed_mean / .data$full_mean - 1
    )
  )

partition_comparison <- analytical_partition %>%
  dplyr::left_join(
    fixed_partition_summary %>%
      dplyr::select(
        weighting,
        q,
        component,
        fixed_mean = mean,
        fixed_median = median,
        fixed_lower_95 = lower_95,
        fixed_upper_95 = upper_95
      ),
    by = c("weighting", "q", "component")
  ) %>%
  dplyr::left_join(
    full_partition_summary,
    by = c("weighting", "q", "component")
  ) %>%
  dplyr::mutate(
    analytical_minus_fixed_mean = .data$analytical_estimate - .data$fixed_mean,
    analytical_percent_error = 100 * (
      .data$analytical_estimate / .data$fixed_mean - 1
    ),
    analytical_inside_fixed_95 =
      .data$analytical_estimate >= .data$fixed_lower_95 &
      .data$analytical_estimate <= .data$fixed_upper_95,
    fixed_minus_full_mean = .data$fixed_mean - .data$full_mean,
    fixed_percent_difference_from_full = 100 * (
      .data$fixed_mean / .data$full_mean - 1
    )
  )

approximation_error_summary <- dplyr::bind_rows(
  monthly_comparison %>%
    dplyr::group_by(.data$q, .data$q_label) %>%
    dplyr::summarise(
      comparison_level = "Monthly",
      weighting = NA_character_,
      n_comparisons = dplyr::n(),
      mean_signed_percent_error = mean(
        .data$analytical_percent_error,
        na.rm = TRUE
      ),
      mean_absolute_percent_error = mean(
        abs(.data$analytical_percent_error),
        na.rm = TRUE
      ),
      maximum_absolute_percent_error = max(
        abs(.data$analytical_percent_error),
        na.rm = TRUE
      ),
      proportion_inside_fixed_mc_95 = mean(
        .data$analytical_inside_fixed_95,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),
  partition_comparison %>%
    dplyr::group_by(.data$weighting, .data$q, .data$q_label) %>%
    dplyr::summarise(
      comparison_level = "Partition",
      n_comparisons = dplyr::n(),
      mean_signed_percent_error = mean(
        .data$analytical_percent_error,
        na.rm = TRUE
      ),
      mean_absolute_percent_error = mean(
        abs(.data$analytical_percent_error),
        na.rm = TRUE
      ),
      maximum_absolute_percent_error = max(
        abs(.data$analytical_percent_error),
        na.rm = TRUE
      ),
      proportion_inside_fixed_mc_95 = mean(
        .data$analytical_inside_fixed_95,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)

quality_checks <- tibble::tibble(
  check = c(
    "Species-month expectation rows",
    "Monthly analytical rows",
    "Analytical partition rows",
    "Fixed-probability simulations",
    "Missing fixed monthly values",
    "Missing fixed partition values",
    "Monthly comparison rows",
    "Partition comparison rows"
  ),
  value = c(
    nrow(species_month_expectations),
    nrow(monthly_plugin),
    nrow(analytical_partition),
    dplyr::n_distinct(fixed_monthly_draws$simulation_id),
    sum(!is.finite(fixed_monthly_draws$diversity)),
    sum(!is.finite(fixed_partition_draws$diversity)),
    nrow(monthly_comparison),
    nrow(partition_comparison)
  ),
  expected = c(
    length(month_levels) * length(species_levels),
    length(month_levels) * length(config$q_values),
    2L * length(config$q_values) * 3L,
    config$n_fixed_probability_simulations,
    0,
    0,
    length(month_levels) * length(config$q_values),
    2L * length(config$q_values) * 3L
  )
)

for (row_index in seq_len(nrow(quality_checks))) {
  assert_expected(
    quality_checks$value[[row_index]],
    quality_checks$expected[[row_index]],
    quality_checks$check[[row_index]]
  )
}

model_dir <- ensure_dir(
  silwood_result_dir(paths, "09_analytical_benchmark", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "09_analytical_benchmark", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "09_analytical_benchmark", "figures")
)

readr::write_csv(
  species_month_expectations,
  file.path(table_dir, "analytical_species_month_expectations.csv"),
  na = ""
)
readr::write_csv(
  monthly_plugin,
  file.path(table_dir, "analytical_monthly_diversity.csv"),
  na = ""
)
readr::write_csv(
  analytical_partition,
  file.path(table_dir, "analytical_partition.csv"),
  na = ""
)
readr::write_csv(
  fixed_monthly_summary,
  file.path(table_dir, "fixed_probability_mc_monthly_summary.csv"),
  na = ""
)
readr::write_csv(
  fixed_partition_summary,
  file.path(table_dir, "fixed_probability_mc_partition_summary.csv"),
  na = ""
)
readr::write_csv(
  monthly_comparison,
  file.path(table_dir, "analytical_vs_monte_carlo_monthly.csv"),
  na = ""
)
readr::write_csv(
  partition_comparison,
  file.path(table_dir, "analytical_vs_monte_carlo_partition.csv"),
  na = ""
)
readr::write_csv(
  approximation_error_summary,
  file.path(table_dir, "analytical_approximation_error_summary.csv"),
  na = ""
)
readr::write_csv(
  quality_checks,
  file.path(table_dir, "analytical_benchmark_quality_checks.csv")
)

saveRDS(
  list(
    species_month_expectations = species_month_expectations,
    expected_count_matrix = expected_count_matrix,
    presence_probability_matrix = presence_probability_matrix,
    monthly_analytical = monthly_plugin,
    partition_analytical = analytical_partition,
    fixed_monthly_draws = fixed_monthly_draws,
    fixed_partition_draws = fixed_partition_draws,
    fixed_monthly_summary = fixed_monthly_summary,
    fixed_partition_summary = fixed_partition_summary,
    monthly_comparison = monthly_comparison,
    partition_comparison = partition_comparison,
    approximation_error_summary = approximation_error_summary,
    config = config
  ),
  file.path(model_dir, "analytical_benchmark_results.rds"),
  compress = "xz"
)

manuscript_table_dir <- ensure_dir(
  silwood_result_dir(paths, "09_analytical_benchmark", "tables")
)
analytical_publication <- partition_comparison %>%
  dplyr::filter(.data$weighting == "equal_month") %>%
  dplyr::transmute(
    Diversity_order = .data$q_label,
    Component = .data$component_label,
    Analytical_or_plugin = round(.data$analytical_estimate, 3),
    Fixed_MC_mean = round(.data$fixed_mean, 3),
    Fixed_MC_lower_95 = round(.data$fixed_lower_95, 3),
    Fixed_MC_upper_95 = round(.data$fixed_upper_95, 3),
    Full_MC_mean = round(.data$full_mean, 3),
    Full_MC_lower_95 = round(.data$full_lower_95, 3),
    Full_MC_upper_95 = round(.data$full_upper_95, 3),
    Analytical_percent_error = round(.data$analytical_percent_error, 2)
  )
readr::write_csv(
  analytical_publication,
  file.path(manuscript_table_dir, "Table_analytical_benchmark.csv"),
  na = ""
)

monthly_identity_data <- monthly_comparison %>%
  dplyr::mutate(
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    ),
    absolute_error = abs(.data$analytical_percent_error)
  )

monthly_identity <- ggplot2::ggplot(
  monthly_identity_data,
  ggplot2::aes(
    x = .data$fixed_mean,
    y = .data$analytical_estimate,
    colour = .data$q_label
  )
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$fixed_lower_95, xmax = .data$fixed_upper_95),
    orientation = "y",
    width = 0,
    alpha = 0.5,
    linewidth = 0.5
  ) +
  ggplot2::geom_point(size = 2.5, alpha = 0.85) +
  ggplot2::facet_wrap(~q_label, scales = "free", nrow = 1) +
  ggplot2::scale_colour_manual(values = palette$q) +
  ggplot2::labs(
    title = "Analytical estimates versus fixed-probability Monte Carlo",
    subtitle = "Horizontal intervals show 95% fixed-probability Monte Carlo simulation intervals",
    x = "Fixed-probability Monte Carlo mean",
    y = "Analytical or plug-in estimate",
    colour = "Diversity order",
    caption = "The dashed line is exact agreement. q = 1 and q = 2 are plug-in approximations based on expected candidate counts."
  ) +
  theme_silwood(9) +
  ggplot2::theme(legend.position = "none")

save_figure(
  monthly_identity,
  file.path(figure_dir, "figure_12_analytical_vs_fixed_mc_monthly"),
  width = 12.5,
  height = 5.8
)

partition_plot_base <- partition_comparison %>%
  dplyr::filter(.data$weighting == "equal_month") %>%
  dplyr::select(
    q,
    q_label,
    component,
    component_label,
    analytical_estimate,
    fixed_median,
    fixed_lower_95,
    fixed_upper_95,
    full_median,
    full_lower_95,
    full_upper_95
  )

# Build the three plotting approaches explicitly.  This avoids a data-mask
# failure after pivot_longer(), where analytical_estimate no longer exists as
# a column once it has been gathered into the generic `estimate` column.
partition_plot_data <- dplyr::bind_rows(
  partition_plot_base %>%
    dplyr::transmute(
      q,
      q_label,
      component,
      component_label,
      approach = "analytical_estimate",
      estimate = .data$analytical_estimate,
      lower_95 = .data$analytical_estimate,
      upper_95 = .data$analytical_estimate
    ),
  partition_plot_base %>%
    dplyr::transmute(
      q,
      q_label,
      component,
      component_label,
      approach = "fixed_median",
      estimate = .data$fixed_median,
      lower_95 = .data$fixed_lower_95,
      upper_95 = .data$fixed_upper_95
    ),
  partition_plot_base %>%
    dplyr::transmute(
      q,
      q_label,
      component,
      component_label,
      approach = "full_median",
      estimate = .data$full_median,
      lower_95 = .data$full_lower_95,
      upper_95 = .data$full_upper_95
    )
) %>%
  dplyr::mutate(
    approach_label = dplyr::recode(
      .data$approach,
      analytical_estimate = "Analytical / plug-in",
      fixed_median = "Fixed-probability MC",
      full_median = "Full MC"
    ),
    approach_label = factor(
      .data$approach_label,
      levels = c("Analytical / plug-in", "Fixed-probability MC", "Full MC")
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

partition_comparison_figure <- ggplot2::ggplot(
  partition_plot_data,
  ggplot2::aes(
    x = .data$estimate,
    y = .data$approach_label,
    colour = .data$approach_label
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$lower_95, xmax = .data$upper_95),
    orientation = "y",
    width = 0,
    linewidth = 0.75
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$component_label),
    cols = ggplot2::vars(.data$q_label),
    scales = "free_x"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Analytical / plug-in" = "#E69F00",
      "Fixed-probability MC" = "#009E73",
      "Full MC" = "#0072B2"
    )
  ) +
  ggplot2::labs(
    title = "Analytical benchmark, fixed-probability MC and full MC",
    subtitle = "Equal-month alpha, gamma and multiplicative beta diversity",
    x = "Effective diversity",
    y = NULL,
    colour = "Approach",
    caption = "Analytical values are point estimates. Fixed-probability MC adds candidate-state uncertainty; full MC also adds calibration-parameter uncertainty."
  ) +
  theme_silwood(8.5) +
  ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))

save_figure(
  partition_comparison_figure,
  file.path(figure_dir, "figure_13_analytical_partition_comparison"),
  width = 13.5,
  height = 8.5
)

error_heatmap <- dplyr::bind_rows(
  monthly_comparison %>%
    dplyr::group_by(.data$q, .data$q_label) %>%
    dplyr::summarise(
      comparison = "Monthly",
      mean_absolute_percent_error = mean(
        abs(.data$analytical_percent_error),
        na.rm = TRUE
      ),
      .groups = "drop"
    ),
  partition_comparison %>%
    dplyr::filter(.data$weighting == "equal_month") %>%
    dplyr::transmute(
      q,
      q_label,
      comparison = paste0("Partition: ", .data$component_label),
      mean_absolute_percent_error = abs(.data$analytical_percent_error)
    )
) %>%
  dplyr::mutate(
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    ),
    comparison = factor(
      .data$comparison,
      levels = c("Monthly", "Partition: Alpha", "Partition: Gamma", "Partition: Beta")
    ),
    label = sprintf("%.2f%%", .data$mean_absolute_percent_error)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data$q_label,
      y = .data$comparison,
      fill = .data$mean_absolute_percent_error
    )
  ) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = .data$label), fontface = "bold") +
  ggplot2::scale_fill_gradient(
    low = "#F7FBFF",
    high = "#B2182B",
    labels = scales::label_number(suffix = "%", accuracy = 0.1)
  ) +
  ggplot2::labs(
    title = "Absolute error of the analytical / plug-in benchmark",
    subtitle = "Error is measured relative to the fixed-probability Monte Carlo mean",
    x = NULL,
    y = NULL,
    fill = "Absolute error"
  ) +
  theme_silwood(9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")

save_figure(
  error_heatmap,
  file.path(figure_dir, "figure_S06_analytical_approximation_error"),
  width = 9.5,
  height = 5.3
)

message("Analytical benchmark completed.")
message("Tables: ", table_dir)
message("Figures: ", figure_dir)

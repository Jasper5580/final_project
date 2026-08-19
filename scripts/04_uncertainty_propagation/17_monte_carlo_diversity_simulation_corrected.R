# ============================================================
# 17_monte_carlo_diversity_simulation_corrected.R
#
# Output 1: full Monte Carlo uncertainty propagation.
#
# Propagated uncertainty:
#   1. fixed-effect covariance;
#   2. conditional predicted-species random-intercept uncertainty;
#   3. candidate-level Bernoulli classification uncertainty.
#
# Scenarios:
#   Primary_M3, Primary_M2, Audio_only.
#
# The Primary M3 species-by-month count array is retained so the
# higher-coverage analysis can reuse identical simulation draws.
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
model_bundle_file <- file.path(
  paths$results_root,
  "02_calibration",
  "models",
  "calibration_models.rds"
)

if (!file.exists(calibrated_file)) {
  stop("Calibrated candidate RDS not found: ", calibrated_file)
}
if (!file.exists(model_bundle_file)) {
  stop("Calibration model bundle not found: ", model_bundle_file)
}

calibrated <- readRDS(calibrated_file)
model_bundle <- readRDS(model_bundle_file)

required_columns <- c(
  "id",
  "event_id",
  "year_month",
  "predicted_species",
  "logit_score",
  "actual_multi_event",
  "p_primary",
  "p_primary_m2",
  "p_audio_only"
)
assert_columns(calibrated, required_columns, "calibrated candidates")

calibrated <- calibrated %>%
  dplyr::mutate(
    id = as.character(.data$id),
    event_id = as.character(.data$event_id),
    year_month = as.character(.data$year_month),
    predicted_species = as.character(.data$predicted_species),
    logit_score = as.numeric(.data$logit_score),
    actual_multi_event = as.integer(.data$actual_multi_event)
  ) %>%
  dplyr::arrange(.data$year_month, .data$predicted_species, .data$id)

assert_expected(nrow(calibrated), 65255, "Monte Carlo candidate rows")

primary_m2 <- extract_fitted_model(model_bundle, "Primary__M2")
primary_m3 <- extract_fitted_model(model_bundle, "Primary__M3")
audio_m2 <- extract_fitted_model(model_bundle, "Audio_only__M2")
audio_m3 <- extract_fitted_model(model_bundle, "Audio_only__M3")

primary_m2_fixed <- prepare_fixed_effect_information(primary_m2)
primary_m3_fixed <- prepare_fixed_effect_information(primary_m3)
primary_m3_random <- prepare_random_intercept_information(primary_m3)
audio_m2_fixed <- prepare_fixed_effect_information(audio_m2)
audio_m3_fixed <- prepare_fixed_effect_information(audio_m3)
audio_m3_random <- prepare_random_intercept_information(audio_m3)

aggregation <- build_candidate_aggregation(calibrated)

scenario_order <- c("Primary_M3", "Primary_M2", "Audio_only")
scenario_labels <- c(
  Primary_M3 = "Primary M3: species-aware",
  Primary_M2 = "Primary M2: pooled",
  Audio_only = "Audio-only sensitivity"
)

weighting_order <- c("equal_month", "abundance_weighted")

probability_draw_function <- function(scenario) {
  if (scenario == "Primary_M3") {
    return(
      draw_m3_hybrid_probabilities(
        primary_m3_fixed,
        primary_m2_fixed,
        primary_m3_random,
        calibrated
      )
    )
  }

  if (scenario == "Primary_M2") {
    return(draw_m2_probabilities(primary_m2_fixed, calibrated))
  }

  if (scenario == "Audio_only") {
    return(
      draw_m3_hybrid_probabilities(
        audio_m3_fixed,
        audio_m2_fixed,
        audio_m3_random,
        calibrated
      )
    )
  }

  stop("Unknown probability scenario: ", scenario)
}

set.seed(config$seed)

n_scenarios <- length(scenario_order)
monthly_draws <- vector(
  "list",
  config$n_simulations * n_scenarios
)
partition_draws <- vector(
  "list",
  config$n_simulations * n_scenarios
)
scenario_draw_summary <- vector(
  "list",
  config$n_simulations * n_scenarios
)

primary_count_array <- array(
  0L,
  dim = c(
    config$n_simulations,
    aggregation$n_months,
    aggregation$n_species
  ),
  dimnames = list(
    simulation = as.character(seq_len(config$n_simulations)),
    year_month = aggregation$month_levels,
    predicted_species = aggregation$species_levels
  )
)

output_index <- 1L

message(
  "Running full Monte Carlo propagation: ",
  config$n_simulations,
  " simulations x ",
  n_scenarios,
  " probability scenarios."
)

for (scenario in scenario_order) {
  scenario_value <- as.character(scenario)
  scenario_label_value <- unname(scenario_labels[[scenario_value]])
  message("Scenario: ", scenario_label_value)

  for (simulation_id in seq_len(config$n_simulations)) {
    simulation_id_value <- as.integer(simulation_id)
    if (
      simulation_id == 1L ||
        simulation_id %% max(1L, floor(config$n_simulations / 10L)) == 0L
    ) {
      message("  simulation ", simulation_id, "/", config$n_simulations)
    }

    probability_draw <- probability_draw_function(scenario_value)
    candidate_state <- stats::rbinom(
      n = nrow(calibrated),
      size = 1L,
      prob = probability_draw
    )

    count_matrix <- aggregate_candidate_states(
      candidate_state,
      aggregation
    )

    if (scenario_value == "Primary_M3") {
      primary_count_array[simulation_id, , ] <- count_matrix
    }

    monthly_draws[[output_index]] <- monthly_diversity_from_matrix(
      count_matrix,
      aggregation$month_levels,
      config$q_values
    ) %>%
      dplyr::mutate(
        scenario = .env$scenario_value,
        scenario_label = .env$scenario_label_value,
        simulation_id = .env$simulation_id_value,
        q_label = q_label(.data$q)
      )

    partition_draws[[output_index]] <- purrr::map_dfr(
      weighting_order,
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
                scenario = .env$scenario_value,
                scenario_label = .env$scenario_label_value,
                simulation_id = .env$simulation_id_value,
                weighting = .env$weighting,
                q = .env$q_value,
                q_label = q_label(.env$q_value),
                component_label = component_label(.data$component)
              )
          }
        )
      }
    )

    scenario_draw_summary[[output_index]] <- tibble::tibble(
      scenario = scenario_value,
      scenario_label = scenario_label_value,
      simulation_id = simulation_id_value,
      mean_candidate_probability = mean(probability_draw),
      expected_correct_candidates = sum(probability_draw),
      accepted_candidates = sum(candidate_state),
      accepted_species_month_cells = sum(count_matrix > 0)
    )

    output_index <- output_index + 1L
  }
}

monthly_draws <- dplyr::bind_rows(monthly_draws)
partition_draws <- dplyr::bind_rows(partition_draws)
scenario_draw_summary <- dplyr::bind_rows(scenario_draw_summary)

monthly_summary <- monthly_draws %>%
  dplyr::group_by(
    .data$scenario,
    .data$scenario_label,
    .data$year_month,
    .data$q,
    .data$q_label
  ) %>%
  dplyr::group_modify(
    ~summarise_distribution(.x$diversity)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    factor(.data$scenario, levels = scenario_order),
    .data$year_month,
    .data$q
  )

partition_summary <- partition_draws %>%
  dplyr::group_by(
    .data$scenario,
    .data$scenario_label,
    .data$weighting,
    .data$q,
    .data$q_label,
    .data$component,
    .data$component_label
  ) %>%
  dplyr::group_modify(
    ~summarise_distribution(.x$diversity)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    factor(.data$scenario, levels = scenario_order),
    factor(.data$weighting, levels = weighting_order),
    .data$q,
    factor(.data$component, levels = c("alpha", "gamma", "beta"))
  )

scenario_summary <- scenario_draw_summary %>%
  dplyr::group_by(.data$scenario, .data$scenario_label) %>%
  dplyr::summarise(
    simulations = dplyr::n(),
    mean_of_mean_candidate_probability = mean(.data$mean_candidate_probability),
    mean_expected_correct_candidates = mean(.data$expected_correct_candidates),
    mean_accepted_candidates = mean(.data$accepted_candidates),
    sd_accepted_candidates = stats::sd(.data$accepted_candidates),
    lower_accepted_candidates = stats::quantile(
      .data$accepted_candidates,
      0.025,
      names = FALSE,
      type = 8
    ),
    upper_accepted_candidates = stats::quantile(
      .data$accepted_candidates,
      0.975,
      names = FALSE,
      type = 8
    ),
    .groups = "drop"
  )

quality_checks <- tibble::tibble(
  check = c(
    "Simulations per scenario",
    "Monthly draw rows",
    "Partition draw rows",
    "Primary count-array simulations",
    "Primary count-array months",
    "Primary count-array species",
    "Missing monthly diversity values",
    "Missing partition diversity values"
  ),
  value = c(
    min(scenario_summary$simulations),
    nrow(monthly_draws),
    nrow(partition_draws),
    dim(primary_count_array)[1],
    dim(primary_count_array)[2],
    dim(primary_count_array)[3],
    sum(!is.finite(monthly_draws$diversity)),
    sum(!is.finite(partition_draws$diversity))
  ),
  expected = c(
    config$n_simulations,
    config$n_simulations * n_scenarios * aggregation$n_months *
      length(config$q_values),
    config$n_simulations * n_scenarios * length(weighting_order) *
      length(config$q_values) * 3L,
    config$n_simulations,
    aggregation$n_months,
    aggregation$n_species,
    0,
    0
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
  silwood_result_dir(paths, "05_monte_carlo", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "05_monte_carlo", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "05_monte_carlo", "figures")
)

results_object <- list(
  monthly_draws = monthly_draws,
  partition_draws = partition_draws,
  scenario_draw_summary = scenario_draw_summary,
  monthly_summary = monthly_summary,
  partition_summary = partition_summary,
  scenario_summary = scenario_summary,
  primary_count_array = primary_count_array,
  month_levels = aggregation$month_levels,
  species_levels = aggregation$species_levels,
  config = config,
  interpretation = list(
    analysis_unit = "event-species candidate",
    count_unit = "retained candidate record, not individual bird",
    random_effect_uncertainty = "conditional random-intercept approximation",
    random_effect_variance_uncertainty = "not explicitly propagated"
  )
)

saveRDS(
  results_object,
  file.path(model_dir, "monte_carlo_results.rds"),
  compress = "xz"
)
readr::write_csv(
  monthly_summary,
  file.path(table_dir, "monte_carlo_monthly_diversity_summary.csv"),
  na = ""
)
readr::write_csv(
  partition_summary,
  file.path(table_dir, "monte_carlo_partition_summary.csv"),
  na = ""
)
readr::write_csv(
  scenario_summary,
  file.path(table_dir, "monte_carlo_scenario_summary.csv"),
  na = ""
)
readr::write_csv(
  quality_checks,
  file.path(table_dir, "monte_carlo_quality_checks.csv")
)

manuscript_table_dir <- ensure_dir(
  silwood_result_dir(paths, "05_monte_carlo", "tables")
)
primary_partition_publication <- partition_summary %>%
  dplyr::filter(
    .data$scenario == "Primary_M3",
    .data$weighting == "equal_month"
  ) %>%
  dplyr::transmute(
    Diversity_order = .data$q_label,
    Component = .data$component_label,
    Median = round(.data$median, 3),
    Lower_95 = round(.data$lower_95, 3),
    Upper_95 = round(.data$upper_95, 3)
  )
readr::write_csv(
  primary_partition_publication,
  file.path(manuscript_table_dir, "Table_primary_M3_partition.csv"),
  na = ""
)

monthly_plot_data <- monthly_summary %>%
  dplyr::mutate(
    month_date = as.Date(paste0(.data$year_month, "-01"))
  )

# Insert calendar months with no represented detections as explicit NA rows.
# This prevents ggplot from drawing a continuous line across long periods for
# which no community was observed.
full_month_sequence <- seq(
  min(monthly_plot_data$month_date, na.rm = TRUE),
  max(monthly_plot_data$month_date, na.rm = TRUE),
  by = "month"
)

monthly_plot_data <- monthly_plot_data %>%
  dplyr::group_by(.data$scenario, .data$scenario_label, .data$q, .data$q_label) %>%
  tidyr::complete(month_date = full_month_sequence) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(scenario_labels[scenario_order])
    ),
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    )
  )

monthly_figure <- ggplot2::ggplot(
  monthly_plot_data,
  ggplot2::aes(
    x = .data$month_date,
    y = .data$median,
    colour = .data$scenario_label,
    fill = .data$scenario_label,
    group = .data$scenario_label
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_point(size = 1.5, alpha = 0.75) +
  ggplot2::facet_wrap(~q_label, scales = "free_y", ncol = 1) +
  ggplot2::scale_colour_manual(values = palette$scenario) +
  ggplot2::scale_fill_manual(values = palette$scenario) +
  ggplot2::scale_x_date(
    date_breaks = "3 months",
    labels = function(x) english_month_year_labels(x, line_break = TRUE),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::labs(
    title = "Monthly diversity after classification-uncertainty propagation",
    subtitle = "Lines show Monte Carlo medians; ribbons show 95% simulation intervals",
    x = NULL,
    y = "Effective species number",
    colour = "Probability scenario",
    fill = "Probability scenario",
    caption = "Candidate-record counts are acoustic classification units, not individual-bird abundance."
  ) +
  theme_silwood(10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom"
  )

save_figure(
  monthly_figure,
  file.path(figure_dir, "figure_05_monthly_monte_carlo_diversity"),
  width = 12,
  height = 10.5
)

partition_plot_data <- partition_summary %>%
  dplyr::filter(.data$weighting == "equal_month") %>%
  dplyr::mutate(
    scenario_label = factor(
      .data$scenario_label,
      levels = rev(unname(scenario_labels[scenario_order]))
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

partition_figure <- ggplot2::ggplot(
  partition_plot_data,
  ggplot2::aes(
    x = .data$median,
    y = .data$scenario_label,
    colour = .data$scenario_label
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
  ggplot2::scale_colour_manual(values = palette$scenario) +
  ggplot2::labs(
    title = "Equal-month alpha, gamma and multiplicative beta diversity",
    subtitle = "Points show medians; horizontal intervals show 95% Monte Carlo simulation intervals",
    x = "Effective diversity",
    y = NULL,
    colour = "Probability scenario",
    caption = "Beta is gamma divided by alpha within each simulation."
  ) +
  theme_silwood(9) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text.y = ggplot2::element_text(angle = 0)
  )

save_figure(
  partition_figure,
  file.path(figure_dir, "figure_06_monte_carlo_partition"),
  width = 13,
  height = 8.5
)

message("Full Monte Carlo propagation completed.")
message("Simulations per scenario: ", config$n_simulations)
message("Result object: ", file.path(model_dir, "monte_carlo_results.rds"))

# ============================================================
# 16_apply_calibration_to_full_data.R
#
# Purpose:
# Apply fitted calibration models to every cleaned event-species
# candidate and save point-estimate correctness probabilities.
#
# Primary rule:
#   known M3 species -> conditional M3 probability
#   species without an M3 effect -> pooled M2 fallback
#
# Audio-only uses the same hybrid rule. Primary M2 is retained
# as a pooled probability sensitivity scenario. M4 probabilities
# are exported for model-structure inspection only.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()

clean_data <- read_rds_or_csv(
  paths$clean_rds,
  paths$clean_csv,
  "cleaned BirdNET detections"
)

model_bundle_file <- file.path(
  paths$results_root,
  "02_calibration",
  "models",
  "calibration_models.rds"
)
if (!file.exists(model_bundle_file)) {
  stop("Calibration model bundle not found: ", model_bundle_file)
}
model_bundle <- readRDS(model_bundle_file)

selection_file <- file.path(
  paths$results_root,
  "03_model_validation",
  "tables",
  "model_selection_decision.csv"
)
if (!file.exists(selection_file)) {
  stop("Model-selection decision not found: ", selection_file)
}
selection_decision <- readr::read_csv(
  selection_file,
  show_col_types = FALSE
)
if (
  nrow(selection_decision) != 1L ||
    selection_decision$selected_primary_model[[1]] != "M3"
) {
  stop(
    "The downstream propagation code is configured for M3, but the ",
    "prespecified validation rule did not select M3. Review the M4 ",
    "diagnostics and update the propagation model deliberately before ",
    "continuing."
  )
}

validation <- read_rds_or_csv(
  file.path(
    paths$processed_root,
    "validation",
    "validation_analysis_ready.rds"
  ),
  file.path(
    paths$processed_root,
    "validation",
    "validation_analysis_ready.csv"
  ),
  "analysis-ready validation data"
)

required_clean_columns <- c(
  "id",
  "event_id",
  "audio_id",
  "species",
  "birdnet_score",
  "year_month"
)
assert_columns(clean_data, required_clean_columns, "cleaned detections")

calibrated <- clean_data %>%
  dplyr::mutate(
    id = as.character(.data$id),
    event_id = as.character(.data$event_id),
    audio_id = as.character(.data$audio_id),
    predicted_species = as.character(.data$species),
    birdnet_score = as.numeric(.data$birdnet_score),
    year_month = as.character(.data$year_month)
  ) %>%
  dplyr::group_by(.data$event_id) %>%
  dplyr::mutate(
    n_species_predicted = dplyr::n_distinct(.data$predicted_species),
    event_candidate_species = paste(
      sort(unique(.data$predicted_species)),
      collapse = " | "
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    actual_multi_event = as.integer(.data$n_species_predicted > 1L),
    score_bounded = bound_score(.data$birdnet_score),
    logit_score = stats::qlogis(.data$score_bounded)
  )

assert_expected(nrow(calibrated), 65255, "Calibrated candidate rows")
assert_expected(
  dplyr::n_distinct(calibrated$event_id),
  63985,
  "Calibrated acoustic events"
)
assert_expected(
  dplyr::n_distinct(calibrated$predicted_species),
  229,
  "Calibrated predicted species"
)

primary_m2 <- extract_fitted_model(model_bundle, "Primary__M2")
primary_m3 <- extract_fitted_model(model_bundle, "Primary__M3")
primary_m4 <- extract_fitted_model(
  model_bundle,
  "Primary__M4",
  required = FALSE
)
audio_m2 <- extract_fitted_model(model_bundle, "Audio_only__M2")
audio_m3 <- extract_fitted_model(model_bundle, "Audio_only__M3")

predict_pooled <- function(model, data) {
  clip_probability(
    stats::predict(model, newdata = data, type = "response")
  )
}

predict_conditional <- function(model, data) {
  random_table <- lme4::ranef(model)$predicted_species
  known_species <- rownames(random_table)
  all_species <- sort(unique(as.character(data$predicted_species)))

  prediction_data <- data %>%
    dplyr::mutate(
      predicted_species = factor(
        as.character(.data$predicted_species),
        levels = union(known_species, all_species)
      )
    )

  probability <- clip_probability(
    stats::predict(
      model,
      newdata = prediction_data,
      type = "response",
      allow.new.levels = TRUE,
      re.form = NULL
    )
  )

  list(
    probability = probability,
    known_species = known_species,
    species_known = as.integer(
      as.character(data$predicted_species) %in% known_species
    )
  )
}

primary_m2_probability <- predict_pooled(primary_m2, calibrated)
primary_m3_prediction <- predict_conditional(primary_m3, calibrated)
audio_m2_probability <- predict_pooled(audio_m2, calibrated)
audio_m3_prediction <- predict_conditional(audio_m3, calibrated)

primary_probability <- dplyr::if_else(
  primary_m3_prediction$species_known == 1L,
  primary_m3_prediction$probability,
  primary_m2_probability
)

audio_probability <- dplyr::if_else(
  audio_m3_prediction$species_known == 1L,
  audio_m3_prediction$probability,
  audio_m2_probability
)

if (!is.null(primary_m4)) {
  primary_m4_prediction <- predict_conditional(primary_m4, calibrated)
  primary_m4_probability <- primary_m4_prediction$probability
  primary_m4_species_known <- primary_m4_prediction$species_known
} else {
  primary_m4_probability <- rep(NA_real_, nrow(calibrated))
  primary_m4_species_known <- rep(0L, nrow(calibrated))
}

species_validation_summary <- validation %>%
  dplyr::filter(.data$analysis_primary_include == 1L) %>%
  dplyr::group_by(.data$predicted_species) %>%
  dplyr::summarise(
    species_validation_n = dplyr::n(),
    species_validation_events = dplyr::n_distinct(.data$event_id),
    species_correct_n = sum(.data$prediction_correct == 1L),
    species_incorrect_n = sum(.data$prediction_correct == 0L),
    species_confirmation_rate = mean(.data$prediction_correct),
    species_audio_only_n = sum(.data$validation_basis == "audio_based"),
    species_ecological_prior_n = sum(
      .data$validation_basis == "audio_plus_ecological_prior"
    ),
    .groups = "drop"
  )

calibrated <- calibrated %>%
  dplyr::mutate(
    species_in_primary_m3 = primary_m3_prediction$species_known,
    p_primary_m2 = primary_m2_probability,
    p_primary_m3_conditional = primary_m3_prediction$probability,
    p_primary = primary_probability,
    primary_model_used = dplyr::if_else(
      .data$species_in_primary_m3 == 1L,
      "Primary M3 conditional species effect",
      "Primary M2 pooled fallback"
    ),
    species_in_primary_m4 = primary_m4_species_known,
    p_primary_m4_conditional = primary_m4_probability,
    species_in_audio_m3 = audio_m3_prediction$species_known,
    p_audio_m2 = audio_m2_probability,
    p_audio_m3_conditional = audio_m3_prediction$probability,
    p_audio_only = audio_probability,
    audio_model_used = dplyr::if_else(
      .data$species_in_audio_m3 == 1L,
      "Audio-only M3 conditional species effect",
      "Audio-only M2 pooled fallback"
    )
  ) %>%
  dplyr::left_join(
    species_validation_summary,
    by = "predicted_species"
  )

probability_columns <- c(
  "p_primary_m2",
  "p_primary_m3_conditional",
  "p_primary",
  "p_audio_m2",
  "p_audio_m3_conditional",
  "p_audio_only"
)

invalid_probability_n <- sum(
  !is.finite(as.matrix(calibrated[, probability_columns])) |
    as.matrix(calibrated[, probability_columns]) <= 0 |
    as.matrix(calibrated[, probability_columns]) >= 1
)
if (invalid_probability_n > 0L) {
  stop("Invalid calibrated probabilities found: ", invalid_probability_n)
}

application_summary <- tibble::tibble(
  scenario = c(
    "Primary M3 hybrid",
    "Primary M2 pooled",
    "Audio-only hybrid"
  ),
  probability_column = c("p_primary", "p_primary_m2", "p_audio_only"),
  n_records = nrow(calibrated),
  n_species = dplyr::n_distinct(calibrated$predicted_species),
  expected_correct_candidates = c(
    sum(calibrated$p_primary),
    sum(calibrated$p_primary_m2),
    sum(calibrated$p_audio_only)
  ),
  mean_probability = c(
    mean(calibrated$p_primary),
    mean(calibrated$p_primary_m2),
    mean(calibrated$p_audio_only)
  ),
  median_probability = c(
    stats::median(calibrated$p_primary),
    stats::median(calibrated$p_primary_m2),
    stats::median(calibrated$p_audio_only)
  ),
  m2_fallback_records = c(
    sum(calibrated$species_in_primary_m3 == 0L),
    nrow(calibrated),
    sum(calibrated$species_in_audio_m3 == 0L)
  ),
  m2_fallback_species = c(
    dplyr::n_distinct(
      calibrated$predicted_species[
        calibrated$species_in_primary_m3 == 0L
      ]
    ),
    dplyr::n_distinct(calibrated$predicted_species),
    dplyr::n_distinct(
      calibrated$predicted_species[
        calibrated$species_in_audio_m3 == 0L
      ]
    )
  )
)

species_probability_summary <- calibrated %>%
  dplyr::group_by(.data$predicted_species) %>%
  dplyr::summarise(
    n_candidates = dplyr::n(),
    mean_birdnet_score = mean(.data$birdnet_score),
    mean_p_primary = mean(.data$p_primary),
    median_p_primary = stats::median(.data$p_primary),
    expected_correct_primary = sum(.data$p_primary),
    mean_p_primary_m2 = mean(.data$p_primary_m2),
    expected_correct_primary_m2 = sum(.data$p_primary_m2),
    mean_p_audio_only = mean(.data$p_audio_only),
    expected_correct_audio_only = sum(.data$p_audio_only),
    species_validation_n = dplyr::first(.data$species_validation_n),
    species_confirmation_rate = dplyr::first(
      .data$species_confirmation_rate
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    expected_difference_m2_minus_m3 =
      .data$expected_correct_primary_m2 - .data$expected_correct_primary
  ) %>%
  dplyr::arrange(dplyr::desc(abs(.data$expected_difference_m2_minus_m3)))

processed_dir <- ensure_dir(
  file.path(paths$processed_root, "calibration")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "04_probability_assignment", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "04_probability_assignment", "figures")
)

full_rds <- file.path(processed_dir, "calibrated_candidates.rds")
minimal_csv <- file.path(processed_dir, "calibrated_candidates_minimal.csv")

saveRDS(calibrated, full_rds)

minimal_output <- calibrated %>%
  dplyr::select(
    id,
    event_id,
    audio_id,
    year_month,
    predicted_species,
    birdnet_score,
    logit_score,
    n_species_predicted,
    actual_multi_event,
    p_primary,
    p_primary_m2,
    p_audio_only,
    p_primary_m4_conditional,
    primary_model_used,
    audio_model_used
  )
readr::write_csv(minimal_output, minimal_csv)
readr::write_csv(
  application_summary,
  file.path(table_dir, "probability_application_summary.csv")
)
readr::write_csv(
  species_probability_summary,
  file.path(table_dir, "species_probability_summary.csv")
)

plot_data <- calibrated %>%
  dplyr::select(
    id,
    actual_multi_event,
    Primary_M3 = p_primary,
    Primary_M2 = p_primary_m2,
    Audio_only = p_audio_only
  ) %>%
  tidyr::pivot_longer(
    cols = c(Primary_M3, Primary_M2, Audio_only),
    names_to = "scenario",
    values_to = "probability"
  ) %>%
  dplyr::mutate(
    scenario_label = dplyr::recode(
      .data$scenario,
      Primary_M3 = "Primary M3: species-aware",
      Primary_M2 = "Primary M2: pooled",
      Audio_only = "Audio-only sensitivity"
    ),
    event_type = dplyr::if_else(
      .data$actual_multi_event == 1L,
      "Multi-candidate event",
      "Single-candidate event"
    )
  )

probability_density_panel <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$probability,
    colour = .data$scenario_label,
    fill = .data$scenario_label
  )
) +
  ggplot2::geom_density(alpha = 0.12, linewidth = 0.9, adjust = 1.1) +
  ggplot2::facet_wrap(~event_type, nrow = 1) +
  ggplot2::scale_colour_manual(values = palette$scenario) +
  ggplot2::scale_fill_manual(values = palette$scenario) +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1)
  ) +
  ggplot2::labs(
    title = "A. Candidate-level calibrated probability distributions",
    x = "Estimated probability that the candidate label is correct",
    y = "Density",
    colour = "Scenario",
    fill = "Scenario"
  ) +
  theme_silwood(9)

species_difference_panel <- species_probability_summary %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::mutate(
    predicted_species = forcats::fct_reorder(
      .data$predicted_species,
      .data$expected_difference_m2_minus_m3
    ),
    direction = dplyr::if_else(
      .data$expected_difference_m2_minus_m3 >= 0,
      "M2 higher",
      "M3 higher"
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data$expected_difference_m2_minus_m3,
      y = .data$predicted_species,
      fill = .data$direction
    )
  ) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey25") +
  ggplot2::scale_fill_manual(
    values = c("M2 higher" = "#D55E00", "M3 higher" = "#0072B2")
  ) +
  ggplot2::labs(
    title = "B. Species most affected by pooled versus species-aware calibration",
    x = "Expected correct candidate count: M2 minus M3",
    y = NULL,
    fill = NULL
  ) +
  theme_silwood(9)

probability_overview <- probability_density_panel /
  species_difference_panel +
  patchwork::plot_annotation(
    title = "Application of calibrated correctness probabilities",
    subtitle = "Primary M2 provides a pooled sensitivity scenario; M3 adds conditional species effects",
    caption = "Probabilities refer to candidate-label correctness. They are not occupancy or individual-abundance probabilities."
  )

save_figure(
  probability_overview,
  file.path(figure_dir, "figure_04_probability_assignment"),
  width = 12,
  height = 9
)

message("Probability assignment completed.")
message("Full calibrated RDS: ", full_rds)
message("Minimal calibrated CSV: ", minimal_csv)

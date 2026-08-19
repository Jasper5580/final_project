# ============================================================
# 13_fit_final_calibration_models.R
#
# Purpose:
# Fit the streamlined calibration models and create detailed
# model-structure tables and publication-ready calibration plots.
#
# Primary models:
#   M1: score only
#   M2: score + multi-candidate event
#   M3: M2 + predicted-species random intercept
#   M4: M3 + predicted-species random score slope (sensitivity)
#
# Audio-only sensitivity models:
#   M2 and M3 only
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()

validation_rds <- file.path(
  paths$processed_root,
  "validation",
  "validation_analysis_ready.rds"
)
validation_csv <- file.path(
  paths$processed_root,
  "validation",
  "validation_analysis_ready.csv"
)

validation <- read_rds_or_csv(
  validation_rds,
  validation_csv,
  "analysis-ready validation data"
)

required_columns <- c(
  "validation_id",
  "event_id",
  "predicted_species",
  "birdnet_score",
  "score_bin",
  "logit_score",
  "actual_multi_event",
  "prediction_correct",
  "analysis_primary_include",
  "analysis_audio_only_sensitivity_include",
  "species_validation_n",
  "species_confirmation_rate"
)
assert_columns(validation, required_columns, "validation data")

prepare_model_data <- function(data) {
  data %>%
    dplyr::filter(
      .data$prediction_correct %in% c(0L, 1L),
      is.finite(.data$logit_score),
      .data$actual_multi_event %in% c(0L, 1L),
      !is.na(.data$predicted_species),
      !is.na(.data$event_id)
    ) %>%
    dplyr::mutate(
      predicted_species = factor(.data$predicted_species)
    )
}

primary_data <- validation %>%
  dplyr::filter(.data$analysis_primary_include == 1L) %>%
  prepare_model_data()

audio_only_data <- validation %>%
  dplyr::filter(.data$analysis_audio_only_sensitivity_include == 1L) %>%
  prepare_model_data()

assert_expected(nrow(primary_data), 1147, "Primary validation rows")
assert_expected(nrow(audio_only_data), 1016, "Audio-only validation rows")
assert_expected(
  dplyr::n_distinct(primary_data$predicted_species),
  229,
  "Primary validated species"
)

message("Fitting Primary M1-M4...")
fit_primary_m1 <- fit_calibration_model("M1", primary_data)
fit_primary_m2 <- fit_calibration_model("M2", primary_data)
fit_primary_m3 <- fit_calibration_model("M3", primary_data)
fit_primary_m4 <- fit_calibration_model("M4", primary_data)

message("Fitting Audio-only M2-M3...")
fit_audio_m2 <- fit_calibration_model("M2", audio_only_data)
fit_audio_m3 <- fit_calibration_model("M3", audio_only_data)

if (!fit_primary_m3$fit_ok) {
  stop("Primary M3 failed: ", fit_primary_m3$error)
}
if (!fit_audio_m3$fit_ok) {
  stop("Audio-only M3 failed: ", fit_audio_m3$error)
}

fit_results <- list(
  Primary__M1 = fit_primary_m1,
  Primary__M2 = fit_primary_m2,
  Primary__M3 = fit_primary_m3,
  Primary__M4 = fit_primary_m4,
  Audio_only__M2 = fit_audio_m2,
  Audio_only__M3 = fit_audio_m3
)

model_registry <- tibble::tribble(
  ~analysis_set, ~model_id, ~fit_key, ~model_label, ~role,
  "Primary", "M1", "Primary__M1", "M1: Score", "Baseline",
  "Primary", "M2", "Primary__M2", "M2: Score + event", "Pooled candidate",
  "Primary", "M3", "Primary__M3", "M3: Species-aware", "Primary candidate",
  "Primary", "M4", "Primary__M4", "M4: Random score slope", "Structural sensitivity",
  "Audio_only", "M2", "Audio_only__M2", "M2: Score + event", "Sensitivity",
  "Audio_only", "M3", "Audio_only__M3", "M3: Species-aware", "Sensitivity"
)

analysis_set_summary <- dplyr::bind_rows(
  primary_data %>%
    dplyr::summarise(
      analysis_set = "Primary",
      n_records = dplyr::n(),
      n_events = dplyr::n_distinct(.data$event_id),
      n_species = dplyr::n_distinct(.data$predicted_species),
      confirmation_rate = mean(.data$prediction_correct),
      multi_candidate_proportion = mean(.data$actual_multi_event)
    ),
  audio_only_data %>%
    dplyr::summarise(
      analysis_set = "Audio_only",
      n_records = dplyr::n(),
      n_events = dplyr::n_distinct(.data$event_id),
      n_species = dplyr::n_distinct(.data$predicted_species),
      confirmation_rate = mean(.data$prediction_correct),
      multi_candidate_proportion = mean(.data$actual_multi_event)
    )
)

summarise_model <- function(
    analysis_set,
    model_id,
    fit_key,
    model_label,
    role
) {
  fit_result <- fit_results[[fit_key]]
  data_used <- if (analysis_set == "Primary") primary_data else audio_only_data

  if (!fit_result$fit_ok) {
    return(
      tibble::tibble(
        analysis_set = analysis_set,
        model_id = model_id,
        model_label = model_label,
        role = role,
        fit_ok = FALSE,
        n_records = nrow(data_used),
        n_events = dplyr::n_distinct(data_used$event_id),
        n_species = dplyr::n_distinct(data_used$predicted_species),
        n_parameters = NA_real_,
        aic = NA_real_,
        bic = NA_real_,
        log_likelihood = NA_real_,
        singular_fit = NA,
        convergence_message = NA_character_,
        max_absolute_gradient = NA_real_,
        species_intercept_sd = NA_real_,
        species_score_slope_sd = NA_real_,
        warnings = paste(fit_result$warnings, collapse = " | "),
        error = fit_result$error
      )
    )
  }

  model <- fit_result$model
  n_parameters <- if (inherits(model, "merMod")) {
    as.numeric(attr(stats::logLik(model), "df"))
  } else {
    length(stats::coef(model))
  }

  tibble::tibble(
    analysis_set = analysis_set,
    model_id = model_id,
    model_label = model_label,
    role = role,
    fit_ok = TRUE,
    n_records = nrow(data_used),
    n_events = dplyr::n_distinct(data_used$event_id),
    n_species = dplyr::n_distinct(data_used$predicted_species),
    n_parameters = n_parameters,
    aic = stats::AIC(model),
    bic = stats::BIC(model),
    log_likelihood = as.numeric(stats::logLik(model)),
    singular_fit = if (inherits(model, "merMod")) {
      lme4::isSingular(model, tol = 1e-4)
    } else {
      NA
    },
    convergence_message = extract_convergence_message(model),
    max_absolute_gradient = extract_max_gradient(model),
    species_intercept_sd = extract_random_sd(model, "(Intercept)"),
    species_score_slope_sd = extract_random_sd(model, "logit_score"),
    warnings = if (length(fit_result$warnings) == 0L) {
      NA_character_
    } else {
      paste(fit_result$warnings, collapse = " | ")
    },
    error = NA_character_
  )
}

model_comparison <- purrr::pmap_dfr(model_registry, summarise_model) %>%
  dplyr::group_by(.data$analysis_set) %>%
  dplyr::mutate(
    delta_aic = dplyr::if_else(
      is.finite(.data$aic),
      .data$aic - min(.data$aic, na.rm = TRUE),
      NA_real_
    )
  ) %>%
  dplyr::ungroup()

extract_fixed_effects <- function(
    analysis_set,
    model_id,
    fit_key,
    model_label,
    role
) {
  fit_result <- fit_results[[fit_key]]
  if (!fit_result$fit_ok) return(tibble::tibble())
  coefficient_matrix <- summary(fit_result$model)$coefficients
  tibble::tibble(
    analysis_set = analysis_set,
    model_id = model_id,
    model_label = model_label,
    term = rownames(coefficient_matrix),
    estimate = coefficient_matrix[, 1],
    std_error = coefficient_matrix[, 2],
    statistic = coefficient_matrix[, 3],
    p_value = coefficient_matrix[, 4]
  ) %>%
    dplyr::mutate(
      conf_low = .data$estimate - 1.96 * .data$std_error,
      conf_high = .data$estimate + 1.96 * .data$std_error,
      odds_ratio = exp(.data$estimate),
      odds_ratio_low = exp(.data$conf_low),
      odds_ratio_high = exp(.data$conf_high)
    )
}

fixed_effects <- purrr::pmap_dfr(model_registry, extract_fixed_effects)

extract_species_effects <- function(model, model_id, validation_data) {
  random_table <- lme4::ranef(model, condVar = TRUE)$predicted_species
  intercept_se <- extract_ranef_conditional_se(
    random_table,
    column_index = 1L
  )
  slope_se <- if (ncol(random_table) >= 2L) {
    extract_ranef_conditional_se(
      random_table,
      column_index = 2L
    )
  } else {
    rep(NA_real_, nrow(random_table))
  }

  output <- tibble::tibble(
    predicted_species = rownames(random_table),
    random_intercept = as.numeric(random_table[, "(Intercept)"]),
    random_intercept_se = as.numeric(intercept_se)
  )

  if ("logit_score" %in% colnames(random_table)) {
    output <- output %>%
      dplyr::mutate(
        random_score_slope = as.numeric(random_table[, "logit_score"]),
        random_score_slope_se = as.numeric(slope_se)
      )
  } else {
    output <- output %>%
      dplyr::mutate(
        random_score_slope = NA_real_,
        random_score_slope_se = NA_real_
      )
  }

  validation_summary <- validation_data %>%
    dplyr::mutate(predicted_species = as.character(.data$predicted_species)) %>%
    dplyr::group_by(.data$predicted_species) %>%
    dplyr::summarise(
      validation_n = dplyr::n(),
      validation_events = dplyr::n_distinct(.data$event_id),
      confirmation_rate = mean(.data$prediction_correct),
      score_min = min(.data$birdnet_score),
      score_max = max(.data$birdnet_score),
      .groups = "drop"
    )

  output %>%
    dplyr::left_join(validation_summary, by = "predicted_species") %>%
    dplyr::mutate(model_id = model_id)
}

primary_m3 <- fit_primary_m3$model
species_effects_m3 <- extract_species_effects(
  primary_m3,
  "M3",
  primary_data
)

species_effects_m4 <- if (fit_primary_m4$fit_ok) {
  extract_species_effects(
    fit_primary_m4$model,
    "M4",
    primary_data
  )
} else {
  tibble::tibble()
}

model_dir <- ensure_dir(
  silwood_result_dir(paths, "02_calibration", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "02_calibration", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "02_calibration", "figures")
)

model_bundle <- list(
  fitted_models = fit_results,
  model_registry = model_registry,
  analysis_set_summary = analysis_set_summary,
  model_comparison = model_comparison,
  fixed_effects = fixed_effects,
  species_effects_m3 = species_effects_m3,
  species_effects_m4 = species_effects_m4,
  primary_model_for_propagation = "M3",
  audio_only_model_for_propagation = "M3",
  config = config
)

saveRDS(
  model_bundle,
  file.path(model_dir, "calibration_models.rds")
)
readr::write_csv(
  analysis_set_summary,
  file.path(table_dir, "calibration_analysis_set_summary.csv")
)
readr::write_csv(
  model_comparison,
  file.path(table_dir, "calibration_model_comparison.csv"),
  na = ""
)
readr::write_csv(
  fixed_effects,
  file.path(table_dir, "calibration_fixed_effects.csv"),
  na = ""
)
readr::write_csv(
  species_effects_m3,
  file.path(table_dir, "calibration_species_effects_m3.csv"),
  na = ""
)
readr::write_csv(
  species_effects_m4,
  file.path(table_dir, "calibration_species_effects_m4.csv"),
  na = ""
)

manuscript_table_dir <- ensure_dir(
  silwood_result_dir(paths, "02_calibration", "tables")
)
calibration_model_publication <- model_comparison %>%
  dplyr::filter(.data$analysis_set == "Primary") %>%
  dplyr::transmute(
    Model = .data$model_label,
    Role = .data$role,
    Records = .data$n_records,
    Events = .data$n_events,
    Species = .data$n_species,
    AIC = round(.data$aic, 2),
    Delta_AIC = round(.data$delta_aic, 2),
    Species_intercept_SD = round(.data$species_intercept_sd, 3),
    Species_score_slope_SD = round(.data$species_score_slope_sd, 3),
    Singular_fit = .data$singular_fit,
    Maximum_gradient = signif(.data$max_absolute_gradient, 3)
  )
readr::write_csv(
  calibration_model_publication,
  file.path(manuscript_table_dir, "Table_calibration_model_structure.csv"),
  na = ""
)

score_grid <- tidyr::expand_grid(
  birdnet_score = seq(0.45, 0.999, length.out = 300),
  actual_multi_event = c(0L, 1L)
) %>%
  dplyr::mutate(
    logit_score = stats::qlogis(bound_score(.data$birdnet_score)),
    event_type = dplyr::if_else(
      .data$actual_multi_event == 1L,
      "Multi-candidate event",
      "Single-candidate event"
    )
  )

primary_model_keys <- c(
  M1 = "Primary__M1",
  M2 = "Primary__M2",
  M3 = "Primary__M3",
  M4 = "Primary__M4"
)

curve_predictions <- purrr::imap_dfr(
  primary_model_keys,
  function(fit_key, model_id) {
    fit_result <- fit_results[[fit_key]]
    if (!fit_result$fit_ok) return(tibble::tibble())
    fixed_prediction <- fixed_effect_prediction(
      fit_result$model,
      score_grid
    )
    dplyr::bind_cols(score_grid, fixed_prediction) %>%
      dplyr::mutate(
        model_id = model_id,
        model_label = dplyr::recode(
          model_id,
          M1 = "M1: Score",
          M2 = "M2: Score + event",
          M3 = "M3: Species-aware",
          M4 = "M4: Random score slope"
        )
      )
  }
)

observed_bins <- primary_data %>%
  dplyr::mutate(
    event_type = dplyr::if_else(
      .data$actual_multi_event == 1L,
      "Multi-candidate event",
      "Single-candidate event"
    )
  ) %>%
  dplyr::group_by(.data$event_type, .data$score_bin) %>%
  dplyr::summarise(
    mean_score = mean(.data$birdnet_score),
    confirmation_rate = mean(.data$prediction_correct),
    n_records = dplyr::n(),
    .groups = "drop"
  )

calibration_panel <- ggplot2::ggplot(
  curve_predictions,
  ggplot2::aes(
    x = .data$birdnet_score,
    y = .data$probability,
    colour = .data$model_label,
    fill = .data$model_label
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
    alpha = 0.10,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(
    data = observed_bins,
    ggplot2::aes(
      x = .data$mean_score,
      y = .data$confirmation_rate,
      size = .data$n_records
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.6
  ) +
  ggplot2::facet_wrap(~event_type, nrow = 1) +
  ggplot2::scale_colour_manual(values = palette$model) +
  ggplot2::scale_fill_manual(values = palette$model) +
  ggplot2::scale_size_continuous(range = c(2, 6)) +
  ggplot2::scale_x_continuous(
    breaks = c(0.45, 0.60, 0.75, 0.90, 1.00),
    limits = c(0.45, 1)
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1)
  ) +
  ggplot2::labs(
    title = "A. Validation-calibrated correctness curves",
    subtitle = "Bands show fixed-effect 95% intervals; points show observed score-bin rates",
    x = "BirdNET confidence score",
    y = "Probability that the candidate label is correct",
    colour = "Model",
    fill = "Model",
    size = "Validated n"
  ) +
  theme_silwood(10) +
  ggplot2::theme(legend.position = "bottom")

species_effects_m3_plot <- species_effects_m3 %>%
  dplyr::mutate(
    lower = .data$random_intercept - 1.96 * .data$random_intercept_se,
    upper = .data$random_intercept + 1.96 * .data$random_intercept_se,
    direction = dplyr::if_else(
      .data$random_intercept >= 0,
      "Above population baseline",
      "Below population baseline"
    )
  ) %>%
  dplyr::arrange(.data$random_intercept)

focus_n <- min(12L, max(1L, floor(nrow(species_effects_m3_plot) / 2L)))
species_focus <- dplyr::bind_rows(
  dplyr::slice_head(species_effects_m3_plot, n = focus_n),
  dplyr::slice_tail(species_effects_m3_plot, n = focus_n)
) %>%
  dplyr::distinct(.data$predicted_species, .keep_all = TRUE) %>%
  dplyr::mutate(
    predicted_species = forcats::fct_reorder(
      .data$predicted_species,
      .data$random_intercept
    )
  )

random_intercept_panel <- ggplot2::ggplot(
  species_focus,
  ggplot2::aes(
    x = .data$random_intercept,
    y = .data$predicted_species,
    colour = .data$direction,
    size = .data$validation_n
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
    orientation = "y",
    width = 0,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_colour_manual(
    values = c(
      "Above population baseline" = "#0072B2",
      "Below population baseline" = "#D55E00"
    )
  ) +
  ggplot2::scale_size_continuous(range = c(2, 5)) +
  ggplot2::labs(
    title = "B. Largest M3 species baseline deviations",
    x = "Conditional random intercept (log-odds scale)",
    y = NULL,
    colour = NULL,
    size = "Validated n"
  ) +
  theme_silwood(9) +
  ggplot2::theme(legend.position = "bottom")

if (nrow(species_effects_m4) > 0L) {
  random_slope_panel <- species_effects_m4 %>%
    dplyr::filter(is.finite(.data$random_score_slope)) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$random_intercept,
        y = .data$random_score_slope,
        size = .data$validation_n,
        colour = .data$confirmation_rate
      )
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_colour_gradient2(
      low = "#D55E00",
      mid = "#F7F7F7",
      high = "#0072B2",
      midpoint = 0.5,
      labels = scales::label_percent(accuracy = 1)
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 6)) +
    ggplot2::labs(
      title = "C. M4 species intercept and score-slope deviations",
      x = "Random intercept deviation",
      y = "Random score-slope deviation",
      colour = "Observed\nconfirmation",
      size = "Validated n"
    ) +
    theme_silwood(9)
} else {
  random_slope_panel <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "M4 did not fit; no random-slope panel available."
    ) +
    ggplot2::theme_void()
}

calibration_overview <- calibration_panel /
  (random_intercept_panel | random_slope_panel) +
  patchwork::plot_annotation(
    title = "BirdNET calibration model structure",
    subtitle = "M3 captures species baseline heterogeneity; M4 tests additional species variation in score slopes",
    caption = "M4 is a structural sensitivity model. Downstream propagation uses M3 unless M4 improves stable out-of-fold prediction."
  )

save_figure(
  calibration_overview,
  file.path(figure_dir, "figure_01_calibration_model_structure"),
  width = 13,
  height = 10
)

# Save plot components so that the final manuscript Figure 2 can be rebuilt
# after cross-validation without refitting any calibration model.
saveRDS(
  list(
    calibration_panel = calibration_panel,
    random_intercept_panel = random_intercept_panel,
    random_slope_panel = random_slope_panel
  ),
  file.path(model_dir, "calibration_plot_components.rds")
)

full_caterpillar <- species_effects_m3 %>%
  dplyr::mutate(
    lower = .data$random_intercept - 1.96 * .data$random_intercept_se,
    upper = .data$random_intercept + 1.96 * .data$random_intercept_se,
    predicted_species = forcats::fct_reorder(
      .data$predicted_species,
      .data$random_intercept
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data$random_intercept,
      y = .data$predicted_species,
      colour = .data$confirmation_rate
    )
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
    orientation = "y",
    width = 0,
    linewidth = 0.35,
    alpha = 0.75
  ) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_gradient2(
    low = "#D55E00",
    mid = "#F7F7F7",
    high = "#0072B2",
    midpoint = 0.5,
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "M3 conditional predicted-species random intercepts",
    x = "Random intercept (log-odds scale)",
    y = NULL,
    colour = "Observed\nconfirmation"
  ) +
  theme_silwood(8) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 5),
    legend.position = "right"
  )

save_figure(
  full_caterpillar,
  file.path(figure_dir, "figure_S01_m3_all_species_caterpillar"),
  width = 10,
  height = 22
)

message("Calibration model fitting completed.")
message("Model bundle: ", file.path(model_dir, "calibration_models.rds"))
message("Tables: ", table_dir)
message("Figures: ", figure_dir)

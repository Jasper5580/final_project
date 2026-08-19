# ============================================================
# 14_event_grouped_cross_validation.R
#
# Purpose:
# Compare Primary M1-M4 using repeated event-grouped
# cross-validation. All candidate labels from one acoustic event
# remain in the same fold, preventing event-level information
# leakage.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()

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

primary_data <- validation %>%
  dplyr::filter(.data$analysis_primary_include == 1L) %>%
  dplyr::transmute(
    validation_id = as.character(.data$validation_id),
    event_id = as.character(.data$event_id),
    predicted_species = as.character(.data$predicted_species),
    prediction_correct = as.integer(.data$prediction_correct),
    birdnet_score = as.numeric(.data$birdnet_score),
    logit_score = as.numeric(.data$logit_score),
    actual_multi_event = as.integer(.data$actual_multi_event)
  ) %>%
  dplyr::arrange(.data$event_id, .data$validation_id)

assert_expected(nrow(primary_data), 1147, "Primary CV records")
assert_unique(primary_data$validation_id, "Primary CV validation IDs")

model_ids <- c("M1", "M2", "M3", "M4")
model_labels <- c(
  M1 = "M1: Score",
  M2 = "M2: Score + event",
  M3 = "M3: Species-aware",
  M4 = "M4: Random score slope"
)

predict_fold_model <- function(model, model_id, test_data) {
  prediction <- tryCatch(
    {
      if (model_id %in% c("M3", "M4")) {
        stats::predict(
          model,
          newdata = test_data,
          type = "response",
          allow.new.levels = TRUE
        )
      } else {
        stats::predict(model, newdata = test_data, type = "response")
      }
    },
    error = function(e) e
  )

  if (inherits(prediction, "error")) {
    return(
      list(
        ok = FALSE,
        probability = rep(NA_real_, nrow(test_data)),
        error = conditionMessage(prediction)
      )
    )
  }

  list(
    ok = TRUE,
    probability = clip_probability(prediction),
    error = NA_character_
  )
}

run_cv_fold <- function(repeat_id, fold_id, fold_map) {
  test_events <- fold_map %>%
    dplyr::filter(.data$fold_id == .env$fold_id) %>%
    dplyr::pull(event_id)

  training_data <- primary_data %>%
    dplyr::filter(!.data$event_id %in% test_events)
  test_data <- primary_data %>%
    dplyr::filter(.data$event_id %in% test_events)

  if (
    nrow(training_data) == 0L ||
      nrow(test_data) == 0L ||
      dplyr::n_distinct(training_data$prediction_correct) < 2L
  ) {
    stop("Invalid train/test split at repeat ", repeat_id, ", fold ", fold_id)
  }

  training_species <- sort(unique(training_data$predicted_species))
  test_species <- sort(unique(test_data$predicted_species))
  unseen_species <- setdiff(test_species, training_species)

  training_data <- training_data %>%
    dplyr::mutate(
      predicted_species = factor(
        .data$predicted_species,
        levels = training_species
      )
    )

  test_data <- test_data %>%
    dplyr::mutate(
      predicted_species = factor(
        .data$predicted_species,
        levels = union(training_species, test_species)
      )
    )

  fold_predictions <- list()
  fold_diagnostics <- list()

  for (model_id in model_ids) {
    fit_result <- fit_calibration_model(model_id, training_data)

    if (!fit_result$fit_ok) {
      fold_diagnostics[[model_id]] <- tibble::tibble(
        repeat_id = repeat_id,
        fold_id = fold_id,
        model_id = .env$model_id,
        model_label = .env$model_labels[[.env$model_id]],
        fit_ok = FALSE,
        predict_ok = FALSE,
        n_training_records = nrow(training_data),
        n_test_records = nrow(test_data),
        n_training_events = dplyr::n_distinct(training_data$event_id),
        n_test_events = dplyr::n_distinct(test_data$event_id),
        n_unseen_test_species = length(unseen_species),
        n_unseen_test_records = sum(
          as.character(test_data$predicted_species) %in% unseen_species
        ),
        singular_fit = NA,
        convergence_message = NA_character_,
        max_absolute_gradient = NA_real_,
        warning = paste(fit_result$warnings, collapse = " | "),
        fit_error = fit_result$error,
        prediction_error = NA_character_
      )
      next
    }

    prediction_result <- predict_fold_model(
      fit_result$model,
      model_id,
      test_data
    )

    fold_diagnostics[[model_id]] <- tibble::tibble(
      repeat_id = repeat_id,
      fold_id = fold_id,
      model_id = .env$model_id,
      model_label = .env$model_labels[[.env$model_id]],
      fit_ok = TRUE,
      predict_ok = prediction_result$ok,
      n_training_records = nrow(training_data),
      n_test_records = nrow(test_data),
      n_training_events = dplyr::n_distinct(training_data$event_id),
      n_test_events = dplyr::n_distinct(test_data$event_id),
      n_unseen_test_species = length(unseen_species),
      n_unseen_test_records = sum(
        as.character(test_data$predicted_species) %in% unseen_species
      ),
      singular_fit = if (inherits(fit_result$model, "merMod")) {
        lme4::isSingular(fit_result$model, tol = 1e-4)
      } else {
        NA
      },
      convergence_message = extract_convergence_message(fit_result$model),
      max_absolute_gradient = extract_max_gradient(fit_result$model),
      warning = if (length(fit_result$warnings) == 0L) {
        NA_character_
      } else {
        paste(fit_result$warnings, collapse = " | ")
      },
      fit_error = NA_character_,
      prediction_error = prediction_result$error
    )

    if (!prediction_result$ok) next

    fold_predictions[[model_id]] <- test_data %>%
      dplyr::transmute(
        repeat_id = repeat_id,
        fold_id = fold_id,
        model_id = .env$model_id,
        model_label = .env$model_labels[[.env$model_id]],
        validation_id,
        event_id,
        predicted_species = as.character(.data$predicted_species),
        species_seen_in_training = as.integer(
          !as.character(.data$predicted_species) %in% unseen_species
        ),
        prediction_correct,
        birdnet_score,
        logit_score,
        actual_multi_event,
        predicted_probability = prediction_result$probability
      )
  }

  list(
    predictions = dplyr::bind_rows(fold_predictions),
    diagnostics = dplyr::bind_rows(fold_diagnostics)
  )
}

all_predictions <- list()
all_diagnostics <- list()
result_index <- 1L

message(
  "Running ",
  config$n_folds,
  "-fold event-grouped CV repeated ",
  config$n_repeats,
  " times."
)

for (repeat_id in seq_len(config$n_repeats)) {
  fold_map <- make_event_folds(
    primary_data,
    k = config$n_folds,
    seed_value = config$seed + repeat_id
  )

  for (fold_id in seq_len(config$n_folds)) {
    message("CV repeat ", repeat_id, "/", config$n_repeats,
            "; fold ", fold_id, "/", config$n_folds)
    fold_result <- run_cv_fold(repeat_id, fold_id, fold_map)
    all_predictions[[result_index]] <- fold_result$predictions
    all_diagnostics[[result_index]] <- fold_result$diagnostics
    result_index <- result_index + 1L
  }
}

cv_predictions <- dplyr::bind_rows(all_predictions)
cv_diagnostics <- dplyr::bind_rows(all_diagnostics)

if (nrow(cv_predictions) == 0L) stop("No CV predictions were produced.")

repeat_metrics <- cv_predictions %>%
  dplyr::group_by(.data$repeat_id, .data$model_id, .data$model_label) %>%
  dplyr::summarise(
    n_predictions = dplyr::n(),
    n_validation_records = dplyr::n_distinct(.data$validation_id),
    n_events = dplyr::n_distinct(.data$event_id),
    observed_confirmation_rate = mean(.data$prediction_correct),
    mean_predicted_probability = mean(.data$predicted_probability),
    brier_score = calculate_brier(
      .data$prediction_correct,
      .data$predicted_probability
    ),
    log_loss = calculate_log_loss(
      .data$prediction_correct,
      .data$predicted_probability
    ),
    calibration_intercept = calculate_calibration_intercept(
      .data$prediction_correct,
      .data$predicted_probability
    ),
    calibration_slope = calculate_calibration_slope(
      .data$prediction_correct,
      .data$predicted_probability
    ),
    unseen_species_record_proportion = mean(
      .data$species_seen_in_training == 0L
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    complete_repeat = .data$n_validation_records == nrow(primary_data)
  )

fold_status <- cv_diagnostics %>%
  dplyr::group_by(.data$model_id, .data$model_label) %>%
  dplyr::summarise(
    expected_folds = config$n_folds * config$n_repeats,
    successful_fit_folds = sum(.data$fit_ok),
    successful_prediction_folds = sum(.data$fit_ok & .data$predict_ok),
    failed_fit_folds = sum(!.data$fit_ok),
    failed_prediction_folds = sum(.data$fit_ok & !.data$predict_ok),
    singular_folds = sum(.data$singular_fit %in% TRUE, na.rm = TRUE),
    maximum_gradient_across_folds = if (
      all(is.na(.data$max_absolute_gradient))
    ) {
      NA_real_
    } else {
      max(.data$max_absolute_gradient, na.rm = TRUE)
    },
    .groups = "drop"
  )

model_summary <- repeat_metrics %>%
  dplyr::filter(.data$complete_repeat) %>%
  dplyr::group_by(.data$model_id, .data$model_label) %>%
  dplyr::summarise(
    complete_repeats = dplyr::n(),
    mean_brier_score = mean(.data$brier_score),
    se_brier_score = standard_error(.data$brier_score),
    mean_log_loss = mean(.data$log_loss),
    se_log_loss = standard_error(.data$log_loss),
    mean_calibration_intercept = mean(
      .data$calibration_intercept,
      na.rm = TRUE
    ),
    se_calibration_intercept = standard_error(.data$calibration_intercept),
    mean_calibration_slope = mean(.data$calibration_slope, na.rm = TRUE),
    se_calibration_slope = standard_error(.data$calibration_slope),
    mean_unseen_species_record_proportion = mean(
      .data$unseen_species_record_proportion
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(fold_status, by = c("model_id", "model_label")) %>%
  dplyr::mutate(
    brier_rank = dplyr::min_rank(.data$mean_brier_score),
    log_loss_rank = dplyr::min_rank(.data$mean_log_loss)
  ) %>%
  dplyr::arrange(.data$brier_rank, .data$model_id)

m3_summary <- model_summary %>% dplyr::filter(.data$model_id == "M3")
m4_summary <- model_summary %>% dplyr::filter(.data$model_id == "M4")

m4_stable <- nrow(m4_summary) == 1L &&
  m4_summary$failed_fit_folds[[1]] == 0L &&
  m4_summary$singular_folds[[1]] == 0L &&
  is.finite(m4_summary$maximum_gradient_across_folds[[1]]) &&
  m4_summary$maximum_gradient_across_folds[[1]] <= 0.01

m4_improved <- nrow(m3_summary) == 1L &&
  nrow(m4_summary) == 1L &&
  m4_summary$mean_brier_score[[1]] < m3_summary$mean_brier_score[[1]] &&
  m4_summary$mean_log_loss[[1]] < m3_summary$mean_log_loss[[1]]

selected_model <- if (config$mode == "quick") {
  "M3"
} else if (m4_stable && m4_improved) {
  "M4"
} else {
  "M3"
}

selection_decision <- tibble::tibble(
  selected_primary_model = selected_model,
  m4_stable = m4_stable,
  m4_improved_brier_and_log_loss = m4_improved,
  decision_rule = if (config$mode == "quick") {
    paste0(
      "Quick mode is a smoke test only. M3 is forced for downstream ",
      "pipeline testing; formal M3/M4 selection is performed only in final mode."
    )
  } else {
    paste0(
      "M4 replaces M3 only if all folds fit without singularity, ",
      "maximum gradient <= 0.01, and both mean Brier score and ",
      "mean log loss improve."
    )
  },
  interpretation = if (selected_model == "M3") {
    "M3 retained for downstream propagation; M4 remains a structural sensitivity model."
  } else {
    "M4 passed the prespecified rule and should be reviewed before downstream propagation."
  }
)

model_dir <- ensure_dir(
  silwood_result_dir(paths, "03_model_validation", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "03_model_validation", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "03_model_validation", "figures")
)

readr::write_csv(
  model_summary,
  file.path(table_dir, "event_cv_model_summary.csv"),
  na = ""
)
readr::write_csv(
  repeat_metrics,
  file.path(table_dir, "event_cv_repeat_metrics.csv"),
  na = ""
)
readr::write_csv(
  cv_diagnostics,
  file.path(table_dir, "event_cv_fold_diagnostics.csv"),
  na = ""
)
readr::write_csv(
  cv_predictions,
  file.path(table_dir, "event_cv_predictions.csv"),
  na = ""
)
readr::write_csv(
  selection_decision,
  file.path(table_dir, "model_selection_decision.csv")
)
saveRDS(
  list(
    predictions = cv_predictions,
    diagnostics = cv_diagnostics,
    repeat_metrics = repeat_metrics,
    model_summary = model_summary,
    selection_decision = selection_decision,
    config = config
  ),
  file.path(model_dir, "event_grouped_cv_results.rds")
)

manuscript_table_dir <- ensure_dir(
  silwood_result_dir(paths, "03_model_validation", "tables")
)
event_cv_publication <- model_summary %>%
  dplyr::transmute(
    Model = .data$model_label,
    Mean_Brier = round(.data$mean_brier_score, 4),
    SE_Brier = round(.data$se_brier_score, 4),
    Mean_log_loss = round(.data$mean_log_loss, 4),
    SE_log_loss = round(.data$se_log_loss, 4),
    Calibration_intercept = round(.data$mean_calibration_intercept, 3),
    Calibration_slope = round(.data$mean_calibration_slope, 3),
    Singular_folds = .data$singular_folds,
    Maximum_gradient = signif(.data$maximum_gradient_across_folds, 3),
    Selected_for_propagation = .data$model_id ==
      selection_decision$selected_primary_model[[1]]
  )
readr::write_csv(
  event_cv_publication,
  file.path(manuscript_table_dir, "Table_event_grouped_model_validation.csv"),
  na = ""
)

model_level_order <- unname(model_labels[model_ids])
repeat_plot_data <- repeat_metrics %>%
  dplyr::mutate(
    model_label = factor(.data$model_label, levels = model_level_order)
  )

brier_panel <- ggplot2::ggplot(
  repeat_plot_data,
  ggplot2::aes(
    x = .data$model_label,
    y = .data$brier_score,
    colour = .data$model_label,
    fill = .data$model_label
  )
) +
  ggplot2::geom_boxplot(width = 0.55, alpha = 0.18, outlier.shape = NA) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    size = 2.4,
    alpha = 0.85
  ) +
  ggplot2::scale_colour_manual(values = palette$model) +
  ggplot2::scale_fill_manual(values = palette$model) +
  ggplot2::labs(
    title = "A. Brier score",
    x = NULL,
    y = "Out-of-fold Brier score"
  ) +
  theme_silwood(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 22, hjust = 1),
    legend.position = "none"
  )

logloss_panel <- ggplot2::ggplot(
  repeat_plot_data,
  ggplot2::aes(
    x = .data$model_label,
    y = .data$log_loss,
    colour = .data$model_label,
    fill = .data$model_label
  )
) +
  ggplot2::geom_boxplot(width = 0.55, alpha = 0.18, outlier.shape = NA) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    size = 2.4,
    alpha = 0.85
  ) +
  ggplot2::scale_colour_manual(values = palette$model) +
  ggplot2::scale_fill_manual(values = palette$model) +
  ggplot2::labs(
    title = "B. Log loss",
    x = NULL,
    y = "Out-of-fold log loss"
  ) +
  theme_silwood(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 22, hjust = 1),
    legend.position = "none"
  )

calibration_intercept_panel <- ggplot2::ggplot(
  repeat_plot_data,
  ggplot2::aes(
    x = .data$model_label,
    y = .data$calibration_intercept,
    colour = .data$model_label
  )
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08),
    size = 2.4,
    alpha = 0.85
  ) +
  ggplot2::scale_colour_manual(values = palette$model) +
  ggplot2::labs(
    title = "C. Calibration intercept",
    subtitle = "Ideal value = 0",
    x = NULL,
    y = "Intercept"
  ) +
  theme_silwood(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 22, hjust = 1),
    legend.position = "none"
  )

calibration_slope_panel <- ggplot2::ggplot(
  repeat_plot_data,
  ggplot2::aes(
    x = .data$model_label,
    y = .data$calibration_slope,
    colour = .data$model_label
  )
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08),
    size = 2.4,
    alpha = 0.85
  ) +
  ggplot2::scale_colour_manual(values = palette$model) +
  ggplot2::labs(
    title = "D. Calibration slope",
    subtitle = "Ideal value = 1",
    x = NULL,
    y = "Slope"
  ) +
  theme_silwood(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 22, hjust = 1),
    legend.position = "none"
  )

cv_overview <- (brier_panel | logloss_panel) /
  (calibration_intercept_panel | calibration_slope_panel) +
  patchwork::plot_annotation(
    title = "Repeated event-grouped validation of calibration models",
    subtitle = paste0(
      config$n_folds,
      " folds x ",
      config$n_repeats,
      " repeats; each point represents one complete repeat"
    ),
    caption = paste0(
      "All candidate labels from one acoustic event remained in the same fold. ",
      selection_decision$interpretation[[1]]
    )
  )

save_figure(
  cv_overview,
  file.path(figure_dir, "figure_02_event_grouped_cv_overview"),
  width = 12,
  height = 9
)

summary_plot_data <- model_summary %>%
  dplyr::select(
    model_id,
    model_label,
    mean_brier_score,
    se_brier_score,
    mean_log_loss,
    se_log_loss
  ) %>%
  tidyr::pivot_longer(
    cols = c(mean_brier_score, mean_log_loss),
    names_to = "metric",
    values_to = "estimate"
  ) %>%
  dplyr::mutate(
    standard_error = dplyr::if_else(
      .data$metric == "mean_brier_score",
      .data$se_brier_score,
      .data$se_log_loss
    ),
    metric_label = dplyr::recode(
      .data$metric,
      mean_brier_score = "Brier score",
      mean_log_loss = "Log loss"
    ),
    model_label = factor(.data$model_label, levels = model_level_order)
  )

cv_summary_figure <- ggplot2::ggplot(
  summary_plot_data,
  ggplot2::aes(
    x = .data$model_label,
    y = .data$estimate,
    colour = .data$model_label
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = .data$estimate - .data$standard_error,
      ymax = .data$estimate + .data$standard_error
    ),
    width = 0.12,
    linewidth = 0.75
  ) +
  ggplot2::geom_point(size = 3.3) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", .data$estimate)),
    nudge_y = 0.01,
    size = 3.1,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~metric_label, scales = "free_y", nrow = 1) +
  ggplot2::scale_colour_manual(values = palette$model) +
  ggplot2::labs(
    title = "Mean out-of-fold prediction error",
    subtitle = "Points show means across repeated folds; error bars show one standard error",
    x = NULL,
    y = "Prediction error",
    colour = "Model"
  ) +
  theme_silwood(10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

save_figure(
  cv_summary_figure,
  file.path(figure_dir, "figure_03_event_grouped_cv_summary"),
  width = 10.5,
  height = 5.8
)

message("Event-grouped cross-validation completed.")
message("Selected downstream model: ", selected_model)
message("Tables: ", table_dir)
message("Figures: ", figure_dir)

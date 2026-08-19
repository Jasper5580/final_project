# ============================================================
# 15_leave_one_species_out_validation.R
#
# Optional supplementary analysis.
# Event-safe leave-one-species-out validation tests prediction
# for a completely unseen predicted species. The held-out
# species and every event containing that species are removed
# from the training data.
#
# This script is not required by the primary pipeline. Set
# SILWOOD_RUN_LOSO=true before running run_all.R to include it.
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
  )

held_species <- sort(unique(primary_data$predicted_species))
model_ids <- c("M1", "M2", "M3")
model_labels <- c(
  M1 = "M1: Score",
  M2 = "M2: Score + event",
  M3 = "M3: Population fixed effects"
)

loso_predictions <- vector("list", length(held_species) * length(model_ids))
loso_diagnostics <- vector("list", length(held_species) * length(model_ids))
output_index <- 1L

for (species_name in held_species) {
  test_data <- primary_data %>%
    dplyr::filter(.data$predicted_species == species_name)

  held_event_ids <- unique(test_data$event_id)
  training_data <- primary_data %>%
    dplyr::filter(!.data$event_id %in% held_event_ids)

  training_species <- sort(unique(training_data$predicted_species))
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
        levels = c(training_species, species_name)
      )
    )

  for (model_id in model_ids) {
    fit_result <- fit_calibration_model(model_id, training_data)

    if (!fit_result$fit_ok) {
      loso_diagnostics[[output_index]] <- tibble::tibble(
        held_species = species_name,
        model_id = .env$model_id,
        model_label = .env$model_labels[[.env$model_id]],
        fit_ok = FALSE,
        predict_ok = FALSE,
        n_test_records = nrow(test_data),
        n_training_records = nrow(training_data),
        n_removed_events = length(held_event_ids),
        singular_fit = NA,
        max_absolute_gradient = NA_real_,
        fit_error = fit_result$error,
        prediction_error = NA_character_
      )
      output_index <- output_index + 1L
      next
    }

    prediction <- tryCatch(
      {
        if (model_id == "M3") {
          stats::predict(
            fit_result$model,
            newdata = test_data,
            type = "response",
            re.form = NA,
            allow.new.levels = TRUE
          )
        } else {
          stats::predict(
            fit_result$model,
            newdata = test_data,
            type = "response"
          )
        }
      },
      error = function(e) e
    )

    predict_ok <- !inherits(prediction, "error")
    prediction_error <- if (predict_ok) NA_character_ else conditionMessage(prediction)

    loso_diagnostics[[output_index]] <- tibble::tibble(
      held_species = species_name,
      model_id = .env$model_id,
      model_label = .env$model_labels[[.env$model_id]],
      fit_ok = TRUE,
      predict_ok = predict_ok,
      n_test_records = nrow(test_data),
      n_training_records = nrow(training_data),
      n_removed_events = length(held_event_ids),
      singular_fit = if (inherits(fit_result$model, "merMod")) {
        lme4::isSingular(fit_result$model, tol = 1e-4)
      } else {
        NA
      },
      max_absolute_gradient = extract_max_gradient(fit_result$model),
      fit_error = NA_character_,
      prediction_error = prediction_error
    )

    if (predict_ok) {
      loso_predictions[[output_index]] <- test_data %>%
        dplyr::transmute(
          held_species = species_name,
          model_id = .env$model_id,
          model_label = .env$model_labels[[.env$model_id]],
          validation_id,
          event_id,
          prediction_correct,
          predicted_probability = clip_probability(prediction)
        )
    }

    output_index <- output_index + 1L
  }
}

prediction_table <- dplyr::bind_rows(loso_predictions)
diagnostic_table <- dplyr::bind_rows(loso_diagnostics)

species_metrics <- prediction_table %>%
  dplyr::group_by(.data$held_species, .data$model_id, .data$model_label) %>%
  dplyr::summarise(
    validation_n = dplyr::n(),
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
    .groups = "drop"
  )

record_weighted_summary <- prediction_table %>%
  dplyr::group_by(.data$model_id, .data$model_label) %>%
  dplyr::summarise(
    aggregation = "Record weighted",
    n_records = dplyr::n(),
    n_species = dplyr::n_distinct(.data$held_species),
    brier_score = calculate_brier(
      .data$prediction_correct,
      .data$predicted_probability
    ),
    log_loss = calculate_log_loss(
      .data$prediction_correct,
      .data$predicted_probability
    ),
    .groups = "drop"
  )

species_balanced_summary <- purrr::map_dfr(
  c(1L, 4L, 8L),
  function(minimum_n) {
    species_metrics %>%
      dplyr::filter(.data$validation_n >= minimum_n) %>%
      dplyr::group_by(.data$model_id, .data$model_label) %>%
      dplyr::summarise(
        aggregation = paste0("Species balanced, n >= ", minimum_n),
        n_records = sum(.data$validation_n),
        n_species = dplyr::n(),
        brier_score = mean(.data$brier_score),
        log_loss = mean(.data$log_loss),
        .groups = "drop"
      )
  }
)

summary_table <- dplyr::bind_rows(
  record_weighted_summary,
  species_balanced_summary
)

model_dir <- ensure_dir(
  silwood_result_dir(paths, "99_optional_loso", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "99_optional_loso", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "99_optional_loso", "figures")
)

readr::write_csv(
  prediction_table,
  file.path(table_dir, "loso_predictions.csv"),
  na = ""
)
readr::write_csv(
  diagnostic_table,
  file.path(table_dir, "loso_diagnostics.csv"),
  na = ""
)
readr::write_csv(
  species_metrics,
  file.path(table_dir, "loso_species_metrics.csv"),
  na = ""
)
readr::write_csv(
  summary_table,
  file.path(table_dir, "loso_summary.csv"),
  na = ""
)
saveRDS(
  list(
    predictions = prediction_table,
    diagnostics = diagnostic_table,
    species_metrics = species_metrics,
    summary = summary_table
  ),
  file.path(model_dir, "loso_results.rds")
)

plot_data <- summary_table %>%
  tidyr::pivot_longer(
    cols = c(brier_score, log_loss),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      .data$metric,
      brier_score = "Brier score",
      log_loss = "Log loss"
    ),
    model_label = factor(
      .data$model_label,
      levels = unname(model_labels[model_ids])
    )
  )

loso_figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$aggregation,
    y = .data$value,
    colour = .data$model_label,
    group = .data$model_label
  )
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
  ggplot2::scale_colour_manual(
    values = c(
      "M1: Score" = palette$model[["M1: Score"]],
      "M2: Score + event" = palette$model[["M2: Score + event"]],
      "M3: Population fixed effects" = palette$model[["M3: Species-aware"]]
    )
  ) +
  ggplot2::labs(
    title = "Supplementary leave-one-species-out generalisation",
    subtitle = "M3 predictions for unseen species use population-level fixed effects",
    x = NULL,
    y = "Prediction error",
    colour = "Model"
  ) +
  theme_silwood(9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

save_figure(
  loso_figure,
  file.path(figure_dir, "figure_S02_loso_generalisation"),
  width = 11,
  height = 5.8
)

message("Optional LOSO analysis completed.")

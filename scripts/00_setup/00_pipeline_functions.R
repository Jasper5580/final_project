# ============================================================
# 00_pipeline_functions.R
# Shared settings, validation helpers, model helpers, diversity
# functions, Monte Carlo helpers, and publication-style themes.
#
# This file is sourced by the complete final pipeline. It is not intended to
# be run independently.
# ============================================================

silwood_required_packages <- function() {
  c(
    "tidyverse",
    "here",
    "lme4",
    "Matrix",
    "MASS",
    "patchwork",
    "lubridate",
    "scales"
  )
}

load_silwood_packages <- function() {
  required_packages <- silwood_required_packages()
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Install the following R packages before running the pipeline: ",
      paste(missing_packages, collapse = ", "),
      "\nRun:\ninstall.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "),
      "))"
    )
  }

  suppressPackageStartupMessages({
    library(tidyverse)
    library(here)
    library(lme4)
    library(Matrix)
    library(patchwork)
    library(lubridate)
    library(scales)
  })

  invisible(TRUE)
}

silwood_config <- function() {
  mode <- tolower(Sys.getenv("SILWOOD_MODE", unset = "final"))
  if (!mode %in% c("quick", "final")) {
    stop("SILWOOD_MODE must be either 'quick' or 'final'.")
  }

  default_simulations <- if (mode == "quick") 20L else 1000L
  requested_simulations <- suppressWarnings(
    as.integer(Sys.getenv("SILWOOD_N_SIMULATIONS", unset = ""))
  )
  n_simulations <- if (
    is.na(requested_simulations) || requested_simulations < 1L
  ) {
    default_simulations
  } else {
    requested_simulations
  }

  run_loso <- tolower(Sys.getenv("SILWOOD_RUN_LOSO", unset = "false")) %in%
    c("true", "1", "yes", "y")

  default_folds <- if (mode == "quick") 3L else 5L
  default_repeats <- if (mode == "quick") 1L else 5L

  requested_folds <- suppressWarnings(
    as.integer(Sys.getenv("SILWOOD_N_FOLDS", unset = ""))
  )
  requested_repeats <- suppressWarnings(
    as.integer(Sys.getenv("SILWOOD_N_REPEATS", unset = ""))
  )

  n_folds <- if (is.na(requested_folds) || requested_folds < 2L) {
    default_folds
  } else {
    requested_folds
  }

  n_repeats <- if (is.na(requested_repeats) || requested_repeats < 1L) {
    default_repeats
  } else {
    requested_repeats
  }

  default_output_root <- "results"

  list(
    mode = mode,
    n_simulations = n_simulations,
    n_fixed_probability_simulations = n_simulations,
    n_folds = n_folds,
    n_repeats = n_repeats,
    seed = 20260729L,
    run_loso = run_loso,
    output_root = Sys.getenv(
      "SILWOOD_OUTPUT_ROOT",
      unset = default_output_root
    ),
    processed_root = Sys.getenv(
      "SILWOOD_PROCESSED_ROOT",
      unset = file.path("data", "03_analysis_ready")
    ),
    q_values = c(0, 1, 2),
    thresholds = c(0.50, 0.70, 0.90)
  )
}

silwood_paths <- function(config = silwood_config()) {
  results_root <- here::here(config$output_root)
  list(
    validation_input = here::here(
      "data", "02_manual_validation", "final",
      "validation_sample_silwood_combined_final.csv"
    ),
    clean_rds = here::here(
      "data", "01_cleaned", "birdnet_detections_clean.rds"
    ),
    clean_csv = here::here(
      "data", "01_cleaned", "birdnet_detections_clean.csv"
    ),
    raw_detection_csv = here::here(
      "data", "00_raw", "detections-export-proj_silwood_park(in).csv"
    ),
    processed_root = here::here(config$processed_root),
    results_root = results_root,
    output_root = results_root,
    final_output_root = file.path(results_root, "10_final_outputs"),
    log_root = here::here("log"),
    static_asset_root = here::here("data", "04_static_assets")
  )
}

silwood_result_dir <- function(paths, step, kind = NULL) {
  output <- file.path(paths$results_root, step)
  if (!is.null(kind)) output <- file.path(output, kind)
  ensure_dir(output)
}

silwood_final_dir <- function(paths, kind = NULL) {
  output <- paths$final_output_root
  if (!is.null(kind)) output <- file.path(output, kind)
  ensure_dir(output)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

assert_columns <- function(data, required_columns, data_label = "data") {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      data_label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  invisible(TRUE)
}

assert_unique <- function(x, label) {
  duplicated_n <- sum(duplicated(x))
  if (duplicated_n > 0L) {
    stop(label, " contains ", duplicated_n, " duplicate values.")
  }
  invisible(TRUE)
}

assert_expected <- function(actual, expected, label, tolerance = 0) {
  if (is.numeric(actual) && is.numeric(expected)) {
    valid <- isTRUE(all.equal(actual, expected, tolerance = tolerance))
  } else {
    valid <- identical(actual, expected)
  }
  if (!valid) {
    stop(
      label,
      " failed the regression check. Expected ",
      expected,
      ", observed ",
      actual,
      "."
    )
  }
  invisible(TRUE)
}

read_rds_or_csv <- function(rds_path, csv_path, label = "input") {
  if (file.exists(rds_path)) {
    message(label, " loaded from RDS: ", rds_path)
    return(readRDS(rds_path))
  }
  if (file.exists(csv_path)) {
    message(label, " loaded from CSV: ", csv_path)
    return(
      readr::read_csv(
        csv_path,
        na = c("", "NA", "N/A"),
        show_col_types = FALSE
      )
    )
  }
  stop(
    "Could not find ", label, ".\nExpected either:\n",
    rds_path, "\n", csv_path
  )
}

clip_probability <- function(probability, lower = 1e-8, upper = 1 - 1e-8) {
  pmin(pmax(as.numeric(probability), lower), upper)
}

bound_score <- function(score, lower = 0.001, upper = 0.999) {
  pmin(pmax(as.numeric(score), lower), upper)
}

as_binary_integer <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  if (is.numeric(x)) {
    return(as.integer(x))
  }
  normalised <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    normalised %in% c("1", "true", "yes", "y") ~ 1L,
    normalised %in% c("0", "false", "no", "n") ~ 0L,
    TRUE ~ NA_integer_
  )
}

safe_fit <- function(expression) {
  warning_messages <- character()
  result <- withCallingHandlers(
    tryCatch(
      eval.parent(substitute(expression)),
      error = function(error_object) error_object
    ),
    warning = function(warning_object) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(warning_object)
      )
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(result, "error")) {
    return(
      list(
        fit_ok = FALSE,
        model = NULL,
        warnings = unique(warning_messages),
        error = conditionMessage(result)
      )
    )
  }

  list(
    fit_ok = TRUE,
    model = result,
    warnings = unique(warning_messages),
    error = NA_character_
  )
}

glmm_control_silwood <- function() {
  lme4::glmerControl(
    optimizer = c("bobyqa", "Nelder_Mead"),
    optCtrl = list(maxfun = 200000),
    calc.derivs = TRUE
  )
}

fit_calibration_model <- function(model_id, data) {
  model_id <- as.character(model_id)
  control <- glmm_control_silwood()

  if (model_id == "M1") {
    return(
      safe_fit(
        glm(
          prediction_correct ~ logit_score,
          data = data,
          family = binomial(link = "logit")
        )
      )
    )
  }

  if (model_id == "M2") {
    return(
      safe_fit(
        glm(
          prediction_correct ~ logit_score + actual_multi_event,
          data = data,
          family = binomial(link = "logit")
        )
      )
    )
  }

  if (model_id == "M3") {
    return(
      safe_fit(
        lme4::glmer(
          prediction_correct ~
            logit_score +
            actual_multi_event +
            (1 | predicted_species),
          data = data,
          family = binomial(link = "logit"),
          control = control,
          nAGQ = 1
        )
      )
    )
  }

  if (model_id == "M4") {
    return(
      safe_fit(
        lme4::glmer(
          prediction_correct ~
            logit_score +
            actual_multi_event +
            (1 + logit_score || predicted_species),
          data = data,
          family = binomial(link = "logit"),
          control = control,
          nAGQ = 1
        )
      )
    )
  }

  stop("Unknown calibration model ID: ", model_id)
}

extract_fitted_model <- function(model_bundle, model_key, required = TRUE) {
  fitted_models <- model_bundle$fitted_models
  if (is.null(fitted_models)) {
    if (required) stop("The model bundle has no 'fitted_models' element.")
    return(NULL)
  }

  model_entry <- fitted_models[[model_key]]
  if (is.null(model_entry)) {
    if (required) stop("Model not found in bundle: ", model_key)
    return(NULL)
  }

  if (inherits(model_entry, c("glm", "merMod"))) {
    return(model_entry)
  }

  if (
    is.list(model_entry) &&
      !is.null(model_entry$model) &&
      inherits(model_entry$model, c("glm", "merMod"))
  ) {
    return(model_entry$model)
  }

  if (required) {
    stop("Model entry does not contain a fitted model: ", model_key)
  }
  NULL
}

extract_convergence_message <- function(model) {
  if (!inherits(model, "merMod")) return("not_applicable")
  messages <- model@optinfo$conv$lme4$messages
  if (is.null(messages)) return("none")
  paste(messages, collapse = " | ")
}

extract_max_gradient <- function(model) {
  if (!inherits(model, "merMod")) return(NA_real_)
  gradient <- model@optinfo$derivs$gradient
  if (is.null(gradient)) return(NA_real_)
  max(abs(gradient))
}

extract_random_sd <- function(model, term = "(Intercept)") {
  if (!inherits(model, "merMod")) return(NA_real_)
  variance_table <- as.data.frame(lme4::VarCorr(model))
  candidate <- variance_table %>%
    dplyr::filter(
      stringr::str_detect(.data$grp, "^predicted_species"),
      .data$var1 == term,
      is.na(.data$var2)
    ) %>%
    dplyr::slice_head(n = 1)
  if (nrow(candidate) == 0L) return(NA_real_)
  as.numeric(candidate$sdcor[[1]])
}

fixed_effect_vector <- function(model) {
  if (inherits(model, "merMod")) lme4::fixef(model) else stats::coef(model)
}

fixed_effect_covariance <- function(model) {
  as.matrix(stats::vcov(model))
}

fixed_effect_prediction <- function(model, new_data, level = 0.95) {
  coefficients <- fixed_effect_vector(model)
  covariance <- fixed_effect_covariance(model)
  coefficient_names <- names(coefficients)

  design <- matrix(
    0,
    nrow = nrow(new_data),
    ncol = length(coefficient_names),
    dimnames = list(NULL, coefficient_names)
  )

  if ("(Intercept)" %in% coefficient_names) {
    design[, "(Intercept)"] <- 1
  }
  if ("logit_score" %in% coefficient_names) {
    design[, "logit_score"] <- new_data$logit_score
  }
  if ("actual_multi_event" %in% coefficient_names) {
    design[, "actual_multi_event"] <- new_data$actual_multi_event
  }

  eta <- as.numeric(design %*% coefficients)
  variance_eta <- rowSums((design %*% covariance) * design)
  se_eta <- sqrt(pmax(variance_eta, 0))
  critical <- stats::qnorm(1 - (1 - level) / 2)

  tibble::tibble(
    eta = eta,
    se_eta = se_eta,
    probability = stats::plogis(eta),
    lower = stats::plogis(eta - critical * se_eta),
    upper = stats::plogis(eta + critical * se_eta)
  )
}

predict_model_probability <- function(model, model_id, new_data) {
  if (model_id %in% c("M3", "M4")) {
    return(
      clip_probability(
        stats::predict(
          model,
          newdata = new_data,
          type = "response",
          allow.new.levels = TRUE
        )
      )
    )
  }
  clip_probability(
    stats::predict(model, newdata = new_data, type = "response")
  )
}

make_event_folds <- function(data, k, seed_value) {
  event_table <- data %>%
    dplyr::group_by(.data$event_id) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      event_outcome_class = dplyr::case_when(
        all(.data$prediction_correct == 1L) ~ "all_correct",
        all(.data$prediction_correct == 0L) ~ "all_incorrect",
        TRUE ~ "mixed_outcome"
      ),
      .groups = "drop"
    )

  set.seed(seed_value)
  event_table %>%
    dplyr::group_by(.data$event_outcome_class) %>%
    dplyr::mutate(random_value = stats::runif(dplyr::n())) %>%
    dplyr::arrange(
      .data$event_outcome_class,
      .data$random_value,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(fold_id = ((dplyr::row_number() - 1L) %% k) + 1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(event_id, fold_id)
}

calculate_brier <- function(observed, predicted) {
  mean((observed - predicted)^2)
}

calculate_log_loss <- function(observed, predicted) {
  predicted <- clip_probability(predicted)
  -mean(
    observed * log(predicted) +
      (1 - observed) * log(1 - predicted)
  )
}

calculate_calibration_intercept <- function(observed, predicted) {
  if (dplyr::n_distinct(observed) < 2L) return(NA_real_)
  linear_predictor <- stats::qlogis(clip_probability(predicted))
  fitted <- tryCatch(
    suppressWarnings(
      stats::glm(
        observed ~ 1,
        offset = linear_predictor,
        family = stats::binomial(link = "logit")
      )
    ),
    error = function(e) NULL
  )
  if (is.null(fitted)) return(NA_real_)
  as.numeric(stats::coef(fitted)[[1]])
}

calculate_calibration_slope <- function(observed, predicted) {
  if (dplyr::n_distinct(observed) < 2L) return(NA_real_)
  linear_predictor <- stats::qlogis(clip_probability(predicted))
  fitted <- tryCatch(
    suppressWarnings(
      stats::glm(
        observed ~ linear_predictor,
        family = stats::binomial(link = "logit")
      )
    ),
    error = function(e) NULL
  )
  if (is.null(fitted)) return(NA_real_)
  value <- stats::coef(fitted)[["linear_predictor"]]
  if (is.null(value) || !is.finite(value)) return(NA_real_)
  as.numeric(value)
}

standard_error <- function(x) {
  valid <- x[is.finite(x)]
  if (length(valid) <= 1L) return(NA_real_)
  stats::sd(valid) / sqrt(length(valid))
}

prepare_fixed_effect_information <- function(model) {
  list(
    mean = fixed_effect_vector(model),
    covariance = fixed_effect_covariance(model)
  )
}

draw_fixed_effects <- function(fixed_information) {
  covariance <- as.matrix(fixed_information$covariance)
  covariance <- (covariance + t(covariance)) / 2

  draw <- tryCatch(
    suppressWarnings(
      MASS::mvrnorm(
        n = 1L,
        mu = fixed_information$mean,
        Sigma = covariance,
        tol = 1e-8
      )
    ),
    error = function(first_error) {
      adjusted_covariance <- as.matrix(
        Matrix::nearPD(
          covariance,
          corr = FALSE,
          keepDiag = TRUE
        )$mat
      )
      tryCatch(
        MASS::mvrnorm(
          n = 1L,
          mu = fixed_information$mean,
          Sigma = adjusted_covariance,
          tol = 1e-8
        ),
        error = function(second_error) {
          stop(
            "Could not draw calibration fixed effects from the fitted ",
            "variance-covariance matrix. Initial error: ",
            conditionMessage(first_error),
            "; near-positive-definite retry: ",
            conditionMessage(second_error)
          )
        }
      )
    }
  )

  draw <- as.numeric(draw)
  names(draw) <- names(fixed_information$mean)
  draw
}

extract_ranef_conditional_se <- function(random_table, column_index = 1L) {
  conditional_variance <- attr(random_table, "postVar")
  n_levels <- nrow(random_table)

  if (is.null(conditional_variance)) {
    return(rep(NA_real_, n_levels))
  }

  # lme4 returns either one 3-D array or, when the same grouping
  # factor appears in multiple independent random-effect terms, a
  # list of 3-D arrays. Handle both representations explicitly.
  if (is.list(conditional_variance)) {
    if (length(conditional_variance) < column_index) {
      return(rep(NA_real_, n_levels))
    }
    variance_array <- conditional_variance[[column_index]]
    if (length(dim(variance_array)) != 3L) {
      return(rep(NA_real_, n_levels))
    }
    diagonal_index <- if (dim(variance_array)[1] >= column_index) {
      column_index
    } else {
      1L
    }
    return(
      sqrt(pmax(variance_array[diagonal_index, diagonal_index, ], 0))
    )
  }

  if (length(dim(conditional_variance)) != 3L) {
    return(rep(NA_real_, n_levels))
  }
  if (dim(conditional_variance)[1] < column_index) {
    return(rep(NA_real_, n_levels))
  }

  sqrt(
    pmax(
      conditional_variance[column_index, column_index, ],
      0
    )
  )
}

prepare_random_intercept_information <- function(model) {
  if (!inherits(model, "merMod")) {
    stop("Random-intercept information requires a mixed model.")
  }
  random_table <- lme4::ranef(model, condVar = TRUE)$predicted_species
  random_se <- extract_ranef_conditional_se(
    random_table,
    column_index = 1L
  )
  random_se[!is.finite(random_se)] <- 0
  tibble::tibble(
    predicted_species = rownames(random_table),
    mean = as.numeric(random_table[, "(Intercept)"]),
    se = as.numeric(random_se)
  )
}

draw_species_random_intercepts <- function(random_information) {
  values <- stats::rnorm(
    nrow(random_information),
    mean = random_information$mean,
    sd = random_information$se
  )
  names(values) <- random_information$predicted_species
  values
}

linear_predictor_from_fixed <- function(coefficient_draw, data) {
  eta <- rep(0, nrow(data))
  if ("(Intercept)" %in% names(coefficient_draw)) {
    eta <- eta + coefficient_draw[["(Intercept)"]]
  }
  if ("logit_score" %in% names(coefficient_draw)) {
    eta <- eta + coefficient_draw[["logit_score"]] * data$logit_score
  }
  if ("actual_multi_event" %in% names(coefficient_draw)) {
    eta <- eta +
      coefficient_draw[["actual_multi_event"]] * data$actual_multi_event
  }
  eta
}

draw_m2_probabilities <- function(fixed_information, data) {
  fixed_draw <- draw_fixed_effects(fixed_information)
  clip_probability(
    stats::plogis(linear_predictor_from_fixed(fixed_draw, data))
  )
}

draw_m3_hybrid_probabilities <- function(
    m3_fixed_information,
    m2_fixed_information,
    random_information,
    data
) {
  m3_draw <- draw_fixed_effects(m3_fixed_information)
  m2_draw <- draw_fixed_effects(m2_fixed_information)
  random_draw <- draw_species_random_intercepts(random_information)
  species_names <- as.character(data$predicted_species)
  known <- species_names %in% names(random_draw)

  m3_eta <- linear_predictor_from_fixed(m3_draw, data)
  m3_eta[known] <- m3_eta[known] + random_draw[species_names[known]]
  m2_eta <- linear_predictor_from_fixed(m2_draw, data)
  eta <- ifelse(known, m3_eta, m2_eta)
  clip_probability(stats::plogis(eta))
}

hill_number_from_counts <- function(counts, q) {
  counts <- as.numeric(counts)
  counts <- counts[is.finite(counts) & counts > 0]
  if (length(counts) == 0L) return(NA_real_)
  relative <- counts / sum(counts)
  if (q == 0) return(as.numeric(length(relative)))
  if (q == 1) return(exp(-sum(relative * log(relative))))
  (sum(relative^q))^(1 / (1 - q))
}

monthly_diversity_from_matrix <- function(
    count_matrix,
    month_names,
    q_values = c(0, 1, 2)
) {
  purrr::map_dfr(
    seq_len(nrow(count_matrix)),
    function(month_index) {
      purrr::map_dfr(
        q_values,
        function(q_value) {
          tibble::tibble(
            year_month = month_names[[month_index]],
            q = q_value,
            diversity = hill_number_from_counts(
              count_matrix[month_index, ],
              q_value
            ),
            accepted_event_species = sum(count_matrix[month_index, ])
          )
        }
      )
    }
  )
}

partition_hill_diversity <- function(count_matrix, q, weighting) {
  month_totals <- rowSums(count_matrix)
  valid_months <- month_totals > 0
  if (!any(valid_months)) {
    return(
      tibble::tibble(
        alpha = NA_real_,
        gamma = NA_real_,
        beta = NA_real_,
        n_months_used = 0L
      )
    )
  }

  counts <- count_matrix[valid_months, , drop = FALSE]
  totals <- month_totals[valid_months]
  within_relative <- sweep(counts, 1, totals, "/")

  if (weighting == "equal_month") {
    month_weights <- rep(1 / nrow(within_relative), nrow(within_relative))
  } else if (weighting == "abundance_weighted") {
    month_weights <- totals / sum(totals)
  } else {
    stop("Unknown weighting: ", weighting)
  }

  pooled_relative <- colSums(within_relative * month_weights)
  pooled_relative <- pooled_relative[pooled_relative > 0]

  if (q == 0) {
    alpha <- sum(month_weights * rowSums(counts > 0))
    gamma <- length(pooled_relative)
  } else if (q == 1) {
    within_entropy <- -rowSums(
      ifelse(
        within_relative > 0,
        within_relative * log(within_relative),
        0
      )
    )
    alpha <- exp(sum(month_weights * within_entropy))
    gamma <- exp(-sum(pooled_relative * log(pooled_relative)))
  } else {
    alpha <- (
      sum(month_weights * rowSums(within_relative^q))
    )^(1 / (1 - q))
    gamma <- (sum(pooled_relative^q))^(1 / (1 - q))
  }

  beta <- if (is.finite(alpha) && alpha > 0) gamma / alpha else NA_real_
  tibble::tibble(
    alpha = alpha,
    gamma = gamma,
    beta = beta,
    n_months_used = sum(valid_months)
  )
}

summarise_distribution <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(
      tibble::tibble(
        mean = NA_real_,
        sd = NA_real_,
        median = NA_real_,
        lower_95 = NA_real_,
        upper_95 = NA_real_
      )
    )
  }
  tibble::tibble(
    mean = mean(values),
    sd = stats::sd(values),
    median = stats::median(values),
    lower_95 = as.numeric(
      stats::quantile(values, 0.025, names = FALSE, type = 8)
    ),
    upper_95 = as.numeric(
      stats::quantile(values, 0.975, names = FALSE, type = 8)
    )
  )
}

build_candidate_aggregation <- function(data) {
  month_levels <- sort(unique(as.character(data$year_month)))
  species_levels <- sort(unique(as.character(data$predicted_species)))
  month_index <- match(as.character(data$year_month), month_levels)
  species_index <- match(as.character(data$predicted_species), species_levels)
  n_months <- length(month_levels)
  n_species <- length(species_levels)
  cell_index <- (month_index - 1L) * n_species + species_index

  aggregation_matrix <- Matrix::sparseMatrix(
    i = cell_index,
    j = seq_len(nrow(data)),
    x = 1L,
    dims = c(n_months * n_species, nrow(data))
  )

  list(
    aggregation_matrix = aggregation_matrix,
    month_levels = month_levels,
    species_levels = species_levels,
    n_months = n_months,
    n_species = n_species
  )
}

aggregate_candidate_states <- function(candidate_states, aggregation) {
  count_vector <- as.integer(
    aggregation$aggregation_matrix %*% as.integer(candidate_states)
  )
  matrix(
    count_vector,
    nrow = aggregation$n_months,
    ncol = aggregation$n_species,
    byrow = TRUE,
    dimnames = list(
      aggregation$month_levels,
      aggregation$species_levels
    )
  )
}

build_species_month_count_matrix <- function(
    data,
    month_levels = NULL,
    species_levels = NULL
) {
  if (is.null(month_levels)) {
    month_levels <- sort(unique(as.character(data$year_month)))
  }
  if (is.null(species_levels)) {
    species_levels <- sort(unique(as.character(data$predicted_species)))
  }

  counts <- data %>%
    dplyr::count(.data$year_month, .data$predicted_species, name = "n") %>%
    tidyr::complete(
      year_month = month_levels,
      predicted_species = species_levels,
      fill = list(n = 0L)
    ) %>%
    dplyr::mutate(
      year_month = factor(.data$year_month, levels = month_levels),
      predicted_species = factor(
        .data$predicted_species,
        levels = species_levels
      )
    ) %>%
    dplyr::arrange(.data$year_month, .data$predicted_species)

  matrix(
    counts$n,
    nrow = length(month_levels),
    ncol = length(species_levels),
    byrow = TRUE,
    dimnames = list(month_levels, species_levels)
  )
}

q_label <- function(q) {
  dplyr::case_when(
    q == 0 ~ "q = 0: Richness",
    q == 1 ~ "q = 1: Exponential Shannon",
    q == 2 ~ "q = 2: Inverse Simpson",
    TRUE ~ paste0("q = ", q)
  )
}

component_label <- function(component) {
  dplyr::recode(
    as.character(component),
    alpha = "Alpha",
    gamma = "Gamma",
    beta = "Beta",
    .default = as.character(component)
  )
}


english_month_year_labels <- function(x, line_break = TRUE) {
  # Locale-independent English month labels for manuscript figures.
  # ggplot2::scale_x_date(date_labels = "%b") follows LC_TIME and can
  # therefore produce non-English month names on the same analysis.
  x <- as.Date(x)
  month_index <- suppressWarnings(as.integer(format(x, "%m")))
  year_label <- format(x, "%Y")
  separator <- if (isTRUE(line_break)) "\n" else " "
  output <- rep(NA_character_, length(x))
  valid <- !is.na(x) & !is.na(month_index)
  output[valid] <- paste0(month.abb[month_index[valid]], separator, year_label[valid])
  output
}

truncate_identifier <- function(x, prefix = 10L, suffix = 5L) {
  x <- as.character(x)
  vapply(
    x,
    function(value) {
      if (is.na(value) || !nzchar(value)) return(NA_character_)
      n <- nchar(value)
      if (n <= prefix + suffix + 1L) return(value)
      paste0(substr(value, 1L, prefix), "...", substr(value, n - suffix + 1L, n))
    },
    FUN.VALUE = character(1)
  )
}

copy_if_exists <- function(from, to) {
  if (!file.exists(from)) return(FALSE)
  ensure_dir(dirname(to))
  isTRUE(file.copy(from, to, overwrite = TRUE))
}

silwood_palette <- function() {
  list(
    q = c(
      "q = 0: Richness" = "#0072B2",
      "q = 1: Exponential Shannon" = "#009E73",
      "q = 2: Inverse Simpson" = "#D55E00"
    ),
    model = c(
      "M1: Score" = "#999999",
      "M2: Score + event" = "#E69F00",
      "M3: Species-aware" = "#0072B2",
      "M4: Random score slope" = "#CC79A7"
    ),
    scenario = c(
      "Primary M3: species-aware" = "#0072B2",
      "Primary M2: pooled" = "#E69F00",
      "Audio-only sensitivity" = "#009E73"
    ),
    method = c(
      "Primary M3: species-aware" = "#0072B2",
      "Primary M2: pooled" = "#56B4E9",
      "Audio-only sensitivity" = "#009E73",
      "Naive: all candidates" = "#D55E00",
      "Hard top-1" = "#CC79A7",
      "Threshold >= 0.50" = "#E69F00",
      "Threshold >= 0.70" = "#A6761D",
      "Threshold >= 0.90" = "#666666"
    ),
    component = c(
      "Alpha" = "#0072B2",
      "Gamma" = "#D55E00",
      "Beta" = "#009E73"
    ),
    event = c(
      "Single-candidate event" = "#4C78A8",
      "Multi-candidate event" = "#F58518"
    )
  )
}

theme_silwood <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = base_size - 2,
        colour = "grey30"
      ),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

save_figure <- function(plot, file_stem, width, height, dpi = 400) {
  ensure_dir(dirname(file_stem))
  png_file <- paste0(file_stem, ".png")
  pdf_file <- paste0(file_stem, ".pdf")

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      device = ragg::agg_png
    )
  } else {
    ggplot2::ggsave(
      png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }

  pdf_device <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    "pdf"
  }
  ggplot2::ggsave(
    pdf_file,
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    device = pdf_device
  )

  invisible(c(png_file, pdf_file))
}

parse_time_safely <- function(x, timezone = "UTC") {
  if (inherits(x, "POSIXct")) return(x)
  parsed <- suppressWarnings(
    lubridate::ymd_hms(x, quiet = TRUE, tz = timezone)
  )
  if (all(is.na(parsed))) stop("Timestamp parsing failed.")
  parsed
}

format_p_value <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "p = NA",
    p < 0.001 ~ "p < 0.001",
    TRUE ~ paste0("p = ", formatC(p, digits = 3, format = "f"))
  )
}

# ============================================================
# 27_collect_final_outputs.R
#
# Purpose:
# Collect the exact figures and tables used in the thesis into one
# clearly indexed folder. Analytical quantities are read from the
# step-specific results; no models are refitted here.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)

main_figure_dir <- silwood_final_dir(paths, "main_figures")
main_table_dir <- silwood_final_dir(paths, "main_tables")

# Figure 1 is a static conceptual workflow, not a statistical plot.
figure1_source <- file.path(
  paths$static_asset_root,
  "Figure_01_validation_calibrated_framework.png"
)
if (!file.exists(figure1_source)) stop("Static Figure 1 not found: ", figure1_source)
file.copy(
  figure1_source,
  file.path(main_figure_dir, basename(figure1_source)),
  overwrite = TRUE
)

# Table 1: notation.
notation <- tibble::tribble(
  ~Symbol, ~Definition,
  "r", "Event-species candidate record.",
  "j", "Sampled month.",
  "s", "Predicted species.",
  "b", "Monte Carlo simulation.",
  "q", "Hill diversity order.",
  "Y_r", "Manually validated correctness outcome for candidate r.",
  "x_r", "Original BirdNET confidence score.",
  "x_r*", "Logit-transformed BirdNET confidence score.",
  "m_r", "Indicator that candidate r came from a multi-candidate acoustic event.",
  "p_r", "Candidate-level correctness probability for candidate r.",
  "Z_r^(b)", "Simulated correctness state of candidate r in simulation b.",
  "R_js", "Set of records predicted as species s in month j.",
  "n_js^(b)", "Retained candidate count for species s in month j.",
  "pi_js^(b)", "Relative candidate-record frequency.",
  "w_j", "Weight assigned to month j.",
  "J", "Number of sampled months included in a partition."
)
readr::write_csv(notation, file.path(main_table_dir, "Table_01_mathematical_notation.csv"))

# Table 2: model structure + event-grouped validation.
model_fit <- readr::read_csv(
  file.path(paths$results_root, "02_calibration", "tables", "calibration_model_comparison.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(.data$analysis_set == "Primary", .data$model_id %in% c("M1","M2","M3","M4"))

cv <- readr::read_csv(
  file.path(paths$results_root, "03_model_validation", "tables", "event_cv_model_summary.csv"),
  show_col_types = FALSE
)

roles <- tibble::tribble(
  ~model_id, ~Model_role,
  "M1", "Score-only baseline",
  "M2", "Pooled probability sensitivity / fallback",
  "M3", "Selected for downstream propagation",
  "M4", "Random-slope structural sensitivity"
)

table2 <- model_fit %>%
  dplyr::left_join(cv, by = c("model_id", "model_label")) %>%
  dplyr::left_join(roles, by = "model_id") %>%
  dplyr::transmute(
    Model = .data$model_label,
    Model_role = .data$Model_role,
    AIC = round(.data$aic, 2),
    Delta_AIC = round(.data$delta_aic, 2),
    Species_intercept_SD = round(.data$species_intercept_sd, 3),
    Species_score_slope_SD = round(.data$species_score_slope_sd, 3),
    Mean_Brier = round(.data$mean_brier_score, 4),
    SE_Brier = round(.data$se_brier_score, 4),
    Mean_log_loss = round(.data$mean_log_loss, 4),
    SE_log_loss = round(.data$se_log_loss, 4),
    Singular_CV_folds = .data$singular_folds,
    Maximum_CV_gradient = .data$maximum_gradient_across_folds
  ) %>%
  dplyr::arrange(match(.data$Model, c("M1: Score","M2: Score + event","M3: Species-aware","M4: Random score slope")))
readr::write_csv(table2, file.path(main_table_dir, "Table_02_calibration_model_selection.csv"), na = "")

# Table 3: Primary species-aware equal-month partition.
partition <- readr::read_csv(
  file.path(paths$results_root, "05_monte_carlo", "tables", "monte_carlo_partition_summary.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(.data$scenario == "Primary_M3", .data$weighting == "equal_month") %>%
  dplyr::mutate(
    component_order = match(.data$component, c("alpha","beta","gamma")),
    q_order = match(.data$q, c(0,1,2))
  ) %>%
  dplyr::arrange(.data$q_order, .data$component_order) %>%
  dplyr::transmute(
    Diversity_order = .data$q_label,
    Component = .data$component_label,
    Median = round(.data$median, 3),
    Lower_95_simulation_interval = round(.data$lower_95, 3),
    Upper_95_simulation_interval = round(.data$upper_95, 3)
  )
readr::write_csv(partition, file.path(main_table_dir, "Table_03_primary_species_aware_diversity.csv"), na = "")

# Table 4: analytical benchmark.
analytical <- readr::read_csv(
  file.path(paths$results_root, "09_analytical_benchmark", "tables", "analytical_vs_monte_carlo_partition.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(.data$weighting == "equal_month") %>%
  dplyr::mutate(
    component_order = match(.data$component, c("alpha","beta","gamma")),
    q_order = match(.data$q, c(0,1,2))
  ) %>%
  dplyr::arrange(.data$q_order, .data$component_order) %>%
  dplyr::transmute(
    Diversity_order = .data$q_label,
    Component = .data$component_label,
    Analytical_or_plugin = round(.data$analytical_estimate, 3),
    Fixed_MC_mean = round(.data$fixed_mean, 3),
    Fixed_MC_lower_95 = round(.data$fixed_lower_95, 3),
    Fixed_MC_upper_95 = round(.data$fixed_upper_95, 3),
    Full_MC_median = round(.data$full_median, 3),
    Full_MC_lower_95 = round(.data$full_lower_95, 3),
    Full_MC_upper_95 = round(.data$full_upper_95, 3),
    Analytical_percent_error_vs_fixed_MC = round(.data$analytical_percent_error, 2)
  )
readr::write_csv(analytical, file.path(main_table_dir, "Table_04_analytical_benchmark.csv"), na = "")

# Indexes for rapid navigation.
figure_files <- list.files(paths$final_output_root, pattern = "\\.(pdf|png)$", recursive = TRUE, full.names = TRUE)
table_files <- list.files(paths$final_output_root, pattern = "\\.(csv|tex)$", recursive = TRUE, full.names = TRUE)

make_index <- function(files, type) {
  tibble::tibble(
    type = type,
    relative_path = gsub("\\\\", "/", substring(files, nchar(paths$final_output_root) + 2L)),
    file_name = basename(files),
    size_bytes = file.info(files)$size
  ) %>% dplyr::arrange(.data$relative_path)
}

readr::write_csv(
  dplyr::bind_rows(make_index(figure_files, "figure"), make_index(table_files, "table")),
  file.path(paths$final_output_root, "FINAL_OUTPUT_INDEX.csv")
)

message("Final thesis outputs collected in: ", paths$final_output_root)

# ============================================================
# run_all.R
# One-command launcher for the final, reproducible Silwood pipeline.
#
# From the project root:
#   Rscript run_all.R
#
# Default final settings:
#   1,000 Monte Carlo simulations
#   5-fold event-grouped cross-validation repeated 5 times
#   frozen final manual-validation decisions
#   no LOSO analysis (not used in the final thesis)
#
# Optional environment variables:
#   SILWOOD_MODE=quick|final
#   SILWOOD_N_SIMULATIONS=1000
#   SILWOOD_N_FOLDS=5
#   SILWOOD_N_REPEATS=5
#   SILWOOD_REBUILD_VALIDATION_SAMPLING=true|false
#   SILWOOD_RUN_LOSO=true|false
#   SILWOOD_BUILD_REPORTING=true|false
#   SILWOOD_START_AT=1
#   SILWOOD_END_AT=27
#   SILWOOD_CLEAN_RESULTS=true|false
# ============================================================

if (!file.exists(file.path("scripts", "00_setup", "00_pipeline_functions.R"))) {
  stop("run_all.R must be executed from the final_project root.")
}

source(file.path("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()
options(readr.show_progress = FALSE)

config <- silwood_config()
paths <- silwood_paths(config)
ensure_dir(paths$results_root)
ensure_dir(paths$processed_root)
ensure_dir(paths$log_root)

is_true <- function(name, default = "false") {
  tolower(Sys.getenv(name, unset = default)) %in% c("true", "1", "yes", "y")
}
parse_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = "")))
  if (is.na(value)) default else value
}

rebuild_sampling <- is_true("SILWOOD_REBUILD_VALIDATION_SAMPLING")
run_loso <- is_true("SILWOOD_RUN_LOSO")
build_reporting <- is_true("SILWOOD_BUILD_REPORTING", "true")
clean_results <- is_true("SILWOOD_CLEAN_RESULTS")
start_at <- parse_integer("SILWOOD_START_AT", 1L)
end_at <- parse_integer("SILWOOD_END_AT", 27L)

if (start_at > end_at) stop("SILWOOD_START_AT cannot exceed SILWOOD_END_AT.")

if (clean_results && start_at <= 1L) {
  generated_steps <- file.path(paths$results_root, sprintf("%02d", 0:10))
  generated_steps <- c(
    file.path(paths$results_root, "00_data_preparation"),
    file.path(paths$results_root, "01_validation"),
    file.path(paths$results_root, "02_calibration"),
    file.path(paths$results_root, "03_model_validation"),
    file.path(paths$results_root, "04_probability_assignment"),
    file.path(paths$results_root, "05_monte_carlo"),
    file.path(paths$results_root, "06_method_comparison"),
    file.path(paths$results_root, "07_sampling_coverage"),
    file.path(paths$results_root, "08_coverage_robustness"),
    file.path(paths$results_root, "09_analytical_benchmark"),
    file.path(paths$results_root, "10_final_outputs")
  )
  unlink(generated_steps, recursive = TRUE, force = TRUE)
}

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_log_dir <- ensure_dir(file.path(paths$log_root, paste0("run_", run_id)))
master_log <- file.path(run_log_dir, "pipeline_log.txt")
manifest_file <- file.path(run_log_dir, "pipeline_manifest.csv")

log_message <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(line, "\n")
  cat(line, "\n", file = master_log, append = TRUE)
}

steps <- tibble::tribble(
  ~step, ~name, ~script, ~default_run,
  1L, "Clean raw BirdNET export", "scripts/01_data_preparation/01_clean_birdnet_data.R", TRUE,
  2L, "Rebuild Round-1 sample design", "scripts/02_manual_validation/02_generate_round1_validation_sample.R", rebuild_sampling,
  10L, "Build Round-2 species sampling frame", "scripts/02_manual_validation/10_build_species_sampling_frame.R", rebuild_sampling,
  11L, "Generate Round-2 screening sample", "scripts/02_manual_validation/11_generate_round2_screening_sample.R", rebuild_sampling,
  12L, "Prepare frozen final validation data", "scripts/02_manual_validation/12_prepare_final_validation_data.R", TRUE,
  13L, "Fit calibration models", "scripts/03_calibration/13_fit_final_calibration_models.R", TRUE,
  14L, "Run event-grouped cross-validation", "scripts/03_calibration/14_event_grouped_cross_validation.R", TRUE,
  15L, "Optional leave-one-species-out validation", "scripts/99_optional_diagnostics/15_leave_one_species_out_validation.R", run_loso,
  16L, "Assign calibrated probabilities", "scripts/03_calibration/16_apply_calibration_to_full_data.R", TRUE,
  17L, "Run full Monte Carlo propagation", "scripts/04_uncertainty_propagation/17_monte_carlo_diversity_simulation_corrected.R", TRUE,
  18L, "Compare identification workflows", "scripts/05_method_comparison/18_compare_diversity_methods_final.R", TRUE,
  19L, "Calculate sampling-coverage proxies", "scripts/06_sampling_coverage/19_sampling_effort_proxy_sensitivity_fixed5min.R", TRUE,
  20L, "Run higher-coverage robustness analysis", "scripts/06_sampling_coverage/20_low_coverage_month_sensitivity.R", TRUE,
  21L, "Build analytical benchmark", "scripts/07_analytical_benchmark/21_analytical_vs_monte_carlo.R", TRUE,
  22L, "Rebuild manuscript Figure 2", "scripts/08_reporting/22_rebuild_manuscript_figure2.R", build_reporting,
  23L, "Rebuild manuscript Figure 3", "scripts/08_reporting/23_rebuild_manuscript_figure3.R", build_reporting,
  24L, "Rebuild manuscript Figure 4", "scripts/08_reporting/24_rebuild_manuscript_figure4.R", build_reporting,
  25L, "Rebuild manuscript Figure 5", "scripts/08_reporting/25_rebuild_manuscript_figure5.R", build_reporting,
  26L, "Rebuild supplementary figures and tables", "scripts/08_reporting/26_rebuild_supplementary_outputs.R", build_reporting,
  27L, "Collect final thesis outputs", "scripts/08_reporting/27_collect_final_outputs.R", build_reporting
) %>%
  dplyr::filter(.data$default_run, .data$step >= start_at, .data$step <= end_at)

if (nrow(steps) == 0L) stop("No pipeline steps selected.")

manifest <- vector("list", nrow(steps))
log_message("Silwood final pipeline started. Mode=", config$mode,
            "; simulations=", config$n_simulations,
            "; CV=", config$n_folds, "x", config$n_repeats,
            "; run_id=", run_id)

for (i in seq_len(nrow(steps))) {
  step <- steps$step[[i]]
  name <- steps$name[[i]]
  script <- steps$script[[i]]
  if (!file.exists(script)) stop("Missing pipeline script: ", script)

  step_log <- file.path(run_log_dir, sprintf("step_%02d_%s.log", step, gsub("[^A-Za-z0-9]+", "_", tolower(name))))
  started <- Sys.time()
  log_message("START step ", step, ": ", name)

  output_depth <- sink.number(type = "output")
  message_connection <- sink.number(type = "message")
  con <- file(step_log, open = "wt")
  sink(con, type = "output", split = TRUE)
  sink(con, type = "message")

  result <- tryCatch(
    {
      source(script, local = new.env(parent = globalenv()), echo = FALSE, chdir = FALSE)
      list(ok = TRUE, error = NA_character_)
    },
    error = function(e) list(ok = FALSE, error = conditionMessage(e))
  )

  if (sink.number(type = "message") != message_connection) sink(type = "message")
  while (sink.number(type = "output") > output_depth) sink(type = "output")
  close(con)

  completed <- Sys.time()
  manifest[[i]] <- tibble::tibble(
    step = step,
    name = name,
    script = script,
    log_file = step_log,
    started = as.character(started),
    completed = as.character(completed),
    elapsed_seconds = as.numeric(difftime(completed, started, units = "secs")),
    status = if (result$ok) "completed" else "failed",
    error = result$error
  )
  readr::write_csv(dplyr::bind_rows(manifest), manifest_file, na = "")

  if (!result$ok) {
    log_message("FAILED step ", step, ": ", result$error)
    stop("Pipeline failed at step ", step, ". See ", step_log)
  }
  log_message("DONE step ", step, " in ", round(manifest[[i]]$elapsed_seconds, 1), " seconds")
}

session_file <- file.path(run_log_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)

log_message("Pipeline completed successfully. Final outputs: ", paths$final_output_root)

# Stable pointers to the latest successful run.
file.copy(master_log, file.path(paths$log_root, "latest_pipeline_log.txt"), overwrite = TRUE)
file.copy(manifest_file, file.path(paths$log_root, "latest_pipeline_manifest.csv"), overwrite = TRUE)
file.copy(session_file, file.path(paths$log_root, "latest_sessionInfo.txt"), overwrite = TRUE)

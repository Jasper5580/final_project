# ============================================================
# 12_prepare_final_validation_data.R
#
# Purpose:
# Prepare the frozen, manually reviewed validation dataset for
# the streamlined calibration pipeline.
#
# This script does not repeat manual listening. It reads the
# final combined validation file, checks the recorded decisions,
# derives model variables, and saves analysis-ready data.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)

input_file <- paths$validation_input
processed_dir <- ensure_dir(
  file.path(paths$processed_root, "validation")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "01_validation", "tables")
)

if (!file.exists(input_file)) {
  stop("Final combined validation file not found: ", input_file)
}

validation_raw <- readr::read_csv(
  input_file,
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

required_columns <- c(
  "validation_id",
  "validation_round",
  "id",
  "event_id",
  "predicted_species",
  "birdnet_score",
  "score_bin",
  "n_species_predicted",
  "actual_multi_event",
  "prediction_correct",
  "uncertain",
  "multiple_species",
  "validation_basis"
)
assert_columns(validation_raw, required_columns, "combined validation data")

validation <- validation_raw %>%
  dplyr::mutate(
    validation_id = as.character(.data$validation_id),
    id = as.character(.data$id),
    event_id = as.character(.data$event_id),
    audio_id = as.character(.data$audio_id),
    predicted_species = as.character(.data$predicted_species),
    validation_round = as.integer(.data$validation_round),
    birdnet_score = as.numeric(.data$birdnet_score),
    n_species_predicted = as.integer(.data$n_species_predicted),
    actual_multi_event = as_binary_integer(.data$actual_multi_event),
    prediction_correct = as_binary_integer(.data$prediction_correct),
    uncertain = as_binary_integer(.data$uncertain),
    multiple_species = as_binary_integer(.data$multiple_species),
    validation_basis = as.character(.data$validation_basis)
  )

assert_unique(validation$id, "Validation detection IDs")
assert_unique(validation$validation_id, "Validation IDs")

valid_basis <- c("audio_based", "audio_plus_ecological_prior")
invalid_basis_n <- sum(
  is.na(validation$validation_basis) |
    !validation$validation_basis %in% valid_basis
)
if (invalid_basis_n > 0L) {
  stop("Invalid validation_basis values found in ", invalid_basis_n, " rows.")
}

invalid_outcome_n <- sum(
  is.na(validation$prediction_correct) |
    !validation$prediction_correct %in% c(0L, 1L)
)
if (invalid_outcome_n > 0L) {
  stop("Invalid prediction_correct values found in ", invalid_outcome_n, " rows.")
}

invalid_uncertain_n <- sum(
  is.na(validation$uncertain) |
    !validation$uncertain %in% c(0L, 1L)
)
if (invalid_uncertain_n > 0L) {
  stop("Invalid uncertain values found in ", invalid_uncertain_n, " rows.")
}

event_field_mismatch_n <- sum(
  validation$actual_multi_event !=
    as.integer(validation$n_species_predicted > 1L),
  na.rm = TRUE
)
if (event_field_mismatch_n > 0L) {
  stop(
    "actual_multi_event disagrees with n_species_predicted for ",
    event_field_mismatch_n,
    " records."
  )
}

validation <- validation %>%
  dplyr::mutate(
    score_bounded = bound_score(.data$birdnet_score),
    logit_score = stats::qlogis(.data$score_bounded),
    analysis_primary_include = as.integer(
      .data$uncertain == 0L &
        .data$prediction_correct %in% c(0L, 1L)
    ),
    analysis_audio_only_sensitivity_include = as.integer(
      .data$analysis_primary_include == 1L &
        .data$validation_basis == "audio_based"
    ),
    analysis_ecological_prior_include = as.integer(
      .data$analysis_primary_include == 1L &
        .data$validation_basis == "audio_plus_ecological_prior"
    ),
    human_multiple_species_heard = .data$multiple_species
  )

primary_data <- validation %>%
  dplyr::filter(.data$analysis_primary_include == 1L)

species_summary <- primary_data %>%
  dplyr::group_by(.data$predicted_species) %>%
  dplyr::summarise(
    species_validation_n = dplyr::n(),
    species_validation_events = dplyr::n_distinct(.data$event_id),
    species_correct_n = sum(.data$prediction_correct == 1L),
    species_incorrect_n = sum(.data$prediction_correct == 0L),
    species_confirmation_rate = mean(.data$prediction_correct),
    species_mean_score = mean(.data$birdnet_score),
    species_min_score = min(.data$birdnet_score),
    species_max_score = max(.data$birdnet_score),
    species_score_bin_n = dplyr::n_distinct(
      .data$score_bin[!is.na(.data$score_bin)]
    ),
    species_multi_event_n = sum(.data$actual_multi_event == 1L),
    species_audio_only_n = sum(.data$validation_basis == "audio_based"),
    species_ecological_prior_n = sum(
      .data$validation_basis == "audio_plus_ecological_prior"
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    validation_depth = dplyr::case_when(
      .data$species_validation_n <= 3L ~ "n = 1-3",
      .data$species_validation_n <= 7L ~ "n = 4-7",
      .data$species_validation_n <= 19L ~ "n = 8-19",
      TRUE ~ "n >= 20"
    ),
    random_slope_support = dplyr::case_when(
      .data$species_validation_n >= 20L &
        .data$species_score_bin_n >= 3L ~ "reasonable",
      .data$species_validation_n >= 8L &
        .data$species_score_bin_n >= 3L ~ "limited",
      TRUE ~ "insufficient"
    )
  )

replace_columns <- setdiff(names(species_summary), "predicted_species")
validation_analysis_ready <- validation %>%
  dplyr::select(-dplyr::any_of(replace_columns)) %>%
  dplyr::left_join(species_summary, by = "predicted_species") %>%
  dplyr::arrange(
    .data$validation_round,
    .data$predicted_species,
    .data$validation_id
  )

quality_checks <- tibble::tibble(
  check = c(
    "Total validation records",
    "Primary analysis records",
    "Audio-only sensitivity records",
    "Ecological-prior records",
    "Predicted species labels",
    "Unique acoustic events",
    "Round 1 records",
    "Round 2 records",
    "Duplicate detection IDs",
    "Duplicate validation IDs",
    "Event-field mismatches"
  ),
  value = c(
    nrow(validation_analysis_ready),
    sum(validation_analysis_ready$analysis_primary_include == 1L),
    sum(
      validation_analysis_ready$analysis_audio_only_sensitivity_include == 1L
    ),
    sum(validation_analysis_ready$analysis_ecological_prior_include == 1L),
    dplyr::n_distinct(validation_analysis_ready$predicted_species),
    dplyr::n_distinct(validation_analysis_ready$event_id),
    sum(validation_analysis_ready$validation_round == 1L),
    sum(validation_analysis_ready$validation_round == 2L),
    sum(duplicated(validation_analysis_ready$id)),
    sum(duplicated(validation_analysis_ready$validation_id)),
    event_field_mismatch_n
  ),
  expected = c(
    1147, 1147, 1016, 131, 229, 1138, 439, 708, 0, 0, 0
  )
)

for (row_index in seq_len(nrow(quality_checks))) {
  assert_expected(
    quality_checks$value[[row_index]],
    quality_checks$expected[[row_index]],
    quality_checks$check[[row_index]]
  )
}

round_summary <- validation_analysis_ready %>%
  dplyr::group_by(.data$validation_round) %>%
  dplyr::summarise(
    n_records = dplyr::n(),
    n_species = dplyr::n_distinct(.data$predicted_species),
    n_events = dplyr::n_distinct(.data$event_id),
    score_min = min(.data$birdnet_score),
    score_max = max(.data$birdnet_score),
    .groups = "drop"
  )

basis_summary <- validation_analysis_ready %>%
  dplyr::count(.data$validation_basis, name = "n_records") %>%
  dplyr::mutate(proportion = .data$n_records / sum(.data$n_records))

example_audio_correct <- validation_analysis_ready %>%
  dplyr::filter(
    .data$validation_basis == "audio_based",
    .data$prediction_correct == 1L,
    .data$actual_multi_event == 0L
  ) %>%
  dplyr::slice_head(n = 1)

example_audio_incorrect_multi <- validation_analysis_ready %>%
  dplyr::filter(
    .data$validation_basis == "audio_based",
    .data$prediction_correct == 0L,
    .data$actual_multi_event == 1L
  ) %>%
  dplyr::slice_head(n = 1)

example_ecological_prior <- validation_analysis_ready %>%
  dplyr::filter(
    .data$validation_basis == "audio_plus_ecological_prior",
    .data$prediction_correct == 0L
  ) %>%
  dplyr::arrange(
    dplyr::desc(.data$predicted_species == "Red-billed Chough")
  ) %>%
  dplyr::slice_head(n = 1)

methods_examples <- dplyr::bind_rows(
  example_audio_correct,
  example_audio_incorrect_multi,
  example_ecological_prior
) %>%
  dplyr::transmute(
    validation_id,
    validation_round,
    predicted_species,
    birdnet_score,
    actual_multi_event,
    prediction_correct,
    validation_basis,
    analysis_primary_include,
    analysis_audio_only_sensitivity_include
  )

analysis_ready_csv <- file.path(
  processed_dir,
  "validation_analysis_ready.csv"
)
analysis_ready_rds <- file.path(
  processed_dir,
  "validation_analysis_ready.rds"
)

readr::write_csv(validation_analysis_ready, analysis_ready_csv, na = "")
saveRDS(validation_analysis_ready, analysis_ready_rds)
readr::write_csv(
  quality_checks,
  file.path(table_dir, "validation_quality_checks.csv")
)
readr::write_csv(
  round_summary,
  file.path(table_dir, "validation_round_summary.csv")
)
readr::write_csv(
  basis_summary,
  file.path(table_dir, "validation_basis_summary.csv")
)
readr::write_csv(
  species_summary,
  file.path(table_dir, "validation_species_summary.csv")
)
readr::write_csv(
  methods_examples,
  file.path(table_dir, "methods_table_validation_examples.csv")
)

# ------------------------------------------------------------
# Methods Table 1: reproducible raw-to-cleaned illustration
# ------------------------------------------------------------
# This table is descriptive only.  It selects one single-candidate
# event and two candidate rows from one multi-candidate event, then
# joins them back to the original raw export by detection ID.
# Selection is deterministic so that the Methods example is regenerated
# automatically rather than copied manually into the thesis.
clean_for_methods <- read_rds_or_csv(
  paths$clean_rds,
  paths$clean_csv,
  "cleaned BirdNET detections for Methods Table 1"
)

raw_for_methods <- readr::read_csv(
  paths$raw_detection_csv,
  show_col_types = FALSE,
  progress = FALSE
)

assert_columns(
  clean_for_methods,
  c(
    "id", "event_id", "audio_id", "analysis", "species",
    "birdnet_score", "year_month", "start_secs", "end_secs",
    "duration_secs"
  ),
  "cleaned detections for Methods Table 1"
)
assert_columns(
  raw_for_methods,
  c(
    "id", "analysis", "tags", "confidence", "audio_id",
    "start_secs", "end_secs", "detected_time"
  ),
  "raw detections for Methods Table 1"
)

method_event_structure <- clean_for_methods %>%
  dplyr::count(.data$event_id, name = "candidate_labels_in_event")

clean_method_pool <- clean_for_methods %>%
  dplyr::transmute(
    id = as.character(.data$id),
    event_id = as.character(.data$event_id),
    audio_id = as.character(.data$audio_id),
    analysis = as.character(.data$analysis),
    predicted_species = as.character(.data$species),
    birdnet_score = as.numeric(.data$birdnet_score),
    year_month = as.character(.data$year_month),
    start_secs = as.numeric(.data$start_secs),
    end_secs = as.numeric(.data$end_secs),
    duration_secs = as.numeric(.data$duration_secs)
  ) %>%
  dplyr::left_join(method_event_structure, by = "event_id")

single_example <- clean_method_pool %>%
  dplyr::filter(.data$candidate_labels_in_event == 1L) %>%
  dplyr::mutate(distance_from_target = abs(.data$birdnet_score - 0.75)) %>%
  dplyr::arrange(.data$distance_from_target, .data$id) %>%
  dplyr::slice_head(n = 1L) %>%
  dplyr::select(-distance_from_target)

multi_event_choice <- clean_method_pool %>%
  dplyr::group_by(.data$event_id) %>%
  dplyr::summarise(
    candidate_labels_in_event = dplyr::n(),
    distinct_species = dplyr::n_distinct(.data$predicted_species),
    mean_score = mean(.data$birdnet_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    .data$candidate_labels_in_event >= 2L,
    .data$distinct_species >= 2L
  ) %>%
  dplyr::mutate(distance_from_target = abs(.data$mean_score - 0.70)) %>%
  dplyr::arrange(
    .data$candidate_labels_in_event,
    .data$distance_from_target,
    .data$event_id
  ) %>%
  dplyr::slice_head(n = 1L) %>%
  dplyr::pull(event_id)

if (length(multi_event_choice) != 1L) {
  stop("Could not identify a multi-candidate event for Methods Table 1.")
}

multi_examples <- clean_method_pool %>%
  dplyr::filter(.data$event_id == .env$multi_event_choice) %>%
  dplyr::arrange(dplyr::desc(.data$birdnet_score), .data$predicted_species, .data$id) %>%
  dplyr::slice_head(n = 2L)

methods_table1_clean_full <- dplyr::bind_rows(
  single_example,
  multi_examples
) %>%
  dplyr::mutate(example_order = dplyr::row_number())

assert_expected(
  nrow(methods_table1_clean_full),
  3L,
  "Methods Table 1 cleaned example rows"
)

methods_table1_raw_full <- raw_for_methods %>%
  dplyr::mutate(id = as.character(.data$id)) %>%
  dplyr::semi_join(
    methods_table1_clean_full %>% dplyr::select(id),
    by = "id"
  ) %>%
  dplyr::left_join(
    methods_table1_clean_full %>%
      dplyr::select(id, example_order),
    by = "id"
  ) %>%
  dplyr::arrange(.data$example_order)

assert_expected(
  nrow(methods_table1_raw_full),
  3L,
  "Methods Table 1 raw example rows"
)

methods_table1_panel_a <- methods_table1_raw_full %>%
  dplyr::transmute(
    Example = .data$example_order,
    Detection_ID = truncate_identifier(.data$id),
    Analysis = as.character(.data$analysis),
    Raw_species_tag = as.character(.data$tags),
    Confidence = round(as.numeric(.data$confidence), 6),
    Audio_ID = truncate_identifier(.data$audio_id),
    Start_s = as.numeric(.data$start_secs),
    End_s = as.numeric(.data$end_secs),
    Detection_time = as.character(.data$detected_time)
  )

methods_table1_panel_b <- methods_table1_clean_full %>%
  dplyr::transmute(
    Example = .data$example_order,
    Detection_ID = truncate_identifier(.data$id),
    Event_ID = truncate_identifier(.data$event_id, prefix = 14L, suffix = 8L),
    Predicted_species = .data$predicted_species,
    BirdNET_score = round(.data$birdnet_score, 6),
    Year_month = .data$year_month,
    Duration_s = .data$duration_secs,
    Candidate_labels_in_event = .data$candidate_labels_in_event
  )

methods_table1_audit <- methods_table1_clean_full %>%
  dplyr::left_join(
    methods_table1_raw_full %>%
      dplyr::select(
        id,
        raw_analysis = analysis,
        raw_species_tag = tags,
        raw_confidence = confidence,
        raw_detected_time = detected_time
      ),
    by = "id"
  )

readr::write_csv(
  methods_table1_panel_a,
  file.path(table_dir, "methods_table1_panelA_raw.csv"),
  na = ""
)
readr::write_csv(
  methods_table1_panel_b,
  file.path(table_dir, "methods_table1_panelB_cleaned.csv"),
  na = ""
)
readr::write_csv(
  methods_table1_audit,
  file.path(table_dir, "methods_table1_raw_to_cleaned_audit.csv"),
  na = ""
)

message("Validation preparation completed.")
message("Analysis-ready RDS: ", analysis_ready_rds)
message("Validation tables: ", table_dir)

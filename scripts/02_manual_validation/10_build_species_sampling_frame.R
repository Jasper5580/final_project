# ============================================================
# 10_build_species_sampling_frame.R
#
# Reproducible manual-validation sampling-design audit.
# This script does not recreate human listening decisions.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)

# ------------------------------------------------------------
# 1. Input and output paths
# ------------------------------------------------------------

clean_rds <- paths$clean_rds

clean_csv <- paths$clean_csv

round1_candidates <- c(
  here::here("data", "02_manual_validation", "round1_completed",
             "validation_sample_silwood_round1_completed.csv")
)

output_dir <- here::here("data", "02_manual_validation", "round2_design")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Read complete cleaned detections
# ------------------------------------------------------------

if (file.exists(clean_rds)) {
  detections <- readRDS(clean_rds)
} else if (file.exists(clean_csv)) {
  detections <- read_csv(
    clean_csv,
    show_col_types = FALSE
  )
} else {
  stop(
    "Could not find birdnet_detections_clean.rds or .csv."
  )
}

required_detection_columns <- c(
  "id",
  "event_id",
  "species",
  "birdnet_score",
  "year_month"
)

missing_detection_columns <- setdiff(
  required_detection_columns,
  names(detections)
)

if (length(missing_detection_columns) > 0) {
  stop(
    "Missing detection columns: ",
    paste(missing_detection_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. Read round-1 manual validation
# ------------------------------------------------------------

round1_files <- round1_candidates[
  file.exists(round1_candidates)
]

if (length(round1_files) == 0) {
  stop(
    paste0(
      "No round-1 validation file found. Expected either:\n",
      "data/02_manual_validation/round1_completed/",
      "validation_sample_silwood_round1.csv or ",
      "validation_sample_silwood(4).csv\n",
      "or data/02_manual_validation/final/validation_sample_silwood_combined_final.csv"
    )
  )
}

round1_file <- round1_files[[1]]

round1 <- read_csv(
  round1_file,
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

required_round1_columns <- c(
  "id",
  "predicted_species",
  "birdnet_score",
  "prediction_correct",
  "uncertain"
)

missing_round1_columns <- setdiff(
  required_round1_columns,
  names(round1)
)

if (length(missing_round1_columns) > 0) {
  stop(
    "Missing round-1 columns: ",
    paste(missing_round1_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 4. Standardise and derive event information
# ------------------------------------------------------------

score_bin_levels <- c(
  "0.45-0.60",
  "0.60-0.75",
  "0.75-0.90",
  "0.90-1.00"
)

detections <- detections %>%
  mutate(
    id = as.character(id),
    event_id = as.character(event_id),
    predicted_species = as.character(species),
    birdnet_score = as.numeric(birdnet_score),
    score_bin = case_when(
      birdnet_score >= 0.45 & birdnet_score < 0.60 ~
        "0.45-0.60",
      birdnet_score >= 0.60 & birdnet_score < 0.75 ~
        "0.60-0.75",
      birdnet_score >= 0.75 & birdnet_score < 0.90 ~
        "0.75-0.90",
      birdnet_score >= 0.90 & birdnet_score <= 1.00 ~
        "0.90-1.00",
      TRUE ~ NA_character_
    ),
    score_bin = factor(
      score_bin,
      levels = score_bin_levels,
      ordered = TRUE
    )
  )

event_lookup <- detections %>%
  group_by(event_id) %>%
  summarise(
    n_species_predicted = n_distinct(
      predicted_species,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

detections <- detections %>%
  left_join(
    event_lookup,
    by = "event_id"
  ) %>%
  mutate(
    actual_multi_event = as.integer(
      n_species_predicted > 1
    )
  )

round1 <- round1 %>%
  mutate(
    id = as.character(id),
    predicted_species = as.character(predicted_species),
    prediction_correct = as.integer(prediction_correct),
    uncertain = as.integer(uncertain),
    birdnet_score = as.numeric(birdnet_score)
  )

# ------------------------------------------------------------
# 5. Full-data species summary
# ------------------------------------------------------------

full_species_summary <- detections %>%
  group_by(predicted_species) %>%
  summarise(
    total_detections = n(),
    n_events = n_distinct(event_id),
    mean_score = mean(
      birdnet_score,
      na.rm = TRUE
    ),
    median_score = median(
      birdnet_score,
      na.rm = TRUE
    ),
    min_score = min(
      birdnet_score,
      na.rm = TRUE
    ),
    max_score = max(
      birdnet_score,
      na.rm = TRUE
    ),
    n_months = n_distinct(
      year_month,
      na.rm = TRUE
    ),
    proportion_multi_event = mean(
      actual_multi_event == 1,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

score_bin_counts <- detections %>%
  count(
    predicted_species,
    score_bin,
    name = "n_in_score_bin"
  ) %>%
  mutate(
    score_bin = as.character(score_bin)
  ) %>%
  pivot_wider(
    names_from = score_bin,
    values_from = n_in_score_bin,
    values_fill = 0,
    names_prefix = "n_"
  )

# ------------------------------------------------------------
# 6. Existing validation summary
# ------------------------------------------------------------

round1_species_summary <- round1 %>%
  group_by(predicted_species) %>%
  summarise(
    current_validation_n = n_distinct(id),
    confident_validation_n = sum(
      uncertain == 0 &
        prediction_correct %in% c(0L, 1L),
      na.rm = TRUE
    ),
    n_confirmed = sum(
      uncertain == 0 &
        prediction_correct == 1L,
      na.rm = TRUE
    ),
    current_confirmation_rate = if_else(
      confident_validation_n > 0,
      n_confirmed / confident_validation_n,
      NA_real_
    ),
    n_high_score_errors = sum(
      uncertain == 0 &
        prediction_correct == 0L &
        birdnet_score >= 0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. Apply broad-screening targets
# ------------------------------------------------------------

species_sampling_frame <- full_species_summary %>%
  left_join(
    score_bin_counts,
    by = "predicted_species"
  ) %>%
  left_join(
    round1_species_summary,
    by = "predicted_species"
  ) %>%
  mutate(
    across(
      starts_with("n_"),
      ~ replace_na(.x, 0L)
    ),
    current_validation_n = replace_na(
      current_validation_n,
      0L
    ),
    confident_validation_n = replace_na(
      confident_validation_n,
      0L
    ),
    n_confirmed = replace_na(
      n_confirmed,
      0L
    ),
    n_high_score_errors = replace_na(
      n_high_score_errors,
      0L
    ),
    previously_validated = current_validation_n > 0,
    screening_target_n = case_when(
      previously_validated ~ current_validation_n,
      total_detections <= 3 ~ total_detections,
      total_detections <= 20 ~ 4L,
      total_detections <= 100 ~ 6L,
      TRUE ~ 10L
    ),
    additional_screening_n = case_when(
      previously_validated ~ 0L,
      TRUE ~ screening_target_n
    ),
    screening_priority = case_when(
      previously_validated ~
        "already_validated_round1",
      total_detections <= 3 ~
        "unvalidated_complete_census",
      total_detections <= 20 ~
        "unvalidated_screening_n4",
      total_detections <= 100 ~
        "unvalidated_screening_n6",
      TRUE ~
        "unvalidated_screening_n10"
    ),
    systematic_error_followup_needed = case_when(
      n_high_score_errors > 0 ~
        "yes_round1_high_score_error",
      confident_validation_n >= 3 &
        current_confirmation_rate < 0.75 ~
        "yes_round1_low_confirmation",
      TRUE ~
        "not_assessed_or_no_signal"
    )
  ) %>%
  arrange(
    previously_validated,
    desc(total_detections),
    predicted_species
  )

# ------------------------------------------------------------
# 8. Quality checks
# ------------------------------------------------------------

if (n_distinct(detections$predicted_species) != 229) {
  warning(
    "Expected 229 species labels, but found ",
    n_distinct(detections$predicted_species),
    "."
  )
}

if (
  sum(species_sampling_frame$additional_screening_n) <= 0
) {
  stop(
    "No additional screening records were requested."
  )
}

# ------------------------------------------------------------
# 9. Save sampling frame
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "species_sampling_frame.csv"
)

write_csv(
  species_sampling_frame,
  output_file,
  na = ""
)

message("")
message("Species sampling frame created successfully.")
message("Output: ", output_file)
message(
  "Species labels in full data: ",
  n_distinct(detections$predicted_species)
)
message(
  "Species represented in round 1: ",
  sum(species_sampling_frame$previously_validated)
)
message(
  "Species requiring broad screening: ",
  sum(species_sampling_frame$additional_screening_n > 0)
)
message(
  "Planned screening records: ",
  sum(species_sampling_frame$additional_screening_n)
)

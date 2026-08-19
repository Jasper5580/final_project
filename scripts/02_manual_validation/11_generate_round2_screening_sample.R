# ============================================================
# 11_generate_round2_screening_sample.R
#
# Reproducible manual-validation sampling-design audit.
# This script does not recreate human listening decisions.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)

sampling_seed_text <- "20260728"

score_bin_levels <- c(
  "0.45-0.60",
  "0.60-0.75",
  "0.75-0.90",
  "0.90-1.00"
)

# When the target is smaller than the number of represented score strata,
# low and high score bins are prioritised before intermediate bins.
score_bin_priority <- c(
  "0.45-0.60",
  "0.90-1.00",
  "0.60-0.75",
  "0.75-0.90"
)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

clean_rds <- paths$clean_rds

clean_csv <- paths$clean_csv

round1_candidates <- c(
  here::here("data", "02_manual_validation", "round1_completed",
             "validation_sample_silwood_round1_completed.csv")
)

output_dir <- here::here("data", "02_manual_validation", "round2_design")

sampling_frame_file <- file.path(
  output_dir,
  "species_sampling_frame.csv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(sampling_frame_file)) {
  stop(
    paste0(
      "species_sampling_frame.csv was not found. ",
      "Run script 10 first."
    )
  )
}

# ------------------------------------------------------------
# 2. Read inputs
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

round1_files <- round1_candidates[
  file.exists(round1_candidates)
]

if (length(round1_files) == 0) {
  stop(
    "No round-1 validation file was found."
  )
}

round1 <- read_csv(
  round1_files[[1]],
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)

sampling_frame <- read_csv(
  sampling_frame_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 3. Standardise detection data
# ------------------------------------------------------------

required_columns <- c(
  "id",
  "event_id",
  "audio_id",
  "species",
  "birdnet_score",
  "detected_time",
  "date",
  "year_month",
  "start_secs",
  "end_secs",
  "duration_secs",
  "clip_link",
  "audio_link"
)

missing_columns <- setdiff(
  required_columns,
  names(detections)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing detection columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

detections <- detections %>%
  mutate(
    id = as.character(id),
    event_id = as.character(event_id),
    audio_id = as.character(audio_id),
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

# ------------------------------------------------------------
# 4. Add event-level candidate information
# ------------------------------------------------------------

event_lookup <- detections %>%
  group_by(event_id) %>%
  summarise(
    n_species_predicted = n_distinct(
      predicted_species,
      na.rm = TRUE
    ),
    event_candidate_species = paste(
      sort(unique(predicted_species)),
      collapse = "; "
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

round1_ids <- round1 %>%
  transmute(
    id = as.character(id)
  ) %>%
  distinct(id) %>%
  pull(id)

# ------------------------------------------------------------
# 5. Allocation helper
# ------------------------------------------------------------

allocate_score_bins <- function(
    candidate_data,
    target_n
) {

  count_table <- candidate_data %>%
    count(
      score_bin,
      name = "eligible_n",
      .drop = FALSE
    ) %>%
    mutate(
      score_bin = as.character(score_bin)
    ) %>%
    filter(
      !is.na(score_bin),
      eligible_n > 0
    )

  counts <- setNames(
    rep(0L, length(score_bin_levels)),
    score_bin_levels
  )

  counts[count_table$score_bin] <-
    count_table$eligible_n

  allocation <- setNames(
    rep(0L, length(score_bin_levels)),
    score_bin_levels
  )

  available_bins <- score_bin_priority[
    counts[score_bin_priority] > 0
  ]

  # First pass: represent each available bin, with low and high
  # score bins given first priority.
  for (bin_name in available_bins) {
    if (sum(allocation) >= target_n) {
      break
    }

    allocation[bin_name] <-
      allocation[bin_name] + 1L
  }

  # Fill remaining positions from the strata with the greatest
  # number of unused candidates.
  while (sum(allocation) < target_n) {

    possible_bins <- score_bin_levels[
      allocation[score_bin_levels] <
        counts[score_bin_levels]
    ]

    if (length(possible_bins) == 0) {
      break
    }

    remaining_capacity <-
      counts[possible_bins] -
      allocation[possible_bins]

    maximum_capacity <- max(
      remaining_capacity
    )

    tied_bins <- possible_bins[
      remaining_capacity == maximum_capacity
    ]

    selected_bin <- score_bin_priority[
      score_bin_priority %in% tied_bins
    ][[1]]

    allocation[selected_bin] <-
      allocation[selected_bin] + 1L
  }

  tibble(
    score_bin = names(allocation),
    selected_n = as.integer(allocation),
    eligible_n = as.integer(counts)
  ) %>%
    filter(selected_n > 0)
}

# ------------------------------------------------------------
# 6. Reproducible hash ordering
# ------------------------------------------------------------

make_sampling_hash <- function(id_value) {
  digest::digest(
    paste0(
      sampling_seed_text,
      "|",
      id_value
    ),
    algo = "md5",
    serialize = FALSE
  )
}

# ------------------------------------------------------------
# 7. Select all unvalidated species
# ------------------------------------------------------------

target_species_table <- sampling_frame %>%
  filter(
    current_validation_n == 0,
    additional_screening_n > 0
  ) %>%
  select(
    predicted_species,
    total_detections,
    screening_target_n,
    additional_screening_n,
    screening_priority
  )

sampling_candidates <- detections %>%
  filter(
    !id %in% round1_ids
  ) %>%
  inner_join(
    target_species_table,
    by = "predicted_species"
  )

species_list <- sort(
  unique(target_species_table$predicted_species)
)

selected_records <- map_dfr(
  species_list,
  function(focal_species) {

    species_candidates <- sampling_candidates %>%
      filter(
        predicted_species == focal_species
      )

    target_n <- unique(
      species_candidates$additional_screening_n
    )

    if (length(target_n) != 1) {
      stop(
        "Invalid target for species: ",
        focal_species
      )
    }

    allocation <- allocate_score_bins(
      candidate_data = species_candidates,
      target_n = target_n
    )

    map_dfr(
      seq_len(nrow(allocation)),
      function(row_number) {

        bin_name <- allocation$score_bin[[row_number]]
        selected_n <- allocation$selected_n[[row_number]]
        eligible_n <- allocation$eligible_n[[row_number]]

        species_candidates %>%
          filter(
            as.character(score_bin) == bin_name
          ) %>%
          mutate(
            sampling_hash = vapply(
              id,
              make_sampling_hash,
              FUN.VALUE = character(1)
            )
          ) %>%
          arrange(
            sampling_hash,
            id
          ) %>%
          slice_head(
            n = selected_n
          ) %>%
          mutate(
            target_species = focal_species,
            sampling_stage = "round2_screening",
            sampling_reason = screening_priority,
            sampling_stratum = paste(
              focal_species,
              bin_name,
              sep = " | "
            ),
            eligible_n = eligible_n,
            selected_n = selected_n,
            selection_probability =
              selected_n / eligible_n,
            sampling_weight =
              eligible_n / selected_n,
            species_total_detections =
              total_detections,
            species_screening_target_n =
              screening_target_n
          )
      }
    )
  }
)

# ------------------------------------------------------------
# 8. Quality checks
# ------------------------------------------------------------

planned_n <- sum(
  target_species_table$additional_screening_n
)

if (nrow(selected_records) != planned_n) {
  stop(
    "Expected ", planned_n,
    " rows but selected ",
    nrow(selected_records),
    "."
  )
}

if (anyDuplicated(selected_records$id) > 0) {
  stop(
    "The round-2 sample contains duplicate detection IDs."
  )
}

if (any(selected_records$id %in% round1_ids)) {
  stop(
    "The round-2 sample overlaps with round 1."
  )
}

if (
  n_distinct(selected_records$predicted_species) !=
    nrow(target_species_table)
) {
  stop(
    "Not every unvalidated species was represented."
  )
}

# ------------------------------------------------------------
# 9. Prepare manual-validation worksheet
# ------------------------------------------------------------

validation_sample <- selected_records %>%
  arrange(
    desc(species_total_detections),
    target_species,
    score_bin,
    sampling_hash
  ) %>%
  mutate(
    validation_id = paste0(
      "R2S_",
      stringr::str_pad(
        row_number(),
        width = 4,
        pad = "0"
      )
    ),
    validation_round = 2L,
    validation_type =
      "round2_species_screening",
    prediction_correct = NA_integer_,
    uncertain = NA_integer_,
    multiple_species = NA_integer_,
    human_species = NA_character_,
    audio_quality = NA_character_,
    validator = NA_character_,
    validation_basis = NA_character_,
    notes = NA_character_,
    score_bin = as.character(score_bin)
  ) %>%
  select(
    validation_id,
    validation_round,
    validation_type,
    sampling_stage,
    sampling_reason,
    target_species,
    sampling_stratum,
    eligible_n,
    selected_n,
    selection_probability,
    sampling_weight,
    species_total_detections,
    species_screening_target_n,
    id,
    event_id,
    audio_id,
    predicted_species,
    birdnet_score,
    score_bin,
    detected_time,
    date,
    year_month,
    start_secs,
    end_secs,
    duration_secs,
    clip_link,
    audio_link,
    n_species_predicted,
    event_candidate_species,
    actual_multi_event,
    prediction_correct,
    uncertain,
    multiple_species,
    human_species,
    audio_quality,
    validator,
    validation_basis,
    notes
  )

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

sample_file <- file.path(
  output_dir,
  "validation_round2_screening_to_validate.csv"
)

write_csv(
  validation_sample,
  sample_file,
  na = ""
)

summary_table <- bind_rows(
  tibble(
    measure = c(
      "full_detection_rows",
      "full_species_labels",
      "full_acoustic_events",
      "round1_validation_rows",
      "round1_species_covered",
      "unvalidated_species_before_round2",
      "round2_screening_rows",
      "round2_species_covered",
      "species_coverage_after_round2",
      "round1_overlap_rows",
      "duplicate_round2_detection_ids",
      "round2_multi_prediction_rows",
      "round2_high_score_rows_ge_0_9"
    ),
    value = c(
      nrow(detections),
      n_distinct(detections$predicted_species),
      n_distinct(detections$event_id),
      nrow(round1),
      n_distinct(round1$predicted_species),
      nrow(target_species_table),
      nrow(validation_sample),
      n_distinct(validation_sample$predicted_species),
      n_distinct(round1$predicted_species) +
        n_distinct(validation_sample$predicted_species),
      sum(validation_sample$id %in% round1_ids),
      sum(duplicated(validation_sample$id)),
      sum(validation_sample$actual_multi_event == 1),
      sum(validation_sample$birdnet_score >= 0.90)
    )
  ),
  validation_sample %>%
    count(
      sampling_reason,
      name = "value"
    ) %>%
    transmute(
      measure = paste0(
        "rows_",
        sampling_reason
      ),
      value
    ),
  validation_sample %>%
    count(
      score_bin,
      name = "value"
    ) %>%
    transmute(
      measure = paste0(
        "rows_score_bin_",
        score_bin
      ),
      value
    )
)

summary_file <- file.path(
  output_dir,
  "round2_screening_summary.csv"
)

write_csv(
  summary_table,
  summary_file
)

message("")
message("Round-2 screening sample generated successfully.")
message("Sample: ", sample_file)
message("Summary: ", summary_file)
message("Rows selected: ", nrow(validation_sample))
message(
  "Species represented: ",
  n_distinct(validation_sample$predicted_species)
)
message(
  "Total species coverage after round 2: ",
  n_distinct(round1$predicted_species) +
    n_distinct(validation_sample$predicted_species)
)

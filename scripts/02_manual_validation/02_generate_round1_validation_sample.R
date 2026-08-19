# ============================================================
# 02_generate_round1_validation_sample.R
#
# Purpose:
# Reproduce the Round-1 stratified sampling design. The output is an
# unreviewed audit sample and does not replace the completed manual
# validation records supplied in data/02_manual_validation/round1_completed.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)

birdnet <- readr::read_csv(paths$clean_csv, show_col_types = FALSE)
output_dir <- ensure_dir(here::here("data", "02_manual_validation", "round1_design"))
table_dir <- silwood_result_dir(paths, "01_validation", "tables")

birdnet_valid_pool <- birdnet %>%
  dplyr::mutate(
    score_bin = dplyr::case_when(
      .data$birdnet_score >= 0.45 & .data$birdnet_score < 0.60 ~ "0.45-0.60",
      .data$birdnet_score >= 0.60 & .data$birdnet_score < 0.75 ~ "0.60-0.75",
      .data$birdnet_score >= 0.75 & .data$birdnet_score < 0.90 ~ "0.75-0.90",
      .data$birdnet_score >= 0.90 & .data$birdnet_score <= 1.00 ~ "0.90-1.00",
      TRUE ~ NA_character_
    ),
    year_week = paste0(
      lubridate::isoyear(.data$date), "-W",
      stringr::str_pad(lubridate::isoweek(.data$date), 2, pad = "0")
    )
  ) %>%
  dplyr::filter(!is.na(.data$score_bin), !is.na(.data$clip_link) | !is.na(.data$audio_link))

top_species <- birdnet_valid_pool %>%
  dplyr::count(.data$species, sort = TRUE) %>%
  dplyr::slice_head(n = 10L) %>%
  dplyr::pull(.data$species)

multi_event_ids <- birdnet_valid_pool %>%
  dplyr::count(.data$event_id) %>%
  dplyr::filter(.data$n > 1L) %>%
  dplyr::pull(.data$event_id)

set.seed(123L)

sample_common <- birdnet_valid_pool %>%
  dplyr::filter(.data$species %in% top_species) %>%
  dplyr::group_by(.data$species, .data$score_bin) %>%
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(8L, nrow(.x)))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(validation_type = "common_species_score_bin")

sample_multi <- birdnet_valid_pool %>%
  dplyr::filter(.data$event_id %in% multi_event_ids) %>%
  dplyr::group_by(.data$score_bin) %>%
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(15L, nrow(.x)))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(validation_type = "multi_species_event")

sample_less_common <- birdnet_valid_pool %>%
  dplyr::filter(!.data$species %in% top_species) %>%
  dplyr::group_by(.data$score_bin) %>%
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(15L, nrow(.x)))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(validation_type = "less_common_species")

validation_sample <- dplyr::bind_rows(sample_common, sample_multi, sample_less_common) %>%
  dplyr::distinct(.data$id, .keep_all = TRUE) %>%
  dplyr::arrange(.data$species, .data$score_bin, dplyr::desc(.data$birdnet_score)) %>%
  dplyr::mutate(
    validation_id = dplyr::row_number(),
    predicted_species = .data$species,
    prediction_correct = NA_integer_,
    uncertain = NA_integer_,
    multiple_species = NA_integer_,
    human_species = NA_character_,
    audio_quality = NA_character_,
    validator = NA_character_,
    validation_basis = NA_character_,
    notes = NA_character_
  ) %>%
  dplyr::select(
    .data$validation_id, .data$validation_type, .data$id, .data$event_id,
    .data$audio_id, .data$predicted_species, .data$birdnet_score,
    .data$score_bin, .data$detected_time, .data$date, .data$year_month,
    .data$year_week, .data$start_secs, .data$end_secs, .data$duration_secs,
    .data$clip_link, .data$audio_link, .data$prediction_correct,
    .data$uncertain, .data$multiple_species, .data$human_species,
    .data$audio_quality, .data$validator, .data$validation_basis, .data$notes
  )

summary <- validation_sample %>%
  dplyr::count(.data$validation_type, .data$score_bin, name = "n_records")

readr::write_csv(
  validation_sample,
  file.path(output_dir, "round1_sample_to_validate.csv"),
  na = ""
)
readr::write_csv(summary, file.path(table_dir, "round1_sampling_design_summary.csv"), na = "")

message("Round-1 design sample written to: ", output_dir)
message("This file contains no manual decisions and is not used by the final pipeline.")

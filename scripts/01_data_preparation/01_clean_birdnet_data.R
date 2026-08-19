# ============================================================
# 01_clean_birdnet_data.R
#
# Purpose:
# Clean the raw Silwood BirdNET export and create the event-species
# candidate dataset used by all later steps.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)

clean_dir <- ensure_dir(here::here("data", "01_cleaned"))
table_dir <- silwood_result_dir(paths, "00_data_preparation", "tables")

raw <- readr::read_csv(
  paths$raw_detection_csv,
  show_col_types = FALSE,
  progress = FALSE
)

required_raw_columns <- c(
  "id", "analysis", "tags", "start_secs", "end_secs", "confidence",
  "detected_time", "audio_id", "site", "recorder", "latitude",
  "longitude", "clip_link", "audio_link"
)
assert_columns(raw, required_raw_columns, "raw BirdNET export")

analysis_summary <- raw %>%
  dplyr::count(.data$analysis, sort = TRUE, name = "n_rows")

confidence_summary <- raw %>%
  dplyr::group_by(.data$analysis) %>%
  dplyr::summarise(
    n = dplyr::n(),
    min_confidence = min(.data$confidence, na.rm = TRUE),
    median_confidence = stats::median(.data$confidence, na.rm = TRUE),
    mean_confidence = mean(.data$confidence, na.rm = TRUE),
    max_confidence = max(.data$confidence, na.rm = TRUE),
    .groups = "drop"
  )

birdnet_clean <- raw %>%
  dplyr::filter(
    .data$analysis == "birdnet-lite",
    !is.na(.data$tags),
    stringr::str_squish(.data$tags) != "",
    !is.na(.data$confidence),
    .data$confidence >= 0,
    .data$confidence <= 1
  ) %>%
  dplyr::mutate(
    species = stringr::str_squish(.data$tags),
    detected_time = lubridate::ymd_hms(
      .data$detected_time,
      quiet = TRUE,
      tz = "UTC"
    ),
    date = as.Date(.data$detected_time),
    year = lubridate::year(.data$detected_time),
    month = lubridate::month(.data$detected_time),
    month_name = month.abb[.data$month],
    year_month = format(.data$date, "%Y-%m"),
    season = dplyr::case_when(
      .data$month %in% c(12L, 1L, 2L) ~ "Winter",
      .data$month %in% c(3L, 4L, 5L) ~ "Spring",
      .data$month %in% c(6L, 7L, 8L) ~ "Summer",
      .data$month %in% c(9L, 10L, 11L) ~ "Autumn",
      TRUE ~ NA_character_
    ),
    duration_secs = .data$end_secs - .data$start_secs,
    event_id = paste(.data$audio_id, .data$start_secs, .data$end_secs, sep = "_"),
    birdnet_score = .data$confidence
  ) %>%
  dplyr::select(
    .data$id, .data$event_id, .data$audio_id, .data$analysis,
    .data$species, .data$birdnet_score, .data$detected_time, .data$date,
    .data$year, .data$month, .data$month_name, .data$year_month,
    .data$season, .data$start_secs, .data$end_secs, .data$duration_secs,
    .data$site, .data$recorder, .data$latitude, .data$longitude,
    .data$clip_link, .data$audio_link, dplyr::everything()
  )

# Regression checks against the frozen thesis dataset.
assert_expected(nrow(raw), 79398L, "Raw detection rows")
assert_expected(nrow(birdnet_clean), 65255L, "Cleaned candidate rows")
assert_expected(dplyr::n_distinct(birdnet_clean$event_id), 63985L, "Acoustic events")
assert_expected(dplyr::n_distinct(birdnet_clean$species), 229L, "Predicted species")
assert_expected(dplyr::n_distinct(birdnet_clean$year_month), 23L, "Represented months")

# Event-level summary.
event_summary <- birdnet_clean %>%
  dplyr::group_by(
    .data$event_id, .data$audio_id, .data$start_secs, .data$end_secs,
    .data$detected_time, .data$year_month, .data$season
  ) %>%
  dplyr::summarise(
    n_species_predicted = dplyr::n_distinct(.data$species),
    n_rows = dplyr::n(),
    max_score = max(.data$birdnet_score, na.rm = TRUE),
    species_candidates = paste(
      .data$species[order(.data$birdnet_score, decreasing = TRUE)],
      round(.data$birdnet_score[order(.data$birdnet_score, decreasing = TRUE)], 3),
      sep = " = ",
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$n_species_predicted), dplyr::desc(.data$max_score))

multi_species_events <- event_summary %>%
  dplyr::filter(.data$n_species_predicted > 1L)

species_month_counts_naive <- birdnet_clean %>%
  dplyr::count(.data$year_month, .data$season, .data$species, name = "n_detections") %>%
  dplyr::arrange(.data$year_month, dplyr::desc(.data$n_detections))

species_month_counts_hard <- birdnet_clean %>%
  dplyr::group_by(.data$event_id) %>%
  dplyr::slice_max(order_by = .data$birdnet_score, n = 1L, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::count(.data$year_month, .data$season, .data$species, name = "n_detections") %>%
  dplyr::arrange(.data$year_month, dplyr::desc(.data$n_detections))

threshold_counts <- purrr::map_dfr(c(0.50, 0.70, 0.90), function(threshold) {
  birdnet_clean %>%
    dplyr::filter(.data$birdnet_score >= threshold) %>%
    dplyr::count(.data$year_month, .data$season, .data$species, name = "n_detections") %>%
    dplyr::mutate(threshold = threshold)
})

monthly_summary <- birdnet_clean %>%
  dplyr::group_by(.data$year_month, .data$season) %>%
  dplyr::summarise(
    n_detections = dplyr::n(),
    n_events = dplyr::n_distinct(.data$event_id),
    n_audio_files = dplyr::n_distinct(.data$audio_id),
    n_species_naive = dplyr::n_distinct(.data$species),
    mean_score = mean(.data$birdnet_score, na.rm = TRUE),
    median_score = stats::median(.data$birdnet_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$year_month)

time_range <- birdnet_clean %>%
  dplyr::summarise(
    first_detection = min(.data$detected_time, na.rm = TRUE),
    last_detection = max(.data$detected_time, na.rm = TRUE),
    n_detections = dplyr::n(),
    n_events = dplyr::n_distinct(.data$event_id),
    n_species = dplyr::n_distinct(.data$species),
    n_months = dplyr::n_distinct(.data$year_month),
    n_audio_files = dplyr::n_distinct(.data$audio_id),
    n_sites = dplyr::n_distinct(.data$site),
    n_recorders = dplyr::n_distinct(.data$recorder)
  )

top_species <- birdnet_clean %>%
  dplyr::count(.data$species, sort = TRUE, name = "n_detections")

species_confidence_summary <- birdnet_clean %>%
  dplyr::group_by(.data$species) %>%
  dplyr::summarise(
    n_detections = dplyr::n(),
    mean_score = mean(.data$birdnet_score, na.rm = TRUE),
    median_score = stats::median(.data$birdnet_score, na.rm = TRUE),
    min_score = min(.data$birdnet_score, na.rm = TRUE),
    max_score = max(.data$birdnet_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$n_detections))

duration_summary <- birdnet_clean %>%
  dplyr::summarise(
    min_duration = min(.data$duration_secs, na.rm = TRUE),
    median_duration = stats::median(.data$duration_secs, na.rm = TRUE),
    mean_duration = mean(.data$duration_secs, na.rm = TRUE),
    max_duration = max(.data$duration_secs, na.rm = TRUE),
    n_non_positive_duration = sum(.data$duration_secs <= 0, na.rm = TRUE)
  )

monthly_detection_summary <- birdnet_clean %>%
  dplyr::count(.data$year_month, name = "n_detections") %>%
  dplyr::arrange(.data$year_month)

quality_checks <- tibble::tibble(
  check = c("raw_rows", "clean_rows", "events", "species", "months"),
  observed = c(
    nrow(raw), nrow(birdnet_clean), dplyr::n_distinct(birdnet_clean$event_id),
    dplyr::n_distinct(birdnet_clean$species),
    dplyr::n_distinct(birdnet_clean$year_month)
  ),
  expected = c(79398L, 65255L, 63985L, 229L, 23L),
  passed = observed == expected
)

readr::write_csv(birdnet_clean, file.path(clean_dir, "birdnet_detections_clean.csv"), na = "")
saveRDS(birdnet_clean, file.path(clean_dir, "birdnet_detections_clean.rds"))
readr::write_csv(event_summary, file.path(clean_dir, "acoustic_event_summary.csv"), na = "")
saveRDS(event_summary, file.path(clean_dir, "acoustic_event_summary.rds"))
readr::write_csv(multi_species_events, file.path(clean_dir, "multi_species_events.csv"), na = "")
readr::write_csv(species_month_counts_naive, file.path(clean_dir, "species_month_counts_naive.csv"), na = "")
readr::write_csv(species_month_counts_hard, file.path(clean_dir, "species_month_counts_hard_assignment.csv"), na = "")
readr::write_csv(threshold_counts, file.path(clean_dir, "species_month_counts_thresholds.csv"), na = "")
readr::write_csv(monthly_summary, file.path(clean_dir, "monthly_summary.csv"), na = "")

readr::write_csv(analysis_summary, file.path(table_dir, "analysis_summary.csv"), na = "")
readr::write_csv(confidence_summary, file.path(table_dir, "confidence_summary_by_analysis.csv"), na = "")
readr::write_csv(time_range, file.path(table_dir, "time_range_summary.csv"), na = "")
readr::write_csv(top_species, file.path(table_dir, "top_species.csv"), na = "")
readr::write_csv(species_confidence_summary, file.path(table_dir, "species_confidence_summary.csv"), na = "")
readr::write_csv(duration_summary, file.path(table_dir, "duration_summary.csv"), na = "")
readr::write_csv(monthly_detection_summary, file.path(table_dir, "monthly_detection_summary.csv"), na = "")
readr::write_csv(quality_checks, file.path(table_dir, "data_cleaning_quality_checks.csv"), na = "")

message("Data cleaning completed: ", nrow(birdnet_clean), " event-species candidates.")
message("Cleaned data: ", clean_dir)
message("Cleaning tables: ", table_dir)

# ============================================================
# 19_sampling_effort_proxy_sensitivity_fixed5min.R
#
# Purpose:
# Calculate monthly represented-audio coverage proxies from the
# complete raw detection export and examine their association
# with Primary M3 monthly diversity.
#
# represented_audio_hours counts only unique five-minute source
# recordings appearing in the export. It is not complete recorder
# uptime because files with no exported detections may be absent.
# ============================================================

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()

if (!file.exists(paths$raw_detection_csv)) {
  stop("Raw detection export not found: ", paths$raw_detection_csv)
}

mc_file <- file.path(
  paths$results_root,
  "05_monte_carlo",
  "models",
  "monte_carlo_results.rds"
)
if (!file.exists(mc_file)) {
  stop("Monte Carlo result object not found: ", mc_file)
}

raw_data <- readr::read_csv(
  paths$raw_detection_csv,
  na = c("", "NA", "N/A"),
  show_col_types = FALSE
)
mc_results <- readRDS(mc_file)

required_raw_columns <- c(
  "id",
  "audio_id",
  "start_secs",
  "end_secs",
  "detected_time",
  "upload_time",
  "recorder"
)
assert_columns(raw_data, required_raw_columns, "raw detection export")

raw_standardised <- raw_data %>%
  dplyr::transmute(
    detection_id = as.character(.data$id),
    audio_id = as.character(.data$audio_id),
    recorder = as.character(.data$recorder),
    start_secs = as.numeric(.data$start_secs),
    end_secs = as.numeric(.data$end_secs),
    detected_time = parse_time_safely(.data$detected_time),
    upload_time = parse_time_safely(.data$upload_time)
  ) %>%
  dplyr::mutate(
    recording_start = .data$detected_time - lubridate::seconds(.data$start_secs),
    recording_date = as.Date(.data$recording_start),
    recording_hour = lubridate::floor_date(.data$recording_start, unit = "hour"),
    year_month = format(.data$recording_start, "%Y-%m"),
    upload_start_difference_seconds = abs(
      as.numeric(
        difftime(
          .data$recording_start,
          .data$upload_time,
          units = "secs"
        )
      )
    )
  )

assert_unique(raw_standardised$detection_id, "Raw detection IDs")

invalid_offset_n <- sum(
  is.na(raw_standardised$start_secs) |
    is.na(raw_standardised$end_secs) |
    raw_standardised$start_secs < 0 |
    raw_standardised$end_secs < raw_standardised$start_secs
)
if (invalid_offset_n > 0L) {
  stop("Invalid raw start/end offsets found: ", invalid_offset_n)
}

audio_table <- raw_standardised %>%
  dplyr::arrange(.data$audio_id, .data$recording_start) %>%
  dplyr::group_by(.data$audio_id) %>%
  dplyr::summarise(
    recorder = dplyr::first(.data$recorder),
    recording_start = dplyr::first(.data$recording_start),
    recording_date = dplyr::first(.data$recording_date),
    recording_hour = dplyr::first(.data$recording_hour),
    year_month = dplyr::first(.data$year_month),
    detections_in_export = dplyr::n(),
    within_audio_start_range_seconds = as.numeric(
      difftime(
        max(.data$recording_start),
        min(.data$recording_start),
        units = "secs"
      )
    ),
    .groups = "drop"
  )

coverage <- audio_table %>%
  dplyr::group_by(.data$year_month) %>%
  dplyr::summarise(
    unique_audio_ids = dplyr::n_distinct(.data$audio_id),
    represented_audio_hours = .data$unique_audio_ids * 5 / 60,
    active_detection_days = dplyr::n_distinct(.data$recording_date),
    active_detection_hours = dplyr::n_distinct(.data$recording_hour),
    first_recording = min(.data$recording_start),
    last_recording = max(.data$recording_start),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$year_month)

primary_monthly <- mc_results$monthly_summary %>%
  dplyr::filter(.data$scenario == "Primary_M3") %>%
  dplyr::select(
    year_month,
    q,
    q_label,
    diversity_median = median,
    diversity_lower_95 = lower_95,
    diversity_upper_95 = upper_95
  )

monthly_with_coverage <- primary_monthly %>%
  dplyr::left_join(coverage, by = "year_month")

proxy_names <- c(
  "represented_audio_hours",
  "active_detection_days",
  "active_detection_hours"
)

correlation_table <- tidyr::expand_grid(
  q_value = config$q_values,
  proxy_name = proxy_names
) %>%
  purrr::pmap_dfr(
    function(q_value, proxy_name) {
      working <- monthly_with_coverage %>%
        dplyr::filter(.data$q == .env$q_value) %>%
        dplyr::filter(
          is.finite(.data$diversity_median),
          is.finite(.data[[proxy_name]])
        )

      if (nrow(working) < 3L) {
        result <- list(n = nrow(working), rho = NA_real_, p = NA_real_)
      } else {
        test <- suppressWarnings(
          stats::cor.test(
            working$diversity_median,
            working[[proxy_name]],
            method = "spearman",
            exact = FALSE
          )
        )
        result <- list(
          n = nrow(working),
          rho = unname(test$estimate),
          p = test$p.value
        )
      }

      tibble::tibble(
        q = .env$q_value,
        q_label = q_label(.env$q_value),
        proxy = proxy_name,
        proxy_label = dplyr::recode(
          proxy_name,
          represented_audio_hours = "Represented audio hours",
          active_detection_days = "Active detection days",
          active_detection_hours = "Active detection hours"
        ),
        n_months = result$n,
        spearman_rho = result$rho,
        p_value = result$p
      )
    }
  )

represented_cutoff <- stats::quantile(
  coverage$represented_audio_hours,
  0.25,
  type = 8,
  names = FALSE
)
days_cutoff <- stats::quantile(
  coverage$active_detection_days,
  0.25,
  type = 8,
  names = FALSE
)

coverage <- coverage %>%
  dplyr::mutate(
    low_represented_audio = .data$represented_audio_hours <= represented_cutoff,
    low_active_days = .data$active_detection_days <= days_cutoff,
    low_coverage = .data$low_represented_audio | .data$low_active_days,
    coverage_set = dplyr::if_else(
      .data$low_coverage,
      "Low represented coverage",
      "Higher represented coverage"
    )
  )

quality_checks <- tibble::tibble(
  check = c(
    "Raw detection rows",
    "Unique source audio IDs",
    "Months with represented audio",
    "Months with Primary M3 diversity",
    "Months classified as low coverage",
    "Upload-time match within 5 seconds"
  ),
  value = c(
    nrow(raw_standardised),
    nrow(audio_table),
    nrow(coverage),
    dplyr::n_distinct(primary_monthly$year_month),
    sum(coverage$low_coverage),
    mean(
      raw_standardised$upload_start_difference_seconds <= 5,
      na.rm = TRUE
    )
  )
)

assert_expected(nrow(raw_standardised), 79398, "Raw detection rows")
assert_expected(nrow(coverage), 23, "Coverage months")
assert_expected(
  dplyr::n_distinct(primary_monthly$year_month),
  23,
  "Primary M3 diversity months"
)

model_dir <- ensure_dir(
  silwood_result_dir(paths, "07_sampling_coverage", "models")
)
table_dir <- ensure_dir(
  silwood_result_dir(paths, "07_sampling_coverage", "tables")
)
figure_dir <- ensure_dir(
  silwood_result_dir(paths, "07_sampling_coverage", "figures")
)

readr::write_csv(
  coverage,
  file.path(table_dir, "monthly_sampling_coverage.csv"),
  na = ""
)
readr::write_csv(
  monthly_with_coverage,
  file.path(table_dir, "monthly_diversity_with_coverage.csv"),
  na = ""
)
readr::write_csv(
  correlation_table,
  file.path(table_dir, "diversity_coverage_correlations.csv"),
  na = ""
)
readr::write_csv(
  coverage %>% dplyr::filter(.data$low_coverage),
  file.path(table_dir, "low_coverage_months.csv"),
  na = ""
)
readr::write_csv(
  quality_checks,
  file.path(table_dir, "sampling_coverage_quality_checks.csv"),
  na = ""
)
saveRDS(
  list(
    audio_table = audio_table,
    coverage = coverage,
    monthly_with_coverage = monthly_with_coverage,
    correlations = correlation_table,
    represented_cutoff = represented_cutoff,
    days_cutoff = days_cutoff
  ),
  file.path(model_dir, "sampling_coverage_results.rds"),
  compress = "xz"
)

manuscript_table_dir <- ensure_dir(
  silwood_result_dir(paths, "07_sampling_coverage", "tables")
)
coverage_correlation_publication <- correlation_table %>%
  dplyr::transmute(
    Diversity_order = .data$q_label,
    Coverage_proxy = .data$proxy_label,
    Months = .data$n_months,
    Spearman_rho = round(.data$spearman_rho, 3),
    P_value = signif(.data$p_value, 3)
  )
readr::write_csv(
  coverage_correlation_publication,
  file.path(manuscript_table_dir, "Table_diversity_coverage_correlations.csv"),
  na = ""
)

coverage_long <- coverage %>%
  dplyr::select(
    year_month,
    coverage_set,
    represented_audio_hours,
    active_detection_days,
    active_detection_hours
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      represented_audio_hours,
      active_detection_days,
      active_detection_hours
    ),
    names_to = "proxy",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    month_date = as.Date(paste0(.data$year_month, "-01")),
    proxy_label = dplyr::recode(
      .data$proxy,
      represented_audio_hours = "Represented audio hours",
      active_detection_days = "Active detection days",
      active_detection_hours = "Active detection hours"
    )
  )

coverage_figure <- ggplot2::ggplot(
  coverage_long,
  ggplot2::aes(
    x = .data$month_date,
    y = .data$value,
    group = 1
  )
) +
  ggplot2::geom_line(linewidth = 0.8, colour = "#4C78A8") +
  ggplot2::geom_point(
    ggplot2::aes(fill = .data$coverage_set),
    shape = 21,
    size = 3,
    colour = "white",
    stroke = 0.7
  ) +
  ggplot2::facet_wrap(~proxy_label, scales = "free_y", ncol = 1) +
  ggplot2::scale_fill_manual(
    values = c(
      "Higher represented coverage" = "#009E73",
      "Low represented coverage" = "#D55E00"
    )
  ) +
  ggplot2::scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
  ggplot2::labs(
    title = "Monthly recording coverage represented in the detection export",
    x = NULL,
    y = NULL,
    fill = "Coverage classification",
    caption = "These proxies exclude source recordings that produced no exported detections and therefore do not equal complete recorder uptime."
  ) +
  theme_silwood(10)

save_figure(
  coverage_figure,
  file.path(figure_dir, "figure_09_monthly_sampling_coverage"),
  width = 11.5,
  height = 9
)

scatter_data <- monthly_with_coverage %>%
  dplyr::mutate(
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    ),
    coverage_set = dplyr::if_else(
      .data$year_month %in% coverage$year_month[coverage$low_coverage],
      "Low represented coverage",
      "Higher represented coverage"
    )
  )

scatter_annotations <- correlation_table %>%
  dplyr::filter(.data$proxy == "represented_audio_hours") %>%
  dplyr::mutate(
    annotation = paste0(
      "Spearman rho = ",
      sprintf("%.2f", .data$spearman_rho),
      "\n",
      format_p_value(.data$p_value)
    )
  )

coverage_scatter <- ggplot2::ggplot(
  scatter_data,
  ggplot2::aes(
    x = .data$represented_audio_hours,
    y = .data$diversity_median,
    fill = .data$coverage_set
  )
) +
  ggplot2::geom_smooth(
    ggplot2::aes(
      x = .data$represented_audio_hours,
      y = .data$diversity_median,
      group = 1
    ),
    method = "loess",
    se = TRUE,
    colour = "grey35",
    fill = "grey80",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(shape = 21, size = 3.2, colour = "white", stroke = 0.7) +
  ggplot2::geom_text(
    data = scatter_annotations,
    ggplot2::aes(x = -Inf, y = Inf, label = .data$annotation),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.2,
    size = 3.2
  ) +
  ggplot2::facet_wrap(~q_label, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(
    values = c(
      "Higher represented coverage" = "#009E73",
      "Low represented coverage" = "#D55E00"
    )
  ) +
  ggplot2::labs(
    title = "Primary M3 monthly diversity and represented audio hours",
    x = "Represented audio hours",
    y = "Monthly diversity median",
    fill = "Coverage classification"
  ) +
  theme_silwood(9)

save_figure(
  coverage_scatter,
  file.path(figure_dir, "figure_10_diversity_vs_represented_audio"),
  width = 12.5,
  height = 5.8
)

correlation_heatmap <- correlation_table %>%
  dplyr::mutate(
    q_label = factor(
      .data$q_label,
      levels = vapply(config$q_values, q_label, character(1))
    ),
    proxy_label = factor(
      .data$proxy_label,
      levels = c(
        "Represented audio hours",
        "Active detection hours",
        "Active detection days"
      )
    ),
    label = sprintf("%.2f", .data$spearman_rho)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data$proxy_label,
      y = .data$q_label,
      fill = .data$spearman_rho
    )
  ) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = .data$label), fontface = "bold") +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Spearman correlations between diversity and coverage proxies",
    x = NULL,
    y = NULL,
    fill = "Spearman rho"
  ) +
  theme_silwood(9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
    legend.position = "right"
  )

save_figure(
  correlation_heatmap,
  file.path(figure_dir, "figure_S04_coverage_correlation_heatmap"),
  width = 9.5,
  height = 4.8
)

message("Sampling-coverage analysis completed.")
message("Low-coverage months: ", sum(coverage$low_coverage))

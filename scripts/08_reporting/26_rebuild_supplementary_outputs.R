# ============================================================
# 26_rebuild_supplementary_outputs.R
#
# Purpose:
# Rebuild all six supplementary manuscript figures and all six
# supplementary tables from the saved final pipeline outputs.
#
# This script changes presentation only. It does not refit models,
# rerun Monte Carlo simulations, or alter any numerical result.
#
# Main presentation changes:
#   * explicit solid axes or panel borders for all plotted panels;
#   * explanatory sentences removed from inside figures and moved
#     to the LaTeX captions;
#   * alpha -> beta -> gamma ordering in all partition displays;
#   * neutral wording based on differences from Primary M3 rather
#     than treating Primary M3 as an independently known truth;
#   * existing supplementary figure and CSV filenames overwritten;
#   * complete LaTeX table fragments regenerated in
#     results/10_final_outputs/supplementary_latex/.
#
# Run after scripts 12-21, or after the complete final pipeline.
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()

figure_dir <- ensure_dir(
  silwood_final_dir(paths, "supplementary_figures")
)

table_csv_dir <- ensure_dir(
  silwood_final_dir(paths, "supplementary_tables")
)

table_tex_dir <- ensure_dir(
  silwood_final_dir(paths, "supplementary_latex")
)


# ------------------------------------------------------------
# 2. Shared helpers
# ------------------------------------------------------------

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

read_rds_or_csv <- function(rds_path, csv_path, label) {
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  if (file.exists(csv_path)) {
    return(readr::read_csv(csv_path, show_col_types = FALSE))
  }
  stop(
    label,
    " not found. Expected either: ",
    rds_path,
    " or ",
    csv_path
  )
}

supp_axis_theme <- function(
    base_size = 9,
    legend_position = "bottom",
    panel_border = FALSE
) {
  border_element <- if (isTRUE(panel_border)) {
    ggplot2::element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.45
    )
  } else {
    ggplot2::element_blank()
  }

  theme_silwood(base_size) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(
        colour = "black",
        linewidth = 0.50
      ),
      axis.line.y = ggplot2::element_line(
        colour = "black",
        linewidth = 0.50
      ),
      axis.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = 0.45
      ),
      axis.ticks.length = grid::unit(2.5, "pt"),
      panel.border = border_element,
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.position = legend_position
    )
}

supp_heatmap_theme <- function(base_size = 9) {
  theme_silwood(base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.45
      ),
      axis.line.x = ggplot2::element_line(
        colour = "black",
        linewidth = 0.45
      ),
      axis.line.y = ggplot2::element_line(
        colour = "black",
        linewidth = 0.45
      ),
      axis.ticks = ggplot2::element_line(
        colour = "black",
        linewidth = 0.40
      ),
      axis.ticks.length = grid::unit(2.5, "pt"),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      legend.position = "right"
    )
}

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}
latex_row_break <- "\\\\"


fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "--",
    formatC(as.numeric(x), format = "f", digits = digits)
  )
}

fmt_int <- function(x) {
  ifelse(is.na(x), "--", formatC(as.integer(x), format = "d"))
}

fmt_percent_value <- function(x, digits = 2, signed = FALSE) {
  ifelse(
    is.na(x),
    "--",
    paste0(
      ifelse(signed & as.numeric(x) > 0, "+", ""),
      formatC(as.numeric(x), format = "f", digits = digits),
      "\\%"
    )
  )
}

fmt_p_value <- function(x) {
  vapply(
    as.numeric(x),
    function(value) {
      if (!is.finite(value)) return("--")
      if (value < 0.001) {
        exponent <- floor(log10(value))
        coefficient <- value / (10^exponent)
        return(
          paste0(
            "$",
            formatC(coefficient, format = "f", digits = 2),
            "\\times 10^{",
            exponent,
            "}$"
          )
        )
      }
      formatC(value, format = "f", digits = 3)
    },
    FUN.VALUE = character(1)
  )
}

write_longtable_fragment <- function(
    file,
    table_number,
    caption,
    label,
    column_spec,
    header_cells,
    body_rows,
    landscape = TRUE,
    font_command = "\\scriptsize",
    tabcolsep = "2.5pt",
    arraystretch = "1.05"
) {
  n_columns <- length(header_cells)
  header_line <- paste(header_cells, collapse = " & ")

  lines <- c(
    "\\clearpage",
    if (isTRUE(landscape)) "\\begin{landscape}" else character(),
    font_command,
    paste0("\\setlength{\\tabcolsep}{", tabcolsep, "}"),
    paste0("\\renewcommand{\\arraystretch}{", arraystretch, "}"),
    paste0("\\begin{longtable}{", column_spec, "}"),
    paste0("\\caption{", caption, "}\\label{", label, "}", latex_row_break),
    "\\toprule",
    paste0(header_line, " \\\\"),
    "\\midrule",
    "\\endfirsthead",
    paste0(
      "\\multicolumn{",
      n_columns,
      "}{l}{\\textit{Table ",
      table_number,
      " continued from previous page}}\\\\"
    ),
    "\\toprule",
    paste0(header_line, " \\\\"),
    "\\midrule",
    "\\endhead",
    "\\midrule",
    paste0(
      "\\multicolumn{",
      n_columns,
      "}{r}{\\textit{Continued on next page}}\\\\"
    ),
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot",
    body_rows,
    "\\end{longtable}",
    if (isTRUE(landscape)) "\\end{landscape}" else character()
  )

  writeLines(lines, file, useBytes = TRUE)
  invisible(file)
}

write_small_table_fragment <- function(
    file,
    caption,
    label,
    column_spec,
    header_cells,
    body_rows,
    width = "0.92\\textwidth"
) {
  lines <- c(
    "\\clearpage",
    "\\begin{table}[p]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\small",
    paste0("\\begin{tabular}{", column_spec, "}"),
    "\\toprule",
    paste0(paste(header_cells, collapse = " & "), " \\\\"),
    "\\midrule",
    body_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(lines, file, useBytes = TRUE)
  invisible(file)
}


# ------------------------------------------------------------
# 3. Shared ordering
# ------------------------------------------------------------

q_order <- c(
  "q = 0: Richness",
  "q = 1: Exponential Shannon",
  "q = 2: Inverse Simpson"
)

component_order <- c("Alpha", "Beta", "Gamma")

probability_scenario_order <- c(
  "Primary M3: species-aware",
  "Primary M2: pooled",
  "Audio-only sensitivity"
)

workflow_order <- c(
  "Primary M2: pooled",
  "Audio-only sensitivity",
  "Naive: all candidates",
  "Hard top-1",
  "Threshold >= 0.50",
  "Threshold >= 0.70",
  "Threshold >= 0.90"
)

deterministic_order <- c(
  "Naive: all candidates",
  "Hard top-1",
  "Threshold >= 0.50",
  "Threshold >= 0.70",
  "Threshold >= 0.90"
)


# ============================================================
# SUPPLEMENTARY FIGURES
# ============================================================


# ------------------------------------------------------------
# Figure S1: all Primary M3 random intercepts
# ------------------------------------------------------------

species_effects_m3 <- read_required_csv(
  file.path(
    paths$results_root,
    "02_calibration",
    "tables",
    "calibration_species_effects_m3.csv"
  ),
  "Primary M3 species effects"
)

assert_columns(
  species_effects_m3,
  c(
    "predicted_species",
    "random_intercept",
    "random_intercept_se",
    "confirmation_rate"
  ),
  "Primary M3 species effects"
)

figure_s1_data <- species_effects_m3 %>%
  dplyr::mutate(
    lower = .data$random_intercept - 1.96 * .data$random_intercept_se,
    upper = .data$random_intercept + 1.96 * .data$random_intercept_se,
    predicted_species = forcats::fct_reorder(
      .data$predicted_species,
      .data$random_intercept
    )
  )

figure_s1 <- ggplot2::ggplot(
  figure_s1_data,
  ggplot2::aes(
    x = .data$random_intercept,
    y = .data$predicted_species,
    colour = .data$confirmation_rate
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.55
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = .data$lower,
      xmax = .data$upper
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.38,
    alpha = 0.78
  ) +
  ggplot2::geom_point(size = 1.65) +
  ggplot2::scale_colour_gradient2(
    low = "#D55E00",
    mid = "#F7F7F7",
    high = "#0072B2",
    midpoint = 0.5,
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Primary M3 conditional predicted-species random intercepts",
    x = "Conditional random intercept (log-odds scale)",
    y = "Predicted species",
    colour = "Observed\nconfirmation"
  ) +
  supp_axis_theme(8, legend_position = "right") +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 4.7),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0)
  )

save_figure(
  figure_s1,
  file.path(
    figure_dir,
    "Figure_S01_M3_all_species_random_intercepts"
  ),
  width = 10,
  height = 22,
  dpi = 400
)


# ------------------------------------------------------------
# Figure S2: M4 random-slope sensitivity
# ------------------------------------------------------------

species_effects_m4 <- read_required_csv(
  file.path(
    paths$results_root,
    "02_calibration",
    "tables",
    "calibration_species_effects_m4.csv"
  ),
  "Primary M4 species effects"
)

assert_columns(
  species_effects_m4,
  c(
    "random_intercept",
    "random_score_slope",
    "validation_n",
    "confirmation_rate"
  ),
  "Primary M4 species effects"
)

figure_s2_data <- species_effects_m4 %>%
  dplyr::filter(
    is.finite(.data$random_intercept),
    is.finite(.data$random_score_slope)
  )

figure_s2 <- ggplot2::ggplot(
  figure_s2_data,
  ggplot2::aes(
    x = .data$random_intercept,
    y = .data$random_score_slope,
    size = .data$validation_n,
    colour = .data$confirmation_rate
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.55
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.55
  ) +
  ggplot2::geom_point(alpha = 0.80) +
  ggplot2::scale_colour_gradient2(
    low = "#D55E00",
    mid = "#F7F7F7",
    high = "#0072B2",
    midpoint = 0.5,
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::scale_size_continuous(range = c(1.8, 6.5)) +
  ggplot2::labs(
    title = "Primary M4 random-slope sensitivity",
    x = "Conditional random-intercept deviation",
    y = "Conditional random score-slope deviation",
    colour = "Observed\nconfirmation",
    size = "Validated n"
  ) +
  supp_axis_theme(9, legend_position = "bottom")

save_figure(
  figure_s2,
  file.path(
    figure_dir,
    "Figure_S02_M4_random_slope_sensitivity"
  ),
  width = 10.5,
  height = 6.5,
  dpi = 400
)


# ------------------------------------------------------------
# Figure S3: application of calibrated probabilities
# ------------------------------------------------------------

calibrated_candidates <- read_rds_or_csv(
  file.path(
    paths$processed_root,
    "calibration",
    "calibrated_candidates.rds"
  ),
  file.path(
    paths$processed_root,
    "calibration",
    "calibrated_candidates_minimal.csv"
  ),
  "calibrated event-species candidates"
)

assert_columns(
  calibrated_candidates,
  c(
    "actual_multi_event",
    "p_primary",
    "p_primary_m2",
    "p_audio_only"
  ),
  "calibrated event-species candidates"
)

species_probability_summary <- read_required_csv(
  file.path(
    paths$results_root,
    "04_probability_assignment",
    "tables",
    "species_probability_summary.csv"
  ),
  "species probability summary"
)

probability_plot_data <- calibrated_candidates %>%
  dplyr::select(
    actual_multi_event,
    Primary_M3 = p_primary,
    Primary_M2 = p_primary_m2,
    Audio_only = p_audio_only
  ) %>%
  tidyr::pivot_longer(
    cols = c(Primary_M3, Primary_M2, Audio_only),
    names_to = "scenario",
    values_to = "probability"
  ) %>%
  dplyr::mutate(
    scenario_label = dplyr::recode(
      .data$scenario,
      Primary_M3 = "Primary M3: species-aware",
      Primary_M2 = "Primary M2: pooled",
      Audio_only = "Audio-only sensitivity"
    ),
    scenario_label = factor(
      .data$scenario_label,
      levels = probability_scenario_order
    ),
    event_type = dplyr::if_else(
      .data$actual_multi_event == 1L,
      "Multi-candidate event",
      "Single-candidate event"
    )
  )

make_probability_density_panel <- function(event_value, panel_title) {
  probability_plot_data %>%
    dplyr::filter(.data$event_type == event_value) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$probability,
        colour = .data$scenario_label,
        fill = .data$scenario_label
      )
    ) +
    ggplot2::geom_density(
      alpha = 0.12,
      linewidth = 0.90,
      adjust = 1.1
    ) +
    ggplot2::scale_colour_manual(
      values = palette$scenario,
      breaks = probability_scenario_order,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = palette$scenario,
      breaks = probability_scenario_order,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      title = panel_title,
      x = "Candidate-level correctness probability",
      y = "Density",
      colour = "Probability scenario",
      fill = "Probability scenario"
    ) +
    supp_axis_theme(9, legend_position = "bottom")
}

panel_s3_a <- make_probability_density_panel(
  "Multi-candidate event",
  "A. Multi-candidate events"
)

panel_s3_b <- make_probability_density_panel(
  "Single-candidate event",
  "B. Single-candidate events"
)

species_difference_data <- species_probability_summary %>%
  dplyr::slice_max(
    order_by = abs(.data$expected_difference_m2_minus_m3),
    n = 25,
    with_ties = FALSE
  ) %>%
  dplyr::mutate(
    predicted_species = forcats::fct_reorder(
      .data$predicted_species,
      .data$expected_difference_m2_minus_m3
    ),
    direction = dplyr::if_else(
      .data$expected_difference_m2_minus_m3 >= 0,
      "Primary M2 higher",
      "Primary M3 higher"
    )
  )

panel_s3_c <- ggplot2::ggplot(
  species_difference_data,
  ggplot2::aes(
    x = .data$expected_difference_m2_minus_m3,
    y = .data$predicted_species,
    fill = .data$direction
  )
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_vline(
    xintercept = 0,
    colour = "black",
    linewidth = 0.50
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Primary M2 higher" = "#D55E00",
      "Primary M3 higher" = "#0072B2"
    )
  ) +
  ggplot2::labs(
    title = "C. Predicted species with the largest pooled-versus-species-aware differences",
    x = "Expected correct candidate-record total: Primary M2 minus Primary M3",
    y = "Predicted species",
    fill = NULL
  ) +
  supp_axis_theme(8.5, legend_position = "bottom") +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))

figure_s3 <- (
  (panel_s3_a | panel_s3_b) /
    panel_s3_c
) +
  patchwork::plot_layout(
    guides = "collect",
    heights = c(0.90, 1.10)
  ) +
  patchwork::plot_annotation(
    title = "Application of calibrated candidate-level correctness probabilities",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12.5,
        hjust = 0
      )
    )
  ) &
  ggplot2::theme(legend.position = "bottom")

save_figure(
  figure_s3,
  file.path(
    figure_dir,
    "Figure_S03_probability_assignment"
  ),
  width = 12,
  height = 9.5,
  dpi = 400
)


# ------------------------------------------------------------
# Figure S4: coverage correlations
# ------------------------------------------------------------

coverage_correlations <- read_required_csv(
  file.path(
    paths$results_root,
    "07_sampling_coverage",
    "tables",
    "diversity_coverage_correlations.csv"
  ),
  "sampling-coverage correlations"
)

figure_s4_data <- coverage_correlations %>%
  dplyr::mutate(
    q_label = factor(
      .data$q_label,
      levels = rev(q_order)
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
  )

figure_s4 <- ggplot2::ggplot(
  figure_s4_data,
  ggplot2::aes(
    x = .data$proxy_label,
    y = .data$q_label,
    fill = .data$spearman_rho
  )
) +
  ggplot2::geom_tile(colour = "black", linewidth = 0.45) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$label),
    fontface = "bold"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  ggplot2::labs(
    title = "Spearman correlations between monthly diversity and coverage proxies",
    x = "Sampling-coverage proxy",
    y = "Hill diversity order",
    fill = "Spearman rho"
  ) +
  ggplot2::coord_fixed(ratio = 0.70) +
  supp_heatmap_theme(9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)
  )

save_figure(
  figure_s4,
  file.path(
    figure_dir,
    "Figure_S04_sampling_coverage_correlations"
  ),
  width = 9.5,
  height = 5.2,
  dpi = 400
)


# ------------------------------------------------------------
# Figure S5: higher-coverage robustness
# ------------------------------------------------------------

coverage_robustness <- read_required_csv(
  file.path(
    paths$results_root,
    "08_coverage_robustness",
    "tables",
    "method_bias_robustness.csv"
  ),
  "low-coverage robustness results"
)

figure_s5_data <- coverage_robustness %>%
  dplyr::select(
    method_label,
    q_label,
    component_label,
    all_percentage_difference,
    high_percentage_difference
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      all_percentage_difference,
      high_percentage_difference
    ),
    names_to = "month_set",
    values_to = "percentage_difference"
  ) %>%
  dplyr::mutate(
    month_set = dplyr::recode(
      .data$month_set,
      all_percentage_difference = "All represented months",
      high_percentage_difference = "Higher-coverage months"
    ),
    month_set = factor(
      .data$month_set,
      levels = c(
        "All represented months",
        "Higher-coverage months"
      )
    ),
    method_label = factor(
      .data$method_label,
      levels = rev(deterministic_order)
    ),
    component_label = factor(
      .data$component_label,
      levels = component_order
    ),
    q_label = factor(
      .data$q_label,
      levels = q_order
    )
  )

figure_s5 <- ggplot2::ggplot(
  figure_s5_data,
  ggplot2::aes(
    x = .data$percentage_difference,
    y = .data$method_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.55
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      group = interaction(
        .data$method_label,
        .data$q_label,
        .data$component_label
      )
    ),
    colour = "grey65",
    linewidth = 0.70
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = .data$month_set),
    size = 2.7
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$component_label),
    cols = ggplot2::vars(.data$q_label),
    scales = "free_x"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "All represented months" = "#0072B2",
      "Higher-coverage months" = "#D55E00"
    )
  ) +
  ggplot2::labs(
    title = "Robustness of deterministic-workflow differences after excluding low-coverage months",
    x = "Percentage difference from Primary M3",
    y = "Identification workflow",
    colour = "Month set"
  ) +
  supp_axis_theme(
    8.5,
    legend_position = "bottom",
    panel_border = TRUE
  ) +
  ggplot2::theme(
    strip.text.y = ggplot2::element_text(angle = 0),
    axis.line = ggplot2::element_blank()
  )

save_figure(
  figure_s5,
  file.path(
    figure_dir,
    "Figure_S05_low_coverage_robustness"
  ),
  width = 14,
  height = 9.5,
  dpi = 400
)


# ------------------------------------------------------------
# Figure S6: analytical approximation error
# ------------------------------------------------------------

analytical_monthly <- read_required_csv(
  file.path(
    paths$results_root,
    "09_analytical_benchmark",
    "tables",
    "analytical_vs_monte_carlo_monthly.csv"
  ),
  "monthly analytical benchmark comparison"
)

analytical_partition <- read_required_csv(
  file.path(
    paths$results_root,
    "09_analytical_benchmark",
    "tables",
    "analytical_vs_monte_carlo_partition.csv"
  ),
  "partition analytical benchmark comparison"
)

figure_s6_data <- dplyr::bind_rows(
  analytical_monthly %>%
    dplyr::group_by(.data$q, .data$q_label) %>%
    dplyr::summarise(
      comparison = "Monthly",
      absolute_percent_error = mean(
        abs(.data$analytical_percent_error),
        na.rm = TRUE
      ),
      .groups = "drop"
    ),
  analytical_partition %>%
    dplyr::filter(.data$weighting == "equal_month") %>%
    dplyr::transmute(
      q,
      q_label,
      comparison = paste0("Partition: ", .data$component_label),
      absolute_percent_error = abs(.data$analytical_percent_error)
    )
) %>%
  dplyr::mutate(
    q_label = factor(
      .data$q_label,
      levels = q_order
    ),
    # Discrete y axes are drawn bottom-to-top. These levels produce
    # top-to-bottom: Monthly, Alpha, Beta, Gamma.
    comparison = factor(
      .data$comparison,
      levels = c(
        "Partition: Gamma",
        "Partition: Beta",
        "Partition: Alpha",
        "Monthly"
      )
    ),
    label = sprintf("%.2f%%", .data$absolute_percent_error)
  )

figure_s6 <- ggplot2::ggplot(
  figure_s6_data,
  ggplot2::aes(
    x = .data$q_label,
    y = .data$comparison,
    fill = .data$absolute_percent_error
  )
) +
  ggplot2::geom_tile(colour = "black", linewidth = 0.45) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$label),
    fontface = "bold"
  ) +
  ggplot2::scale_fill_gradient(
    low = "#F7FBFF",
    high = "#B2182B",
    labels = scales::label_number(suffix = "%", accuracy = 0.1)
  ) +
  ggplot2::labs(
    title = "Absolute error of the analytical / plug-in benchmark",
    x = "Hill diversity order",
    y = "Comparison level",
    fill = "Absolute percentage error"
  ) +
  supp_heatmap_theme(9)

save_figure(
  figure_s6,
  file.path(
    figure_dir,
    "Figure_S06_analytical_approximation_error"
  ),
  width = 9.5,
  height = 5.5,
  dpi = 400
)


# ============================================================
# SUPPLEMENTARY TABLE DATA AND LATEX FRAGMENTS
# ============================================================


# ------------------------------------------------------------
# Table S1: validation coverage by predicted species
# ------------------------------------------------------------

table_s1 <- read_required_csv(
  file.path(
    paths$results_root,
    "01_validation",
    "tables",
    "validation_species_summary.csv"
  ),
  "validation species summary"
) %>%
  dplyr::arrange(.data$predicted_species)

readr::write_csv(
  table_s1,
  file.path(table_csv_dir, "Table_S01_validation_by_species.csv")
)

rows_s1 <- vapply(
  seq_len(nrow(table_s1)),
  function(i) {
    row <- table_s1[i, ]
    cells <- c(
      latex_escape(row$predicted_species),
      fmt_int(row$species_validation_n),
      fmt_int(row$species_validation_events),
      fmt_int(row$species_correct_n),
      fmt_int(row$species_incorrect_n),
      fmt_num(row$species_confirmation_rate, 2),
      fmt_num(row$species_mean_score, 3),
      fmt_num(row$species_min_score, 3),
      fmt_num(row$species_max_score, 3),
      fmt_int(row$species_score_bin_n),
      fmt_int(row$species_multi_event_n),
      fmt_int(row$species_audio_only_n),
      fmt_int(row$species_ecological_prior_n),
      latex_escape(row$validation_depth),
      latex_escape(row$random_slope_support)
    )
    paste0(paste(cells, collapse = " & "), " \\\\")
  },
  FUN.VALUE = character(1)
)

write_longtable_fragment(
  file.path(table_tex_dir, "Table_S01_validation_by_species.tex"),
  table_number = "S1",
  caption = paste0(
    "Predicted-species-level coverage of the final manual-validation dataset. ",
    "Each row represents one of the 229 predicted species in the Primary ",
    "calibration analysis. Columns report validation depth, observed binary ",
    "outcomes, BirdNET confidence-score coverage, multi-candidate event ",
    "representation and validation-evidence basis. Confirmation proportions ",
    "describe the deliberately stratified validation records and are not ",
    "unbiased estimates of complete-dataset species-specific BirdNET accuracy. ",
    "The random-slope-support category is a diagnostic of available validation ",
    "information rather than a formal hypothesis test."
  ),
  label = "tab:supp-validation",
  column_spec = "L{3.6cm}*{12}{C{0.90cm}}L{1.45cm}L{1.65cm}",
  header_cells = c(
    "Predicted species",
    "$n$",
    "Evts.",
    "Corr.",
    "Incorr.",
    "Conf.",
    "Mean",
    "Min",
    "Max",
    "Bins",
    "Multi",
    "Audio",
    "Eco.",
    "Depth",
    "Slope"
  ),
  body_rows = rows_s1,
  landscape = TRUE,
  font_command = "\\fontsize{7.0}{8.0}\\selectfont",
  tabcolsep = "1.25pt",
  arraystretch = "1.10"
)


# ------------------------------------------------------------
# Table S2: calibration fixed effects
# ------------------------------------------------------------

table_s2 <- read_required_csv(
  file.path(
    paths$results_root,
    "02_calibration",
    "tables",
    "calibration_fixed_effects.csv"
  ),
  "calibration fixed effects"
) %>%
  dplyr::mutate(
    analysis_set = factor(
      .data$analysis_set,
      levels = c("Primary", "Audio_only")
    ),
    model_id = factor(
      .data$model_id,
      levels = c("M1", "M2", "M3", "M4")
    ),
    term = factor(
      .data$term,
      levels = c("(Intercept)", "logit_score", "actual_multi_event")
    )
  ) %>%
  dplyr::arrange(.data$analysis_set, .data$model_id, .data$term) %>%
  dplyr::mutate(
    analysis_set = dplyr::recode(
      as.character(.data$analysis_set),
      Audio_only = "Audio-only"
    ),
    model_id = as.character(.data$model_id),
    term = as.character(.data$term)
  )

readr::write_csv(
  table_s2,
  file.path(table_csv_dir, "Table_S02_calibration_fixed_effects.csv")
)

term_display <- function(term) {
  dplyr::recode(
    as.character(term),
    `(Intercept)` = "Intercept",
    logit_score = "logit score",
    actual_multi_event = "multi-candidate event",
    .default = as.character(term)
  )
}

rows_s2 <- vapply(
  seq_len(nrow(table_s2)),
  function(i) {
    row <- table_s2[i, ]
    cells <- c(
      latex_escape(row$analysis_set),
      latex_escape(row$model_id),
      latex_escape(row$model_label),
      latex_escape(term_display(row$term)),
      fmt_num(row$estimate, 3),
      fmt_num(row$std_error, 3),
      fmt_num(row$statistic, 2),
      fmt_p_value(row$p_value),
      fmt_num(row$conf_low, 3),
      fmt_num(row$conf_high, 3),
      fmt_num(row$odds_ratio, 3),
      fmt_num(row$odds_ratio_low, 3),
      fmt_num(row$odds_ratio_high, 3)
    )
    paste0(paste(cells, collapse = " & "), " \\\\")
  },
  FUN.VALUE = character(1)
)

write_longtable_fragment(
  file.path(table_tex_dir, "Table_S02_calibration_fixed_effects.tex"),
  table_number = "S2",
  caption = paste0(
    "Fixed-effect estimates from the fitted calibration models. Estimates ",
    "are reported on the log-odds scale for Primary M1--M4 and for the ",
    "Audio-only M2--M3 sensitivity models. The table gives model-based ",
    "standard errors, test statistics, $p$-values, 95\\% intervals and ",
    "corresponding odds ratios. The multi-candidate event coefficient compares ",
    "event-species candidates from acoustic events containing more than one ",
    "BirdNET candidate label with those from single-candidate events. Random-",
    "effect standard deviations are reported separately."
  ),
  label = "tab:supp-fixed-effects",
  column_spec = "L{1.5cm}L{1.0cm}L{2.6cm}L{2.4cm}*{9}{C{1.05cm}}",
  header_cells = c(
    "Analysis",
    "Model",
    "Model label",
    "Term",
    "Est.",
    "SE",
    "Stat.",
    "$p$",
    "CI low",
    "CI high",
    "OR",
    "OR low",
    "OR high"
  ),
  body_rows = rows_s2,
  landscape = TRUE,
  font_command = "\\scriptsize",
  tabcolsep = "2.4pt",
  arraystretch = "1.02"
)


# ------------------------------------------------------------
# Table S3: full equal-month workflow comparison
# ------------------------------------------------------------

table_s3 <- read_required_csv(
  file.path(
    paths$results_root,
    "06_method_comparison",
    "tables",
    "method_difference_from_m3_partition.csv"
  ),
  "full workflow comparison"
) %>%
  dplyr::filter(
    .data$weighting == "equal_month",
    .data$method_label %in% workflow_order
  ) %>%
  dplyr::mutate(
    method_label = factor(.data$method_label, levels = workflow_order),
    q_label = factor(.data$q_label, levels = q_order),
    component_label = factor(
      .data$component_label,
      levels = component_order
    )
  ) %>%
  dplyr::arrange(
    .data$method_label,
    .data$q_label,
    .data$component_label
  ) %>%
  dplyr::mutate(
    method_label = as.character(.data$method_label),
    q_label = as.character(.data$q_label),
    component_label = as.character(.data$component_label)
  )

readr::write_csv(
  table_s3,
  file.path(table_csv_dir, "Table_S03_full_workflow_comparison.csv")
)

rows_s3 <- vapply(
  seq_len(nrow(table_s3)),
  function(i) {
    row <- table_s3[i, ]
    cells <- c(
      latex_escape(row$method_label),
      latex_escape(row$q_label),
      latex_escape(row$component_label),
      fmt_num(row$estimate, 3),
      fmt_num(row$lower_95, 3),
      fmt_num(row$upper_95, 3),
      fmt_num(row$reference_estimate, 3),
      fmt_num(row$absolute_difference, 3),
      fmt_num(row$percentage_difference, 2)
    )
    paste0(paste(cells, collapse = " & "), " \\\\")
  },
  FUN.VALUE = character(1)
)

write_longtable_fragment(
  file.path(table_tex_dir, "Table_S03_full_workflow_comparison.tex"),
  table_number = "S3",
  caption = paste0(
    "Complete equal-month comparison of alternative probability and ",
    "deterministic identification workflows with Primary M3. Within each ",
    "diversity order, components are presented as alpha diversity, ",
    "multiplicative beta diversity and gamma diversity. Probabilistic ",
    "scenarios are summarised from 1,000 Monte Carlo simulations, whereas ",
    "deterministic workflows provide calculated point estimates. Absolute and ",
    "percentage differences are measured relative to Primary M3. Primary M2 ",
    "and the Audio-only sensitivity analysis are probability-model ",
    "sensitivities rather than deterministic identification rules."
  ),
  label = "tab:supp-workflows",
  column_spec = "L{3.4cm}L{2.5cm}L{1.6cm}*{6}{C{1.55cm}}",
  header_cells = c(
    "Method",
    "Diversity order",
    "Component",
    "Estimate",
    "Lower 95\\%",
    "Upper 95\\%",
    "M3 ref.",
    "Abs. diff.",
    "\\% diff."
  ),
  body_rows = rows_s3,
  landscape = TRUE,
  font_command = "\\scriptsize",
  tabcolsep = "3.0pt",
  arraystretch = "0.98"
)


# ------------------------------------------------------------
# Table S4: coverage correlations
# ------------------------------------------------------------

table_s4 <- coverage_correlations %>%
  dplyr::mutate(
    q_label = factor(.data$q_label, levels = q_order),
    proxy_label = factor(
      .data$proxy_label,
      levels = c(
        "Represented audio hours",
        "Active detection hours",
        "Active detection days"
      )
    )
  ) %>%
  dplyr::arrange(.data$q_label, .data$proxy_label) %>%
  dplyr::mutate(
    q_label = as.character(.data$q_label),
    proxy_label = as.character(.data$proxy_label)
  )

readr::write_csv(
  table_s4,
  file.path(table_csv_dir, "Table_S04_coverage_correlations.csv")
)

rows_s4 <- vapply(
  seq_len(nrow(table_s4)),
  function(i) {
    row <- table_s4[i, ]
    cells <- c(
      latex_escape(row$q_label),
      latex_escape(row$proxy_label),
      fmt_int(row$n_months),
      fmt_num(row$spearman_rho, 3),
      fmt_p_value(row$p_value)
    )
    paste0(paste(cells, collapse = " & "), " \\\\")
  },
  FUN.VALUE = character(1)
)

write_small_table_fragment(
  file.path(table_tex_dir, "Table_S04_coverage_correlations.tex"),
  caption = paste0(
    "Spearman rank correlations between Primary M3 monthly Hill diversity ",
    "and sampling-coverage proxies across the 23 represented months. ",
    "Represented audio hours include only five-minute source recordings ",
    "appearing in the raw detection export and therefore do not reconstruct ",
    "complete recorder uptime. These analyses are diagnostic associations and ",
    "should not be interpreted as causal effects of sampling coverage on ",
    "ecological diversity."
  ),
  label = "tab:supp-coverage",
  column_spec = "L{3.8cm}L{4.2cm}C{1.6cm}C{1.8cm}C{2.3cm}",
  header_cells = c(
    "Diversity order",
    "Coverage proxy",
    "Months",
    "Spearman $\\rho$",
    "$p$"
  ),
  body_rows = rows_s4
)


# ------------------------------------------------------------
# Table S5: low-coverage robustness
# ------------------------------------------------------------

table_s5 <- coverage_robustness %>%
  dplyr::filter(.data$method_label %in% deterministic_order) %>%
  dplyr::mutate(
    method_label = factor(.data$method_label, levels = deterministic_order),
    q_label = factor(.data$q_label, levels = q_order),
    component_label = factor(
      .data$component_label,
      levels = component_order
    )
  ) %>%
  dplyr::arrange(
    .data$method_label,
    .data$q_label,
    .data$component_label
  ) %>%
  dplyr::transmute(
    method = .data$method,
    method_label = as.character(.data$method_label),
    q = .data$q,
    q_label = as.character(.data$q_label),
    component = .data$component,
    component_label = as.character(.data$component_label),
    all_estimate = .data$all_estimate,
    all_reference = .data$all_reference,
    all_percentage_difference = .data$all_percentage_difference,
    high_estimate = .data$high_estimate,
    high_reference = .data$high_reference,
    high_percentage_difference = .data$high_percentage_difference,
    absolute_percentage_difference_change_points =
      .data$absolute_bias_change_points,
    direction_stable = .data$direction_stable
  )

readr::write_csv(
  table_s5,
  file.path(table_csv_dir, "Table_S05_low_coverage_robustness.csv")
)

rows_s5 <- vapply(
  seq_len(nrow(table_s5)),
  function(i) {
    row <- table_s5[i, ]
    cells <- c(
      latex_escape(row$method_label),
      latex_escape(row$q_label),
      latex_escape(row$component_label),
      fmt_num(row$all_estimate, 3),
      fmt_num(row$all_reference, 3),
      fmt_num(row$all_percentage_difference, 2),
      fmt_num(row$high_estimate, 3),
      fmt_num(row$high_reference, 3),
      fmt_num(row$high_percentage_difference, 2),
      fmt_num(row$absolute_percentage_difference_change_points, 2),
      ifelse(isTRUE(row$direction_stable), "Yes", "No")
    )
    paste0(paste(cells, collapse = " & "), " \\\\")
  },
  FUN.VALUE = character(1)
)

write_longtable_fragment(
  file.path(table_tex_dir, "Table_S05_low_coverage_robustness.tex"),
  table_number = "S5",
  caption = paste0(
    "Robustness of deterministic-workflow differences from Primary M3 after ",
    "exclusion of low-coverage months. Within each diversity order, ",
    "components are presented as alpha diversity, multiplicative beta ",
    "diversity and gamma diversity. Results compare all 23 represented months ",
    "with the 14 higher-coverage months. Direction stable indicates whether ",
    "the sign of the difference from Primary M3 was unchanged; the absolute ",
    "change column reports the change in percentage difference, measured in ",
    "percentage points. The analysis evaluates robustness of method ",
    "comparisons rather than correction of absolute diversity for incomplete ",
    "recording coverage."
  ),
  label = "tab:supp-coverage-robustness",
  column_spec = "L{3.0cm}L{2.4cm}L{1.5cm}*{6}{C{1.25cm}}C{1.65cm}C{1.35cm}",
  header_cells = c(
    "Method",
    "Order",
    "Component",
    "All est.",
    "All M3",
    "All \\% diff.",
    "High est.",
    "High M3",
    "High \\% diff.",
    "\\makecell{Abs.\\\\change\\\\(p.p.)}",
    "\\makecell{Stable\\\\sign}"
  ),
  body_rows = rows_s5,
  landscape = TRUE,
  font_command = "\\scriptsize",
  tabcolsep = "2.2pt",
  arraystretch = "0.98"
)


# ------------------------------------------------------------
# Table S6: abundance-weighted partition sensitivity
# ------------------------------------------------------------

table_s6 <- read_required_csv(
  file.path(
    paths$results_root,
    "05_monte_carlo",
    "tables",
    "monte_carlo_partition_summary.csv"
  ),
  "Monte Carlo partition summary"
) %>%
  dplyr::filter(
    .data$weighting == "abundance_weighted",
    .data$scenario_label %in% probability_scenario_order
  ) %>%
  dplyr::mutate(
    scenario_label = factor(
      .data$scenario_label,
      levels = probability_scenario_order
    ),
    q_label = factor(.data$q_label, levels = q_order),
    component_label = factor(
      .data$component_label,
      levels = component_order
    )
  ) %>%
  dplyr::arrange(
    .data$scenario_label,
    .data$q_label,
    .data$component_label
  ) %>%
  dplyr::mutate(
    scenario_label = as.character(.data$scenario_label),
    q_label = as.character(.data$q_label),
    component_label = as.character(.data$component_label)
  )

readr::write_csv(
  table_s6,
  file.path(table_csv_dir, "Table_S06_abundance_weighted_partitions.csv")
)

rows_s6 <- vapply(
  seq_len(nrow(table_s6)),
  function(i) {
    row <- table_s6[i, ]
    cells <- c(
      latex_escape(row$scenario_label),
      latex_escape(row$q_label),
      latex_escape(row$component_label),
      fmt_num(row$mean, 3),
      fmt_num(row$sd, 3),
      fmt_num(row$median, 3),
      fmt_num(row$lower_95, 3),
      fmt_num(row$upper_95, 3)
    )
    paste0(paste(cells, collapse = " & "), " \\\\")
  },
  FUN.VALUE = character(1)
)

write_longtable_fragment(
  file.path(table_tex_dir, "Table_S06_abundance_weighted_partitions.tex"),
  table_number = "S6",
  caption = paste0(
    "Abundance-weighted sensitivity analysis of temporal Hill-diversity ",
    "partitioning for Primary M3, Primary M2 and the Audio-only sensitivity ",
    "analysis. Within each diversity order, components are presented as alpha ",
    "diversity, multiplicative beta diversity and gamma diversity. Month ",
    "weights were proportional to simulation-specific retained candidate-record ",
    "totals. The table reports Monte Carlo means, standard deviations, medians ",
    "and 2.5\\% and 97.5\\% quantiles from 1,000 simulations. Candidate-record ",
    "totals do not represent numbers of individual birds; equal-month weighting ",
    "remains the primary temporal partition."
  ),
  label = "tab:supp-abundance-weighted",
  column_spec = "L{4.0cm}L{2.5cm}L{1.6cm}*{5}{C{1.65cm}}",
  header_cells = c(
    "Scenario",
    "Diversity order",
    "Component",
    "Mean",
    "SD",
    "Median",
    "Lower 95\\%",
    "Upper 95\\%"
  ),
  body_rows = rows_s6,
  landscape = TRUE,
  font_command = "\\scriptsize",
  tabcolsep = "3.5pt",
  arraystretch = "1.00"
)


# ------------------------------------------------------------
# 4. Completion message
# ------------------------------------------------------------

message("Supplementary outputs rebuilt successfully.")
message("Figures overwritten in: ", figure_dir)
message("Supplementary CSV tables overwritten in: ", table_csv_dir)
message("Supplementary LaTeX fragments written to: ", table_tex_dir)

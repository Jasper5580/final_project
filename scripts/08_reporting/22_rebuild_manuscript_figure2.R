# ============================================================
# 22_rebuild_manuscript_figure2.R
#
# Purpose:
# Rebuild manuscript Figure 2 after model fitting and
# event-grouped cross-validation.
#
# This script does NOT refit any model.
# It uses the saved outputs from scripts 13 and 14 and
# overwrites:
#
# results/10_final_outputs/main_figures/
# Figure_02_calibration_and_validation.pdf
# Figure_02_calibration_and_validation.png
#
# Final panel structure:
#   A. Multi-candidate events
#   B. Single-candidate events
#   C. Largest Primary M3 species baseline deviations
#   D. Event-grouped out-of-fold prediction error
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)
palette <- silwood_palette()


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

calibration_component_file <- file.path(
  paths$results_root,
  "02_calibration",
  "models",
  "calibration_plot_components.rds"
)

cv_summary_file <- file.path(
  paths$results_root,
  "03_model_validation",
  "tables",
  "event_cv_model_summary.csv"
)

if (!file.exists(calibration_component_file)) {
  stop(
    "Calibration plot components not found: ",
    calibration_component_file,
    "\nRun script 13 before rebuilding Figure 2."
  )
}

if (!file.exists(cv_summary_file)) {
  stop(
    "Event-grouped CV summary not found: ",
    cv_summary_file,
    "\nRun script 14 before rebuilding Figure 2."
  )
}


# ------------------------------------------------------------
# 3. Load saved calibration components
# ------------------------------------------------------------

calibration_components <- readRDS(calibration_component_file)

if (
  !"calibration_panel" %in% names(calibration_components) ||
  !"random_intercept_panel" %in% names(calibration_components)
) {
  stop(
    "The saved calibration plot object does not contain the ",
    "required calibration_panel and random_intercept_panel."
  )
}

# The original calibration_panel was produced in script 13.
# Its main data contain the fitted curves.
curve_predictions <- calibration_components$calibration_panel$data

# The third layer of the original plot is the explicitly supplied
# observed score-bin data used for the white circles.
observed_bins <- calibration_components$calibration_panel$layers[[3]]$data

# The saved random-intercept plot contains the 12 most positive
# and 12 most negative Primary M3 conditional random intercepts.
species_focus <- calibration_components$random_intercept_panel$data


# ------------------------------------------------------------
# 4. Common plot settings
# ------------------------------------------------------------

model_order <- c(
  "M1: Score",
  "M2: Score + event",
  "M3: Species-aware",
  "M4: Random score slope"
)

curve_predictions <- curve_predictions %>%
  dplyr::mutate(
    model_label = factor(
      .data$model_label,
      levels = model_order
    )
  )

observed_bins <- observed_bins %>%
  dplyr::mutate(
    event_type = as.character(.data$event_type)
  )

# Add explicit solid x- and y-axis lines to every panel.
axis_theme <- function(base_size = 10) {
  theme_silwood(base_size) +
    ggplot2::theme(
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
      panel.border = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 1,
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank()
    )
}


# ------------------------------------------------------------
# 5. Panels A and B:
#    validation-calibrated correctness curves
# ------------------------------------------------------------

make_calibration_panel <- function(
    event_label,
    panel_title
) {
  
  curve_data <- curve_predictions %>%
    dplyr::filter(.data$event_type == event_label)
  
  observed_data <- observed_bins %>%
    dplyr::filter(.data$event_type == event_label)
  
  ggplot2::ggplot(
    curve_data,
    ggplot2::aes(
      x = .data$birdnet_score,
      y = .data$probability,
      colour = .data$model_label,
      fill = .data$model_label
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower,
        ymax = .data$upper
      ),
      alpha = 0.10,
      colour = NA
    ) +
    ggplot2::geom_line(
      linewidth = 0.95
    ) +
    ggplot2::geom_point(
      data = observed_data,
      ggplot2::aes(
        x = .data$mean_score,
        y = .data$confirmation_rate,
        size = .data$n_records
      ),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      colour = "black",
      stroke = 0.65
    ) +
    ggplot2::scale_colour_manual(
      values = palette$model,
      breaks = model_order,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = palette$model,
      breaks = model_order,
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      range = c(2.2, 6),
      name = "Validated n"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(
        0.45,
        0.60,
        0.75,
        0.90,
        1.00
      ),
      limits = c(0.45, 1.00),
      expand = ggplot2::expansion(
        mult = c(0, 0.01)
      )
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(
        0,
        0.25,
        0.50,
        0.75,
        1.00
      ),
      labels = scales::label_percent(
        accuracy = 1
      ),
      limits = c(0, 1),
      expand = ggplot2::expansion(
        mult = c(0, 0.02)
      )
    ) +
    ggplot2::labs(
      title = panel_title,
      x = "BirdNET confidence score",
      y = "Probability candidate label is correct",
      colour = "Model",
      fill = "Model"
    ) +
    axis_theme(10)
}


panel_A <- make_calibration_panel(
  event_label = "Multi-candidate event",
  panel_title = "A. Multi-candidate events"
)

panel_B <- make_calibration_panel(
  event_label = "Single-candidate event",
  panel_title = "B. Single-candidate events"
)


# ------------------------------------------------------------
# 6. Panel C:
#    Primary M3 species random-intercept deviations
# ------------------------------------------------------------

species_focus <- species_focus %>%
  dplyr::mutate(
    predicted_species = forcats::fct_reorder(
      as.character(.data$predicted_species),
      .data$random_intercept
    )
  )

panel_C <- ggplot2::ggplot(
  species_focus,
  ggplot2::aes(
    x = .data$random_intercept,
    y = .data$predicted_species,
    colour = .data$direction,
    size = .data$validation_n
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.55
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      xmin = .data$lower,
      xmax = .data$upper
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    alpha = 0.90
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Above population baseline" = "#0072B2",
      "Below population baseline" = "#D55E00"
    )
  ) +
  ggplot2::scale_size_continuous(
    range = c(2, 5)
  ) +
  ggplot2::labs(
    title = "C. Largest Primary M3 species baseline deviations",
    x = "Conditional random intercept (log-odds scale)",
    y = "Predicted species",
    colour = NULL,
    size = "Validated n"
  ) +
  axis_theme(9) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      size = 7
    ),
    legend.position = "bottom"
  )


# ------------------------------------------------------------
# 7. Panel D:
#    event-grouped out-of-fold prediction error
# ------------------------------------------------------------

cv_summary <- readr::read_csv(
  cv_summary_file,
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    model_label = factor(
      .data$model_label,
      levels = model_order
    )
  )

model_axis_labels <- c(
  "M1: Score" = "M1:\nScore",
  "M2: Score + event" = "M2:\nScore + event",
  "M3: Species-aware" = "M3:\nSpecies-aware",
  "M4: Random score slope" = "M4:\nRandom slope"
)


make_cv_metric_panel <- function(
    estimate_column,
    se_column,
    metric_title,
    y_label
) {
  
  plot_data <- cv_summary %>%
    dplyr::transmute(
      model_label = .data$model_label,
      estimate = .data[[estimate_column]],
      standard_error = .data[[se_column]]
    )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$model_label,
      y = .data$estimate,
      colour = .data$model_label
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$estimate -
          .data$standard_error,
        ymax = .data$estimate +
          .data$standard_error
      ),
      width = 0.12,
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      size = 3.2
    ) +
    ggplot2::scale_colour_manual(
      values = palette$model,
      breaks = model_order,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      labels = model_axis_labels
    ) +
    ggplot2::labs(
      title = metric_title,
      x = "Calibration model",
      y = y_label
    ) +
    axis_theme(8.5) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        angle = 25,
        hjust = 1,
        vjust = 1
      )
    )
}


panel_D_brier <- make_cv_metric_panel(
  estimate_column = "mean_brier_score",
  se_column = "se_brier_score",
  metric_title = "Brier score",
  y_label = "Mean Brier score"
)

panel_D_logloss <- make_cv_metric_panel(
  estimate_column = "mean_log_loss",
  se_column = "se_log_loss",
  metric_title = "Log loss",
  y_label = "Mean log loss"
)


# Build the D heading as an actual ggplot panel rather than a nested
# plot_annotation(). Nested patchwork annotations can be dropped when the
# object is inserted into a larger patchwork, which is why the "D." label
# was missing in the previous manuscript figure.
panel_D_header <- ggplot2::ggplot() +
  ggplot2::annotate(
    "text",
    x = 0,
    y = 0.5,
    label = "D. Event-grouped out-of-fold prediction error",
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 3.4
  ) +
  ggplot2::xlim(0, 1) +
  ggplot2::ylim(0, 1) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 0,
      r = 0,
      b = 1,
      l = 2
    )
  )

panel_D_metrics <- (
  panel_D_brier |
    panel_D_logloss
) +
  patchwork::plot_layout(
    widths = c(1, 1)
  )

panel_D <- (
  panel_D_header /
    panel_D_metrics
) +
  patchwork::plot_layout(
    heights = c(0.09, 0.91)
  )


# ------------------------------------------------------------
# 8. Combine four manuscript panels
# ------------------------------------------------------------

top_row <- (
  panel_A |
    panel_B
) +
  patchwork::plot_layout(
    guides = "collect",
    widths = c(1, 1)
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

bottom_row <- (
  panel_C |
    panel_D
) +
  patchwork::plot_layout(
    widths = c(1.08, 0.92)
  )

figure_02 <- (
  top_row /
    bottom_row
) +
  patchwork::plot_layout(
    heights = c(0.92, 1.08)
  )


# ------------------------------------------------------------
# 9. Save directly over the manuscript Figure 2
# ------------------------------------------------------------

manuscript_figure_dir <- silwood_final_dir(paths, "main_figures")

output_stem <- file.path(
  manuscript_figure_dir,
  "Figure_02_calibration_and_validation"
)

save_figure(
  figure_02,
  output_stem,
  width = 13.5,
  height = 10.2,
  dpi = 400
)

message(
  "Rebuilt manuscript Figure 2:"
)

message(
  paste0(
    output_stem,
    ".pdf"
  )
)

message(
  paste0(
    output_stem,
    ".png"
  )
)
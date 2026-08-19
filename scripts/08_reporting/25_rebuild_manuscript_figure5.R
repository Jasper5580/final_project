# ============================================================
# 25_rebuild_manuscript_figure5.R
#
# Purpose:
# Rebuild manuscript Figure 5 from the saved analytical-vs-Monte-Carlo
# monthly comparison without rerunning the analytical framework or
# fixed-probability Monte Carlo simulations.
#
# Presentation changes only:
#   * retain the original monthly analytical / plug-in estimates;
#   * retain the original fixed-probability Monte Carlo means and 95% intervals;
#   * give each q panel explicit solid x- and y-axis lines;
#   * use the same numerical scale on both axes within each panel so the
#     one-to-one agreement line is visually meaningful;
#   * remove explanatory subtitle and caption text from inside the figure;
#   * preserve q-panel order: q = 0, q = 1, q = 2;
#   * overwrite the manuscript Figure 5 PDF and PNG.
#
# Note:
# Figure 5 is a monthly q-order benchmark and does not contain alpha, beta,
# or gamma partition components. The alpha -> beta -> gamma presentation
# rule applies to Table 6 and partition/supplementary outputs, not this figure.
#
# Run after script 21 / the complete pipeline.
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
# 2. Input and output paths
# ------------------------------------------------------------

monthly_analytical_file <- file.path(
  paths$results_root,
  "09_analytical_benchmark",
  "tables",
  "analytical_vs_monte_carlo_monthly.csv"
)

if (!file.exists(monthly_analytical_file)) {
  stop(
    "Analytical-vs-Monte-Carlo monthly comparison not found: ",
    monthly_analytical_file,
    "\nRun script 21 (or the complete pipeline) before rebuilding Figure 5."
  )
}

manuscript_figure_dir <- silwood_final_dir(paths, "main_figures")

output_stem <- file.path(
  manuscript_figure_dir,
  "Figure_05_output2_analytical_benchmark"
)


# ------------------------------------------------------------
# 3. Read saved analytical benchmark results
# ------------------------------------------------------------

monthly_analytical <- readr::read_csv(
  monthly_analytical_file,
  show_col_types = FALSE
)

required_columns <- c(
  "q",
  "q_label",
  "fixed_mean",
  "fixed_lower_95",
  "fixed_upper_95",
  "analytical_estimate"
)

assert_columns(
  monthly_analytical,
  required_columns,
  "Analytical-vs-Monte-Carlo monthly comparison"
)

q_order <- c(
  "q = 0: Richness",
  "q = 1: Exponential Shannon",
  "q = 2: Inverse Simpson"
)

monthly_analytical <- monthly_analytical %>%
  dplyr::filter(
    .data$q_label %in% q_order
  ) %>%
  dplyr::mutate(
    q_label = factor(
      .data$q_label,
      levels = q_order
    )
  )


# ------------------------------------------------------------
# 4. Theme with explicit axes for every panel
# ------------------------------------------------------------

figure5_panel_theme <- function(base_size = 9.5) {
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
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 0.5,
        hjust = 0.5
      ),
      legend.position = "none"
    )
}


# ------------------------------------------------------------
# 5. Build one explicit-axis panel for each Hill order
# ------------------------------------------------------------

make_benchmark_panel <- function(q_label_value, panel_title) {

  panel_data <- monthly_analytical %>%
    dplyr::filter(
      as.character(.data$q_label) == q_label_value
    )

  if (nrow(panel_data) == 0L) {
    stop("No data found for panel: ", q_label_value)
  }

  # Use the same numerical limits on x and y so that the dashed
  # one-to-one line is visually interpretable as exact agreement.
  common_min <- min(
    c(
      panel_data$fixed_lower_95,
      panel_data$analytical_estimate
    ),
    na.rm = TRUE
  )

  common_max <- max(
    c(
      panel_data$fixed_upper_95,
      panel_data$analytical_estimate
    ),
    na.rm = TRUE
  )

  span <- common_max - common_min

  if (!is.finite(span) || span <= 0) {
    span <- max(abs(common_max), 1)
  }

  padding <- 0.04 * span

  common_limits <- c(
    common_min - padding,
    common_max + padding
  )

  ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = .data$fixed_mean,
      y = .data$analytical_estimate,
      colour = .data$q_label
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey35",
      linewidth = 0.70
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = .data$fixed_lower_95,
        xmax = .data$fixed_upper_95
      ),
      orientation = "y",
      width = 0,
      alpha = 0.60,
      linewidth = 0.60
    ) +
    ggplot2::geom_point(
      size = 2.8,
      alpha = 0.90
    ) +
    ggplot2::scale_colour_manual(
      values = palette$q,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      limits = common_limits,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = common_limits,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_fixed(
      ratio = 1
    ) +
    ggplot2::labs(
      title = panel_title,
      x = "Fixed-probability Monte Carlo mean",
      y = "Analytical or plug-in estimate"
    ) +
    figure5_panel_theme(9.5)
}


panel_q0 <- make_benchmark_panel(
  q_label_value = "q = 0: Richness",
  panel_title = "q = 0: Richness"
)

panel_q1 <- make_benchmark_panel(
  q_label_value = "q = 1: Exponential Shannon",
  panel_title = "q = 1: Exponential Shannon"
)

panel_q2 <- make_benchmark_panel(
  q_label_value = "q = 2: Inverse Simpson",
  panel_title = "q = 2: Inverse Simpson"
)


# ------------------------------------------------------------
# 6. Combine the three panels
# ------------------------------------------------------------

figure_05 <- (
  panel_q0 |
    panel_q1 |
    panel_q2
) +
  patchwork::plot_layout(
    widths = c(1, 1, 1)
  ) +
  patchwork::plot_annotation(
    title =
      "Analytical benchmark versus fixed-probability Monte Carlo",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12.5,
        hjust = 0
      )
    )
  )


# ------------------------------------------------------------
# 7. Save directly over manuscript Figure 5
# ------------------------------------------------------------

save_figure(
  figure_05,
  output_stem,
  width = 12.5,
  height = 5.8,
  dpi = 400
)

message("Rebuilt manuscript Figure 5:")
message(paste0(output_stem, ".pdf"))
message(paste0(output_stem, ".png"))

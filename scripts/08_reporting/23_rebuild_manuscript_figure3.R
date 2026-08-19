# ============================================================
# 23_rebuild_manuscript_figure3.R
#
# Purpose:
# Rebuild manuscript Figure 3 from the saved Monte Carlo monthly
# diversity summary without rerunning the Monte Carlo simulation.
#
# This script changes presentation only:
#   * retains the original Monte Carlo medians and 95% intervals;
#   * retains calendar gaps for months without represented communities;
#   * gives each q panel explicit solid x- and y-axis lines;
#   * removes explanatory subtitle/caption text from inside the figure;
#   * writes directly over the manuscript Figure 3 PNG and PDF.
#
# Run after script 17 (or after the complete pipeline).
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

monthly_summary_file <- file.path(
  paths$results_root,
  "05_monte_carlo",
  "tables",
  "monte_carlo_monthly_diversity_summary.csv"
)

if (!file.exists(monthly_summary_file)) {
  stop(
    "Monthly Monte Carlo summary not found: ",
    monthly_summary_file,
    "\nRun script 17 before rebuilding Figure 3."
  )
}

manuscript_figure_dir <- silwood_final_dir(paths, "main_figures")

output_stem <- file.path(
  manuscript_figure_dir,
  "Figure_03_output1_monthly_Hill_diversity"
)


# ------------------------------------------------------------
# 3. Read saved Monte Carlo summary
# ------------------------------------------------------------

monthly_summary <- readr::read_csv(
  monthly_summary_file,
  show_col_types = FALSE
)

required_columns <- c(
  "scenario",
  "scenario_label",
  "year_month",
  "q",
  "q_label",
  "median",
  "lower_95",
  "upper_95"
)

assert_columns(
  monthly_summary,
  required_columns,
  "Monte Carlo monthly diversity summary"
)

scenario_order <- c(
  "Primary_M3",
  "Primary_M2",
  "Audio_only"
)

scenario_labels <- c(
  Primary_M3 = "Primary M3: species-aware",
  Primary_M2 = "Primary M2: pooled",
  Audio_only = "Audio-only sensitivity"
)

q_order <- c(
  "q = 0: Richness",
  "q = 1: Exponential Shannon",
  "q = 2: Inverse Simpson"
)


# ------------------------------------------------------------
# 4. Prepare monthly plotting data
# ------------------------------------------------------------

monthly_plot_data <- monthly_summary %>%
  dplyr::filter(
    .data$scenario %in% scenario_order,
    .data$q_label %in% q_order
  ) %>%
  dplyr::mutate(
    month_date = as.Date(paste0(.data$year_month, "-01")),
    scenario = factor(
      .data$scenario,
      levels = scenario_order
    ),
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(scenario_labels[scenario_order])
    ),
    q_label = factor(
      .data$q_label,
      levels = q_order
    )
  )

# Insert calendar months with no represented detections as explicit NA rows.
# This preserves the original line gaps rather than connecting across periods
# with no represented acoustic community.
full_month_sequence <- seq(
  min(monthly_plot_data$month_date, na.rm = TRUE),
  max(monthly_plot_data$month_date, na.rm = TRUE),
  by = "month"
)

monthly_plot_data <- monthly_plot_data %>%
  dplyr::group_by(
    .data$scenario,
    .data$scenario_label,
    .data$q,
    .data$q_label
  ) %>%
  tidyr::complete(month_date = full_month_sequence) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    scenario = factor(
      .data$scenario,
      levels = scenario_order
    ),
    scenario_label = factor(
      .data$scenario_label,
      levels = unname(scenario_labels[scenario_order])
    ),
    q_label = factor(
      .data$q_label,
      levels = q_order
    )
  )


# ------------------------------------------------------------
# 5. Theme with explicit axes for every panel
# ------------------------------------------------------------

figure3_panel_theme <- function(base_size = 10) {
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
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        size = base_size - 1
      ),
      legend.position = "bottom"
    )
}


# ------------------------------------------------------------
# 6. Build one explicit-axis panel for each Hill order
# ------------------------------------------------------------

make_monthly_hill_panel <- function(q_value, panel_title) {
  
  panel_data <- monthly_plot_data %>%
    dplyr::filter(.data$q == q_value)
  
  ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = .data$month_date,
      y = .data$median,
      colour = .data$scenario_label,
      fill = .data$scenario_label,
      group = .data$scenario_label
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower_95,
        ymax = .data$upper_95
      ),
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(
      linewidth = 0.85
    ) +
    ggplot2::geom_point(
      size = 1.5,
      alpha = 0.75
    ) +
    ggplot2::scale_colour_manual(
      values = palette$scenario,
      breaks = unname(scenario_labels[scenario_order]),
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = palette$scenario,
      breaks = unname(scenario_labels[scenario_order]),
      drop = FALSE
    ) +
    ggplot2::scale_x_date(
      date_breaks = "3 months",
      labels = function(x) {
        english_month_year_labels(x, line_break = TRUE)
      },
      limits = range(full_month_sequence),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      title = panel_title,
      x = "Month",
      y = "Effective species number",
      colour = "Probability scenario",
      fill = "Probability scenario"
    ) +
    figure3_panel_theme(10)
}


panel_q0 <- make_monthly_hill_panel(
  q_value = 0,
  panel_title = "q = 0: Richness"
)

panel_q1 <- make_monthly_hill_panel(
  q_value = 1,
  panel_title = "q = 1: Exponential Shannon"
)

panel_q2 <- make_monthly_hill_panel(
  q_value = 2,
  panel_title = "q = 2: Inverse Simpson"
)


# ------------------------------------------------------------
# 7. Combine the three panels
# ------------------------------------------------------------

figure_03 <- (
  panel_q0 /
    panel_q1 /
    panel_q2
) +
  patchwork::plot_layout(
    guides = "collect",
    heights = c(1, 1, 1)
  ) +
  patchwork::plot_annotation(
    title =
      "Monthly Hill diversity after classification-uncertainty propagation",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        hjust = 0
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )


# ------------------------------------------------------------
# 8. Save directly over manuscript Figure 3
# ------------------------------------------------------------

save_figure(
  figure_03,
  output_stem,
  width = 12,
  height = 10.5,
  dpi = 400
)

message("Rebuilt manuscript Figure 3:")
message(paste0(output_stem, ".pdf"))
message(paste0(output_stem, ".png"))
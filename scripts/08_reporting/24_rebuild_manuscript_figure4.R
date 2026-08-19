# ============================================================
# 24_rebuild_manuscript_figure4.R
#
# Purpose:
# Rebuild manuscript Figure 4 from the saved method-comparison
# partition results without rerunning any diversity calculations.
#
# Presentation changes only:
#   * order partition components as Alpha -> Beta -> Gamma;
#   * retain the original percentage differences from Primary M3;
#   * retain the original workflow and q-order structure;
#   * add explicit solid x- and y-axis lines;
#   * remove explanatory subtitle/caption text from inside the figure;
#   * overwrite the manuscript Figure 4 PNG and PDF.
#
# Run after script 18 (or after the complete pipeline).
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

source(here::here("scripts", "00_setup", "00_pipeline_functions.R"))
load_silwood_packages()

config <- silwood_config()
paths <- silwood_paths(config)


# ------------------------------------------------------------
# 2. Input and output paths
# ------------------------------------------------------------

partition_difference_file <- file.path(
  paths$results_root,
  "06_method_comparison",
  "tables",
  "method_difference_from_m3_partition.csv"
)

if (!file.exists(partition_difference_file)) {
  stop(
    "Method-comparison partition table not found: ",
    partition_difference_file,
    "\nRun script 18 before rebuilding Figure 4."
  )
}

manuscript_figure_dir <- silwood_final_dir(paths, "main_figures")

output_stem <- file.path(
  manuscript_figure_dir,
  "Figure_04_workflow_difference_from_M3"
)


# ------------------------------------------------------------
# 3. Read saved method-comparison results
# ------------------------------------------------------------

partition_difference <- readr::read_csv(
  partition_difference_file,
  show_col_types = FALSE
)

required_columns <- c(
  "method",
  "method_label",
  "method_class",
  "weighting",
  "q",
  "q_label",
  "component",
  "component_label",
  "percentage_difference"
)

assert_columns(
  partition_difference,
  required_columns,
  "method-comparison partition table"
)


# ------------------------------------------------------------
# 4. Define manuscript ordering
# ------------------------------------------------------------

# Top-to-bottom order in the final figure.
workflow_order_top_to_bottom <- c(
  "Primary M2: pooled",
  "Audio-only sensitivity",
  "Naive: all candidates",
  "Hard top-1",
  "Threshold >= 0.50",
  "Threshold >= 0.70",
  "Threshold >= 0.90"
)

# For a discrete y-axis, the first factor level is plotted at the
# bottom, so reverse the desired top-to-bottom order.
workflow_factor_levels <- rev(workflow_order_top_to_bottom)

# Professor-requested ecological presentation order.
component_order <- c(
  "Alpha",
  "Beta",
  "Gamma"
)

q_order <- c(
  "q = 0: Richness",
  "q = 1: Exponential Shannon",
  "q = 2: Inverse Simpson"
)


# ------------------------------------------------------------
# 5. Prepare heatmap data
# ------------------------------------------------------------

heatmap_data <- partition_difference %>%
  dplyr::filter(
    .data$weighting == "equal_month",
    .data$method != "Primary_M3"
  ) %>%
  dplyr::mutate(
    method_label = factor(
      .data$method_label,
      levels = workflow_factor_levels
    ),
    component_label = factor(
      .data$component_label,
      levels = component_order
    ),
    q_label = factor(
      .data$q_label,
      levels = q_order
    ),
    display_label = sprintf(
      "%+.1f%%",
      .data$percentage_difference
    ),
    label_colour = dplyr::if_else(
      abs(.data$percentage_difference) >= 35,
      "white",
      "black"
    )
  )

heatmap_limit <- max(
  abs(heatmap_data$percentage_difference),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 6. Rebuild manuscript Figure 4
# ------------------------------------------------------------

figure_04 <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = .data$component_label,
    y = .data$method_label,
    fill = .data$percentage_difference
  )
) +
  ggplot2::geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = .data$display_label,
      colour = .data$label_colour
    ),
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(
    ~q_label,
    nrow = 1
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2C5AA0",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -heatmap_limit,
      heatmap_limit
    ),
    oob = scales::squish,
    labels = scales::label_number(
      suffix = "%",
      accuracy = 1
    )
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::labs(
    title = "Alternative identification workflows relative to Primary M3",
    x = "Partition component",
    y = "Identification workflow",
    fill = "Difference from\nPrimary M3 (%)"
  ) +
  theme_silwood(10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
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
      linewidth = 0.40
    ),
    axis.ticks.length = grid::unit(2.5, "pt"),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 8)
    ),
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    legend.position = "right",
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  )


# ------------------------------------------------------------
# 7. Save directly over manuscript Figure 4
# ------------------------------------------------------------

save_figure(
  figure_04,
  output_stem,
  width = 12.5,
  height = 6.8,
  dpi = 400
)

message("Rebuilt manuscript Figure 4:")
message(paste0(output_stem, ".pdf"))
message(paste0(output_stem, ".png"))

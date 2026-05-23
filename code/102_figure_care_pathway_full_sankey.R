#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    ofiles <- vapply(sys.frames(), function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile, character(1))
    ofiles <- stats::na.omit(ofiles)
    if (length(ofiles) > 0) return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = TRUE))
    stop("Could not determine script path. Run with Rscript.")
  }
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)

input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggalluvial)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(readr)
  library(ragg)
  library(scales)
  library(stringr)
  library(svglite)
  library(tidyr)
})

pathways <- readr::read_csv(file.path(input_dir, "pooled_care_pathway_summary.csv"), show_col_types = FALSE) |>
  mutate(n = as.numeric(.data$n), pct = as.numeric(.data$pct)) |>
  arrange(desc(.data$n))

n_total <- sum(pathways$n, na.rm = TRUE)

nice_location <- function(x) {
  dplyr::recode(
    x,
    ed = "ED",
    icu = "ICU",
    procedural = "Procedural",
    ward = "Ward",
    stepdown = "Stepdown",
    other = "Other",
    radiology = "Radiology",
    rehab = "Rehab",
    dialysis = "Dialysis",
    admitted = "Admitted",
    .default = str_to_title(x)
  )
}

location_palette <- c(
  Admitted = "#17212B",
  ED = "#1B6CA8",
  ICU = "#7A3E9D",
  Procedural = "#C85A2E",
  Ward = "#3D8B5A",
  Stepdown = "#C79A2B",
  Other = "#6B7280",
  Radiology = "#8E6C4A",
  Rehab = "#607D8B",
  Dialysis = "#536DFE",
  `No additional stop` = "#D5DEE8",
  `Additional stops` = "#A7B0BA"
)

route_tokens <- str_split(pathways$care_pathway, " -> ", simplify = FALSE)
max_steps_shown <- 5L

pathway_wide <- pathways |>
  mutate(
    path_id = row_number(),
    route_pretty = vapply(route_tokens, function(x) paste(nice_location(x), collapse = " -> "), character(1)),
    first_location = vapply(route_tokens, function(x) nice_location(x[[1]]), character(1)),
    flow_group = if_else(.data$first_location %in% names(location_palette), .data$first_location, "Other"),
    cohort = sprintf("Admitted\n%s", comma(n_total)),
    step1 = vapply(route_tokens, function(x) if (length(x) >= 1) nice_location(x[[1]]) else NA_character_, character(1)),
    step2 = vapply(route_tokens, function(x) if (length(x) >= 2) nice_location(x[[2]]) else "No additional stop", character(1)),
    step3 = vapply(route_tokens, function(x) if (length(x) >= 3) nice_location(x[[3]]) else "No additional stop", character(1)),
    step4 = vapply(route_tokens, function(x) if (length(x) >= 4) nice_location(x[[4]]) else "No additional stop", character(1)),
    step5 = vapply(route_tokens, function(x) if (length(x) >= 5) nice_location(x[[5]]) else "No additional stop", character(1)),
    later = vapply(route_tokens, function(x) if (length(x) > max_steps_shown) "Additional stops" else "No additional stop", character(1)),
    final_icu = "ICU"
  )

axis_cols <- c("cohort", "step1", "step2", "step3", "step4", "step5", "later", "final_icu")
axis_labels <- c("Cohort", "1st", "2nd", "3rd", "4th", "5th", "Later", "ICU")
cohort_label <- unique(pathway_wide$cohort)

alluvial_data <- pathway_wide |>
  ggalluvial::to_lodes_form(
    axes = axis_cols,
    id = "path_id",
    key = "axis",
    value = "stratum"
  ) |>
  mutate(
    axis = factor(.data$axis, levels = axis_cols, labels = axis_labels),
    stratum = as.character(.data$stratum),
    stratum_fill = if_else(.data$stratum %in% names(location_palette), .data$stratum, "Other"),
    flow_group = if_else(.data$flow_group %in% names(location_palette), .data$flow_group, "Other")
  )

stratum_totals <- alluvial_data |>
  group_by(.data$axis, .data$stratum, .data$stratum_fill) |>
  summarise(n = sum(.data$n), .groups = "drop") |>
  group_by(.data$axis) |>
  mutate(
    pct_total = 100 * .data$n / n_total,
    label = case_when(
      .data$stratum == "No additional stop" ~ "",
      .data$pct_total >= 2 ~ sprintf("%s\n%s", .data$stratum, comma(.data$n)),
      TRUE ~ ""
    )
  ) |>
  ungroup()

top_routes <- pathway_wide |>
  slice_head(n = 15) |>
  mutate(
    route_label = sprintf("%02d. %s", row_number(), .data$route_pretty),
    display = sprintf("%s (%.1f%%)", comma(.data$n), .data$pct)
  )

legend_x_max <- max(top_routes$n, na.rm = TRUE) * 1.42

legend_plot <- ggplot(top_routes, aes(y = fct_rev(factor(.data$route_label, levels = rev(.data$route_label))), x = .data$n, fill = .data$flow_group)) +
  geom_col(width = 0.62, alpha = 0.95) +
  geom_text(aes(label = .data$display), hjust = -0.04, size = 3.1, color = "#263238") +
  scale_fill_manual(values = location_palette, guide = "none") +
  scale_x_continuous(labels = label_comma(), limits = c(0, legend_x_max), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Top Pathways",
    subtitle = "All 68 pathways are included in the Sankey; top 15 are listed here",
    x = "Admissions",
    y = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5, color = "#17212B"),
    plot.subtitle = element_text(size = 8.8, color = "#455A64"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7.5, color = "#263238"),
    axis.text.x = element_text(size = 8, color = "#607D8B"),
    axis.title.x = element_text(size = 8.5, color = "#455A64"),
    plot.margin = margin(8, 28, 8, 8)
  )

sankey_plot <- ggplot(
  alluvial_data,
  aes(
    x = .data$axis,
    stratum = .data$stratum,
    alluvium = .data$path_id,
    y = .data$n
  )
) +
  ggalluvial::geom_alluvium(
    aes(fill = .data$flow_group),
    width = 0.12,
    alpha = 0.56,
    knot.pos = 0.38,
    curve_type = "quintic",
    color = "white",
    linewidth = 0.08
  ) +
  ggalluvial::geom_stratum(
    aes(fill = .data$stratum_fill),
    width = 0.18,
    color = "white",
    linewidth = 0.55,
    alpha = 0.98
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(ifelse(count >= 250 & !stratum %in% c("No additional stop", cohort_label), stratum, ""))),
    color = "white",
    fontface = "bold",
    size = 2.7,
    lineheight = 0.92
  ) +
  annotate(
    "label",
    x = 1,
    y = n_total * 0.48,
    label = sprintf("Admitted\n%s", comma(n_total)),
    color = "white",
    fill = "#17212B",
    fontface = "bold",
    size = 2.9,
    lineheight = 0.92,
    linewidth = 0,
    label.padding = unit(0.13, "lines")
  ) +
  scale_fill_manual(values = location_palette, guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.58, 0.24))) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.01, 0.04))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Full Care-Pathway Sankey to First ICU Entry",
    subtitle = sprintf("All 68 observed care pathways represented; flow color indicates first recorded location (n=%s)", comma(n_total)),
    x = NULL,
    y = "Admissions",
    caption = "Columns show observed sequence positions. Pathways longer than five locations retain their early sequence and are marked with an 'Additional stops' later column before ICU."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15.5, color = "#17212B"),
    plot.subtitle = element_text(size = 10, color = "#455A64", margin = margin(b = 10)),
    plot.caption = element_text(size = 8.3, color = "#607D8B", hjust = 0, margin = margin(t = 8)),
    axis.text.x = element_text(face = "bold", color = "#263238", size = 10.5),
    axis.text.y = element_text(color = "#607D8B", size = 8.5),
    axis.title.y = element_text(color = "#455A64", size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 14)
  )

full_figure <- sankey_plot + legend_plot + plot_layout(widths = c(2.15, 1))

save_figure <- function(plot, basename, width = 15.5, height = 8.4) {
  pdf(file.path(output_dir, paste0(basename, ".pdf")), width = width, height = height, onefile = FALSE)
  print(plot)
  dev.off()
  svglite::svglite(file.path(output_dir, paste0(basename, ".svg")), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(output_dir, paste0(basename, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(output_dir, paste0(basename, ".tiff")), width = width, height = height, units = "in", res = 600, background = "white", compression = "lzw")
  print(plot)
  dev.off()
}

save_figure(full_figure, "figure_care_pathway_full_sankey")
readr::write_csv(pathway_wide, file.path(output_dir, "figure_care_pathway_full_sankey_source.csv"))

message("Wrote full care-pathway Sankey to: ", output_dir)

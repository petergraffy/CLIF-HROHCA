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
  library(ggforce)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(ragg)
  library(scales)
  library(stringr)
  library(svglite)
  library(tidyr)
})

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
    .default = str_to_title(x)
  )
}

pathways <- readr::read_csv(file.path(input_dir, "pooled_care_pathway_summary.csv"), show_col_types = FALSE) |>
  mutate(n = as.numeric(.data$n), pct = as.numeric(.data$pct)) |>
  arrange(desc(.data$n))

n_total <- sum(pathways$n, na.rm = TRUE)
route_tokens <- str_split(pathways$care_pathway, " -> ", simplify = FALSE)

top5 <- pathways |>
  slice_head(n = 5) |>
  mutate(
    path_id = row_number(),
    route_tokens = route_tokens[seq_len(5)],
    route_pretty = vapply(.data$route_tokens, function(x) paste(nice_location(x), collapse = " -> "), character(1)),
    flow_group = factor(
      sprintf("%02d. %s", .data$path_id, .data$route_pretty),
      levels = sprintf("%02d. %s", .data$path_id, .data$route_pretty)
    ),
    first_location = vapply(.data$route_tokens, function(x) nice_location(x[[1]]), character(1)),
    intervening_step = vapply(.data$route_tokens, function(x) {
      if (length(x) == 1 && x[[1]] == "icu") return("Already in ICU")
      if (length(x) == 2 && x[[2]] == "icu") return("Direct to ICU")
      nice_location(x[[2]])
    }, character(1)),
    icu_entry = "ICU"
  )

n_top5 <- sum(top5$n, na.rm = TRUE)
top5_pct <- 100 * n_top5 / n_total

flow_palette <- c(
  "01. ED -> ICU" = "#1B6CA8",
  "02. ICU" = "#7A3E9D",
  "03. ED -> Procedural -> ICU" = "#2A83B9",
  "04. Procedural -> ICU" = "#C85A2E",
  "05. ED -> Ward -> ICU" = "#3D8B5A"
)

ring_def <- tibble::tribble(
  ~stage, ~r0, ~r, ~alpha,
  "Cohort", 0.72, 1.05, 0.96,
  "First location", 1.16, 1.49, 0.86,
  "Intervening step", 1.60, 1.93, 0.74,
  "First ICU entry", 2.04, 2.37, 0.88
)

angle_gap <- 0.018
top5_angles <- top5 |>
  mutate(
    frac = .data$n / n_top5,
    angle_start_raw = cumsum(dplyr::lag(.data$frac, default = 0)) * 2 * pi + pi / 2,
    angle_end_raw = cumsum(.data$frac) * 2 * pi + pi / 2,
    start = .data$angle_start_raw + angle_gap,
    end = .data$angle_end_raw - angle_gap,
    mid = (.data$start + .data$end) / 2
  )

ring_segments <- top5_angles |>
  select("path_id", "route_pretty", "flow_group", "n", "pct", "start", "end", "mid", "first_location", "intervening_step", "icu_entry") |>
  tidyr::crossing(ring_def) |>
  mutate(
    label = case_when(
      .data$stage == "Cohort" ~ sprintf("Pathway %d", .data$path_id),
      .data$stage == "First location" ~ .data$first_location,
      .data$stage == "Intervening step" ~ .data$intervening_step,
      TRUE ~ .data$icu_entry
    ),
    label_r = (.data$r0 + .data$r) / 2,
    label_angle = (.data$mid * 180 / pi) %% 360,
    label_rotation = if_else(.data$label_angle > 90 & .data$label_angle < 270, .data$label_angle + 180, .data$label_angle),
    label_color = if_else(.data$stage == "Intervening step" & .data$label %in% c("Direct to ICU", "Already in ICU"), "#263238", "white"),
    show_label = case_when(
      .data$stage == "Cohort" ~ FALSE,
      .data$pct >= 4 ~ TRUE,
      .data$stage == "First location" & .data$pct >= 2 ~ TRUE,
      TRUE ~ FALSE
    )
  )

route_table <- top5 |>
  mutate(
    route_label_short = sprintf("%d. %s", .data$path_id, .data$route_pretty),
    display = sprintf("%s (%.1f%%)", comma(.data$n), .data$pct),
    route_label_short = fct_rev(factor(.data$route_label_short, levels = rev(.data$route_label_short)))
  )

circular_plot <- ggplot() +
  ggforce::geom_arc_bar(
    data = ring_segments,
    aes(
      x0 = 0,
      y0 = 0,
      r0 = .data$r0,
      r = .data$r,
      start = .data$start,
      end = .data$end,
      fill = .data$flow_group,
      alpha = .data$alpha
    ),
    color = "white",
    linewidth = 0.65
  ) +
  annotate(
    "label",
    x = 0,
    y = 0,
    label = sprintf("Top 5\n%s\n%.1f%%", comma(n_top5), top5_pct),
    fill = "#17212B",
    color = "white",
    fontface = "bold",
    size = 3.55,
    lineheight = 0.92,
    linewidth = 0,
    label.padding = unit(0.2, "lines")
  ) +
  annotate(
    "text",
    x = 0,
    y = 2.74,
    label = "Read center outward:\nCohort -> First location -> Intervening step -> First ICU entry",
    color = "#263238",
    fontface = "bold",
    size = 3.3,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = flow_palette, guide = "none") +
  scale_alpha_identity() +
  scale_color_identity() +
  coord_fixed(clip = "off") +
  labs(
    title = "Circular Alluvial of Dominant Care Pathways",
    subtitle = sprintf("Top five pathways account for %s of %s classified admissions (%.1f%%)", comma(n_top5), comma(n_total), top5_pct),
    caption = "Each colored sector is one pathway; concentric rings show cohort, first location, intervening step, and first ICU entry."
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, color = "#17212B"),
    plot.subtitle = element_text(size = 11.2, color = "#455A64", margin = margin(b = 10)),
    plot.caption = element_text(size = 9.5, color = "#607D8B", hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(8, 18, 8, 18)
  )

route_table_plot <- ggplot(route_table, aes(x = .data$n, y = .data$route_label_short, fill = .data$flow_group)) +
  geom_col(width = 0.6, alpha = 0.96) +
  geom_text(aes(label = .data$display), hjust = -0.05, size = 3.8, color = "#263238") +
  scale_fill_manual(values = flow_palette, guide = "none") +
  scale_x_continuous(labels = label_comma(), limits = c(0, max(route_table$n) * 1.36), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Pathway Counts",
    x = "Admissions",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#17212B"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10.5, color = "#263238"),
    axis.text.x = element_text(size = 9.5, color = "#607D8B"),
    axis.title.x = element_text(size = 10, color = "#455A64"),
    plot.margin = margin(8, 28, 8, 8)
  )

figure <- circular_plot + route_table_plot + plot_layout(widths = c(1.18, 1))

save_figure <- function(plot, basename, width = 13.8, height = 7.6) {
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

save_figure(figure, "figure_care_pathway_top5_circular_alluvial")
readr::write_csv(top5, file.path(output_dir, "figure_care_pathway_top5_circular_alluvial_source.csv"))

message("Wrote top-5 circular care-pathway alluvial to: ", output_dir)

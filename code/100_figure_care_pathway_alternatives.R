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
  library(ggtext)
  library(patchwork)
  library(readr)
  library(ragg)
  library(scales)
  library(stringr)
  library(svglite)
  library(tidyr)
})

pathway_path <- file.path(input_dir, "pooled_care_pathway_summary.csv")
timing_path <- file.path(input_dir, "pooled_icu_timing_bins.csv")
timing_summary_path <- file.path(input_dir, "pooled_icu_timing_summary.csv")

pathways_raw <- readr::read_csv(pathway_path, show_col_types = FALSE) |>
  mutate(n = as.numeric(.data$n), pct = as.numeric(.data$pct)) |>
  arrange(desc(.data$n))

timing <- readr::read_csv(timing_path, show_col_types = FALSE) |>
  mutate(n = as.numeric(.data$n), pct = as.numeric(.data$pct))

timing_summary <- readr::read_csv(timing_summary_path, show_col_types = FALSE)
n_total <- sum(pathways_raw$n, na.rm = TRUE)
site_count <- if (file.exists(file.path(input_dir, "care_pathway_export_availability.csv"))) {
  nrow(readr::read_csv(file.path(input_dir, "care_pathway_export_availability.csv"), show_col_types = FALSE))
} else {
  NA_integer_
}

loc_palette <- c(
  ED = "#1B6CA8",
  ICU = "#7A3E9D",
  Procedural = "#C85A2E",
  Ward = "#3D8B5A",
  Stepdown = "#C79A2B",
  Other = "#6B7280",
  Direct = "#6D8FA6",
  Multiple = "#A7B0BA",
  Cohort = "#17212B"
)

nice_location <- function(x) {
  dplyr::recode(
    x,
    ed = "ED",
    icu = "ICU",
    procedural = "Procedural",
    ward = "Ward",
    stepdown = "Stepdown",
    other = "Other",
    radiology = "Other",
    rehab = "Other",
    .default = "Other"
  )
}

split_path <- function(x) str_split(x, " -> ", simplify = FALSE)

classify_intermediate <- function(pathway) {
  tokens <- split_path(pathway)
  vapply(tokens, function(x) {
    x <- x[nzchar(x)]
    if (length(x) == 0L) return("Multiple")
    if (identical(x[[1]], "icu")) return("Direct")
    middle <- x[-c(1, length(x))]
    middle <- middle[middle != "icu"]
    if (length(middle) == 0L) return("Direct")
    if (length(middle) == 1L && middle %in% c("procedural", "ward", "stepdown")) return(nice_location(middle))
    "Multiple"
  }, character(1))
}

format_route <- function(x) {
  tokens <- split_path(x)[[1]]
  paste(nice_location(tokens), collapse = " -> ")
}

pathways <- pathways_raw |>
  mutate(
    path_id = row_number(),
    route = vapply(.data$care_pathway, format_route, character(1)),
    first_raw = str_split_i(.data$care_pathway, " -> ", 1),
    first_location = nice_location(.data$first_raw),
    intermediate = classify_intermediate(.data$care_pathway),
    timing_text = sprintf("%s (%.1f%%)", comma(.data$n), .data$pct)
  )

top_n <- 10L
top_pathways <- pathways |>
  slice_head(n = top_n) |>
  mutate(
    route_rank = sprintf("%02d. %s", row_number(), .data$route),
    route_rank = fct_rev(factor(.data$route_rank, levels = rev(.data$route_rank)))
  )

tile_steps <- top_pathways |>
  select("path_id", "route_rank", "care_pathway", "n", "pct") |>
  mutate(tokens = split_path(.data$care_pathway)) |>
  unnest_longer(.data$tokens, values_to = "location_raw", indices_to = "step") |>
  mutate(
    location = nice_location(.data$location_raw),
    step = as.integer(.data$step)
  )

pathway_tile_plot <- ggplot(tile_steps, aes(x = .data$step, y = .data$route_rank, fill = .data$location)) +
  geom_tile(width = 0.86, height = 0.72, color = "white", linewidth = 0.6) +
  geom_text(aes(label = .data$location), color = "white", fontface = "bold", size = 3.0) +
  geom_text(
    data = top_pathways,
    aes(x = 5.05, y = .data$route_rank, label = sprintf("%s  |  %.1f%%", comma(.data$n), .data$pct)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.2,
    color = "#263238"
  ) +
  scale_fill_manual(values = loc_palette, guide = "none") +
  scale_x_continuous(breaks = 1:5, labels = paste("Step", 1:5), limits = c(0.5, 6.4), expand = c(0, 0)) +
  labs(
    title = "Option A. Ranked Pathway Tile Matrix",
    subtitle = "Top pre-ICU routes preserve sequence without tangled flows",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#17212B", size = 13),
    plot.subtitle = element_text(color = "#455A64", size = 9.5),
    axis.text.y = element_text(color = "#263238", size = 8.7),
    axis.text.x = element_text(color = "#455A64", face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(8, 18, 8, 8)
  )

mekko <- pathways |>
  mutate(
    first_location = fct_relevel(.data$first_location, "ED", "ICU", "Procedural", "Ward", "Stepdown", "Other"),
    intermediate = factor(.data$intermediate, levels = c("Direct", "Procedural", "Ward", "Stepdown", "Multiple"))
  ) |>
  group_by(.data$first_location, .data$intermediate) |>
  summarise(n = sum(.data$n), .groups = "drop") |>
  group_by(.data$first_location) |>
  mutate(first_total = sum(.data$n), y0 = cumsum(lag(.data$n, default = 0)), y1 = cumsum(.data$n)) |>
  ungroup() |>
  arrange(.data$first_location) |>
  mutate(first_location = as.character(.data$first_location))

first_widths <- mekko |>
  distinct(.data$first_location, .data$first_total) |>
  arrange(factor(.data$first_location, levels = c("ED", "ICU", "Procedural", "Ward", "Stepdown", "Other"))) |>
  mutate(x0 = cumsum(lag(.data$first_total, default = 0)), x1 = cumsum(.data$first_total))

mekko <- mekko |>
  left_join(first_widths, by = c("first_location", "first_total")) |>
  mutate(
    xmid = (.data$x0 + .data$x1) / 2,
    ymid = (.data$y0 + .data$y1) / 2,
    first_pct = 100 * .data$first_total / n_total,
    cell_pct = 100 * .data$n / n_total,
    label = if_else(.data$cell_pct >= 2, sprintf("%s\n%.1f%%", .data$intermediate, .data$cell_pct), "")
  )

marimekko_plot <- ggplot(mekko) +
  geom_rect(
    aes(xmin = .data$x0, xmax = .data$x1, ymin = .data$y0 / .data$first_total, ymax = .data$y1 / .data$first_total, fill = .data$intermediate),
    color = "white",
    linewidth = 0.6,
    alpha = 0.96
  ) +
  geom_text(
    aes(x = .data$xmid, y = .data$ymid / .data$first_total, label = .data$label),
    color = "white",
    fontface = "bold",
    size = 3.0,
    lineheight = 0.92
  ) +
  geom_text(
    data = first_widths |> mutate(xmid = (.data$x0 + .data$x1) / 2, pct = 100 * .data$first_total / n_total),
    aes(x = .data$xmid, y = -0.06, label = sprintf("%s\n%.1f%%", .data$first_location, .data$pct)),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold",
    color = "#263238",
    lineheight = 0.9
  ) +
  scale_fill_manual(values = loc_palette, guide = guide_legend(title = "Pre-ICU pathway")) +
  scale_x_continuous(labels = NULL, expand = c(0, 0)) +
  scale_y_continuous(labels = label_percent(), limits = c(-0.13, 1), expand = c(0, 0)) +
  labs(
    title = "Option B. Marimekko First Location by Pre-ICU Pathway",
    subtitle = "Column width shows first location; column height shows pathway mix within that first location",
    x = NULL,
    y = "Within first location"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#17212B", size = 13),
    plot.subtitle = element_text(color = "#455A64", size = 9.5),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "#607D8B"),
    axis.title.y = element_text(color = "#455A64", size = 9),
    plot.margin = margin(8, 8, 8, 8)
  )

sun_levels <- pathways |>
  mutate(
    first_location = if_else(.data$first_location %in% c("ED", "ICU", "Procedural", "Ward", "Stepdown"), .data$first_location, "Other"),
    intermediate = if_else(.data$intermediate %in% c("Direct", "Procedural", "Ward", "Stepdown", "Multiple"), .data$intermediate, "Multiple")
  ) |>
  group_by(.data$first_location, .data$intermediate) |>
  summarise(n = sum(.data$n), .groups = "drop")

sun_first <- sun_levels |>
  group_by(.data$first_location) |>
  summarise(n = sum(.data$n), .groups = "drop") |>
  arrange(desc(.data$n)) |>
  mutate(
    start = 2 * pi * cumsum(lag(.data$n, default = 0)) / sum(.data$n),
    end = 2 * pi * cumsum(.data$n) / sum(.data$n),
    r0 = 0.25,
    r = 0.58,
    label_angle = (.data$start + .data$end) / 2,
    label = if_else(.data$n / n_total >= 0.03, sprintf("%s\n%.1f%%", .data$first_location, 100 * .data$n / n_total), "")
  )

sun_second <- sun_levels |>
  left_join(sun_first |> select("first_location", first_start = "start", first_end = "end", first_n = "n"), by = "first_location") |>
  group_by(.data$first_location) |>
  arrange(desc(.data$n), .by_group = TRUE) |>
  mutate(
    cum0 = cumsum(lag(.data$n, default = 0)),
    cum1 = cumsum(.data$n),
    start = .data$first_start + (.data$first_end - .data$first_start) * .data$cum0 / .data$first_n,
    end = .data$first_start + (.data$first_end - .data$first_start) * .data$cum1 / .data$first_n,
    r0 = 0.61,
    r = 1.0,
    label_angle = (.data$start + .data$end) / 2,
    pct = 100 * .data$n / n_total,
    label = if_else(.data$pct >= 3, sprintf("%s\n%.1f%%", .data$intermediate, .data$pct), "")
  ) |>
  ungroup()

sunburst_plot <- ggplot() +
  ggforce::geom_arc_bar(
    data = sun_first,
    aes(x0 = 0, y0 = 0, r0 = .data$r0, r = .data$r, start = .data$start, end = .data$end, fill = .data$first_location),
    color = "white",
    linewidth = 0.55
  ) +
  ggforce::geom_arc_bar(
    data = sun_second,
    aes(x0 = 0, y0 = 0, r0 = .data$r0, r = .data$r, start = .data$start, end = .data$end, fill = .data$intermediate),
    color = "white",
    linewidth = 0.55,
    alpha = 0.95
  ) +
  geom_text(
    data = sun_first,
    aes(x = 0.43 * sin(.data$label_angle), y = 0.43 * cos(.data$label_angle), label = .data$label),
    color = "white",
    fontface = "bold",
    size = 3.0,
    lineheight = 0.9
  ) +
  geom_text(
    data = sun_second,
    aes(x = 0.82 * sin(.data$label_angle), y = 0.82 * cos(.data$label_angle), label = .data$label),
    color = "white",
    fontface = "bold",
    size = 2.8,
    lineheight = 0.9
  ) +
  annotate("text", x = 0, y = 0, label = sprintf("OHCA ICU\n%s", comma(n_total)), fontface = "bold", size = 4.0, color = "#17212B", lineheight = 0.9) +
  scale_fill_manual(values = loc_palette, guide = "none") +
  coord_fixed(clip = "off") +
  labs(
    title = "Option C. Sunburst Hierarchy",
    subtitle = "Inner ring: first location; outer ring: pre-ICU pathway category"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#17212B", size = 13, hjust = 0),
    plot.subtitle = element_text(color = "#455A64", size = 9.5, hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

timing_compact <- timing |>
  mutate(
    timing_group = case_when(
      .data$icu_entry_timing %in% c("<=0 hours", "0-6 hours") ~ "<=6 h",
      .data$icu_entry_timing == "6-12 hours" ~ "6-12 h",
      .data$icu_entry_timing == "12-24 hours" ~ "12-24 h",
      TRUE ~ ">24 h"
    ),
    timing_group = factor(.data$timing_group, levels = c("<=6 h", "6-12 h", "12-24 h", ">24 h"))
  ) |>
  group_by(.data$timing_group) |>
  summarise(n = sum(.data$n), .groups = "drop") |>
  mutate(
    pct = 100 * .data$n / sum(.data$n),
    ymax = cumsum(.data$n),
    ymin = lag(.data$ymax, default = 0),
    label_pos = (.data$ymin + .data$ymax) / 2,
    label = if_else(.data$pct >= 5, sprintf("%s\n%.1f%%", .data$timing_group, .data$pct), "")
  )

timing_palette <- c("<=6 h" = "#2F6F73", "6-12 h" = "#6D8FA6", "12-24 h" = "#B8A35C", ">24 h" = "#A75D5D")

timing_donut <- ggplot(timing_compact, aes(ymax = .data$ymax, ymin = .data$ymin, xmax = 1, xmin = 0.55, fill = .data$timing_group)) +
  geom_rect(color = "white", linewidth = 0.7) +
  coord_polar(theta = "y") +
  geom_text(
    aes(x = 0.79, y = .data$label_pos, label = .data$label),
    color = "white",
    fontface = "bold",
    size = 3.0,
    lineheight = 0.9
  ) +
  annotate("text", x = 0, y = 0, label = sprintf("Median\n%.1f h", timing_summary$median_hours_approx[1]), fontface = "bold", size = 4.0, color = "#17212B", lineheight = 0.9) +
  scale_fill_manual(values = timing_palette, guide = guide_legend(title = "Time to ICU")) +
  xlim(c(0, 1)) +
  labs(title = "ICU Timing Donut") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#17212B", size = 13, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5)
  )

comparison <- (pathway_tile_plot | marimekko_plot) / (sunburst_plot | timing_donut) +
  plot_annotation(
    title = "Care Pathway Figure Alternatives",
    subtitle = sprintf("Pooled OHCA ICU admissions across %d CLIF sites (n=%s)", site_count, comma(n_total)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, color = "#17212B"),
      plot.subtitle = element_text(size = 10.5, color = "#455A64")
    )
  )

save_figure <- function(plot, basename, width, height) {
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

save_figure(pathway_tile_plot, "figure_care_pathway_option_tile_matrix", width = 8.3, height = 5.0)
save_figure(marimekko_plot, "figure_care_pathway_option_marimekko", width = 7.5, height = 5.0)
save_figure(sunburst_plot, "figure_care_pathway_option_sunburst", width = 6.0, height = 5.6)
save_figure(timing_donut, "figure_icu_timing_donut", width = 4.8, height = 4.6)
save_figure(comparison, "figure_care_pathway_alternatives_comparison", width = 13.2, height = 10.0)

readr::write_csv(pathways, file.path(output_dir, "figure_care_pathway_alternatives_source.csv"))

message("Wrote care pathway alternatives to: ", output_dir)

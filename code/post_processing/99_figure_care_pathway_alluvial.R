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

repo_root <- normalizePath(file.path(dirname(get_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)
input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
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

pathway_path <- file.path(input_dir, "pooled_care_pathway_summary.csv")
timing_path <- file.path(input_dir, "pooled_icu_timing_bins.csv")
timing_summary_path <- file.path(input_dir, "pooled_icu_timing_summary.csv")

if (!file.exists(pathway_path) || !file.exists(timing_path)) {
  stop("Run code/post_processing/98_pool_care_pathways.R before building the alluvial figure.")
}

pathways <- readr::read_csv(pathway_path, show_col_types = FALSE) |>
  mutate(
    n = as.numeric(.data$n),
    pct = as.numeric(.data$pct)
  ) |>
  mutate(pct = 100 * .data$n / sum(.data$n)) |>
  arrange(desc(.data$n))

timing <- readr::read_csv(timing_path, show_col_types = FALSE) |>
  mutate(
    n = as.numeric(.data$n),
    pct = as.numeric(.data$pct),
    icu_entry_timing = factor(
      .data$icu_entry_timing,
      levels = c("<=0 hours", "0-6 hours", "6-12 hours", "12-24 hours", "24-48 hours", "48-72 hours", ">72 hours")
    )
  ) |>
  arrange(.data$icu_entry_timing)

timing_summary <- readr::read_csv(timing_summary_path, show_col_types = FALSE)
n_total <- sum(pathways$n, na.rm = TRUE)
site_count <- if (file.exists(file.path(input_dir, "care_pathway_export_availability.csv"))) {
  nrow(readr::read_csv(file.path(input_dir, "care_pathway_export_availability.csv"), show_col_types = FALSE))
} else {
  NA_integer_
}

location_palette <- c(
  admitted = "#17212B",
  ED = "#1B6CA8",
  ICU = "#7A3E9D",
  Procedural = "#C85A2E",
  Ward = "#3D8B5A",
  Stepdown = "#C79A2B",
  Other = "#6B7280",
  `Direct to ICU` = "#7A3E9D",
  `No pre-ICU location` = "#4C566A",
  `Multiple or other` = "#9CA3AF",
  `First ICU` = "#7A3E9D"
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
    radiology = "Radiology",
    rehab = "Rehab",
    admitted = "Admitted",
    .default = str_to_title(x)
  )
}

classify_intermediate <- function(pathway) {
  tokens <- str_split(pathway, " -> ", simplify = FALSE)
  vapply(tokens, function(x) {
    x <- x[nzchar(x)]
    if (length(x) == 0L) return("Multiple or other")
    if (identical(x[[1]], "icu")) return("No pre-ICU location")
    middle <- x[-c(1, length(x))]
    middle <- middle[middle != "icu"]
    if (length(middle) == 0L) return("Direct to ICU")
    if (length(middle) == 1L && middle %in% c("procedural", "ward", "stepdown")) return(nice_location(middle))
    "Multiple or other"
  }, character(1))
}

alluvial_wide <- pathways |>
  mutate(
    path_id = row_number(),
    first_raw = str_split_i(.data$care_pathway, " -> ", 1),
    first_location = nice_location(.data$first_raw),
    first_location = if_else(.data$first_location %in% c("ED", "ICU", "Procedural", "Ward", "Stepdown"), .data$first_location, "Other"),
    intermediate = classify_intermediate(.data$care_pathway),
    intermediate = factor(
      .data$intermediate,
      levels = c("Direct to ICU", "No pre-ICU location", "Procedural", "Ward", "Stepdown", "Multiple or other")
    ),
    cohort = sprintf("Cohort\n%s", comma(n_total)),
    first_icu = "ICU",
    first_fill = if_else(.data$first_location %in% names(location_palette), .data$first_location, "Other")
  ) |>
  group_by(.data$cohort, .data$first_location, .data$intermediate, .data$first_icu, .data$first_fill) |>
  summarise(n = sum(.data$n), .groups = "drop") |>
  arrange(desc(.data$n)) |>
  mutate(path_id = row_number())

alluvial_data <- alluvial_wide |>
  ggalluvial::to_lodes_form(
    axes = c("cohort", "first_location", "intermediate", "first_icu"),
    id = "path_id",
    key = "axis",
    value = "stratum"
  ) |>
  mutate(
    axis = factor(
      .data$axis,
      levels = c("cohort", "first_location", "intermediate", "first_icu"),
      labels = c("Cohort", "First location", "Pre-ICU pathway", "ICU")
    ),
    stratum = as.character(.data$stratum),
    stratum_fill = if_else(.data$stratum %in% names(location_palette), .data$stratum, "Other")
  )

stratum_totals <- alluvial_data |>
  group_by(.data$axis, .data$stratum) |>
  summarise(n = sum(.data$n), .groups = "drop") |>
  mutate(
    pct = 100 * .data$n / n_total,
    label = if_else(.data$pct >= 2, sprintf("%s\n%s (%.1f%%)", .data$stratum, comma(.data$n), .data$pct), "")
  )

alluvial_plot <- ggplot(
  alluvial_data,
  aes(
    x = .data$axis,
    stratum = .data$stratum,
    alluvium = .data$path_id,
    y = .data$n
  )
) +
  ggalluvial::geom_alluvium(aes(fill = .data$first_fill), width = 0.12, alpha = 0.52, knot.pos = 0.35, curve_type = "quintic", color = NA) +
  ggalluvial::geom_stratum(aes(fill = .data$stratum_fill), width = 0.16, color = "white", linewidth = 0.5, alpha = 0.98) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(ifelse(count >= 500, stratum, ""))),
    size = 2.9,
    color = "white",
    fontface = "bold",
    lineheight = 0.92
  ) +
  scale_fill_manual(values = location_palette, guide = "none") +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.01, 0.04))) +
  labs(
    title = "Care Pathways to First ICU Entry",
    subtitle = sprintf("Collapsed ADT pathways for OHCA ICU admissions across %d CLIF sites (n=%s)", site_count, comma(n_total)),
    x = NULL,
    y = "Admissions",
    caption = "Intermediate locations are collapsed to direct ICU transfer, no pre-ICU location, a single procedural/ward/stepdown stop, or multiple/other stops."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#17212B", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "#455A64", margin = margin(b = 10)),
    plot.caption = element_text(size = 8.5, color = "#607D8B", hjust = 0, margin = margin(t = 8)),
    axis.text.x = element_text(face = "bold", color = "#263238", size = 10),
    axis.title.y = element_text(color = "#455A64", size = 9.5),
    axis.text.y = element_text(color = "#607D8B", size = 8.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 4, 10)
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
    label = if_else(.data$pct >= 8, sprintf("%s\n%s (%.1f%%)", .data$timing_group, comma(.data$n), .data$pct), "")
  )

later_timing_note <- timing_compact |>
  filter(.data$pct < 8) |>
  summarise(note = paste(sprintf("%s: %s (%.1f%%)", .data$timing_group, comma(.data$n), .data$pct), collapse = "  |  ")) |>
  pull(.data$note)

timing_palette <- c("<=6 h" = "#2F6F73", "6-12 h" = "#6D8FA6", "12-24 h" = "#B8A35C", ">24 h" = "#A75D5D")

timing_plot <- ggplot(timing_compact, aes(x = .data$pct, y = "Timing", fill = .data$timing_group)) +
  geom_col(width = 0.42, color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = .data$label),
    position = position_stack(vjust = 0.5),
    size = 3.0,
    lineheight = 0.92,
    color = "white",
    fontface = "bold"
  ) +
  annotate(
    "richtext",
    x = 101,
    y = 1,
    label = sprintf(
      "<b>Median time to ICU</b><br><span style='font-size:15pt;color:#17212B'>%.1f h</span> <span style='color:#607D8B'>(IQR %.1f-%.1f)</span>",
      timing_summary$median_hours_approx[1],
      timing_summary$iqr_low_hours_approx[1],
      timing_summary$iqr_high_hours_approx[1]
    ),
    hjust = 0,
    vjust = 0.5,
    fill = "#F5F7FA",
    color = "#17212B",
    label.color = "#D5DEE8",
    label.r = unit(0.08, "lines"),
    size = 3.0,
    lineheight = 1.0
  ) +
  annotate(
    "text",
    x = 0,
    y = 1.38,
    label = later_timing_note,
    hjust = 0,
    vjust = 0.5,
    size = 2.9,
    color = "#455A64"
  ) +
  scale_fill_manual(values = timing_palette, guide = "none") +
  scale_x_continuous(limits = c(0, 128), expand = expansion(mult = c(0, 0)), labels = label_percent(scale = 1)) +
  labs(
    title = sprintf("Timing to first ICU entry: %.1f%% within 6 h; %.1f%% within 24 h",
                    timing_summary$within_6h_pct[1], timing_summary$within_24h_pct[1]),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, color = "#17212B", margin = margin(b = 4)),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    plot.margin = margin(0, 10, 10, 10)
  )

combined <- alluvial_plot / timing_plot +
  plot_layout(heights = c(1, 0.18)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

source_csv <- alluvial_wide |>
  transmute(
    path_id,
    cohort = .data$cohort,
    first_location = .data$first_fill,
    intermediate = as.character(.data$intermediate),
    first_icu = .data$first_icu,
    n
  )
readr::write_csv(source_csv, file.path(output_dir, "figure_care_pathway_alluvial_source.csv"))

save_figure <- function(plot, basename, width = 12, height = 6.3) {
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

save_figure(combined, "figure_care_pathway_alluvial", width = 10.5, height = 6.4)
message("Wrote care-pathway alluvial figure to: ", output_dir)

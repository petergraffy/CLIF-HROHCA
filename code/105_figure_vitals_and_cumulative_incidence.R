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
box_root <- "/Users/saborpete/Library/CloudStorage/Box-Box/CLIF/Projects/CLIF-Heat-Related-OHCA"
input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

HEAT_DEFINITION <- Sys.getenv("HEAT_DEFINITION", unset = "heat95")
OUTPUT_SUFFIX <- if (identical(HEAT_DEFINITION, "heat95")) "" else paste0("_", HEAT_DEFINITION)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(scales)
  library(svglite)
  library(tidyr)
})

theme_manuscript <- function(base_size = 13) {
  theme_classic(base_family = "Helvetica", base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.35),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = grid::unit(18, "pt"),
      panel.spacing.y = grid::unit(16, "pt"),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = grid::unit(0.10, "cm"),
      axis.title = element_text(face = "bold", color = "#111827", size = 12.5),
      axis.text = element_text(color = "black", size = 10.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, color = "#111827", size = 11.5),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(color = "#111827", size = 11.5),
      plot.title = element_text(face = "bold", color = "#111827", size = 16, margin = margin(b = 8)),
      plot.margin = margin(8, 14, 8, 16)
    )
}

save_figure <- function(plot, basename, width, height) {
  ggsave(file.path(output_dir, paste0(basename, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  svglite::svglite(file.path(output_dir, paste0(basename, ".svg")), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(output_dir, paste0(basename, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(output_dir, paste0(basename, ".tiff")), width = width, height = height, units = "in", res = 600, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

group_palette <- c("HROHCA" = "#B23A2E", "Non-HROHCA" = "#2F5D8C")

vital_labels <- c(
  "heart_rate" = "Heart rate, bpm",
  "map" = "Mean arterial pressure, mm Hg",
  "sbp" = "Systolic blood pressure, mm Hg",
  "respiratory_rate" = "Respiratory rate, breaths/min",
  "spo2" = "SpO2, %",
  "temp_c" = "Temperature, °C"
)

vitals_raw <- read.csv(
  file.path(input_dir, "pooled_heat_related_hourly_vital_median_differences.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

vitals <- vitals_raw |>
  filter(
    .data$heat_definition == .env$HEAT_DEFINITION,
    .data$variable %in% names(vital_labels),
    .data$icu_hour >= 0,
    .data$icu_hour <= 72
  ) |>
  mutate(
    vital = factor(vital_labels[.data$variable], levels = vital_labels),
    significant = .data$approx_p_value < 0.05
  )

vitals_long <- bind_rows(
  vitals |>
    transmute(
      icu_hour = .data$icu_hour,
      vital = .data$vital,
      variable = .data$variable,
      group = "HROHCA",
      median_value = .data$heat_weighted_median,
      n_patients = .data$heat_n_patients,
      n_sites = .data$heat_n_sites
    ),
  vitals |>
    transmute(
      icu_hour = .data$icu_hour,
      vital = .data$vital,
      variable = .data$variable,
      group = "Non-HROHCA",
      median_value = .data$non_heat_weighted_median,
      n_patients = .data$non_heat_n_patients,
      n_sites = .data$non_heat_n_sites
  )
) |>
  mutate(group = factor(.data$group, levels = c("HROHCA", "Non-HROHCA")))

smooth_vital_group <- function(dat) {
  dat |>
    group_by(.data$vital, .data$variable, .data$group) |>
    arrange(.data$icu_hour, .by_group = TRUE) |>
    mutate(
      smoothed_median = stats::filter(.data$median_value, rep(1 / 5, 5), sides = 2),
      smoothed_median = if_else(is.na(as.numeric(.data$smoothed_median)), .data$median_value, as.numeric(.data$smoothed_median))
    ) |>
    ungroup()
}

smooth_vital_difference <- function(dat) {
  dat |>
    group_by(.data$vital, .data$variable) |>
    arrange(.data$icu_hour, .by_group = TRUE) |>
    mutate(
      smoothed_difference = stats::filter(.data$approx_difference, rep(1 / 5, 5), sides = 2),
      smoothed_difference = if_else(is.na(as.numeric(.data$smoothed_difference)), .data$approx_difference, as.numeric(.data$smoothed_difference))
    ) |>
    ungroup()
}

vitals_long_smoothed <- smooth_vital_group(vitals_long)
vitals_smoothed <- smooth_vital_difference(vitals)

make_vital_pair <- function(vital_name, col_index = 1L) {
  group_dat <- filter(vitals_long_smoothed, .data$vital == vital_name)
  diff_dat <- filter(vitals_smoothed, .data$vital == vital_name)
  show_y <- col_index == 1L

  trajectory <- ggplot(group_dat, aes(.data$icu_hour, .data$smoothed_median, color = .data$group)) +
    geom_line(linewidth = 1.05) +
    scale_color_manual(values = group_palette) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0.08, 0.12))) +
    labs(title = as.character(vital_name), x = NULL, y = if (show_y) "Pooled median" else NULL) +
    theme_manuscript() +
    theme(
      plot.title = element_text(face = "bold", color = "#111827", size = 11.5, hjust = 0, margin = margin(b = 4)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 10.8, margin = margin(r = 5)),
      legend.position = "none",
      plot.margin = margin(4, 8, 6, 8)
    )

  difference <- ggplot(diff_dat, aes(.data$icu_hour, .data$smoothed_difference)) +
    geom_hline(yintercept = 0, color = "#6B7280", linewidth = 0.45, linetype = "dashed") +
    geom_line(color = "#17212B", linewidth = 1.0) +
    geom_point(
      data = filter(diff_dat, .data$significant),
      color = "#B23A2E",
      fill = "#B23A2E",
      size = 1.55,
      alpha = 0.95
    ) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(accuracy = 0.1), expand = expansion(mult = c(0.14, 0.16))) +
    labs(x = "ICU hour", y = if (show_y) "Median difference" else NULL) +
    theme_manuscript() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title.y = element_text(size = 10.8, margin = margin(r = 5)),
      plot.margin = margin(8, 8, 16, 8)
    )

  trajectory / patchwork::plot_spacer() / difference +
    patchwork::plot_layout(heights = c(0.58, 0.045, 0.42))
}

vital_pair_names <- levels(vitals_long_smoothed$vital)
vital_pairs <- Map(make_vital_pair, vital_pair_names, rep(1:3, length.out = length(vital_pair_names)))

legend_dat <- data.frame(
  icu_hour = rep(c(0, 1), 2),
  median_value = c(0, 1, 0, 1),
  group = factor(rep(c("HROHCA", "Non-HROHCA"), each = 2), levels = c("HROHCA", "Non-HROHCA"))
)

vital_legend <- ggplot(legend_dat, aes(.data$icu_hour, .data$median_value, color = .data$group)) +
  geom_line(linewidth = 1.05, alpha = 0) +
  scale_color_manual(values = group_palette) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.05))) +
  theme_void(base_family = "Helvetica", base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(color = "#111827", size = 11.5),
    plot.margin = margin(0, 0, 0, 0)
  )

vital_grid <- patchwork::wrap_plots(vital_pairs, ncol = 3)

vital_plot <- vital_legend / vital_grid +
  patchwork::plot_layout(heights = c(0.045, 0.955)) +
  patchwork::plot_annotation(
    title = "Hourly vital sign trajectories after OHCA",
    theme = theme(
      plot.title = element_text(
        family = "Helvetica",
        face = "bold",
        size = 16,
        color = "#111827",
        hjust = 0,
        margin = margin(b = 8)
      )
    )
  )

write.csv(vitals_long, file.path(output_dir, paste0("figure_hrohca_vital_trajectories_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
write.csv(vitals_long_smoothed, file.path(output_dir, paste0("figure_hrohca_vital_trajectories_smoothed_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
write.csv(vitals, file.path(output_dir, paste0("figure_hrohca_vital_trajectory_differences_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
write.csv(vitals_smoothed, file.path(output_dir, paste0("figure_hrohca_vital_trajectory_differences_smoothed_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
save_figure(vital_plot, paste0("figure_hrohca_vital_trajectories", OUTPUT_SUFFIX), width = 11.2, height = 11.4)

site_export_files <- list.files(
  box_root,
  pattern = "^[A-Za-z0-9]+_heat_related_hourly_cumulative_incidence[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)

read_export <- function(path) {
  out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"site_name" %in% names(out)) {
    parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
    idx <- match("CLIF-Heat-Related-OHCA", parts)
    out$site_name <- if (!is.na(idx) && length(parts) > idx) parts[[idx + 1L]] else sub("_.*$", "", basename(path))
  }
  out$source_file <- path
  out
}

cumulative_site <- bind_rows(lapply(site_export_files, read_export)) |>
  mutate(
    n_group = suppressWarnings(as.numeric(.data$n_group)),
    n_events = suppressWarnings(as.numeric(.data$n_events)),
    icu_hour = suppressWarnings(as.numeric(.data$icu_hour))
  ) |>
  filter(.data$heat_definition == .env$HEAT_DEFINITION, .data$icu_hour >= 0, .data$icu_hour <= 72)

cumulative_pooled <- cumulative_site |>
  group_by(.data$heat_definition, .data$heat_related_ohca, .data$event, .data$icu_hour) |>
  summarise(
    k_sites = n_distinct(.data$site_name),
    n_group = sum(.data$n_group, na.rm = TRUE),
    n_events = sum(.data$n_events, na.rm = TRUE),
    cumulative_pct = 100 * .data$n_events / .data$n_group,
    .groups = "drop"
  ) |>
  mutate(
    group = factor(
      if_else(.data$heat_related_ohca == "Heat-related OHCA", "HROHCA", "Non-HROHCA"),
      levels = c("HROHCA", "Non-HROHCA")
    )
  )

event_labels <- c(
  "Hospital death" = "Hospital death",
  "Death or hospice" = "Death or hospice",
  "Discharged alive without hospice" = "Discharged alive without hospice",
  "First IMV" = "First IMV",
  "First vasopressor" = "First vasopressor",
  "First CRRT" = "First CRRT"
)

outcome_events <- c("Hospital death", "Death or hospice", "Discharged alive without hospice")
support_events <- c("First IMV", "First vasopressor", "First CRRT")

plot_cumulative <- function(dat, events, title, basename) {
  plot_dat <- dat |>
    filter(.data$event %in% events) |>
    mutate(event = factor(event_labels[.data$event], levels = event_labels[events]))

  p <- ggplot(plot_dat, aes(.data$icu_hour, .data$cumulative_pct, color = .data$group)) +
    geom_step(linewidth = 1.05) +
    facet_wrap(~event, nrow = 1) +
    scale_color_manual(values = group_palette) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(suffix = "%", accuracy = 1), expand = expansion(mult = c(0.04, 0.10))) +
    labs(
      title = title,
      x = "ICU hour",
      y = "Cumulative incidence"
    ) +
    theme_manuscript()

  save_figure(p, paste0(basename, OUTPUT_SUFFIX), width = 9.6, height = 4.4)
  invisible(p)
}

plot_cumulative_with_difference <- function(dat, events, title, basename) {
  plot_dat <- dat |>
    filter(.data$event %in% events) |>
    mutate(event = factor(event_labels[.data$event], levels = event_labels[events]))

  diff_dat <- plot_dat |>
    select("event", "icu_hour", "group", "n_group", "n_events", "cumulative_pct") |>
    tidyr::pivot_wider(
      names_from = "group",
      values_from = c("n_group", "n_events", "cumulative_pct"),
      names_sep = "__"
    ) |>
    mutate(
      diff_pct = .data$`cumulative_pct__HROHCA` - .data$`cumulative_pct__Non-HROHCA`,
      p_value = mapply(
        function(x_heat, x_non, n_heat, n_non) {
          if (!all(is.finite(c(x_heat, x_non, n_heat, n_non))) || n_heat <= 0 || n_non <= 0) return(NA_real_)
          suppressWarnings(stats::prop.test(x = c(x_heat, x_non), n = c(n_heat, n_non))$p.value)
        },
        .data$`n_events__HROHCA`,
        .data$`n_events__Non-HROHCA`,
        .data$`n_group__HROHCA`,
        .data$`n_group__Non-HROHCA`
      ),
      significant = .data$p_value < 0.05
    )

  incidence_plot <- ggplot(plot_dat, aes(.data$icu_hour, .data$cumulative_pct, color = .data$group)) +
    geom_step(linewidth = 1.05) +
    facet_wrap(~event, nrow = 1) +
    scale_color_manual(values = group_palette) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(suffix = "%", accuracy = 1), expand = expansion(mult = c(0.04, 0.10))) +
    labs(
      title = title,
      x = NULL,
      y = "Cumulative incidence"
    ) +
    theme_manuscript()

  difference_plot <- ggplot(diff_dat, aes(.data$icu_hour, .data$diff_pct)) +
    geom_hline(yintercept = 0, color = "#6B7280", linewidth = 0.45, linetype = "dashed") +
    geom_line(color = "#17212B", linewidth = 1.0) +
    geom_point(
      data = filter(diff_dat, .data$significant),
      color = "#B23A2E",
      fill = "#B23A2E",
      size = 1.6,
      alpha = 0.95
    ) +
    facet_wrap(~event, nrow = 1) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(suffix = " pp", accuracy = 1), expand = expansion(mult = c(0.16, 0.18))) +
    labs(
      x = "ICU hour",
      y = "Difference (pp)\nHROHCA - non-HROHCA"
    ) +
    theme_manuscript() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      strip.text = element_blank()
    )

  p <- incidence_plot / difference_plot +
    patchwork::plot_layout(heights = c(0.62, 0.38))

  write.csv(diff_dat, file.path(output_dir, paste0(basename, "_difference_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
  save_figure(p, paste0(basename, OUTPUT_SUFFIX), width = 9.6, height = 6.3)
  invisible(p)
}

build_outcome_difference <- function(plot_dat) {
  plot_dat |>
    select("event", "icu_hour", "group", "n_group", "n_events", "cumulative_pct") |>
    tidyr::pivot_wider(
      names_from = "group",
      values_from = c("n_group", "n_events", "cumulative_pct"),
      names_sep = "__"
    ) |>
    mutate(
      diff_pct = .data$`cumulative_pct__HROHCA` - .data$`cumulative_pct__Non-HROHCA`,
      p_value = mapply(
        function(x_heat, x_non, n_heat, n_non) {
          if (!all(is.finite(c(x_heat, x_non, n_heat, n_non))) || n_heat <= 0 || n_non <= 0) return(NA_real_)
          suppressWarnings(stats::prop.test(x = c(x_heat, x_non), n = c(n_heat, n_non))$p.value)
        },
        .data$`n_events__HROHCA`,
        .data$`n_events__Non-HROHCA`,
        .data$`n_group__HROHCA`,
        .data$`n_group__Non-HROHCA`
      ),
      significant = .data$p_value < 0.05
    )
}

support_labels <- c(
  "Invasive mechanical ventilation" = "IMV",
  "Vasopressor infusion" = "Vasopressors",
  "CRRT" = "CRRT"
)

support_raw <- read.csv(
  file.path(input_dir, "pooled_heat_related_hourly_support_prevalence_differences.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

support_combined <- support_raw |>
  filter(
    .data$heat_definition == .env$HEAT_DEFINITION,
    .data$variable %in% names(support_labels),
    .data$icu_hour >= 0,
    .data$icu_hour <= 72
  ) |>
  mutate(support = factor(support_labels[.data$variable], levels = support_labels))

support_combined_long <- bind_rows(
  support_combined |>
    transmute(
      icu_hour = .data$icu_hour,
      panel = .data$support,
      group = "HROHCA",
      pct = .data$heat_prevalence_pct
    ),
  support_combined |>
    transmute(
      icu_hour = .data$icu_hour,
      panel = .data$support,
      group = "Non-HROHCA",
      pct = .data$non_heat_prevalence_pct
    )
) |>
  mutate(group = factor(.data$group, levels = c("HROHCA", "Non-HROHCA")))

make_combined_support_outcome_figure <- function() {
  outcome_plot_dat <- cumulative_pooled |>
    filter(.data$event %in% outcome_events) |>
    mutate(event = factor(event_labels[.data$event], levels = event_labels[outcome_events]))

  outcome_diff_dat <- build_outcome_difference(outcome_plot_dat)

  support_plot <- ggplot(support_combined_long, aes(.data$icu_hour, .data$pct, color = .data$group)) +
    geom_line(linewidth = 1.05) +
    facet_wrap(~panel, nrow = 1) +
    scale_color_manual(values = group_palette) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(suffix = "%", accuracy = 1), expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = NULL, y = "Receiving support") +
    theme_manuscript() +
    theme(plot.title = element_blank())

  support_difference_plot <- ggplot(support_combined, aes(.data$icu_hour, .data$absolute_difference_pct)) +
    geom_hline(yintercept = 0, color = "#6B7280", linewidth = 0.45, linetype = "dashed") +
    geom_line(color = "#17212B", linewidth = 1.0) +
    geom_point(
      data = filter(support_combined, .data$p_value < 0.05),
      color = "#B23A2E",
      fill = "#B23A2E",
      size = 1.6,
      alpha = 0.95
    ) +
    facet_wrap(~support, nrow = 1) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(suffix = " pp", accuracy = 1), expand = expansion(mult = c(0.16, 0.18))) +
    labs(x = NULL, y = "Difference (pp)\nHROHCA - non-HROHCA") +
    theme_manuscript() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      strip.text = element_blank()
    )

  outcome_plot <- ggplot(outcome_plot_dat, aes(.data$icu_hour, .data$cumulative_pct, color = .data$group)) +
    geom_step(linewidth = 1.05) +
    facet_wrap(~event, nrow = 1) +
    scale_color_manual(values = group_palette) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(suffix = "%", accuracy = 1), expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = NULL, y = "Cumulative incidence") +
    theme_manuscript() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, color = "#111827", size = 11.2, lineheight = 0.92)
    )

  difference_plot <- ggplot(outcome_diff_dat, aes(.data$icu_hour, .data$diff_pct)) +
    geom_hline(yintercept = 0, color = "#6B7280", linewidth = 0.45, linetype = "dashed") +
    geom_line(color = "#17212B", linewidth = 1.0) +
    geom_point(
      data = filter(outcome_diff_dat, .data$significant),
      color = "#B23A2E",
      fill = "#B23A2E",
      size = 1.6,
      alpha = 0.95
    ) +
    facet_wrap(~event, nrow = 1) +
    scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = label_number(suffix = " pp", accuracy = 1), expand = expansion(mult = c(0.16, 0.18))) +
    labs(x = "ICU hour", y = "Difference (pp)\nHROHCA - non-HROHCA") +
    theme_manuscript() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      strip.text = element_blank()
    )

  combined <- support_plot / support_difference_plot / outcome_plot / difference_plot +
    patchwork::plot_layout(heights = c(0.27, 0.21, 0.29, 0.23), guides = "collect") +
    patchwork::plot_annotation(
      title = "Hourly support trajectories and cumulative outcome incidence after OHCA",
      theme = theme(
        plot.title = element_text(
          family = "Helvetica",
          face = "bold",
          size = 17,
          color = "#111827",
          hjust = 0.5,
          margin = margin(b = 8)
        ),
        legend.position = "top"
      )
    ) &
    theme(legend.position = "top")

  write.csv(support_combined_long, file.path(output_dir, paste0("figure_hrohca_support_outcome_combined_support_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
  write.csv(support_combined, file.path(output_dir, paste0("figure_hrohca_support_outcome_combined_support_difference_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
  write.csv(outcome_diff_dat, file.path(output_dir, paste0("figure_hrohca_support_outcome_combined_difference_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
  save_figure(combined, paste0("figure_hrohca_support_outcome_combined", OUTPUT_SUFFIX), width = 11.6, height = 11.8)
  invisible(combined)
}

write.csv(cumulative_site, file.path(input_dir, paste0("all_sites_heat_related_hourly_cumulative_incidence", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
write.csv(cumulative_pooled, file.path(input_dir, paste0("pooled_heat_related_hourly_cumulative_incidence", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)
write.csv(cumulative_pooled, file.path(output_dir, paste0("figure_hrohca_cumulative_incidence_source", OUTPUT_SUFFIX, ".csv")), row.names = FALSE)

plot_cumulative_with_difference(
  cumulative_pooled,
  outcome_events,
  "Hourly cumulative outcome incidence after OHCA",
  "figure_hrohca_cumulative_outcome_incidence"
)
plot_cumulative(
  cumulative_pooled,
  support_events,
  "Hourly cumulative organ-support initiation after OHCA",
  "figure_hrohca_cumulative_support_incidence"
)
make_combined_support_outcome_figure()

message("Wrote vital sign and cumulative incidence figures to ", output_dir)

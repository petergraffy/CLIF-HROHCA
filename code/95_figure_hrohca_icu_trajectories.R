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
input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
HEAT_DEFINITION <- Sys.getenv("HEAT_DEFINITION", unset = "heat95")
OUTPUT_SUFFIX <- if (HEAT_DEFINITION == "heat95") "" else paste0("_", HEAT_DEFINITION)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(scales)
  library(svglite)
})

theme_trajectory <- function(base_size = 12) {
  theme_classic(base_family = "Helvetica", base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.35),
      panel.spacing.x = grid::unit(18, "pt"),
      panel.spacing.y = grid::unit(16, "pt"),
      axis.title = element_text(face = "bold", color = "#111827", size = 12),
      axis.title.y = element_text(margin = margin(r = 7)),
      axis.text = element_text(color = "black", size = 10.5),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = grid::unit(0.10, "cm"),
      strip.text = element_text(face = "bold", color = "#111827", hjust = 0, size = 11.5),
      strip.background = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(color = "#111827", size = 11.5),
      plot.title = element_text(face = "bold", size = 16, color = "#111827", margin = margin(b = 8)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(8, 14, 8, 16)
    )
}

save_figure <- function(plot, basename, width, height) {
  ggsave(
    file.path(output_dir, paste0(basename, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
  svglite::svglite(file.path(output_dir, paste0(basename, ".svg")), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(output_dir, paste0(basename, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(
    file.path(output_dir, paste0(basename, ".tiff")),
    width = width,
    height = height,
    units = "in",
    res = 600,
    compression = "lzw",
    background = "white"
  )
  print(plot)
  dev.off()
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)))
}

support_path <- file.path(input_dir, "pooled_heat_related_hourly_support_prevalence_differences.csv")
lab_path <- file.path(input_dir, "pooled_heat_related_hourly_lab_median_differences.csv")
crrt_path <- file.path(input_dir, "pooled_heat_related_crrt_window_differences.csv")

support_raw <- read.csv(support_path, stringsAsFactors = FALSE, check.names = FALSE)
lab_raw <- read.csv(lab_path, stringsAsFactors = FALSE, check.names = FALSE)
crrt_raw <- read.csv(crrt_path, stringsAsFactors = FALSE, check.names = FALSE)

support_levels <- c("Invasive mechanical ventilation", "Vasopressor infusion", "CRRT")
support_labels <- c(
  "Invasive mechanical ventilation" = "IMV",
  "Vasopressor infusion" = "Vasopressors",
  "CRRT" = "CRRT"
)

support <- support_raw |>
  filter(.data$heat_definition == .env$HEAT_DEFINITION, .data$variable %in% support_levels) |>
  mutate(
    support = factor(support_labels[.data$variable], levels = support_labels[support_levels]),
    significant = .data$p_value < 0.05
  )

support_long <- bind_rows(
  support |>
    transmute(
      icu_hour = .data$icu_hour,
      support = .data$support,
      group = "HROHCA",
      prevalence_pct = .data$heat_prevalence_pct
    ),
  support |>
    transmute(
      icu_hour = .data$icu_hour,
      support = .data$support,
      group = "Non-HROHCA",
      prevalence_pct = .data$non_heat_prevalence_pct
    )
) |>
  mutate(group = factor(.data$group, levels = c("HROHCA", "Non-HROHCA")))

support_prevalence_plot <- ggplot(support_long, aes(.data$icu_hour, .data$prevalence_pct, color = .data$group)) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~support, nrow = 1) +
  scale_color_manual(values = c("HROHCA" = "#B23A2E", "Non-HROHCA" = "#2F5D8C")) +
  scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(labels = label_number(suffix = "%"), expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "Hourly organ-support trajectories after OHCA",
    x = "ICU hour",
    y = "Receiving support (%)"
  ) +
  theme_trajectory()

support_difference_plot <- ggplot(support, aes(.data$icu_hour, .data$absolute_difference_pct, color = .data$support)) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "#6B7280", linetype = "dashed") +
  geom_line(linewidth = 1.05) +
  geom_point(data = filter(support, .data$significant), size = 1.8, alpha = 0.95) +
  facet_wrap(~support, nrow = 1) +
  scale_color_manual(values = c("IMV" = "#6B7280", "Vasopressors" = "#8A4F0F", "CRRT" = "#3F6B3C"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(labels = label_number(suffix = " pp"), expand = expansion(mult = c(0.12, 0.22))) +
  labs(
    x = "ICU hour",
    y = "Difference (pp)\nHROHCA - non-HROHCA"
  ) +
  theme_trajectory() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.y = element_text(size = 11, margin = margin(r = 10))
  )

support_figure <- support_prevalence_plot / support_difference_plot +
  patchwork::plot_layout(heights = c(0.58, 0.42))

lab_labels <- c(
  "lactate" = "Lactate, mmol/L",
  "bicarbonate" = "Bicarbonate, mEq/L",
  "ph_arterial" = "Arterial pH",
  "creatinine" = "Creatinine, mg/dL",
  "bun" = "BUN, mg/dL",
  "potassium" = "Potassium, mEq/L",
  "magnesium" = "Magnesium, mg/dL",
  "phosphate" = "Phosphate, mg/dL",
  "wbc" = "WBC, K/uL"
)
lab_order <- names(lab_labels)
key_hours <- seq(0, 72, by = 6)

labs <- lab_raw |>
  filter(
    .data$heat_definition == .env$HEAT_DEFINITION,
    .data$variable %in% lab_order,
    .data$icu_hour %in% key_hours,
    .data$approx_k_sites >= 3,
    .data$heat_n_patients >= 5,
    .data$non_heat_n_patients >= 100,
    is.finite(.data$approx_difference),
    is.finite(.data$approx_ci_low),
    is.finite(.data$approx_ci_high)
  ) |>
  mutate(
    lab = factor(lab_labels[.data$variable], levels = lab_labels[lab_order]),
    significant = .data$approx_p_value < 0.05,
    direction = case_when(
      .data$approx_difference > 0 ~ "Higher in HROHCA",
      .data$approx_difference < 0 ~ "Lower in HROHCA",
      TRUE ~ "No difference"
    )
  )

lab_figure <- ggplot(labs, aes(.data$icu_hour, .data$approx_difference)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "#6B7280", linetype = "dashed") +
  geom_ribbon(aes(ymin = .data$approx_ci_low, ymax = .data$approx_ci_high), fill = "#9CA3AF", alpha = 0.20) +
  geom_line(color = "#1F4E5F", linewidth = 1.05) +
  geom_point(aes(fill = .data$significant), shape = 21, size = 2.35, color = "#1F4E5F", stroke = 0.4) +
  facet_wrap(~lab, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("TRUE" = "#B23A2E", "FALSE" = "white"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "Laboratory 6-hourly trajectories after OHCA",
    x = "ICU hour",
    y = "Median difference"
  ) +
  theme_trajectory()

crrt_heat95 <- crrt_raw |>
  filter(.data$heat_definition == .env$HEAT_DEFINITION) |>
  mutate(
    window = factor(.data$window, levels = c("0-24h", "0-72h", "0-168h")),
    label = paste0(sprintf("%+.1f pp", .data$absolute_difference_pct), "\n", fmt_p(.data$p_value))
  )

crrt_figure <- ggplot(crrt_heat95, aes(.data$window, .data$absolute_difference_pct)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "#6B7280", linetype = "dashed") +
  geom_col(width = 0.62, fill = "#3F6B3C") +
  geom_text(aes(label = .data$label), vjust = -0.35, size = 2.7, color = "#111827", family = "Helvetica") +
  scale_y_continuous(labels = label_number(suffix = " pp"), expand = expansion(mult = c(0.04, 0.20))) +
  labs(
    title = "Windowed CRRT use was numerically higher in HROHCA",
    x = NULL,
    y = "Absolute difference\nHROHCA minus non-HROHCA"
  ) +
  theme_trajectory() +
  theme(legend.position = "none")

combined_figure <- support_prevalence_plot / support_difference_plot / lab_figure +
  patchwork::plot_layout(heights = c(0.30, 0.24, 0.46)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

write.csv(support, file.path(output_dir, paste0("figure_hrohca_support_trajectories", OUTPUT_SUFFIX, "_source.csv")), row.names = FALSE)
write.csv(labs, file.path(output_dir, paste0("figure_hrohca_lab_trajectories", OUTPUT_SUFFIX, "_source.csv")), row.names = FALSE)
write.csv(crrt_heat95, file.path(output_dir, paste0("figure_hrohca_crrt_windows", OUTPUT_SUFFIX, "_source.csv")), row.names = FALSE)

save_figure(support_figure, paste0("figure_hrohca_support_trajectories", OUTPUT_SUFFIX), width = 9.2, height = 5.8)
save_figure(lab_figure, paste0("figure_hrohca_lab_trajectories", OUTPUT_SUFFIX), width = 10.4, height = 8.4)
save_figure(crrt_figure, paste0("figure_hrohca_crrt_windows", OUTPUT_SUFFIX), width = 4.8, height = 4.1)
save_figure(combined_figure, paste0("figure_hrohca_icu_trajectories_combined", OUTPUT_SUFFIX), width = 10.4, height = 13.2)

message("Wrote HROHCA ICU trajectory figures to ", output_dir)

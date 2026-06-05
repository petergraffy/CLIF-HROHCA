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
box_root <- "/Users/saborpete/Library/CloudStorage/Box-Box/CLIF/Projects/CLIF-Heat-Related-OHCA"
input_path <- file.path(repo_root, "output", "final", "federated_pooled", "pooled_dlnm_random_effects_curves.csv")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ragg)
  library(scales)
  library(svglite)
})

curve_all <- read.csv(input_path, stringsAsFactors = FALSE) |>
  filter(
    .data$stratum == "Overall",
    .data$model == "sensitivity_mrt_reference",
    .data$reference_type == "mrt"
  ) |>
  mutate(
    tmax_mean_c = as.numeric(.data$tmax_mean_c),
    cumulative_rr = as.numeric(.data$cumulative_rr),
    cumulative_rr_low = as.numeric(.data$cumulative_rr_low),
    cumulative_rr_high = as.numeric(.data$cumulative_rr_high)
  )

max_k_sites <- max(curve_all$k_sites, na.rm = TRUE)

curve <- curve_all |>
  filter(.data$k_sites == max_k_sites) |>
  arrange(.data$tmax_mean_c)

y_upper <- max(3.75, max(curve$cumulative_rr_high, na.rm = TRUE) * 1.05)

threshold_files <- list.files(
  box_root,
  pattern = "^[^_]+_heat_related_ohca_thresholds[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
threshold_files <- threshold_files[grepl("/federated_exports/", threshold_files)]
thresholds <- bind_rows(lapply(threshold_files, read.csv, stringsAsFactors = FALSE, check.names = FALSE)) |>
  filter(.data$heat_definition %in% c("heat90", "heat95")) |>
  mutate(heat_threshold_tmax_c = as.numeric(.data$heat_threshold_tmax_c))

heat90_threshold <- thresholds |>
  filter(.data$heat_definition == "heat90") |>
  summarise(threshold = median(.data$heat_threshold_tmax_c, na.rm = TRUE)) |>
  pull(.data$threshold)

heat95_threshold <- thresholds |>
  filter(.data$heat_definition == "heat95") |>
  summarise(threshold = median(.data$heat_threshold_tmax_c, na.rm = TRUE)) |>
  pull(.data$threshold)

mrt_temp <- curve$tmax_mean_c[which.min(curve$cumulative_rr)]

label_targets <- c(mrt_temp, heat90_threshold, heat95_threshold)
label_points <- do.call(rbind, lapply(label_targets, function(target) {
  curve[which.min(abs(curve$tmax_mean_c - target)), , drop = FALSE]
})) |>
  mutate(
    label_type = c("MRT", "90th percentile", "95th percentile"),
    label_x = c(26.8, 31.0, 35.3),
    label_y = c(0.78, 1.78, 2.22),
    label = sprintf(
      "%s\n%.1f \u00b0C: RR %.2f\n95%% CI %.2f-%.2f",
      .data$label_type,
      .data$tmax_mean_c,
      .data$cumulative_rr,
      .data$cumulative_rr_low,
      .data$cumulative_rr_high
    )
  )

figure1 <- ggplot(curve, aes(x = .data$tmax_mean_c, y = .data$cumulative_rr)) +
  annotate(
    "rect",
    xmin = heat90_threshold,
    xmax = heat95_threshold,
    ymin = -Inf,
    ymax = Inf,
    fill = "#F6D7B0",
    alpha = 0.20
  ) +
  annotate(
    "rect",
    xmin = heat95_threshold,
    xmax = max(curve$tmax_mean_c),
    ymin = -Inf,
    ymax = Inf,
    fill = "#F1B56A",
    alpha = 0.18
  ) +
  geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dashed", color = "#5F6368") +
  geom_vline(xintercept = mrt_temp, linewidth = 0.45, linetype = "dotted", color = "#5F6368") +
  geom_ribbon(
    aes(ymin = .data$cumulative_rr_low, ymax = .data$cumulative_rr_high),
    fill = "#7EAED3",
    alpha = 0.35,
    linewidth = 0
  ) +
  geom_line(linewidth = 1.35, color = "#0B3C5D", lineend = "round") +
  geom_segment(
    data = label_points,
    aes(
      x = .data$label_x,
      y = .data$label_y,
      xend = .data$tmax_mean_c,
      yend = .data$cumulative_rr
    ),
    inherit.aes = FALSE,
    color = "#6B7280",
    linewidth = 0.28
  ) +
  geom_point(data = label_points, size = 2.3, color = "#0B3C5D", fill = "white", shape = 21, stroke = 0.8) +
  geom_label(
    data = label_points,
    aes(x = .data$label_x, y = .data$label_y, label = .data$label),
    inherit.aes = FALSE,
    size = 2.75,
    family = "Helvetica",
    color = "#1F2933",
    fill = "white",
    linewidth = 0.18,
    label.padding = unit(0.13, "lines"),
    lineheight = 1.05
  ) +
  scale_x_continuous(
    breaks = seq(18, 37, by = 2),
    minor_breaks = seq(18, 37, by = 1),
    expand = expansion(mult = c(0.015, 0.035))
  ) +
  scale_y_continuous(
    breaks = c(0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4),
    labels = label_number(accuracy = 0.01),
    limits = c(0.55, y_upper),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Cumulative relative risk of OHCA ICU admission per daily maximum temperature",
    x = "County-level daily maximum temperature, Tmax (\u00b0C)",
    y = "Cumulative relative risk of OHCA"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5, color = "#111827", margin = margin(b = 4)),
    axis.title = element_text(face = "bold", size = 9.8, color = "#111827"),
    axis.text = element_text(size = 8.8, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 12, 10, 10)
  )

write.csv(
  curve,
  file.path(output_dir, "figure1_dlnm_mrt_curve_source.csv"),
  row.names = FALSE
)

write.csv(
  label_points,
  file.path(output_dir, "figure1_dlnm_mrt_curve_label_points.csv"),
  row.names = FALSE
)

ggsave(
  file.path(output_dir, "figure1_dlnm_mrt_curve.pdf"),
  figure1,
  width = 7.6,
  height = 5.0,
  device = cairo_pdf,
  bg = "white"
)

svglite::svglite(file.path(output_dir, "figure1_dlnm_mrt_curve.svg"), width = 7.6, height = 5.0, bg = "white")
print(figure1)
dev.off()

ragg::agg_png(
  file.path(output_dir, "figure1_dlnm_mrt_curve.png"),
  width = 7.6,
  height = 5.0,
  units = "in",
  res = 600,
  background = "white"
)
print(figure1)
dev.off()

ragg::agg_tiff(
  file.path(output_dir, "figure1_dlnm_mrt_curve.tiff"),
  width = 7.6,
  height = 5.0,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(figure1)
dev.off()

message("Wrote Figure 1 DLNM curve to ", output_dir)

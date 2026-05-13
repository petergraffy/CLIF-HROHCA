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

curve <- read.csv(input_path, stringsAsFactors = FALSE) |>
  filter(
    .data$stratum == "Overall",
    .data$model == "sensitivity_mrt_reference",
    .data$reference_type == "mrt",
    .data$k_sites == 8
  ) |>
  mutate(
    tmax_mean_c = as.numeric(.data$tmax_mean_c),
    cumulative_rr = as.numeric(.data$cumulative_rr),
    cumulative_rr_low = as.numeric(.data$cumulative_rr_low),
    cumulative_rr_high = as.numeric(.data$cumulative_rr_high)
  ) |>
  arrange(.data$tmax_mean_c)

label_points <- curve |>
  filter(.data$tmax_mean_c %in% c(30, 32.5, 35)) |>
  mutate(
    label = dplyr::case_when(
      .data$tmax_mean_c == 30 ~ "MRT reference\nRR 1.01",
      .data$tmax_mean_c == 32.5 ~ "32.5 C: RR 1.34\n95% CI 1.06-1.70",
      .data$tmax_mean_c == 35 ~ "35.0 C: RR 1.98\n95% CI 1.43-2.73",
      TRUE ~ ""
    )
  )

caption_text <- paste(
  "Pooled random-effects distributed lag nonlinear model (DLNM).",
  "\nCurve is restricted to the Tmax range where all 8 sites contributed estimates.",
  "\nRibbon shows 95% CI; RR is referenced to site-specific minimum-risk temperature (MRT)."
)

figure1 <- ggplot(curve, aes(x = .data$tmax_mean_c, y = .data$cumulative_rr)) +
  annotate(
    "rect",
    xmin = 32.5,
    xmax = max(curve$tmax_mean_c),
    ymin = -Inf,
    ymax = Inf,
    fill = "#F6D7B0",
    alpha = 0.22
  ) +
  geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dashed", color = "#5F6368") +
  geom_vline(xintercept = 30, linewidth = 0.45, linetype = "dotted", color = "#5F6368") +
  geom_ribbon(
    aes(ymin = .data$cumulative_rr_low, ymax = .data$cumulative_rr_high),
    fill = "#7EAED3",
    alpha = 0.35,
    linewidth = 0
  ) +
  geom_line(linewidth = 1.35, color = "#0B3C5D", lineend = "round") +
  geom_point(data = label_points, size = 2.3, color = "#0B3C5D", fill = "white", shape = 21, stroke = 0.8) +
  ggrepel::geom_label_repel(
    data = label_points,
    aes(label = .data$label),
    size = 2.75,
    family = "Helvetica",
    color = "#1F2933",
    fill = "white",
    label.size = 0.18,
    label.padding = unit(0.13, "lines"),
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "#6B7280",
    segment.size = 0.35,
    seed = 20260513
  ) +
  scale_x_continuous(
    breaks = seq(18, 37, by = 2),
    minor_breaks = seq(18, 37, by = 1),
    expand = expansion(mult = c(0.015, 0.035))
  ) +
  scale_y_continuous(
    breaks = c(0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5),
    labels = label_number(accuracy = 0.01),
    limits = c(0.72, 3.75),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Higher daily maximum temperature was associated with increased OHCA risk",
    subtitle = "Pooled MRT-referenced DLNM across 8 CLIF sites, 2018-2024",
    x = "County-level daily maximum temperature, Tmax (C)",
    y = "Cumulative relative risk of OHCA",
    caption = caption_text
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5, color = "#111827", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9.8, color = "#374151", margin = margin(b = 10)),
    plot.caption = element_text(size = 7.5, color = "#4B5563", hjust = 0, lineheight = 1.15, margin = margin(t = 8)),
    axis.title = element_text(face = "bold", size = 9.8, color = "#111827"),
    axis.text = element_text(size = 8.8, color = "#374151"),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
    panel.grid.minor = element_line(color = "#F3F4F6", linewidth = 0.25),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 12, 10, 10)
  )

write.csv(
  curve,
  file.path(output_dir, "figure1_dlnm_mrt_curve_source.csv"),
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

#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    ofiles <- vapply(sys.frames(), function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile, character(1))
    ofiles <- stats::na.omit(ofiles)
    if (length(ofiles) > 0) return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = TRUE))
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      active_path <- rstudioapi::getActiveDocumentContext()$path
      if (nzchar(active_path)) return(normalizePath(active_path, winslash = "/", mustWork = TRUE))
    }
    stop("Could not determine script path. Run with Rscript or source the script from RStudio.")
  }
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(pooled_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

all_site <- readr::read_csv(file.path(pooled_dir, "all_site_dlnm_curves.csv"), show_col_types = FALSE)
pooled <- readr::read_csv(file.path(pooled_dir, "pooled_dlnm_random_effects_curves.csv"), show_col_types = FALSE)

family_labels <- c(
  count_dlnm = "OHCA ICU admission count DLNM",
  rate_dlnm = "OHCA per ICU admission rate DLNM",
  case_crossover_dlnm = "Time-stratified case-crossover DLNM"
)

model_labels <- c(
  primary_humidity_adjusted = "Median ref, humidity",
  sensitivity_humidity_pollution_adjusted = "Median ref, humidity + pollution",
  sensitivity_mrt_reference = "MRT ref, humidity + pollution",
  rate_median_humidity_adjusted = "Median ref, humidity",
  rate_median_humidity_pollution_adjusted = "Median ref, humidity + pollution",
  rate_mrt_humidity_pollution_adjusted = "MRT ref, humidity + pollution",
  stratified_humidity_adjusted = "Median ref, humidity",
  stratified_humidity_pollution_adjusted = "Median ref, humidity + pollution",
  stratified_mrt_reference = "MRT ref, humidity + pollution",
  rate_stratified_humidity_adjusted = "Median ref, humidity",
  rate_stratified_humidity_pollution_adjusted = "Median ref, humidity + pollution",
  rate_stratified_mrt_humidity_pollution_adjusted = "MRT ref, humidity + pollution"
)

family_colors <- c(
  count_dlnm = "#1f6f8b",
  rate_dlnm = "#8a5a00",
  case_crossover_dlnm = "#7f3c8d"
)

pretty_model <- function(x) {
  out <- model_labels[x]
  missing <- is.na(out)
  out[missing] <- str_replace_all(x[missing], "_", " ")
  out
}

prepare_plot_data <- function(family, overall_only = FALSE) {
  site_dat <- all_site |>
    filter(.data$export_family == family) |>
    mutate(
      model_label = pretty_model(.data$model),
      stratum = factor(.data$stratum, levels = c("Overall", "Male", "Female", "<65", ">=65", "Black", "Non-Black", "Hispanic", "Non-Hispanic"))
    )
  pooled_dat <- pooled |>
    filter(.data$export_family == family) |>
    mutate(
      model_label = pretty_model(.data$model),
      stratum = factor(.data$stratum, levels = c("Overall", "Male", "Female", "<65", ">=65", "Black", "Non-Black", "Hispanic", "Non-Hispanic"))
    )
  if (overall_only) {
    site_dat <- site_dat |> filter(.data$stratum == "Overall")
    pooled_dat <- pooled_dat |> filter(.data$stratum == "Overall")
  }
  list(site = site_dat, pooled = pooled_dat)
}

make_curve_plot <- function(family, overall_only = FALSE) {
  dat <- prepare_plot_data(family, overall_only)
  title_suffix <- if (overall_only) "Overall" else "All strata"
  p <- ggplot() +
    geom_hline(yintercept = 1, linewidth = 0.35, color = "grey55") +
    geom_line(
      data = dat$site,
      aes(.data$tmax_mean_c, .data$cumulative_rr, group = interaction(.data$site_name, .data$stratum, .data$model, .data$reference_type)),
      color = "grey72",
      alpha = 0.65,
      linewidth = 0.55
    ) +
    geom_ribbon(
      data = dat$pooled,
      aes(.data$tmax_mean_c, ymin = .data$rr_low, ymax = .data$rr_high),
      fill = family_colors[[family]],
      alpha = 0.16
    ) +
    geom_line(
      data = dat$pooled,
      aes(.data$tmax_mean_c, .data$rr),
      color = family_colors[[family]],
      linewidth = 1.15
    ) +
    scale_y_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
      labels = c("0.25", "0.5", "1", "2", "4", "8", "16")
    ) +
    labs(
      title = paste(family_labels[[family]], "-", title_suffix),
      subtitle = "Pooled curve in color with 95% CI; individual site curves in light grey",
      x = "Daily maximum temperature (degrees C)",
      y = "Cumulative relative risk"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

  if (overall_only) {
    p + facet_wrap(~ model_label, nrow = 1)
  } else {
    p + facet_grid(stratum ~ model_label)
  }
}

for (family in names(family_labels)) {
  overall_plot <- make_curve_plot(family, overall_only = TRUE)
  all_strata_plot <- make_curve_plot(family, overall_only = FALSE)
  overall_path <- file.path(figure_dir, paste0("pooled_", family, "_overall_curves.png"))
  all_strata_path <- file.path(figure_dir, paste0("pooled_", family, "_all_strata_curves.png"))
  ggsave(overall_path, overall_plot, width = 12, height = 4.4, dpi = 300)
  ggsave(all_strata_path, all_strata_plot, width = 12.5, height = 14, dpi = 300)
}

message("Wrote pooled DLNM curve figures to ", figure_dir)

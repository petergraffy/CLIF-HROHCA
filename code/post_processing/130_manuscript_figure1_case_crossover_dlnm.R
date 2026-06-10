#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

all_site_path <- file.path(pooled_dir, "all_site_dlnm_curves.csv")
pooled_path <- file.path(pooled_dir, "pooled_dlnm_random_effects_curves.csv")
if (!file.exists(all_site_path) || !file.exists(pooled_path)) {
  stop("Missing pooled DLNM curve inputs. Run code/post_processing/90_pool_federated_results.R first.", call. = FALSE)
}

all_site <- readr::read_csv(all_site_path, show_col_types = FALSE)
pooled <- readr::read_csv(pooled_path, show_col_types = FALSE)

target_filter <- function(df) {
  df |>
    dplyr::filter(
      .data$export_family == "case_crossover_dlnm",
      .data$stratum == "Overall",
      .data$model == "primary_humidity_adjusted",
      .data$reference_type == "median"
    )
}

site_curve <- target_filter(all_site) |>
  dplyr::mutate(
    rr = dplyr::coalesce(.data$cumulative_rr, exp(.data$log_rr)),
    site_name = as.character(.data$site_name)
  )

pooled_curve <- target_filter(pooled) |>
  dplyr::arrange(.data$tmax_mean_c)

if (nrow(site_curve) == 0 || nrow(pooled_curve) == 0) {
  stop("Could not find pooled case-crossover median-reference humidity-adjusted DLNM curves.", call. = FALSE)
}

max_sites <- max(pooled_curve$k_sites, na.rm = TRUE)
y_breaks <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4)
x_limits_c <- c(-20, 38)
fahrenheit_breaks <- seq(0, 100, by = 20)
x_breaks_c <- (fahrenheit_breaks - 32) * 5 / 9

figure1 <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 1, color = "grey35", linewidth = 0.45, linetype = "dashed") +
  ggplot2::geom_line(
    data = site_curve,
    ggplot2::aes(x = .data$tmax_mean_c, y = .data$rr, group = .data$site_name),
    color = "grey67",
    linewidth = 0.75,
    alpha = 0.75
  ) +
  ggplot2::geom_ribbon(
    data = pooled_curve,
    ggplot2::aes(x = .data$tmax_mean_c, ymin = .data$rr_low, ymax = .data$rr_high),
    fill = "#0f6b78",
    alpha = 0.18
  ) +
  ggplot2::geom_line(
    data = pooled_curve,
    ggplot2::aes(x = .data$tmax_mean_c, y = .data$rr),
    color = "#0f6b78",
    linewidth = 1.35
  ) +
  ggplot2::scale_x_continuous(
    breaks = x_breaks_c,
    labels = fahrenheit_breaks,
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_log10(
    breaks = y_breaks,
    labels = format(y_breaks, trim = TRUE),
    expand = ggplot2::expansion(mult = c(0.03, 0.08))
  ) +
  ggplot2::coord_cartesian(xlim = x_limits_c, ylim = c(0.25, 4)) +
  ggplot2::labs(
    title = "Time-stratified case-crossover DLNM",
    x = "Daily maximum temperature (°F)",
    y = "Cumulative relative risk"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "grey15"),
    axis.line = ggplot2::element_line(color = "grey15"),
    axis.ticks = ggplot2::element_line(color = "grey15"),
    plot.margin = ggplot2::margin(12, 16, 12, 12)
  )

png_path <- file.path(figure_dir, "figure1_pooled_case_crossover_dlnm_median_humidity.png")
pdf_path <- file.path(figure_dir, "figure1_pooled_case_crossover_dlnm_median_humidity.pdf")
ggplot2::ggsave(png_path, figure1, width = 7.4, height = 5.2, dpi = 600)
ggplot2::ggsave(pdf_path, figure1, width = 7.4, height = 5.2)

message("Wrote Figure 1 to ", png_path, " and ", pdf_path)

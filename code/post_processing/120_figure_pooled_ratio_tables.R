#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(pooled_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_if_exists <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  readr::read_csv(path, show_col_types = FALSE)
}

clean_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("temperature humidity demographics mechanism", "temp + humidity + demographics + mechanism", x)
  x <- gsub("temperature humidity demographics", "temp + humidity + demographics", x)
  x <- gsub("humidity pollution adjusted", "humidity + pollution", x)
  x <- gsub("humidity adjusted", "humidity", x)
  x <- gsub("mrt reference", "MRT reference", x)
  x <- gsub("median", "median", x)
  x
}

clean_outcome <- function(x) {
  dplyr::case_when(
    x == "limited_brain_function" ~ "Limited brain function",
    x == "anoxic_brain_injury" ~ "Death within 72h",
    TRUE ~ clean_label(x)
  )
}

clean_adjustment <- function(x) {
  dplyr::case_when(
    grepl("mechanism", x) ~ "Mechanism adjusted",
    TRUE ~ "Primary"
  )
}

clean_term_family <- function(x) {
  dplyr::case_when(
    x == "temperature_spline_basis" ~ "Tmax",
    x == "humidity_spline_basis" ~ "Humidity",
    TRUE ~ clean_label(x)
  )
}

trim_ratio_plot <- function(df, lower = 0.05, upper = 20) {
  df |>
    dplyr::mutate(
      plot_low = pmax(.data$ratio_low, lower, na.rm = TRUE),
      plot_high = pmin(.data$ratio_high, upper, na.rm = TRUE),
      clipped_low = .data$ratio_low < lower,
      clipped_high = .data$ratio_high > upper
    )
}

save_forest <- function(df, output_file, title, subtitle, x_label = "Pooled ratio", width = 10, height = 6) {
  if (nrow(df) == 0) return(invisible(FALSE))
  df <- df |>
    dplyr::mutate(row_label = factor(.data$row_label, levels = rev(unique(.data$row_label))))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$ratio, y = .data$row_label, color = .data$effect_measure)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$plot_low, xmax = .data$plot_high),
      orientation = "y",
      width = 0.18,
      linewidth = 0.7,
      alpha = 0.9
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c("RR" = "#0f6b78", "OR" = "#8f7a2f", "SHR" = "#8c3b3b")) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey35"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
  invisible(TRUE)
}

dlnm <- read_if_exists(file.path(pooled_dir, "pooled_dlnm_ratio_results.csv"))
if (nrow(dlnm) > 0) {
  dlnm_fig <- dlnm |>
    dplyr::filter(.data$stratum == "Overall") |>
    dplyr::mutate(
      row_label = paste(clean_label(.data$analysis), clean_label(.data$model), clean_label(.data$contrast), sep = " | ")
    ) |>
    trim_ratio_plot()

  save_forest(
    dlnm_fig,
    file.path(figure_dir, "pooled_dlnm_ratio_forest_overall.png"),
    "Pooled DLNM Relative Ratios",
    "Overall pooled random-effects RR estimates. Error bars are 95% CIs.",
    x_label = "Relative risk",
    width = 11,
    height = max(4.8, 0.34 * nrow(dlnm_fig))
  )

  dlnm_strata_fig <- dlnm |>
    dplyr::filter(.data$stratum != "Overall") |>
    dplyr::mutate(
      row_label = paste(clean_label(.data$analysis), clean_label(.data$model), .data$stratum, clean_label(.data$contrast), sep = " | ")
    ) |>
    trim_ratio_plot()

  save_forest(
    dlnm_strata_fig,
    file.path(figure_dir, "pooled_dlnm_ratio_forest_stratified.png"),
    "Pooled Stratified DLNM Relative Ratios",
    "Stratum-specific pooled random-effects RR estimates. Error bars are 95% CIs.",
    x_label = "Relative risk",
    width = 12,
    height = max(6, 0.28 * nrow(dlnm_strata_fig))
  )
}

dlnm_lags <- read_if_exists(file.path(pooled_dir, "pooled_dlnm_lag_window_ratio_results.csv"))
if (nrow(dlnm_lags) > 0) {
  lag_fig <- dlnm_lags |>
    dplyr::filter(.data$stratum == "Overall", grepl("primary|mrt", .data$model)) |>
    dplyr::mutate(
      row_label = paste(clean_label(.data$analysis), clean_label(.data$model), .data$lag_window, sep = " | ")
    ) |>
    trim_ratio_plot(lower = 0.1, upper = 10)

  save_forest(
    lag_fig,
    file.path(figure_dir, "pooled_dlnm_lag_window_ratio_forest_overall.png"),
    "Pooled DLNM Cumulative Lag-window Relative Ratios",
    "Overall pooled cumulative RR estimates by lag window. Error bars are 95% CIs.",
    x_label = "Cumulative relative risk",
    width = 11,
    height = max(6, 0.26 * nrow(lag_fig))
  )
}

dlnm_lag_specific <- read_if_exists(file.path(pooled_dir, "pooled_dlnm_lag_specific_ratio_results.csv"))
if (nrow(dlnm_lag_specific) > 0) {
  lag_specific_fig <- dlnm_lag_specific |>
    dplyr::filter(.data$stratum == "Overall", grepl("primary|mrt", .data$model)) |>
    dplyr::mutate(
      row_label = paste(clean_label(.data$analysis), clean_label(.data$model), .data$lag_window, sep = " | ")
    ) |>
    trim_ratio_plot(lower = 0.1, upper = 10)

  save_forest(
    lag_specific_fig,
    file.path(figure_dir, "pooled_dlnm_lag_specific_ratio_forest_overall.png"),
    "Pooled DLNM Lag-specific Relative Ratios",
    "Overall pooled lag-specific RR estimates. Error bars are 95% CIs.",
    x_label = "Lag-specific relative risk",
    width = 11,
    height = max(6, 0.26 * nrow(lag_specific_fig))
  )
}

phenotype <- read_if_exists(file.path(pooled_dir, "pooled_phenotype_ratio_results.csv"))
if (nrow(phenotype) > 0) {
  phenotype_terms <- phenotype |>
    dplyr::filter(.data$term_family %in% c("temperature_spline_basis", "humidity_spline_basis")) |>
    dplyr::mutate(
      basis = dplyr::case_when(
        grepl("1", .data$term) ~ "basis 1",
        grepl("2", .data$term) ~ "basis 2",
        grepl("3", .data$term) ~ "basis 3",
        TRUE ~ .data$term
      ),
      row_label = paste(clean_outcome(.data$outcome), clean_adjustment(.data$analysis), clean_term_family(.data$term_family), .data$basis, sep = " | ")
    ) |>
    trim_ratio_plot(lower = 0.05, upper = 20)

  save_forest(
    phenotype_terms,
    file.path(figure_dir, "pooled_72h_phenotype_spline_ratio_forest.png"),
    "Pooled 72-hour Phenotype Multinomial Ratios",
    "Pooled ORs for spline-basis terms; curves remain the main nonlinear interpretation.",
    x_label = "Odds ratio",
    width = 12,
    height = max(7, 0.25 * nrow(phenotype_terms))
  )
}

fine_gray <- read_if_exists(file.path(pooled_dir, "pooled_fine_gray_ratio_results.csv"))
if (nrow(fine_gray) > 0) {
  fine_gray_terms <- fine_gray |>
    dplyr::filter(.data$term_family %in% c("temperature_spline_basis", "humidity_spline_basis")) |>
    dplyr::mutate(
      basis = dplyr::case_when(
        grepl("1", .data$term) ~ "basis 1",
        grepl("2", .data$term) ~ "basis 2",
        grepl("3", .data$term) ~ "basis 3",
        TRUE ~ .data$term
      ),
      row_label = paste(clean_adjustment(.data$analysis), clean_term_family(.data$term_family), .data$basis, sep = " | ")
    ) |>
    trim_ratio_plot(lower = 0.05, upper = 20)

  save_forest(
    fine_gray_terms,
    file.path(figure_dir, "pooled_fine_gray_spline_ratio_forest.png"),
    "Pooled Fine-Gray Ratios",
    "Pooled SHRs for spline-basis terms; CIFs remain the main competing-risk interpretation.",
    x_label = "Subdistribution hazard ratio",
    width = 11,
    height = max(5, 0.35 * nrow(fine_gray_terms))
  )
}

message("Wrote pooled ratio figures to ", figure_dir)

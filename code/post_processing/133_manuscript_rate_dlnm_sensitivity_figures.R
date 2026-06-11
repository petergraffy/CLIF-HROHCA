#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
export_dir <- file.path(repo_root, "output", "final", "federated_exports")
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(file_name) {
  path <- file.path(pooled_dir, file_name)
  if (!file.exists(path)) stop("Missing required pooled file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

fmt_rr <- function(rr, low, high) sprintf("%.2f (%.2f, %.2f)", rr, low, high)
fmt_i2 <- function(x) ifelse(is.finite(x), paste0(round(x), "%"), "")

read_bind <- function(pattern) {
  files <- list.files(export_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(data.frame())
  dplyr::bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))
}

pool_der_simonian_laird <- function(est, se) {
  keep <- is.finite(est) & is.finite(se) & se > 0
  est <- est[keep]
  se <- se[keep]
  k <- length(est)
  if (k == 0) return(NULL)
  vi <- pmax(se^2, .Machine$double.eps)
  wi <- 1 / vi
  fixed <- sum(wi * est) / sum(wi)
  q <- sum(wi * (est - fixed)^2)
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (k > 1) max(0, (q - (k - 1)) / c_val) else 0
  w_re <- 1 / (vi + tau2)
  pooled <- sum(w_re * est) / sum(w_re)
  pooled_se <- sqrt(1 / sum(w_re))
  data.frame(
    k_sites = k,
    log_rr = pooled,
    log_rr_se = pooled_se,
    rr = exp(pooled),
    rr_low = exp(pooled - 1.96 * pooled_se),
    rr_high = exp(pooled + 1.96 * pooled_se),
    tau2 = tau2,
    q = q,
    i2 = ifelse(k > 1 & q > (k - 1) & q > 0, 100 * (q - (k - 1)) / q, 0),
    stringsAsFactors = FALSE
  )
}

pool_by <- function(df, group_cols, estimate_col = "log_rr", se_col = "log_rr_se") {
  groups <- unique(df[, group_cols, drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    keep <- rep(TRUE, nrow(df))
    for (col in group_cols) keep <- keep & df[[col]] == g[[col]]
    pooled <- pool_der_simonian_laird(df[[estimate_col]][keep], df[[se_col]][keep])
    if (is.null(pooled)) next
    rows[[i]] <- cbind(g, pooled, stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(rows)
}

all_site_curves <- read_required("all_site_dlnm_curves.csv")
pooled_curves <- read_required("pooled_dlnm_random_effects_curves.csv")
pooled_ratios <- read_required("pooled_dlnm_ratio_results.csv")

rate_curve_site <- all_site_curves |>
  dplyr::filter(
    .data$export_family == "rate_dlnm",
    .data$stratum == "Overall",
    .data$model == "rate_median_humidity_adjusted",
    .data$reference_type == "median"
  ) |>
  dplyr::mutate(rr = dplyr::coalesce(.data$cumulative_rate_ratio, .data$cumulative_rr, exp(.data$log_rr)))

rate_curve_pooled <- pooled_curves |>
  dplyr::filter(
    .data$export_family == "rate_dlnm",
    .data$stratum == "Overall",
    .data$model == "rate_median_humidity_adjusted",
    .data$reference_type == "median"
  ) |>
  dplyr::arrange(.data$tmax_mean_c)

if (nrow(rate_curve_site) == 0 || nrow(rate_curve_pooled) == 0) {
  stop("Could not find pooled rate DLNM median-reference humidity-adjusted curves.", call. = FALSE)
}

readr::write_csv(rate_curve_pooled, file.path(figure_dir, "sensitivity_rate_dlnm_curve_source.csv"))

rate_curve_fig <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = 1, color = "grey35", linewidth = 0.45, linetype = "dashed") +
  ggplot2::geom_line(
    data = rate_curve_site,
    ggplot2::aes(x = .data$tmax_mean_c, y = .data$rr, group = .data$site_name),
    color = "grey67",
    linewidth = 0.75,
    alpha = 0.75
  ) +
  ggplot2::geom_ribbon(
    data = rate_curve_pooled,
    ggplot2::aes(x = .data$tmax_mean_c, ymin = .data$rr_low, ymax = .data$rr_high),
    fill = "#8f7a2f",
    alpha = 0.18
  ) +
  ggplot2::geom_line(
    data = rate_curve_pooled,
    ggplot2::aes(x = .data$tmax_mean_c, y = .data$rr),
    color = "#8f7a2f",
    linewidth = 1.35
  ) +
  ggplot2::scale_x_continuous(breaks = seq(-20, 40, by = 10), expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
  ggplot2::scale_y_log10(
    breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4),
    labels = c("0.25", "0.5", "0.75", "1", "1.5", "2", "3", "4")
  ) +
  ggplot2::coord_cartesian(ylim = c(0.25, 4)) +
  ggplot2::labs(
    title = "Rate DLNM sensitivity",
    subtitle = "Median-reference, humidity-adjusted model; pooled curve with 95% CI and site curves in grey",
    x = "Daily maximum temperature (C)",
    y = "Cumulative relative rate"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "grey30", size = 10.5),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "grey15"),
    plot.margin = ggplot2::margin(12, 16, 12, 12)
  )

ggplot2::ggsave(file.path(figure_dir, "sensitivity_rate_dlnm_curve.png"), rate_curve_fig, width = 7.4, height = 5.2, dpi = 600)
ggplot2::ggsave(file.path(figure_dir, "sensitivity_rate_dlnm_curve.pdf"), rate_curve_fig, width = 7.4, height = 5.2)

overall <- pooled_ratios |>
  dplyr::filter(
    .data$analysis == "rate_dlnm",
    .data$stratum == "Overall",
    .data$model == "rate_median_humidity_adjusted",
    .data$contrast == "median-referenced cumulative effect"
  )
strata <- pooled_ratios |>
  dplyr::filter(
    .data$analysis == "rate_dlnm",
    .data$stratum != "Overall",
    .data$model == "rate_stratified_humidity_adjusted",
    .data$contrast == "median-referenced cumulative effect"
  )

forest <- dplyr::bind_rows(overall, strata) |>
  dplyr::mutate(
    group = dplyr::case_when(
      .data$stratum == "Overall" ~ "Overall",
      .data$stratum %in% c("Male", "Female") ~ "Sex",
      .data$stratum %in% c("<65", ">=65") ~ "Age",
      .data$stratum %in% c("Black", "Non-Black") ~ "Race",
      .data$stratum %in% c("Hispanic", "Non-Hispanic") ~ "Ethnicity",
      TRUE ~ "Other"
    ),
    label = dplyr::case_when(
      .data$stratum == "<65" ~ "<65 years",
      .data$stratum == ">=65" ~ ">=65 years",
      TRUE ~ .data$stratum
    ),
    row_order = dplyr::case_when(
      .data$stratum == "Overall" ~ 1,
      .data$stratum == "Male" ~ 3,
      .data$stratum == "Female" ~ 4,
      .data$stratum == "<65" ~ 6,
      .data$stratum == ">=65" ~ 7,
      .data$stratum == "Black" ~ 9,
      .data$stratum == "Non-Black" ~ 10,
      .data$stratum == "Hispanic" ~ 12,
      .data$stratum == "Non-Hispanic" ~ 13,
      TRUE ~ 99
    ),
    rr_label = fmt_rr(.data$ratio, .data$ratio_low, .data$ratio_high),
    i2_label = fmt_i2(.data$i2),
    group = factor(.data$group, levels = c("Overall", "Sex", "Age", "Race", "Ethnicity"))
  ) |>
  dplyr::arrange(.data$row_order)

forest <- forest |>
  dplyr::mutate(
    y = dplyr::case_when(
      .data$stratum == "Overall" ~ 11,
      .data$stratum == "Male" ~ 9,
      .data$stratum == "Female" ~ 8,
      .data$stratum == "<65" ~ 6,
      .data$stratum == ">=65" ~ 5,
      .data$stratum == "Black" ~ 3,
      .data$stratum == "Non-Black" ~ 2,
      .data$stratum == "Hispanic" ~ 0,
      .data$stratum == "Non-Hispanic" ~ -1,
      TRUE ~ -3
    )
  )
readr::write_csv(forest, file.path(figure_dir, "sensitivity_rate_dlnm_stratified_forest_source.csv"))

group_headers <- data.frame(
  y = c(8.5, 5.5, 2.5, -0.5),
  label = c("Sex", "Age", "Race", "Ethnicity"),
  stringsAsFactors = FALSE
)
x_min <- min(0.55, min(forest$ratio_low, na.rm = TRUE) * 0.85)
x_text <- max(forest$ratio_high, na.rm = TRUE) * 1.25
x_i2 <- max(forest$ratio_high, na.rm = TRUE) * 2.35
x_max <- max(4, max(forest$ratio_high, na.rm = TRUE) * 3.6)
colors <- c("Overall" = "#8f7a2f", "Sex" = "#335c67", "Age" = "#8f7a2f", "Race" = "#9b2f37", "Ethnicity" = "#457B9D")

rate_forest_fig <- ggplot2::ggplot(forest, ggplot2::aes(x = .data$ratio, y = .data$y, color = .data$group)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey35") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$ratio_low, xmax = .data$ratio_high),
    orientation = "y",
    width = 0.18,
    linewidth = 0.85,
    alpha = 0.95
  ) +
  ggplot2::geom_point(size = 3.1) +
  ggplot2::geom_text(ggplot2::aes(x = x_text, label = .data$rr_label), hjust = 0, color = "grey12", size = 3.25, show.legend = FALSE) +
  ggplot2::geom_text(ggplot2::aes(x = x_i2, label = .data$i2_label), hjust = 0, color = "grey25", size = 3.25, show.legend = FALSE) +
  ggplot2::annotate("text", x = x_text, y = 11.65, label = "RR (95% CI)", hjust = 0, fontface = "bold", size = 3.35, color = "grey12") +
  ggplot2::annotate("text", x = x_i2, y = 11.65, label = "I2", hjust = 0, fontface = "bold", size = 3.35, color = "grey12") +
  ggplot2::geom_text(data = group_headers, ggplot2::aes(x = x_min, y = .data$y, label = .data$label), inherit.aes = FALSE, hjust = 0, fontface = "bold", color = "grey12", size = 3.35) +
  ggplot2::scale_color_manual(values = colors, guide = "none") +
  ggplot2::scale_y_continuous(breaks = forest$y, labels = forest$label, limits = c(-1.6, 12), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_x_log10(
    breaks = c(0.75, 1, 1.5, 2, 3),
    labels = c("0.75", "1", "1.5", "2", "3"),
    limits = c(x_min, x_max),
    expand = ggplot2::expansion(mult = c(0.015, 0.02))
  ) +
  ggplot2::labs(
    title = "Rate DLNM stratified sensitivity",
    subtitle = "Median-reference and humidity-adjusted rate models",
    x = "Cumulative relative rate",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "grey30", size = 10.5),
    axis.title.x = ggplot2::element_text(face = "bold"),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(12, 22, 12, 12)
  ) +
  ggplot2::coord_cartesian(clip = "off")

ggplot2::ggsave(file.path(figure_dir, "sensitivity_rate_dlnm_stratified_forest.png"), rate_forest_fig, width = 8.2, height = 5.1, dpi = 600)
ggplot2::ggsave(file.path(figure_dir, "sensitivity_rate_dlnm_stratified_forest.pdf"), rate_forest_fig, width = 8.2, height = 5.1)

rate_lag_specific_site <- read_bind("^[^_]+_rate_lag30_hot_temperature_lag_specific_rr[.]csv$")
rate_cumulative_site <- read_bind("^[^_]+_rate_lag30_hot_temperature_cumulative_rr_by_lag[.]csv$")

if (nrow(rate_lag_specific_site) == 0 || nrow(rate_cumulative_site) == 0) {
  stop("Missing rate lag-30 diagnostic federated exports. Run sites through 04e_rate_lag30_diagnostic_plots.R and 07_export_federated_results.R first.", call. = FALSE)
}

rate_lag_specific <- pool_by(rate_lag_specific_site, c("analysis_period", "lag")) |>
  dplyr::arrange(.data$lag)

rate_lag_windows_all <- pool_by(rate_cumulative_site, c("analysis_period", "lag_start", "lag_end")) |>
  dplyr::arrange(.data$lag_end) |>
  dplyr::mutate(lag_label = paste0("Lag 0-", .data$lag_end))
rate_lag_windows_all$lag_label[rate_lag_windows_all$lag_end == 0] <- "Lag 0"

selected_lag_ends <- c(0, 1, 2, 3, 4, 5, 7, 10, 14, 21, 30)
rate_lag_windows <- rate_lag_windows_all |>
  dplyr::filter(.data$lag_end %in% selected_lag_ends) |>
  dplyr::mutate(
    lag_label = factor(.data$lag_label, levels = rev(.data$lag_label[order(.data$lag_end)])),
    rr_label = sprintf("%.2f (%.2f, %.2f); I2 %.0f%%", .data$rr, .data$rr_low, .data$rr_high, .data$i2)
  )

readr::write_csv(rate_lag_specific, file.path(figure_dir, "sensitivity_rate_dlnm_lag_specific_source.csv"))
readr::write_csv(rate_lag_windows_all, file.path(figure_dir, "sensitivity_rate_dlnm_lag_window_all_windows_source.csv"))
readr::write_csv(rate_lag_windows, file.path(figure_dir, "sensitivity_rate_dlnm_lag_window_forest_source.csv"))

lag_a <- ggplot2::ggplot(rate_lag_specific, ggplot2::aes(x = .data$lag, y = .data$rr)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$rr_low, ymax = .data$rr_high), fill = "#8f7a2f", alpha = 0.22) +
  ggplot2::geom_line(color = "#8f7a2f", linewidth = 1.15) +
  ggplot2::geom_point(color = "#8f7a2f", size = 2) +
  ggplot2::scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +
  ggplot2::scale_y_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.25, 1.5), labels = c("0.8", "0.9", "1.0", "1.1", "1.25", "1.5")) +
  ggplot2::labs(title = "A. Rate DLNM lag-specific associations through 30 days", x = "Lag since exposure day", y = "Lag-specific relative rate") +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5), axis.title = ggplot2::element_text(face = "bold"))

x_text_lag <- max(rate_lag_windows$rr_high, na.rm = TRUE) * 1.18
x_max_lag <- max(rate_lag_windows$rr_high, na.rm = TRUE) * 2.7
lag_b <- ggplot2::ggplot(rate_lag_windows, ggplot2::aes(x = .data$rr, y = .data$lag_label)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$rr_low, xmax = .data$rr_high), orientation = "y", width = 0.18, color = "#8f7a2f", linewidth = 0.85) +
  ggplot2::geom_point(color = "#8f7a2f", size = 2.8) +
  ggplot2::geom_text(ggplot2::aes(x = x_text_lag, label = .data$rr_label), hjust = 0, size = 3.0, color = "grey12") +
  ggplot2::annotate("text", x = x_text_lag, y = length(selected_lag_ends) + 0.65, label = "RR (95% CI); I2", hjust = 0, fontface = "bold", size = 3.05) +
  ggplot2::scale_x_log10(breaks = c(0.75, 1, 1.25, 1.5, 2), labels = c("0.75", "1", "1.25", "1.5", "2"), limits = c(0.7, x_max_lag)) +
  ggplot2::labs(title = "B. Rate DLNM cumulative lag-window relative rates", x = "Cumulative relative rate", y = NULL) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(8, 36, 8, 8)
  ) +
  ggplot2::coord_cartesian(clip = "off")

rate_lag_fig <- lag_a / lag_b + patchwork::plot_layout(heights = c(0.48, 0.52))
ggplot2::ggsave(file.path(figure_dir, "sensitivity_rate_dlnm_lag_specific_and_window_forest.png"), rate_lag_fig, width = 8.2, height = 8.2, dpi = 600)
ggplot2::ggsave(file.path(figure_dir, "sensitivity_rate_dlnm_lag_specific_and_window_forest.pdf"), rate_lag_fig, width = 8.2, height = 8.2)

message("Wrote rate DLNM sensitivity figures to ", figure_dir)

#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
export_dir <- file.path(repo_root, "output", "final", "federated_exports")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_bind <- function(pattern) {
  files <- list.files(export_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(tibble::tibble())
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
  tibble::tibble(
    k_sites = k,
    log_rr = pooled,
    log_rr_se = pooled_se,
    rr = exp(pooled),
    rr_low = exp(pooled - 1.96 * pooled_se),
    rr_high = exp(pooled + 1.96 * pooled_se),
    tau2 = tau2,
    q = q,
    i2 = ifelse(k > 1 & q > (k - 1) & q > 0, 100 * (q - (k - 1)) / q, 0)
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

lag_specific_site <- read_bind("^[^_]+_case_crossover_lag30_hot_temperature_lag_specific_rr[.]csv$")
cumulative_site <- read_bind("^[^_]+_case_crossover_lag30_hot_temperature_cumulative_rr_by_lag[.]csv$")

if (nrow(lag_specific_site) == 0 || nrow(cumulative_site) == 0) {
  stop(
    "Missing case-crossover lag-30 federated exports. Sites need to run code/04f_case_crossover_lag30_diagnostic_plots.R and code/07_export_federated_results.R.",
    call. = FALSE
  )
}

lag_specific <- pool_by(lag_specific_site, c("analysis_period", "lag")) |>
  dplyr::arrange(.data$lag)

lag_windows_all <- pool_by(cumulative_site, c("analysis_period", "lag_start", "lag_end")) |>
  dplyr::arrange(.data$lag_end) |>
  dplyr::mutate(lag_label = paste0("Lag 0-", .data$lag_end))
lag_windows_all$lag_label[lag_windows_all$lag_end == 0] <- "Lag 0"

selected_lag_ends <- c(0, 1, 2, 3, 4, 5, 7, 10, 14, 21, 30)
lag_windows <- lag_windows_all |>
  dplyr::filter(.data$lag_end %in% selected_lag_ends) |>
  dplyr::mutate(
    lag_label = factor(.data$lag_label, levels = rev(.data$lag_label[order(.data$lag_end)])),
    rr_label = sprintf("%.2f (%.2f, %.2f); I2 %.0f%%", .data$rr, .data$rr_low, .data$rr_high, .data$i2)
  )

readr::write_csv(lag_specific, file.path(figure_dir, "figure_case_crossover_lag30_panel_a_source.csv"))
readr::write_csv(lag_windows_all, file.path(figure_dir, "figure_case_crossover_lag30_all_windows_source.csv"))
readr::write_csv(lag_windows, file.path(figure_dir, "figure_case_crossover_lag30_panel_b_source.csv"))

panel_a <- ggplot2::ggplot(lag_specific, ggplot2::aes(x = .data$lag, y = .data$rr)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data$rr_low, ymax = .data$rr_high),
    fill = "#0f6b78",
    alpha = 0.18
  ) +
  ggplot2::geom_line(color = "#0f6b78", linewidth = 1.15) +
  ggplot2::geom_point(color = "#0f6b78", size = 2) +
  ggplot2::scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +
  ggplot2::scale_y_log10(
    breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 1.75, 2),
    labels = c("0.7", "0.8", "0.9", "1.0", "1.1", "1.25", "1.5", "1.75", "2.0")
  ) +
  ggplot2::labs(
    title = "A. Case-crossover DLNM lag-specific associations through 30 days",
    x = "Lag since exposure day",
    y = "Lag-specific relative risk"
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    axis.title = ggplot2::element_text(face = "bold")
  )

x_text <- max(lag_windows$rr_high, na.rm = TRUE) * 1.18
x_max <- max(lag_windows$rr_high, na.rm = TRUE) * 2.7

panel_b <- ggplot2::ggplot(lag_windows, ggplot2::aes(x = .data$rr, y = .data$lag_label)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$rr_low, xmax = .data$rr_high),
    orientation = "y",
    width = 0.18,
    color = "#0f6b78",
    linewidth = 0.85
  ) +
  ggplot2::geom_point(color = "#0f6b78", size = 2.8) +
  ggplot2::geom_text(ggplot2::aes(x = x_text, label = .data$rr_label), hjust = 0, size = 3.0, color = "grey12") +
  ggplot2::annotate("text", x = x_text, y = length(selected_lag_ends) + 0.65, label = "RR (95% CI); I2", hjust = 0, fontface = "bold", size = 3.05) +
  ggplot2::scale_x_log10(
    breaks = c(0.75, 1, 1.25, 1.5, 2, 3),
    labels = c("0.75", "1", "1.25", "1.5", "2", "3")
  ) +
  ggplot2::labs(
    title = "B. Case-crossover DLNM cumulative lag-window relative risks",
    x = "Cumulative relative risk",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(8, 36, 8, 8)
  ) +
  ggplot2::coord_cartesian(xlim = c(0.65, x_max), clip = "off")

combined <- panel_a / panel_b + patchwork::plot_layout(heights = c(0.48, 0.52))

png_path <- file.path(figure_dir, "figure_case_crossover_lag30_relative_risk_ab.png")
pdf_path <- file.path(figure_dir, "figure_case_crossover_lag30_relative_risk_ab.pdf")
ggplot2::ggsave(png_path, combined, width = 8.2, height = 8.2, dpi = 600)
ggplot2::ggsave(pdf_path, combined, width = 8.2, height = 8.2)

message("Wrote case-crossover lag-30 A/B figure to ", png_path, " and ", pdf_path)

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

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)

suppressPackageStartupMessages({
  library(dlnm)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(splines)
})

ohca_counts_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_daily_counts_2018_2024.csv")
daily_exposure_path <- file.path(repo_root, "output", "intermediate", "cohorts", "all_icu", "all_icu_daily_patient_address_tmax_2018_2024.csv")
output_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "lag_diagnostics")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

START_DATE <- as.Date("2018-01-01")
END_DATE <- as.Date("2024-12-31")
WARM_MONTHS <- c(5L, 6L, 7L, 8L, 9L)
ANALYSIS_SEASON <- tolower(trimws(Sys.getenv("OHCA_DLNM_SEASON", unset = "all_year")))
if (!ANALYSIS_SEASON %in% c("all_year", "warm")) {
  stop("OHCA_DLNM_SEASON must be either 'all_year' or 'warm'.", call. = FALSE)
}
ANALYSIS_PERIOD_LABEL <- if (ANALYSIS_SEASON == "warm") "warm_season" else "all_year"
MAX_DIAGNOSTIC_LAG <- 30L
TIME_DF_PER_YEAR <- 4L
VAR_DF <- 4L
LAG_DF <- 5L
RMAX_DF <- 3L

safe_quantile <- function(x, prob) {
  stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7)
}

make_prediction_grid <- function(x) {
  unique(as.numeric(seq(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)), by = 0.5)))
}

finite_unique_n <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  length(unique(x[is.finite(x)]))
}

se_from_ci <- function(rr, low, high) {
  (log(high) - log(low)) / (2 * 1.96)
}

ohca_daily <- read.csv(ohca_counts_path, stringsAsFactors = FALSE)
ohca_daily$admission_date <- as.Date(ohca_daily$admission_date)
all_dates <- data.frame(admission_date = seq(START_DATE, END_DATE, by = "day"))
daily_counts <- merge(all_dates, ohca_daily, by = "admission_date", all.x = TRUE, sort = TRUE)
daily_counts$ohca_admissions[is.na(daily_counts$ohca_admissions)] <- 0L

daily_exposure <- read.csv(daily_exposure_path, stringsAsFactors = FALSE)
daily_exposure$admission_date <- as.Date(daily_exposure$admission_date)

model_df <- merge(daily_counts, daily_exposure, by = "admission_date", all.x = TRUE, sort = TRUE)
model_df$year <- as.integer(format(model_df$admission_date, "%Y"))
model_df$month <- as.integer(format(model_df$admission_date, "%m"))
model_df$dow <- factor(weekdays(model_df$admission_date), levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))
model_df$time_index <- seq_len(nrow(model_df))
if (ANALYSIS_SEASON == "warm") {
  model_df <- model_df[model_df$month %in% WARM_MONTHS, , drop = FALSE]
  model_df$time_index <- seq_len(nrow(model_df))
}

needed <- c(
  "icu_patient_address_mean_tmax_c",
  "icu_patient_address_mean_rmax_pct",
  "icu_patient_address_mean_no2",
  "icu_patient_address_mean_pm25"
)
model_df <- model_df[complete.cases(model_df[, needed, drop = FALSE]), , drop = FALSE]
for (needed_col in needed) {
  model_df <- model_df[is.finite(suppressWarnings(as.numeric(model_df[[needed_col]]))), , drop = FALSE]
}

temp_unique <- finite_unique_n(model_df$icu_patient_address_mean_tmax_c)
temp_var_df <- min(VAR_DF, max(2L, temp_unique - 1L))
cb_temp <- crossbasis(
  model_df$icu_patient_address_mean_tmax_c,
  lag = MAX_DIAGNOSTIC_LAG,
  argvar = list(fun = "ns", df = temp_var_df),
  arglag = list(fun = "ns", df = LAG_DF)
)

time_df <- min(length(unique(model_df$year)) * TIME_DF_PER_YEAR, max(1L, nrow(model_df) - 2L))
formula_spec <- ohca_admissions ~ cb_temp +
  ns(time_index, df = time_df) +
  dow +
  factor(year) +
  ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF) +
  icu_patient_address_mean_no2 +
  icu_patient_address_mean_pm25

fit <- glm(formula_spec, data = model_df, family = quasipoisson(link = "log"))

grid <- make_prediction_grid(model_df$icu_patient_address_mean_tmax_c)
reference_temp <- median(model_df$icu_patient_address_mean_tmax_c, na.rm = TRUE)
hot_temp <- grid[which.min(abs(grid - safe_quantile(model_df$icu_patient_address_mean_tmax_c, 0.95)))]
pred <- crosspred(cb_temp, fit, cen = reference_temp, at = grid, cumul = TRUE)
hot_index <- which.min(abs(grid - hot_temp))

lag_values <- seq_len(ncol(pred$matRRfit)) - 1L
surface <- expand.grid(
  tmax_mean_c = grid,
  lag = lag_values
) |>
  mutate(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    reference_temp_c = reference_temp,
    hot_temp_c = hot_temp,
    rr = as.vector(pred$matRRfit),
    rr_low = as.vector(pred$matRRlow),
    rr_high = as.vector(pred$matRRhigh),
    log_rr = log(.data$rr),
    log_rr_se = se_from_ci(.data$rr, .data$rr_low, .data$rr_high)
  )

lag_specific_hot <- tibble::tibble(
  analysis_period = ANALYSIS_PERIOD_LABEL,
  reference_temp_c = reference_temp,
  hot_temp_c = hot_temp,
  lag = lag_values,
  rr = as.numeric(pred$matRRfit[hot_index, ]),
  rr_low = as.numeric(pred$matRRlow[hot_index, ]),
  rr_high = as.numeric(pred$matRRhigh[hot_index, ])
) |>
  mutate(
    log_rr = log(.data$rr),
    log_rr_se = se_from_ci(.data$rr, .data$rr_low, .data$rr_high)
  )

lag_ends <- seq_len(ncol(pred$cumRRfit)) - 1L
cumulative_hot <- tibble::tibble(
  analysis_period = ANALYSIS_PERIOD_LABEL,
  reference_temp_c = reference_temp,
  hot_temp_c = hot_temp,
  lag_start = 0L,
  lag_end = lag_ends,
  cumulative_rr = as.numeric(pred$cumRRfit[hot_index, ]),
  cumulative_rr_low = as.numeric(pred$cumRRlow[hot_index, ]),
  cumulative_rr_high = as.numeric(pred$cumRRhigh[hot_index, ])
) |>
  mutate(
    log_rr = log(.data$cumulative_rr),
    log_rr_se = se_from_ci(.data$cumulative_rr, .data$cumulative_rr_low, .data$cumulative_rr_high)
  )

temperature_distribution <- model_df |>
  transmute(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    admission_date = .data$admission_date,
    ohca_admissions = .data$ohca_admissions,
    tmax_mean_c = .data$icu_patient_address_mean_tmax_c
  )

diagnostic_summary <- tibble::tibble(
  analysis_period = ANALYSIS_PERIOD_LABEL,
  n_days = nrow(model_df),
  n_ohca = sum(model_df$ohca_admissions, na.rm = TRUE),
  max_lag_days = MAX_DIAGNOSTIC_LAG,
  reference_temp_c = reference_temp,
  hot_temp_c = hot_temp,
  hot_lag0_rr = lag_specific_hot$rr[lag_specific_hot$lag == 0],
  hot_lag5_rr = lag_specific_hot$rr[lag_specific_hot$lag == 5],
  hot_lag30_rr = lag_specific_hot$rr[lag_specific_hot$lag == 30],
  cumulative_lag0_5_rr = cumulative_hot$cumulative_rr[cumulative_hot$lag_end == 5],
  cumulative_lag0_30_rr = cumulative_hot$cumulative_rr[cumulative_hot$lag_end == 30],
  dispersion = summary(fit)$dispersion,
  converged = fit$converged
)

readr::write_csv(surface, file.path(output_dir, "lag30_temperature_lag_rr_surface.csv"))
readr::write_csv(lag_specific_hot, file.path(output_dir, "lag30_hot_temperature_lag_specific_rr.csv"))
readr::write_csv(cumulative_hot, file.path(output_dir, "lag30_hot_temperature_cumulative_rr_by_lag.csv"))
readr::write_csv(temperature_distribution, file.path(output_dir, "lag30_temperature_distribution_source.csv"))
readr::write_csv(diagnostic_summary, file.path(output_dir, "lag30_diagnostic_summary.csv"))

figure_surface <- ggplot(surface, aes(x = .data$tmax_mean_c, y = .data$lag, fill = .data$rr)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(
    data = surface,
    aes(x = .data$tmax_mean_c, y = .data$lag, z = .data$rr),
    inherit.aes = FALSE,
    breaks = c(0.75, 1, 1.25, 1.5, 2),
    color = "white",
    linewidth = 0.25,
    alpha = 0.75
  ) +
  geom_vline(xintercept = reference_temp, linetype = "dotted", color = "#1F2933", linewidth = 0.45) +
  geom_vline(xintercept = hot_temp, linetype = "dashed", color = "#7F1D1D", linewidth = 0.45) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "#1F2933", linewidth = 0.45) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 1,
    limits = c(0.5, 2.5),
    oob = scales::squish,
    name = "RR"
  ) +
  scale_y_continuous(breaks = seq(0, 30, by = 5), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Thirty-day temperature-lag DLNM diagnostic surface",
    subtitle = "Overall OHCA model; dashed horizontal line marks the primary 0-5 day lag window",
    x = "Daily assigned-county Tmax (C)",
    y = "Lag since exposure day"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

figure_lag <- ggplot(lag_specific_hot, aes(x = .data$lag, y = .data$rr)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#7F1D1D") +
  geom_ribbon(aes(ymin = .data$rr_low, ymax = .data$rr_high), fill = "#F4A261", alpha = 0.30) +
  geom_line(color = "#9A3412", linewidth = 1.0) +
  geom_point(size = 1.4, color = "#9A3412") +
  labs(
    title = "Lag-specific RR at hot temperature through 30 days",
    subtitle = sprintf("Hot Tmax %.1f C versus median reference %.1f C", hot_temp, reference_temp),
    x = "Lag since exposure day",
    y = "Lag-specific relative risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

figure_cumulative <- ggplot(cumulative_hot, aes(x = .data$lag_end, y = .data$cumulative_rr)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "#7F1D1D") +
  geom_ribbon(aes(ymin = .data$cumulative_rr_low, ymax = .data$cumulative_rr_high), fill = "#8AB6D6", alpha = 0.30) +
  geom_line(color = "#1F4E79", linewidth = 1.0) +
  geom_point(size = 1.4, color = "#1F4E79") +
  labs(
    title = "Cumulative hot-temperature RR by maximum lag",
    subtitle = sprintf("Hot Tmax %.1f C versus median reference %.1f C", hot_temp, reference_temp),
    x = "Cumulative lag window end",
    y = "Cumulative relative risk"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

figure_temperature <- ggplot(temperature_distribution, aes(x = .data$tmax_mean_c)) +
  geom_histogram(binwidth = 2, fill = "#4B5563", color = "white", boundary = 0) +
  geom_vline(xintercept = reference_temp, linetype = "dotted", color = "#1F2933", linewidth = 0.7) +
  geom_vline(xintercept = hot_temp, linetype = "dashed", color = "#7F1D1D", linewidth = 0.7) +
  labs(
    title = "Temperature distribution used in 30-day lag diagnostic",
    subtitle = "Dotted line is median reference; dashed line is 95th percentile hot contrast",
    x = "Daily assigned-county Tmax (C)",
    y = "Days"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

save_plot <- function(plot, basename, width = 8, height = 5.5) {
  ggsave(file.path(figure_dir, paste0(basename, ".png")), plot, width = width, height = height, dpi = 300)
  ggsave(file.path(figure_dir, paste0(basename, ".pdf")), plot, width = width, height = height)
}

save_plot(figure_surface, "figure_lag30_temperature_lag_rr_surface", width = 8, height = 5.8)
save_plot(figure_lag, "figure_lag30_hot_temperature_lag_specific_rr", width = 7.5, height = 5.0)
save_plot(figure_cumulative, "figure_lag30_hot_temperature_cumulative_rr_by_lag", width = 7.5, height = 5.0)
save_plot(figure_temperature, "figure_lag30_temperature_distribution", width = 7.5, height = 4.5)

print(diagnostic_summary)
message("Wrote 30-day lag diagnostic plots and source files to ", output_dir, " and ", figure_dir)

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
  library(readr)
  library(splines)
  library(survival)
  library(tidyr)
})

ohca_cohort_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
output_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "case_crossover")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

START_DATE <- as.Date("2018-01-01")
END_DATE <- as.Date("2024-12-31")
WARM_MONTHS <- c(5L, 6L, 7L, 8L, 9L)
ANALYSIS_SEASON <- tolower(trimws(Sys.getenv("OHCA_DLNM_SEASON", unset = "all_year")))
if (!ANALYSIS_SEASON %in% c("all_year", "warm")) {
  stop("OHCA_DLNM_SEASON must be either 'all_year' or 'warm'.", call. = FALSE)
}
ANALYSIS_PERIOD_LABEL <- if (ANALYSIS_SEASON == "warm") "warm_season" else "all_year"
MAX_DIAGNOSTIC_LAG <- 30L
VAR_DF <- 4L
LAG_DF <- 5L
RMAX_DF <- 3L

safe_quantile <- function(x, prob) stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7)

finite_unique_n <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  length(unique(x[is.finite(x)]))
}

se_from_ci <- function(rr, low, high) {
  (log(high) - log(low)) / (2 * 1.96)
}

filter_analysis_period_cases <- function(df) {
  if (ANALYSIS_SEASON == "warm") df <- df[df$month %in% WARM_MONTHS, , drop = FALSE]
  df
}

make_referent_dates <- function(cases, all_dates) {
  case_list <- split(cases, cases$hospitalization_id)
  bind_rows(lapply(case_list, function(case_row) {
    case_row <- case_row[1, , drop = FALSE]
    referent <- all_dates |>
      filter(
        .data$year == case_row$year[[1]],
        .data$month == case_row$month[[1]],
        .data$dow == case_row$dow[[1]]
      )
    tibble::tibble(
      match_id = case_row$hospitalization_id[[1]],
      hospitalization_id = case_row$hospitalization_id[[1]],
      patient_id = case_row$patient_id[[1]],
      admission_date = case_row$admission_date[[1]],
      referent_date = referent$referent_date,
      case = as.integer(referent$referent_date == case_row$admission_date[[1]]),
      county_fips = case_row$county_fips[[1]]
    )
  }))
}

add_lagged_daily_exposures <- function(df, daily_exposure, value_col, prefix) {
  out <- df
  for (lag in 0:MAX_DIAGNOSTIC_LAG) {
    join_df <- daily_exposure |>
      transmute(
        county_fips = .data$county_fips,
        referent_date = .data$admission_date + lag,
        value = .data[[value_col]]
      )
    names(join_df)[names(join_df) == "value"] <- paste0(prefix, "_lag", lag)
    out <- out |>
      left_join(join_df, by = c("county_fips", "referent_date"))
  }
  out
}

ohca <- readr::read_csv(
  ohca_cohort_path,
  show_col_types = FALSE,
  col_types = readr::cols(
    .default = readr::col_guess(),
    hospitalization_id = readr::col_character(),
    patient_id = readr::col_character(),
    county_fips = readr::col_character()
  )
) |>
  mutate(
    admission_date = as.Date(.data$admission_date),
    county_fips = normalize_county_fips(.data$county_fips),
    year = as.integer(format(.data$admission_date, "%Y")),
    month = as.integer(format(.data$admission_date, "%m")),
    dow = weekdays(.data$admission_date)
  ) |>
  filter(!is.na(.data$admission_date), !is.na(.data$county_fips)) |>
  filter_analysis_period_cases()

tmax <- read_exposome_daymet(repo_root, "daymet_county_tmax_2018_2024_conus.parquet", "tmax_mean_c")
rmax <- read_exposome_daymet(repo_root, "daymet_county_rmax_2018_2024.parquet", "rmax_mean_pct")

daily_exposure <- full_join(tmax, rmax, by = c("admission_date", "county_fips")) |>
  filter(.data$admission_date >= START_DATE - MAX_DIAGNOSTIC_LAG, .data$admission_date <= END_DATE)

all_dates <- tibble::tibble(
  referent_date = seq(START_DATE, END_DATE, by = "day"),
  year = as.integer(format(.data$referent_date, "%Y")),
  month = as.integer(format(.data$referent_date, "%m")),
  dow = weekdays(.data$referent_date)
)

referents <- make_referent_dates(ohca, all_dates)
referents <- add_lagged_daily_exposures(referents, daily_exposure, "tmax_mean_c", "tmax_mean_c")
referents <- add_lagged_daily_exposures(referents, daily_exposure, "rmax_mean_pct", "rmax_mean_pct")

tmax_cols <- paste0("tmax_mean_c_lag", 0:MAX_DIAGNOSTIC_LAG)
needed <- c(tmax_cols, "rmax_mean_pct_lag0")
model_df <- referents[complete.cases(referents[, needed, drop = FALSE]), , drop = FALSE]
for (needed_col in needed) {
  model_df <- model_df[is.finite(suppressWarnings(as.numeric(model_df[[needed_col]]))), , drop = FALSE]
}
model_df <- model_df |>
  group_by(.data$match_id) |>
  filter(sum(.data$case == 1L, na.rm = TRUE) == 1L, n() >= 2L) |>
  ungroup()

if (sum(model_df$case == 1L, na.rm = TRUE) == 0L) {
  stop("No complete case-crossover referent sets available for lag-30 diagnostic.", call. = FALSE)
}

temp_unique <- finite_unique_n(model_df$tmax_mean_c_lag0)
temp_var_df <- min(VAR_DF, max(2L, temp_unique - 1L))
temp_history <- as.matrix(model_df[, tmax_cols, drop = FALSE])
cb_temp <- crossbasis(
  temp_history,
  lag = c(0L, MAX_DIAGNOSTIC_LAG),
  argvar = list(fun = "ns", df = temp_var_df),
  arglag = list(fun = "ns", df = LAG_DF)
)

rmax_term <- if (finite_unique_n(model_df$rmax_mean_pct_lag0) >= RMAX_DF + 1L) {
  "ns(rmax_mean_pct_lag0, df = RMAX_DF)"
} else if (finite_unique_n(model_df$rmax_mean_pct_lag0) >= 2L) {
  "rmax_mean_pct_lag0"
} else {
  NA_character_
}
rhs <- c("cb_temp", rmax_term, "strata(match_id)")
formula_spec <- as.formula(paste("case ~", paste(rhs[!is.na(rhs)], collapse = " + ")))
environment(formula_spec) <- environment()

fit <- survival::clogit(formula_spec, data = model_df, method = "efron")

grid <- unique(as.numeric(seq(
  floor(min(model_df$tmax_mean_c_lag0, na.rm = TRUE)),
  ceiling(max(model_df$tmax_mean_c_lag0, na.rm = TRUE)),
  by = 0.5
)))
reference_temp <- median(model_df$tmax_mean_c_lag0, na.rm = TRUE)
hot_temp <- grid[which.min(abs(grid - safe_quantile(model_df$tmax_mean_c_lag0, 0.95)))]
pred <- crosspred(cb_temp, fit, cen = reference_temp, at = grid, cumul = TRUE)
hot_index <- which.min(abs(grid - hot_temp))

lag_values <- seq_len(ncol(pred$matRRfit)) - 1L
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

diagnostic_summary <- tibble::tibble(
  analysis_period = ANALYSIS_PERIOD_LABEL,
  n_referent_rows = nrow(model_df),
  n_case_events = sum(model_df$case == 1L, na.rm = TRUE),
  n_match_sets = dplyr::n_distinct(model_df$match_id),
  max_lag_days = MAX_DIAGNOSTIC_LAG,
  reference_temp_c = reference_temp,
  hot_temp_c = hot_temp,
  hot_lag0_rr = lag_specific_hot$rr[lag_specific_hot$lag == 0],
  hot_lag5_rr = lag_specific_hot$rr[lag_specific_hot$lag == 5],
  hot_lag30_rr = lag_specific_hot$rr[lag_specific_hot$lag == 30],
  cumulative_lag0_5_rr = cumulative_hot$cumulative_rr[cumulative_hot$lag_end == 5],
  cumulative_lag0_30_rr = cumulative_hot$cumulative_rr[cumulative_hot$lag_end == 30],
  model_family = "conditional_logistic",
  converged = all(is.finite(stats::coef(fit)))
)

readr::write_csv(diagnostic_summary, file.path(output_dir, "case_crossover_lag30_diagnostic_summary.csv"))
readr::write_csv(lag_specific_hot, file.path(output_dir, "case_crossover_lag30_hot_temperature_lag_specific_rr.csv"))
readr::write_csv(cumulative_hot, file.path(output_dir, "case_crossover_lag30_hot_temperature_cumulative_rr_by_lag.csv"))

print(diagnostic_summary)
message("Wrote time-stratified case-crossover lag-30 diagnostic exports to ", output_dir)

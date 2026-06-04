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
  library(scales)
  library(splines)
})

ohca_counts_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_daily_counts_2018_2024.csv")
ohca_cohort_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
daily_exposure_path <- file.path(repo_root, "output", "intermediate", "cohorts", "all_icu", "all_icu_daily_patient_address_tmax_2018_2024.csv")
output_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "icu_admission_rate")
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
MAX_LAG <- 5L
TIME_DF_PER_YEAR <- 4L
VAR_DF <- 4L
LAG_DF <- 3L
RMAX_DF <- 3L

safe_quantile <- function(x, prob) stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7)

make_prediction_grid <- function(x) {
  unique(as.numeric(seq(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)), by = 0.5)))
}

finite_unique_n <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  length(unique(x[is.finite(x)]))
}

choose_spline_term <- function(df, variable, spline_term, requested_df, linear_term = variable) {
  n_unique <- finite_unique_n(df[[variable]])
  if (n_unique >= requested_df + 1L) return(list(term = spline_term, note = NA_character_))
  if (n_unique >= 2L) return(list(term = linear_term, note = paste0(variable, "_spline_reduced_to_linear")))
  list(term = NA_character_, note = paste0(variable, "_omitted_constant_or_unavailable"))
}

sanitize_extra_terms <- function(df, extra_terms) {
  terms <- character()
  notes <- character()
  for (term in extra_terms) {
    if (grepl("icu_patient_address_mean_rmax_pct", term, fixed = TRUE)) {
      selected <- choose_spline_term(
        df,
        "icu_patient_address_mean_rmax_pct",
        "ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)",
        RMAX_DF
      )
    } else if (grepl("icu_patient_address_mean_no2", term, fixed = TRUE)) {
      selected <- choose_spline_term(df, "icu_patient_address_mean_no2", term, 1L)
    } else if (grepl("icu_patient_address_mean_pm25", term, fixed = TRUE)) {
      selected <- choose_spline_term(df, "icu_patient_address_mean_pm25", term, 1L)
    } else {
      selected <- list(term = term, note = NA_character_)
    }
    if (!is.na(selected$term)) terms <- c(terms, selected$term)
    if (!is.na(selected$note)) notes <- c(notes, selected$note)
  }
  list(terms = terms, notes = unique(notes))
}

choose_time_adjustment <- function(df, requested_df) {
  n_days <- finite_unique_n(df$time_index)
  if (n_days < 14L) return(list(term = NA_character_, df = NA_integer_, note = "time_adjustment_omitted_too_few_days"))
  capped_df <- min(as.integer(requested_df), max(1L, n_days - 2L))
  if (capped_df >= 2L) {
    note <- if (capped_df < requested_df) paste0("time_df_capped_from_", requested_df, "_to_", capped_df) else NA_character_
    return(list(term = "ns(time_index, df = time_df)", df = capped_df, note = note))
  }
  list(term = "time_index", df = NA_integer_, note = "time_spline_reduced_to_linear")
}

build_formula <- function(extra_terms = character(), time_term = "ns(time_index, df = time_df)", include_dow = TRUE, include_year = TRUE) {
  rhs <- c("cb_temp", "offset(log(denominator_icu_admissions))")
  if (!is.na(time_term)) rhs <- c(rhs, time_term)
  if (include_dow) rhs <- c(rhs, "dow")
  if (include_year) rhs <- c(rhs, "factor(year)")
  rhs <- c(rhs, extra_terms)
  as.formula(paste("ohca_admissions ~", paste(rhs, collapse = " + ")))
}

run_rate_dlnm_spec <- function(
  df,
  label = "Overall",
  model,
  extra_terms = c("ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)"),
  reference = c("median", "mrt"),
  time_df_per_year = TIME_DF_PER_YEAR,
  include_dow = TRUE,
  include_year = TRUE
) {
  reference <- match.arg(reference)
  needed <- c("denominator_icu_admissions", "icu_patient_address_mean_tmax_c", "icu_patient_address_mean_rmax_pct")
  if (any(grepl("no2", extra_terms, fixed = TRUE))) needed <- c(needed, "icu_patient_address_mean_no2")
  if (any(grepl("pm25", extra_terms, fixed = TRUE))) needed <- c(needed, "icu_patient_address_mean_pm25")
  df <- df[complete.cases(df[, needed, drop = FALSE]), , drop = FALSE]
  for (needed_col in needed) {
    df <- df[is.finite(suppressWarnings(as.numeric(df[[needed_col]]))), , drop = FALSE]
  }
  df <- df[df$denominator_icu_admissions > 0, , drop = FALSE]
  if (nrow(df) <= MAX_LAG + 14L) stop("Only ", nrow(df), " complete denominator/exposure days after filtering", call. = FALSE)
  if (sum(df$ohca_admissions, na.rm = TRUE) == 0L) stop("No OHCA admissions available", call. = FALSE)
  if (finite_unique_n(df$icu_patient_address_mean_tmax_c) < 2L) stop("Temperature exposure is constant or unavailable", call. = FALSE)

  df$time_index <- seq_len(nrow(df))
  temp_unique <- finite_unique_n(df$icu_patient_address_mean_tmax_c)
  temp_var_df <- min(VAR_DF, max(2L, temp_unique - 1L))
  cb_temp <- crossbasis(
    df$icu_patient_address_mean_tmax_c,
    lag = MAX_LAG,
    argvar = list(fun = "ns", df = temp_var_df),
    arglag = list(fun = "ns", df = LAG_DF)
  )

  extra_info <- sanitize_extra_terms(df, extra_terms)
  requested_time_df <- length(unique(df$year)) * time_df_per_year
  time_info <- choose_time_adjustment(df, requested_time_df)
  model_include_dow <- include_dow && length(unique(stats::na.omit(df$dow))) > 1L
  model_include_year <- include_year && length(unique(stats::na.omit(df$year))) > 1L
  formula_spec <- build_formula(
    extra_info$terms,
    time_term = time_info$term,
    include_dow = model_include_dow,
    include_year = model_include_year
  )
  environment(formula_spec) <- environment()
  environment(formula_spec)$time_df <- time_info$df

  fit <- glm(formula_spec, data = df, family = quasipoisson(link = "log"))
  grid <- make_prediction_grid(df$icu_patient_address_mean_tmax_c)
  initial_center <- median(df$icu_patient_address_mean_tmax_c, na.rm = TRUE)
  pred_initial <- crosspred(cb_temp, fit, cen = initial_center, at = grid)
  center <- if (reference == "median") initial_center else grid[which.min(pred_initial$allRRfit)]
  pred <- crosspred(cb_temp, fit, cen = center, at = grid, cumul = TRUE)
  reduced <- crossreduce(cb_temp, fit, cen = center)
  hot_temp <- grid[which.min(abs(grid - safe_quantile(df$icu_patient_address_mean_tmax_c, 0.95)))]
  hot_index <- which.min(abs(grid - hot_temp))

  result <- data.frame(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    n_days = nrow(df),
    n_ohca = sum(df$ohca_admissions, na.rm = TRUE),
    n_icu_admissions = sum(df$denominator_icu_admissions, na.rm = TRUE),
    denominator_source = unique(df$denominator_source)[[1]],
    crude_ohca_per_100_icu_admissions = 100 * sum(df$ohca_admissions, na.rm = TRUE) / sum(df$denominator_icu_admissions, na.rm = TRUE),
    reference_type = reference,
    reference_temp_c = center,
    hot_temp_c = hot_temp,
    cumulative_rate_ratio = as.numeric(pred$allRRfit[hot_index]),
    cumulative_rate_ratio_low = as.numeric(pred$allRRlow[hot_index]),
    cumulative_rate_ratio_high = as.numeric(pred$allRRhigh[hot_index]),
    log_rate_ratio = log(as.numeric(pred$allRRfit[hot_index])),
    log_rate_ratio_se = (log(as.numeric(pred$allRRhigh[hot_index])) - log(as.numeric(pred$allRRlow[hot_index]))) / (2 * 1.96),
    dispersion = summary(fit)$dispersion,
    model_family = "quasipoisson_rate_offset",
    time_df_per_year = time_df_per_year,
    includes_day_of_week = model_include_dow,
    includes_year_fixed_effect = model_include_year,
    estimable = TRUE,
    converged = fit$converged,
    covariate_notes = paste(extra_info$notes, collapse = "; "),
    stringsAsFactors = FALSE
  )

  curve <- data.frame(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    reference_type = reference,
    reference_temp_c = center,
    tmax_mean_c = grid,
    cumulative_rate_ratio = as.numeric(pred$allRRfit),
    cumulative_rate_ratio_low = as.numeric(pred$allRRlow),
    cumulative_rate_ratio_high = as.numeric(pred$allRRhigh),
    log_rate_ratio = log(as.numeric(pred$allRRfit)),
    log_rate_ratio_se = (log(as.numeric(pred$allRRhigh)) - log(as.numeric(pred$allRRlow))) / (2 * 1.96),
    stringsAsFactors = FALSE
  )

  reduced_coef <- data.frame(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    reference_type = reference,
    reference_temp_c = center,
    coefficient = names(coef(reduced)),
    estimate = as.numeric(coef(reduced)),
    stringsAsFactors = FALSE
  )

  reduced_vcov_matrix <- vcov(reduced)
  reduced_vcov <- as.data.frame(as.table(reduced_vcov_matrix), stringsAsFactors = FALSE)
  names(reduced_vcov) <- c("coefficient_row", "coefficient_col", "covariance")
  reduced_vcov$analysis_period <- ANALYSIS_PERIOD_LABEL
  reduced_vcov$stratum <- label
  reduced_vcov$model <- model
  reduced_vcov$reference_type <- reference
  reduced_vcov$reference_temp_c <- center
  reduced_vcov <- reduced_vcov[, c("analysis_period", "stratum", "model", "reference_type", "reference_temp_c", "coefficient_row", "coefficient_col", "covariance")]

  lag_summary <- NULL
  if (!is.null(pred$cumRRfit) && !is.null(pred$cumRRlow) && !is.null(pred$cumRRhigh)) {
    lag_names <- colnames(pred$cumRRfit)
    if (is.null(lag_names)) lag_names <- paste0("lag", seq_len(ncol(pred$cumRRfit)) - 1L)
    lag_ends <- suppressWarnings(as.integer(gsub("[^0-9]", "", lag_names)))
    if (any(!is.finite(lag_ends))) lag_ends <- seq_len(ncol(pred$cumRRfit)) - 1L
    lag_summary <- data.frame(
      analysis_period = ANALYSIS_PERIOD_LABEL,
      stratum = label,
      model = model,
      reference_type = reference,
      reference_temp_c = center,
      hot_temp_c = hot_temp,
      lag_start = 0L,
      lag_end = lag_ends,
      lag_label = ifelse(lag_ends == 0L, "Lag 0", paste0("Lag 0-", lag_ends)),
      cumulative_rate_ratio = as.numeric(pred$cumRRfit[hot_index, ]),
      cumulative_rate_ratio_low = as.numeric(pred$cumRRlow[hot_index, ]),
      cumulative_rate_ratio_high = as.numeric(pred$cumRRhigh[hot_index, ]),
      stringsAsFactors = FALSE
    )
    lag_summary$log_rate_ratio <- log(lag_summary$cumulative_rate_ratio)
    lag_summary$log_rate_ratio_se <- (
      log(lag_summary$cumulative_rate_ratio_high) - log(lag_summary$cumulative_rate_ratio_low)
    ) / (2 * 1.96)
  }

  lag_specific <- NULL
  if (!is.null(pred$matRRfit) && !is.null(pred$matRRlow) && !is.null(pred$matRRhigh)) {
    lag_names <- colnames(pred$matRRfit)
    if (is.null(lag_names)) lag_names <- paste0("lag", seq_len(ncol(pred$matRRfit)) - 1L)
    lag_values <- suppressWarnings(as.integer(gsub("[^0-9]", "", lag_names)))
    if (any(!is.finite(lag_values))) lag_values <- seq_len(ncol(pred$matRRfit)) - 1L
    lag_specific <- data.frame(
      analysis_period = ANALYSIS_PERIOD_LABEL,
      stratum = label,
      model = model,
      reference_type = reference,
      reference_temp_c = center,
      hot_temp_c = hot_temp,
      lag = lag_values,
      lag_label = paste0("Lag ", lag_values),
      rate_ratio = as.numeric(pred$matRRfit[hot_index, ]),
      rate_ratio_low = as.numeric(pred$matRRlow[hot_index, ]),
      rate_ratio_high = as.numeric(pred$matRRhigh[hot_index, ]),
      stringsAsFactors = FALSE
    )
    lag_specific$log_rate_ratio <- log(lag_specific$rate_ratio)
    lag_specific$log_rate_ratio_se <- (
      log(lag_specific$rate_ratio_high) - log(lag_specific$rate_ratio_low)
    ) / (2 * 1.96)
  }

  list(
    result = result,
    curve = curve,
    reduced_coef = reduced_coef,
    reduced_vcov = reduced_vcov,
    lag_summary = lag_summary,
    lag_specific = lag_specific
  )
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
if ("n_total_icu_admissions" %in% names(model_df)) {
  model_df$denominator_icu_admissions <- model_df$n_total_icu_admissions
  model_df$denominator_source <- "n_total_icu_admissions"
} else {
  model_df$denominator_icu_admissions <- model_df$n_icu_admissions
  model_df$denominator_source <- "n_icu_admissions_with_assigned_county_fallback"
}
if (ANALYSIS_SEASON == "warm") model_df <- model_df[model_df$month %in% WARM_MONTHS, , drop = FALSE]
model_df <- model_df[
  !is.na(model_df$icu_patient_address_mean_tmax_c) &
    !is.na(model_df$denominator_icu_admissions) &
    model_df$denominator_icu_admissions > 0,
  ,
  drop = FALSE
]

cohort <- read.csv(ohca_cohort_path, stringsAsFactors = FALSE)
cohort$admission_date <- as.Date(cohort$admission_date)
cohort$age_group <- ifelse(as.numeric(cohort$age_at_admission) >= 65, ">=65", "<65")
cohort$race_group <- ifelse(is_black_race(cohort$race_category), "Black", "Non-Black")

make_stratum_model_df <- function(subset_df) {
  counts <- aggregate(hospitalization_id ~ admission_date, data = subset_df, FUN = function(x) length(unique(x)))
  names(counts)[2] <- "ohca_admissions"
  out <- model_df
  out <- merge(
    out[, setdiff(names(out), "ohca_admissions"), drop = FALSE],
    counts,
    by = "admission_date",
    all.x = TRUE,
    sort = TRUE
  )
  out$ohca_admissions[is.na(out$ohca_admissions)] <- 0L
  out
}

rate_denominator_summary <- data.frame(
  analysis_period = ANALYSIS_PERIOD_LABEL,
  n_days = nrow(model_df),
  n_ohca = sum(model_df$ohca_admissions, na.rm = TRUE),
  n_icu_admissions = sum(model_df$denominator_icu_admissions, na.rm = TRUE),
  denominator_source = unique(model_df$denominator_source)[[1]],
  crude_ohca_per_100_icu_admissions = 100 * sum(model_df$ohca_admissions, na.rm = TRUE) / sum(model_df$denominator_icu_admissions, na.rm = TRUE),
  median_daily_icu_admissions = median(model_df$denominator_icu_admissions, na.rm = TRUE),
  p25_daily_icu_admissions = quantile(model_df$denominator_icu_admissions, 0.25, na.rm = TRUE, names = FALSE),
  p75_daily_icu_admissions = quantile(model_df$denominator_icu_admissions, 0.75, na.rm = TRUE, names = FALSE),
  min_daily_icu_admissions = min(model_df$denominator_icu_admissions, na.rm = TRUE),
  max_daily_icu_admissions = max(model_df$denominator_icu_admissions, na.rm = TRUE),
  stringsAsFactors = FALSE
)

rate_time_series <- model_df |>
  transmute(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    admission_date = .data$admission_date,
    year = .data$year,
    month = .data$month,
    dow = as.character(.data$dow),
    ohca_admissions = .data$ohca_admissions,
    n_icu_admissions = .data$denominator_icu_admissions,
    denominator_source = .data$denominator_source,
    ohca_rate_per_100_icu_admissions = 100 * .data$ohca_admissions / .data$denominator_icu_admissions,
    tmax_mean_c = .data$icu_patient_address_mean_tmax_c,
    rmax_mean_pct = .data$icu_patient_address_mean_rmax_pct,
    no2_mean = .data$icu_patient_address_mean_no2,
    pm25_mean = .data$icu_patient_address_mean_pm25
  )

fits <- list(
  overall_median_humidity = run_rate_dlnm_spec(
    model_df,
    label = "Overall",
    model = "rate_median_humidity_adjusted",
    reference = "median"
  ),
  overall_median_pollution = run_rate_dlnm_spec(
    model_df,
    label = "Overall",
    model = "rate_median_humidity_pollution_adjusted",
    extra_terms = c("ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)", "icu_patient_address_mean_no2", "icu_patient_address_mean_pm25"),
    reference = "median"
  ),
  overall_mrt_pollution = run_rate_dlnm_spec(
    model_df,
    label = "Overall",
    model = "rate_mrt_humidity_pollution_adjusted",
    extra_terms = c("ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)", "icu_patient_address_mean_no2", "icu_patient_address_mean_pm25"),
    reference = "mrt"
  )
)

strata_specs <- list(
  male = list(label = "Male", df = cohort[is_male(cohort$sex_category), ]),
  female = list(label = "Female", df = cohort[is_female(cohort$sex_category), ]),
  age_lt65 = list(label = "<65", df = cohort[cohort$age_group == "<65", ]),
  age_ge65 = list(label = ">=65", df = cohort[cohort$age_group == ">=65", ]),
  race_black = list(label = "Black", df = cohort[cohort$race_group == "Black", ]),
  race_nonblack = list(label = "Non-Black", df = cohort[cohort$race_group == "Non-Black", ])
)

for (nm in names(strata_specs)) {
  spec <- strata_specs[[nm]]
  stratum_model_df <- make_stratum_model_df(spec$df)
  fits[[paste0(nm, "_median_humidity")]] <- run_rate_dlnm_spec(
    stratum_model_df,
    label = spec$label,
    model = "rate_stratified_humidity_adjusted",
    reference = "median"
  )
  fits[[paste0(nm, "_median_pollution")]] <- run_rate_dlnm_spec(
    stratum_model_df,
    label = spec$label,
    model = "rate_stratified_humidity_pollution_adjusted",
    extra_terms = c("ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)", "icu_patient_address_mean_no2", "icu_patient_address_mean_pm25"),
    reference = "median"
  )
  fits[[paste0(nm, "_mrt_pollution")]] <- run_rate_dlnm_spec(
    stratum_model_df,
    label = spec$label,
    model = "rate_stratified_mrt_humidity_pollution_adjusted",
    extra_terms = c("ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)", "icu_patient_address_mean_no2", "icu_patient_address_mean_pm25"),
    reference = "mrt"
  )
}

results_df <- do.call(rbind, lapply(fits, `[[`, "result"))
curves_df <- do.call(rbind, lapply(fits, `[[`, "curve"))
reduced_coef_df <- do.call(rbind, lapply(fits, `[[`, "reduced_coef"))
reduced_vcov_df <- do.call(rbind, lapply(fits, `[[`, "reduced_vcov"))
lag_summary_df <- do.call(rbind, lapply(fits, `[[`, "lag_summary"))
lag_specific_df <- do.call(rbind, lapply(fits, `[[`, "lag_specific"))

readr::write_csv(results_df, file.path(output_dir, "ohca_icu_admission_rate_dlnm_results.csv"))
readr::write_csv(curves_df, file.path(output_dir, "ohca_icu_admission_rate_dlnm_curves.csv"))
readr::write_csv(reduced_coef_df, file.path(output_dir, "ohca_icu_admission_rate_dlnm_reduced_coefficients.csv"))
readr::write_csv(reduced_vcov_df, file.path(output_dir, "ohca_icu_admission_rate_dlnm_reduced_vcov.csv"))
readr::write_csv(lag_summary_df, file.path(output_dir, "ohca_icu_admission_rate_dlnm_lag_summaries.csv"))
readr::write_csv(lag_specific_df, file.path(output_dir, "ohca_icu_admission_rate_dlnm_lag_specific_summaries.csv"))
readr::write_csv(rate_denominator_summary, file.path(output_dir, "ohca_icu_admission_rate_denominator_summary.csv"))
readr::write_csv(rate_time_series, file.path(output_dir, "ohca_icu_admission_rate_daily_timeseries.csv"))

MOVING_AVERAGE_DAYS <- 7L

rate_smooth <- stats::filter(
  rate_time_series$ohca_rate_per_100_icu_admissions,
  rep(1 / MOVING_AVERAGE_DAYS, MOVING_AVERAGE_DAYS),
  sides = 2
)
tmax_smooth <- stats::filter(
  rate_time_series$tmax_mean_c,
  rep(1 / MOVING_AVERAGE_DAYS, MOVING_AVERAGE_DAYS),
  sides = 2
)
plot_df <- rate_time_series |>
  mutate(
    ohca_rate_smoothed = as.numeric(rate_smooth),
    tmax_smoothed = as.numeric(tmax_smooth),
    ohca_rate_display = pmin(
      .data$ohca_rate_per_100_icu_admissions,
      as.numeric(stats::quantile(.data$ohca_rate_per_100_icu_admissions, 0.99, na.rm = TRUE, names = FALSE))
    )
  )

rate_axis_upper <- max(
  10,
  as.numeric(stats::quantile(plot_df$ohca_rate_per_100_icu_admissions, 0.99, na.rm = TRUE, names = FALSE))
)
rate_axis_range <- c(0, rate_axis_upper)
tmax_range <- range(plot_df$tmax_mean_c, na.rm = TRUE)
scale_factor <- diff(rate_axis_range) / diff(tmax_range)
scale_offset <- rate_axis_range[[1]] - tmax_range[[1]] * scale_factor

timeseries_figure <- ggplot(plot_df, aes(x = .data$admission_date)) +
  geom_col(aes(y = .data$ohca_rate_display), fill = "#D6DEE8", width = 1, alpha = 0.50) +
  geom_line(aes(y = .data$ohca_rate_smoothed), color = "#1F4E79", linewidth = 0.75, na.rm = TRUE) +
  geom_line(aes(y = .data$tmax_smoothed * scale_factor + scale_offset), color = "#9A3412", linewidth = 0.75, na.rm = TRUE) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(
    name = "OHCA per 100 ICU admissions",
    limits = rate_axis_range,
    expand = expansion(mult = c(0, 0.04)),
    sec.axis = sec_axis(~ (.x - scale_offset) / scale_factor, name = "Daily Tmax (C)")
  ) +
  labs(
    title = "Daily OHCA ICU admission rate and temperature",
    subtitle = "Bars show daily OHCA per 100 ICU admissions capped at the 99th percentile; lines are centered 7-day moving averages",
    x = "Admission date"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(color = "#1F4E79", face = "bold"),
    axis.title.y.right = element_text(color = "#9A3412", face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(figure_dir, "figure_ohca_icu_admission_rate_temperature_timeseries.png"),
  timeseries_figure,
  width = 10,
  height = 5.6,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "figure_ohca_icu_admission_rate_temperature_timeseries.pdf"),
  timeseries_figure,
  width = 10,
  height = 5.6
)

print(results_df)
message("Wrote OHCA per ICU admission rate DLNM outputs to ", output_dir)

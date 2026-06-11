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
  library(arrow)
  library(dlnm)
  library(dplyr)
  library(readr)
  library(splines)
  library(stringr)
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
MAX_LAG <- 5L
VAR_DF <- 4L
LAG_DF <- 3L
RMAX_DF <- 3L
MIN_STRATUM_OHCA <- 30L
MIN_STRATUM_EVENT_DAYS <- 10L

safe_quantile <- function(x, prob) stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7)

finite_unique_n <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  length(unique(x[is.finite(x)]))
}

make_prediction_grid <- function(x) {
  unique(as.numeric(seq(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)), by = 0.5)))
}

filter_analysis_period_cases <- function(df) {
  out <- df
  if (ANALYSIS_SEASON == "warm") {
    out <- out[out$month %in% WARM_MONTHS, , drop = FALSE]
  }
  out
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
      county_fips = case_row$county_fips[[1]],
      sex_category = case_row$sex_category[[1]],
      age_group = case_row$age_group[[1]],
      race_group = case_row$race_group[[1]]
    )
  }))
}

add_lagged_daily_exposures <- function(df, daily_exposure, value_col, prefix) {
  out <- df
  for (lag in 0:MAX_LAG) {
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

make_nonestimable_result <- function(df, label, model, reference, reason, omitted_terms = NA_character_) {
  tibble::tibble(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    n_referent_rows = nrow(df),
    n_case_events = if ("case" %in% names(df)) sum(df$case == 1, na.rm = TRUE) else 0L,
    n_match_sets = if ("match_id" %in% names(df)) dplyr::n_distinct(df$match_id) else 0L,
    reference_type = reference,
    reference_temp_c = NA_real_,
    hot_temp_c = NA_real_,
    cumulative_rr = NA_real_,
    cumulative_rr_low = NA_real_,
    cumulative_rr_high = NA_real_,
    log_rr = NA_real_,
    log_rr_se = NA_real_,
    model_family = "conditional_logistic",
    estimable = FALSE,
    converged = FALSE,
    model_type = "not_estimable",
    omitted_terms = omitted_terms,
    skip_reason = reason
  )
}

run_case_crossover_dlnm <- function(
  df,
  label,
  model,
  reference = c("median", "mrt"),
  adjust_pollution = FALSE,
  allow_skip = FALSE,
  min_ohca = 0L,
  min_event_days = 0L,
  return_curve = FALSE
) {
  reference <- match.arg(reference)
  tmax_cols <- paste0("tmax_mean_c_lag", 0:MAX_LAG)
  needed <- c(tmax_cols, "rmax_mean_pct_lag0")
  df <- df[complete.cases(df[, needed, drop = FALSE]), , drop = FALSE]
  for (needed_col in needed) {
    df <- df[is.finite(suppressWarnings(as.numeric(df[[needed_col]]))), , drop = FALSE]
  }
  df <- df |>
    group_by(.data$match_id) |>
    filter(sum(.data$case == 1L, na.rm = TRUE) == 1L, n() >= 2L) |>
    ungroup()

  n_ohca <- sum(df$case == 1L, na.rm = TRUE)
  n_event_days <- n_distinct(df$referent_date[df$case == 1L])
  omitted_terms <- if (adjust_pollution) {
    "icu_patient_address_mean_no2 and icu_patient_address_mean_pm25 omitted: annual county pollution is constant within time-stratified county/year/month/day-of-week matched sets"
  } else {
    NA_character_
  }
  skip_reason <- NULL
  if (nrow(df) == 0L) {
    skip_reason <- "No complete referent rows available"
  } else if (n_ohca == 0L) {
    skip_reason <- "No OHCA case events available"
  } else if (n_ohca < min_ohca) {
    skip_reason <- paste0("Only ", n_ohca, " OHCA case events; minimum for this model is ", min_ohca)
  } else if (n_event_days < min_event_days) {
    skip_reason <- paste0("Only ", n_event_days, " event days; minimum for this model is ", min_event_days)
  } else if (finite_unique_n(df$tmax_mean_c_lag0) < 2L) {
    skip_reason <- "Temperature exposure is constant or unavailable"
  }
  if (!is.null(skip_reason)) {
    if (!allow_skip) stop(skip_reason, call. = FALSE)
    skipped <- make_nonestimable_result(df, label, model, reference, skip_reason, omitted_terms)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL, lag_summary = NULL, lag_specific = NULL, lag_temperature_surface = NULL))
    return(skipped)
  }

  temp_unique <- finite_unique_n(df$tmax_mean_c_lag0)
  temp_var_df <- min(VAR_DF, max(2L, temp_unique - 1L))
  temp_history <- as.matrix(df[, tmax_cols, drop = FALSE])
  cb_temp <- tryCatch(
    crossbasis(
      temp_history,
      lag = c(0L, MAX_LAG),
      argvar = list(fun = "ns", df = temp_var_df),
      arglag = list(fun = "ns", df = LAG_DF)
    ),
    error = function(e) e
  )
  if (inherits(cb_temp, "error")) {
    if (!allow_skip) stop(cb_temp)
    skipped <- make_nonestimable_result(df, label, model, reference, paste0("DLNM crossbasis failed: ", conditionMessage(cb_temp)), omitted_terms)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL, lag_summary = NULL, lag_specific = NULL, lag_temperature_surface = NULL))
    return(skipped)
  }

  rmax_term <- if (finite_unique_n(df$rmax_mean_pct_lag0) >= RMAX_DF + 1L) {
    "ns(rmax_mean_pct_lag0, df = RMAX_DF)"
  } else if (finite_unique_n(df$rmax_mean_pct_lag0) >= 2L) {
    "rmax_mean_pct_lag0"
  } else {
    NA_character_
  }
  rhs <- c("cb_temp", rmax_term, "strata(match_id)")
  formula_spec <- as.formula(paste("case ~", paste(rhs[!is.na(rhs)], collapse = " + ")))
  environment(formula_spec) <- environment()

  fit <- tryCatch(
    survival::clogit(formula_spec, data = df, method = "efron"),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    if (!allow_skip) stop(fit)
    skipped <- make_nonestimable_result(df, label, model, reference, paste0("Conditional logistic DLNM failed: ", conditionMessage(fit)), omitted_terms)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL, lag_summary = NULL, lag_specific = NULL, lag_temperature_surface = NULL))
    return(skipped)
  }

  grid <- make_prediction_grid(df$tmax_mean_c_lag0)
  initial_center <- median(df$tmax_mean_c_lag0, na.rm = TRUE)
  prediction_objects <- tryCatch({
    pred_initial <- crosspred(cb_temp, fit, cen = initial_center, at = grid)
    center <- if (reference == "median") initial_center else grid[which.min(pred_initial$allRRfit)]
    list(
      center = center,
      pred = crosspred(cb_temp, fit, cen = center, at = grid, cumul = TRUE),
      reduced = crossreduce(cb_temp, fit, cen = center)
    )
  }, error = function(e) e)
  if (inherits(prediction_objects, "error")) {
    if (!allow_skip) stop(prediction_objects)
    skipped <- make_nonestimable_result(df, label, model, reference, paste0("DLNM prediction failed: ", conditionMessage(prediction_objects)), omitted_terms)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL, lag_summary = NULL, lag_specific = NULL, lag_temperature_surface = NULL))
    return(skipped)
  }

  center <- prediction_objects$center
  pred <- prediction_objects$pred
  reduced <- prediction_objects$reduced
  hot_temp <- grid[which.min(abs(grid - safe_quantile(df$tmax_mean_c_lag0, 0.95)))]
  hot_index <- which.min(abs(grid - hot_temp))
  log_rr <- log(as.numeric(pred$allRRfit[hot_index]))
  log_rr_se <- (log(as.numeric(pred$allRRhigh[hot_index])) - log(as.numeric(pred$allRRlow[hot_index]))) / (2 * 1.96)

  result <- tibble::tibble(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    n_referent_rows = nrow(df),
    n_case_events = sum(df$case == 1L, na.rm = TRUE),
    n_match_sets = n_distinct(df$match_id),
    reference_type = reference,
    reference_temp_c = center,
    hot_temp_c = hot_temp,
    cumulative_rr = as.numeric(pred$allRRfit[hot_index]),
    cumulative_rr_low = as.numeric(pred$allRRlow[hot_index]),
    cumulative_rr_high = as.numeric(pred$allRRhigh[hot_index]),
    log_rr = log_rr,
    log_rr_se = log_rr_se,
    model_family = "conditional_logistic",
    estimable = TRUE,
    converged = all(is.finite(stats::coef(fit))),
    model_type = "case_crossover_dlnm",
    omitted_terms = omitted_terms,
    skip_reason = NA_character_
  )

  curve <- tibble::tibble(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    reference_type = reference,
    reference_temp_c = center,
    tmax_mean_c = grid,
    cumulative_rr = as.numeric(pred$allRRfit),
    cumulative_rr_low = as.numeric(pred$allRRlow),
    cumulative_rr_high = as.numeric(pred$allRRhigh),
    log_rr = log(as.numeric(pred$allRRfit)),
    log_rr_se = (log(as.numeric(pred$allRRhigh)) - log(as.numeric(pred$allRRlow))) / (2 * 1.96)
  )

  reduced_coef <- tibble::tibble(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    reference_type = reference,
    reference_temp_c = center,
    coefficient = names(coef(reduced)),
    estimate = as.numeric(coef(reduced))
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
    lag_summary <- tibble::tibble(
      analysis_period = ANALYSIS_PERIOD_LABEL,
      stratum = label,
      model = model,
      reference_type = reference,
      reference_temp_c = center,
      hot_temp_c = hot_temp,
      lag_start = 0L,
      lag_end = lag_ends,
      lag_label = ifelse(lag_ends == 0L, "Lag 0", paste0("Lag 0-", lag_ends)),
      cumulative_rr = as.numeric(pred$cumRRfit[hot_index, ]),
      cumulative_rr_low = as.numeric(pred$cumRRlow[hot_index, ]),
      cumulative_rr_high = as.numeric(pred$cumRRhigh[hot_index, ])
    ) |>
      mutate(
        log_rr = log(.data$cumulative_rr),
        log_rr_se = (log(.data$cumulative_rr_high) - log(.data$cumulative_rr_low)) / (2 * 1.96)
      )
  }

  lag_specific <- NULL
  if (!is.null(pred$matRRfit) && !is.null(pred$matRRlow) && !is.null(pred$matRRhigh)) {
    lag_names <- colnames(pred$matRRfit)
    if (is.null(lag_names)) lag_names <- paste0("lag", seq_len(ncol(pred$matRRfit)) - 1L)
    lag_values <- suppressWarnings(as.integer(gsub("[^0-9]", "", lag_names)))
    if (any(!is.finite(lag_values))) lag_values <- seq_len(ncol(pred$matRRfit)) - 1L
    lag_specific <- tibble::tibble(
      analysis_period = ANALYSIS_PERIOD_LABEL,
      stratum = label,
      model = model,
      reference_type = reference,
      reference_temp_c = center,
      hot_temp_c = hot_temp,
      lag = lag_values,
      lag_label = paste0("Lag ", lag_values),
      rr = as.numeric(pred$matRRfit[hot_index, ]),
      rr_low = as.numeric(pred$matRRlow[hot_index, ]),
      rr_high = as.numeric(pred$matRRhigh[hot_index, ])
    ) |>
      mutate(
        log_rr = log(.data$rr),
        log_rr_se = (log(.data$rr_high) - log(.data$rr_low)) / (2 * 1.96)
      )
  }

  lag_temperature_surface <- make_dlnm_lag_temperature_surface(
    pred = pred,
    grid = grid,
    label = label,
    model = model,
    reference = reference,
    center = center,
    hot_temp = hot_temp,
    effect_prefix = "rr",
    extra_cols = list(analysis_period = ANALYSIS_PERIOD_LABEL)
  )

  if (return_curve) {
    return(list(
      result = result,
      curve = curve,
      reduced_coef = reduced_coef,
      reduced_vcov = reduced_vcov,
      lag_summary = lag_summary,
      lag_specific = lag_specific,
      lag_temperature_surface = lag_temperature_surface
    ))
  }

  result
}

ohca <- readr::read_csv(
  ohca_cohort_path,
  show_col_types = FALSE,
  col_types = readr::cols(.default = readr::col_guess(), hospitalization_id = readr::col_character(), patient_id = readr::col_character(), county_fips = readr::col_character())
) |>
  mutate(
    admission_date = as.Date(.data$admission_date),
    county_fips = normalize_county_fips(.data$county_fips),
    year = as.integer(format(.data$admission_date, "%Y")),
    month = as.integer(format(.data$admission_date, "%m")),
    dow = weekdays(.data$admission_date),
    age_group = ifelse(suppressWarnings(as.numeric(.data$age_at_admission)) >= 65, ">=65", "<65"),
    race_group = ifelse(is_black_race(.data$race_category), "Black", "Non-Black")
  ) |>
  filter(!is.na(.data$admission_date), !is.na(.data$county_fips)) |>
  filter_analysis_period_cases()

tmax <- read_exposome_daymet(repo_root, "daymet_county_tmax_2018_2024_conus.parquet", "tmax_mean_c")
rmax <- read_exposome_daymet(repo_root, "daymet_county_rmax_2018_2024.parquet", "rmax_mean_pct")
no2 <- read_exposome_pollution(repo_root, "no2_county_year.csv", "no2_mean")
pm25 <- read_exposome_pollution(repo_root, "pm25_county_year.csv", "pm25_mean")

daily_exposure <- full_join(tmax, rmax, by = c("admission_date", "county_fips")) |>
  mutate(year = as.integer(format(.data$admission_date, "%Y"))) |>
  left_join(no2, by = c("county_fips", "year")) |>
  left_join(pm25, by = c("county_fips", "year")) |>
  filter(.data$admission_date >= START_DATE - MAX_LAG, .data$admission_date <= END_DATE)

all_dates <- tibble::tibble(
  referent_date = seq(START_DATE, END_DATE, by = "day"),
  year = as.integer(format(.data$referent_date, "%Y")),
  month = as.integer(format(.data$referent_date, "%m")),
  dow = weekdays(.data$referent_date)
)

referents <- make_referent_dates(ohca, all_dates)
referents <- add_lagged_daily_exposures(referents, daily_exposure, "tmax_mean_c", "tmax_mean_c")
referents <- add_lagged_daily_exposures(referents, daily_exposure, "rmax_mean_pct", "rmax_mean_pct")
referents <- referents |>
  mutate(year = as.integer(format(.data$referent_date, "%Y"))) |>
  left_join(no2, by = c("county_fips", "year")) |>
  left_join(pm25, by = c("county_fips", "year"))

referent_summary <- referents |>
  group_by(.data$match_id) |>
  summarise(
    n_referent_days = n(),
    has_case_day = sum(.data$case == 1L, na.rm = TRUE) == 1L,
    complete_tmax_history_rows = sum(complete.cases(across(all_of(paste0("tmax_mean_c_lag", 0:MAX_LAG))))),
    .groups = "drop"
  )

results <- list()
curve_rows <- list()
reduced_coef_rows <- list()
reduced_vcov_rows <- list()
lag_summary_rows <- list()
lag_specific_rows <- list()
lag_temperature_surface_rows <- list()

run_and_store <- function(key, df, label, model, reference, adjust_pollution = FALSE, allow_skip = FALSE, min_ohca = 0L, min_event_days = 0L) {
  fit <- run_case_crossover_dlnm(
    df = df,
    label = label,
    model = model,
    reference = reference,
    adjust_pollution = adjust_pollution,
    allow_skip = allow_skip,
    min_ohca = min_ohca,
    min_event_days = min_event_days,
    return_curve = TRUE
  )
  results[[key]] <<- fit$result
  if (!is.null(fit$curve)) curve_rows[[key]] <<- fit$curve
  if (!is.null(fit$reduced_coef)) reduced_coef_rows[[key]] <<- fit$reduced_coef
  if (!is.null(fit$reduced_vcov)) reduced_vcov_rows[[key]] <<- fit$reduced_vcov
  if (!is.null(fit$lag_summary)) lag_summary_rows[[key]] <<- fit$lag_summary
  if (!is.null(fit$lag_specific)) lag_specific_rows[[key]] <<- fit$lag_specific
  if (!is.null(fit$lag_temperature_surface)) lag_temperature_surface_rows[[key]] <<- fit$lag_temperature_surface
  invisible(fit)
}

run_and_store("overall_primary", referents, "Overall", "primary_humidity_adjusted", "median")
run_and_store("overall_pollution", referents, "Overall", "sensitivity_humidity_pollution_adjusted", "median", adjust_pollution = TRUE)
run_and_store("overall_mrt", referents, "Overall", "sensitivity_mrt_reference", "mrt", adjust_pollution = TRUE)

strata_specs <- list(
  male = list(label = "Male", mask = is_male(referents$sex_category)),
  female = list(label = "Female", mask = is_female(referents$sex_category)),
  age_lt65 = list(label = "<65", mask = referents$age_group == "<65"),
  age_ge65 = list(label = ">=65", mask = referents$age_group == ">=65"),
  race_black = list(label = "Black", mask = referents$race_group == "Black"),
  race_nonblack = list(label = "Non-Black", mask = referents$race_group == "Non-Black")
)

for (nm in names(strata_specs)) {
  spec <- strata_specs[[nm]]
  df <- referents[spec$mask %in% TRUE, , drop = FALSE]
  run_and_store(
    paste0(nm, "_primary"),
    df,
    spec$label,
    "stratified_humidity_adjusted",
    "median",
    allow_skip = TRUE,
    min_ohca = MIN_STRATUM_OHCA,
    min_event_days = MIN_STRATUM_EVENT_DAYS
  )
  run_and_store(
    paste0(nm, "_pollution"),
    df,
    spec$label,
    "stratified_humidity_pollution_adjusted",
    "median",
    adjust_pollution = TRUE,
    allow_skip = TRUE,
    min_ohca = MIN_STRATUM_OHCA,
    min_event_days = MIN_STRATUM_EVENT_DAYS
  )
  run_and_store(
    paste0(nm, "_mrt"),
    df,
    spec$label,
    "stratified_mrt_reference",
    "mrt",
    adjust_pollution = TRUE,
    allow_skip = TRUE,
    min_ohca = MIN_STRATUM_OHCA,
    min_event_days = MIN_STRATUM_EVENT_DAYS
  )
}

results_df <- bind_rows(results)
curves_df <- bind_rows(curve_rows)
reduced_coef_df <- bind_rows(reduced_coef_rows)
reduced_vcov_df <- bind_rows(reduced_vcov_rows)
lag_summary_df <- bind_rows(lag_summary_rows)
lag_specific_df <- bind_rows(lag_specific_rows)
lag_temperature_surface_df <- bind_rows(lag_temperature_surface_rows)

readr::write_csv(results_df, file.path(output_dir, "case_crossover_dlnm_results.csv"))
readr::write_csv(curves_df, file.path(output_dir, "case_crossover_dlnm_curves.csv"))
readr::write_csv(reduced_coef_df, file.path(output_dir, "case_crossover_dlnm_reduced_coefficients.csv"))
readr::write_csv(reduced_vcov_df, file.path(output_dir, "case_crossover_dlnm_reduced_vcov.csv"))
readr::write_csv(lag_summary_df, file.path(output_dir, "case_crossover_dlnm_lag_summaries.csv"))
readr::write_csv(lag_specific_df, file.path(output_dir, "case_crossover_dlnm_lag_specific_summaries.csv"))
readr::write_csv(lag_temperature_surface_df, file.path(output_dir, "case_crossover_dlnm_lag_temperature_surface.csv"))
invisible(write_dlnm_lag_temperature_surface_pdf(
  lag_temperature_surface_df,
  file.path(output_dir, "case_crossover_dlnm_lag_temperature_surface_plots.pdf"),
  effect_col = "rr",
  title_prefix = paste("Case-crossover DLNM", ANALYSIS_PERIOD_LABEL)
))
readr::write_csv(referent_summary, file.path(output_dir, "case_crossover_referent_set_summary.csv"))

print(results_df)
message("Wrote time-stratified case-crossover DLNM results to ", output_dir)

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
  library(dplyr)
  library(readr)
  library(splines)
  library(survival)
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
VAR_DF <- 4L
RMAX_DF <- 3L
MIN_STRATUM_OHCA <- 30L
MIN_STRATUM_EVENT_DAYS <- 10L

finite_unique_n <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  length(unique(x[is.finite(x)]))
}

safe_quantile <- function(x, prob) stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7)

make_prediction_grid <- function(x) {
  unique(as.numeric(seq(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)), by = 0.5)))
}

filter_analysis_period_cases <- function(df) {
  if (ANALYSIS_SEASON == "warm") return(df[df$month %in% WARM_MONTHS, , drop = FALSE])
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
      county_fips = case_row$county_fips[[1]],
      sex_category = case_row$sex_category[[1]],
      age_group = case_row$age_group[[1]],
      race_group = case_row$race_group[[1]],
      ethnicity_group = case_row$ethnicity_group[[1]]
    )
  }))
}

make_nonestimable_result <- function(df, label, model, reference, reason, omitted_terms = NA_character_) {
  tibble::tibble(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    n_referent_rows = nrow(df),
    n_case_events = if ("case" %in% names(df)) sum(df$case == 1L, na.rm = TRUE) else 0L,
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
    model_type = "case_crossover_lag0",
    omitted_terms = omitted_terms,
    skip_reason = reason
  )
}

run_case_crossover_lag0 <- function(
  df,
  label,
  model,
  reference = c("median", "mrt"),
  adjust_pollution = FALSE,
  allow_skip = FALSE,
  min_ohca = 0L,
  min_event_days = 0L
) {
  reference <- match.arg(reference)
  needed <- c("tmax_mean_c", "rmax_mean_pct")
  df <- df[complete.cases(df[, needed, drop = FALSE]), , drop = FALSE]
  for (needed_col in needed) df <- df[is.finite(suppressWarnings(as.numeric(df[[needed_col]]))), , drop = FALSE]
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
  } else if (finite_unique_n(df$tmax_mean_c) < 2L) {
    skip_reason <- "Temperature exposure is constant or unavailable"
  }
  if (!is.null(skip_reason)) {
    if (!allow_skip) stop(skip_reason, call. = FALSE)
    return(list(result = make_nonestimable_result(df, label, model, reference, skip_reason, omitted_terms)))
  }

  temp_var_df <- min(VAR_DF, max(2L, finite_unique_n(df$tmax_mean_c) - 1L))
  tmax_basis <- splines::ns(df$tmax_mean_c, df = temp_var_df)
  tmax_cols <- paste0("tmax_basis_", seq_len(ncol(tmax_basis)))
  df[, tmax_cols] <- as.data.frame(tmax_basis)

  rmax_cols <- character()
  if (finite_unique_n(df$rmax_mean_pct) >= RMAX_DF + 1L) {
    rmax_basis <- splines::ns(df$rmax_mean_pct, df = RMAX_DF)
    rmax_cols <- paste0("rmax_basis_", seq_len(ncol(rmax_basis)))
    df[, rmax_cols] <- as.data.frame(rmax_basis)
  } else if (finite_unique_n(df$rmax_mean_pct) >= 2L) {
    rmax_cols <- "rmax_mean_pct"
  }

  formula_spec <- as.formula(paste("case ~", paste(c(tmax_cols, rmax_cols, "strata(match_id)"), collapse = " + ")))
  fit <- tryCatch(survival::clogit(formula_spec, data = df, method = "efron"), error = function(e) e)
  if (inherits(fit, "error")) {
    if (!allow_skip) stop(fit)
    return(list(result = make_nonestimable_result(df, label, model, reference, paste0("Conditional logistic lag-0 model failed: ", conditionMessage(fit)), omitted_terms)))
  }

  grid <- make_prediction_grid(df$tmax_mean_c)
  initial_center <- median(df$tmax_mean_c, na.rm = TRUE)
  coef_vec <- stats::coef(fit)
  vc <- stats::vcov(fit)
  temp_coef <- coef_vec[tmax_cols]
  temp_vcov <- vc[tmax_cols, tmax_cols, drop = FALSE]

  basis_grid <- predict(tmax_basis, newx = grid)
  basis_center_initial <- as.numeric(predict(tmax_basis, newx = initial_center))
  log_rr_initial <- as.numeric((basis_grid - matrix(basis_center_initial, nrow = nrow(basis_grid), ncol = length(basis_center_initial), byrow = TRUE)) %*% temp_coef)
  center <- if (reference == "median") initial_center else grid[which.min(exp(log_rr_initial))]
  basis_center <- as.numeric(predict(tmax_basis, newx = center))
  basis_diff <- basis_grid - matrix(basis_center, nrow = nrow(basis_grid), ncol = length(basis_center), byrow = TRUE)
  log_rr <- as.numeric(basis_diff %*% temp_coef)
  log_rr_se <- sqrt(pmax(rowSums((basis_diff %*% temp_vcov) * basis_diff), 0))
  rr <- exp(log_rr)
  rr_low <- exp(log_rr - 1.96 * log_rr_se)
  rr_high <- exp(log_rr + 1.96 * log_rr_se)

  hot_temp <- grid[which.min(abs(grid - safe_quantile(df$tmax_mean_c, 0.95)))]
  hot_index <- which.min(abs(grid - hot_temp))
  result <- tibble::tibble(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    n_referent_rows = nrow(df),
    n_case_events = n_ohca,
    n_match_sets = n_distinct(df$match_id),
    reference_type = reference,
    reference_temp_c = center,
    hot_temp_c = hot_temp,
    cumulative_rr = rr[hot_index],
    cumulative_rr_low = rr_low[hot_index],
    cumulative_rr_high = rr_high[hot_index],
    log_rr = log_rr[hot_index],
    log_rr_se = log_rr_se[hot_index],
    model_family = "conditional_logistic",
    estimable = TRUE,
    converged = all(is.finite(coef_vec)),
    model_type = "case_crossover_lag0",
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
    cumulative_rr = rr,
    cumulative_rr_low = rr_low,
    cumulative_rr_high = rr_high,
    log_rr = log_rr,
    log_rr_se = log_rr_se
  )

  pred_like <- list(allRRfit = rr, allRRlow = rr_low, allRRhigh = rr_high)
  contrast_summary <- make_dlnm_contrast_summary(
    pred = pred_like,
    grid = grid,
    source_x = df$tmax_mean_c,
    label = label,
    model = model,
    reference = reference,
    center = center,
    effect_prefix = "rr",
    extra_cols = list(analysis_period = ANALYSIS_PERIOD_LABEL, n_case_events = n_ohca)
  )

  coefficient <- tibble::tibble(
    analysis_period = ANALYSIS_PERIOD_LABEL,
    stratum = label,
    model = model,
    reference_type = reference,
    reference_temp_c = center,
    coefficient = names(coef_vec),
    estimate = as.numeric(coef_vec)
  )
  vcov_df <- as.data.frame(as.table(vc), stringsAsFactors = FALSE)
  names(vcov_df) <- c("coefficient_row", "coefficient_col", "covariance")
  vcov_df$analysis_period <- ANALYSIS_PERIOD_LABEL
  vcov_df$stratum <- label
  vcov_df$model <- model
  vcov_df$reference_type <- reference
  vcov_df$reference_temp_c <- center
  vcov_df <- vcov_df[, c("analysis_period", "stratum", "model", "reference_type", "reference_temp_c", "coefficient_row", "coefficient_col", "covariance")]

  list(result = result, curve = curve, contrast_summary = contrast_summary, coefficient = coefficient, vcov = vcov_df)
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
    race_group = ifelse(is_black_race(.data$race_category), "Black", "Non-Black"),
    ethnicity_group = ifelse(
      is_hispanic_ethnicity(.data$ethnicity_category),
      "Hispanic",
      ifelse(is_non_hispanic_ethnicity(.data$ethnicity_category), "Non-Hispanic", NA_character_)
    )
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
  filter(.data$admission_date >= START_DATE, .data$admission_date <= END_DATE)

all_dates <- tibble::tibble(
  referent_date = seq(START_DATE, END_DATE, by = "day"),
  year = as.integer(format(.data$referent_date, "%Y")),
  month = as.integer(format(.data$referent_date, "%m")),
  dow = weekdays(.data$referent_date)
)

referents <- make_referent_dates(ohca, all_dates) |>
  left_join(
    daily_exposure |>
      transmute(
        county_fips = .data$county_fips,
        referent_date = .data$admission_date,
        tmax_mean_c = .data$tmax_mean_c,
        rmax_mean_pct = .data$rmax_mean_pct
      ),
    by = c("county_fips", "referent_date")
  ) |>
  mutate(year = as.integer(format(.data$referent_date, "%Y"))) |>
  left_join(no2, by = c("county_fips", "year")) |>
  left_join(pm25, by = c("county_fips", "year"))

results <- list()
curve_rows <- list()
contrast_rows <- list()
coefficient_rows <- list()
vcov_rows <- list()

run_and_store <- function(key, df, label, model, reference, adjust_pollution = FALSE, allow_skip = FALSE, min_ohca = 0L, min_event_days = 0L) {
  fit <- run_case_crossover_lag0(
    df = df,
    label = label,
    model = model,
    reference = reference,
    adjust_pollution = adjust_pollution,
    allow_skip = allow_skip,
    min_ohca = min_ohca,
    min_event_days = min_event_days
  )
  results[[key]] <<- fit$result
  if (!is.null(fit$curve)) curve_rows[[key]] <<- fit$curve
  if (!is.null(fit$contrast_summary)) contrast_rows[[key]] <<- fit$contrast_summary
  if (!is.null(fit$coefficient)) coefficient_rows[[key]] <<- fit$coefficient
  if (!is.null(fit$vcov)) vcov_rows[[key]] <<- fit$vcov
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
  race_nonblack = list(label = "Non-Black", mask = referents$race_group == "Non-Black"),
  ethnicity_hispanic = list(label = "Hispanic", mask = referents$ethnicity_group == "Hispanic"),
  ethnicity_nonhispanic = list(label = "Non-Hispanic", mask = referents$ethnicity_group == "Non-Hispanic")
)

for (nm in names(strata_specs)) {
  spec <- strata_specs[[nm]]
  df <- referents[spec$mask %in% TRUE, , drop = FALSE]
  run_and_store(paste0(nm, "_primary"), df, spec$label, "stratified_humidity_adjusted", "median", allow_skip = TRUE, min_ohca = MIN_STRATUM_OHCA, min_event_days = MIN_STRATUM_EVENT_DAYS)
  run_and_store(paste0(nm, "_pollution"), df, spec$label, "stratified_humidity_pollution_adjusted", "median", adjust_pollution = TRUE, allow_skip = TRUE, min_ohca = MIN_STRATUM_OHCA, min_event_days = MIN_STRATUM_EVENT_DAYS)
  run_and_store(paste0(nm, "_mrt"), df, spec$label, "stratified_mrt_reference", "mrt", adjust_pollution = TRUE, allow_skip = TRUE, min_ohca = MIN_STRATUM_OHCA, min_event_days = MIN_STRATUM_EVENT_DAYS)
}

results_df <- bind_rows(results)
curves_df <- bind_rows(curve_rows)
contrast_df <- bind_rows(contrast_rows)
coefficient_df <- bind_rows(coefficient_rows)
vcov_df <- bind_rows(vcov_rows)

readr::write_csv(results_df, file.path(output_dir, "case_crossover_lag0_results.csv"))
readr::write_csv(curves_df, file.path(output_dir, "case_crossover_lag0_curves.csv"))
readr::write_csv(contrast_df, file.path(output_dir, "case_crossover_lag0_contrast_summaries.csv"))
readr::write_csv(coefficient_df, file.path(output_dir, "case_crossover_lag0_coefficients.csv"))
readr::write_csv(vcov_df, file.path(output_dir, "case_crossover_lag0_vcov.csv"))

print(results_df)
message("Wrote lag-0 time-stratified case-crossover results to ", output_dir)

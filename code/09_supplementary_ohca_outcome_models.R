#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) stop("Could not determine script path from commandArgs().")
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)

suppressPackageStartupMessages({
  library(arrow)
  library(splines)
})

ohca_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
tmax_path <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
rmax_path <- file.path(repo_root, "exposome", "daymet_county_rmax_2018_2024.parquet")
no2_path <- file.path(repo_root, "exposome", "no2_county_year.csv")
pm25_path <- file.path(repo_root, "exposome", "pm25_county_year.csv")
output_dir <- file.path(repo_root, "output", "final", "ohca_outcomes")
intermediate_dir <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca_outcomes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

WARM_MONTHS <- c(5L, 6L, 7L, 8L, 9L)

normalize_county_fips <- function(x) {
  clean <- gsub("\\.0$", "", as.character(x))
  clean <- gsub("[^0-9]", "", clean)
  clean <- trimws(clean)
  clean <- ifelse(nchar(clean) == 0, NA_character_, clean)
  clean <- ifelse(is.na(clean), NA_character_, sprintf("%05s", clean))
  gsub(" ", "0", clean, fixed = TRUE)
}

safe_pct <- function(x, denom) ifelse(denom > 0, 100 * x / denom, NA_real_)

extract_model_term <- function(fit, term, outcome, exposure_label, n, events) {
  coef_table <- summary(fit)$coefficients
  if (!term %in% rownames(coef_table)) {
    return(data.frame(
      outcome = outcome,
      exposure = exposure_label,
      n = n,
      events = events,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      estimable = FALSE,
      converged = fit$converged,
      note = "Term not present in fitted model.",
      stringsAsFactors = FALSE
    ))
  }
  beta <- coef_table[term, "Estimate"]
  se <- coef_table[term, "Std. Error"]
  stable <- is.finite(beta) && is.finite(se) && abs(beta) < 20 && se < 10
  data.frame(
    outcome = outcome,
    exposure = exposure_label,
    n = n,
    events = events,
    odds_ratio = ifelse(stable, exp(beta), NA_real_),
    ci_low = ifelse(stable, exp(beta - 1.96 * se), NA_real_),
    ci_high = ifelse(stable, exp(beta + 1.96 * se), NA_real_),
    p_value = ifelse(stable, coef_table[term, "Pr(>|z|)"], NA_real_),
    estimable = stable,
    converged = fit$converged,
    note = ifelse(stable, "", "Potential separation or unstable coefficient."),
    stringsAsFactors = FALSE
  )
}

fit_outcome_model <- function(df, outcome, exposure, exposure_label) {
  model_df <- df
  model_df$tmax_per_5c <- model_df$tmax_mean_c / 5
  model_df$race_group <- ifelse(model_df$race_category == "Black or African American", "Black", "Non-Black")
  model_df$sex_group <- ifelse(model_df$sex_category %in% c("Male", "Female"), model_df$sex_category, "Other/Unknown")
  time_df <- max(3L, length(unique(model_df$year)) * 4L)
  covariates <- c(
    exposure,
    "ns(rmax_mean_pct, df = 3)",
    "no2_mean",
    "pm25_mean",
    "age_at_admission",
    "sex_group",
    "race_group",
    "dow",
    "factor(year)",
    "ns(time_index, df = time_df)"
  )
  formula_spec <- as.formula(paste(outcome, "~", paste(covariates, collapse = " + ")))
  environment(formula_spec) <- environment()
  needed <- c(outcome, exposure, "rmax_mean_pct", "no2_mean", "pm25_mean", "age_at_admission", "sex_group", "race_group", "dow", "year", "time_index")
  model_df <- model_df[complete.cases(model_df[, needed, drop = FALSE]), , drop = FALSE]
  events <- sum(model_df[[outcome]] == 1)
  non_events <- sum(model_df[[outcome]] == 0)
  if (length(unique(model_df[[outcome]])) < 2 || min(events, non_events) < 20) {
    return(data.frame(
      outcome = outcome,
      exposure = exposure_label,
      n = nrow(model_df),
      events = events,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      estimable = FALSE,
      converged = FALSE,
      note = "Adjusted model not estimated because one outcome level had fewer than 20 observations.",
      stringsAsFactors = FALSE
    ))
  }
  fit <- suppressWarnings(glm(formula_spec, data = model_df, family = binomial()))
  extract_model_term(fit, exposure, outcome, exposure_label, nrow(model_df), sum(model_df[[outcome]] == 1))
}

fit_continuous_outcome_model <- function(df, outcome, exposure, exposure_label, subset_var = NULL, subset_value = NULL) {
  model_df <- df
  if (!is.null(subset_var)) {
    model_df <- model_df[model_df[[subset_var]] == subset_value, , drop = FALSE]
  }
  model_df$tmax_per_5c <- model_df$tmax_mean_c / 5
  model_df$race_group <- ifelse(model_df$race_category == "Black or African American", "Black", "Non-Black")
  model_df$sex_group <- ifelse(model_df$sex_category %in% c("Male", "Female"), model_df$sex_category, "Other/Unknown")
  model_df$outcome_log1p <- log1p(model_df[[outcome]])
  time_df <- max(3L, length(unique(model_df$year)) * 4L)
  covariates <- c(
    exposure,
    "ns(rmax_mean_pct, df = 3)",
    "no2_mean",
    "pm25_mean",
    "age_at_admission",
    "sex_group",
    "race_group",
    "dow",
    "factor(year)",
    "ns(time_index, df = time_df)"
  )
  formula_spec <- as.formula(paste("outcome_log1p ~", paste(covariates, collapse = " + ")))
  environment(formula_spec) <- environment()
  needed <- c(outcome, exposure, "rmax_mean_pct", "no2_mean", "pm25_mean", "age_at_admission", "sex_group", "race_group", "dow", "year", "time_index")
  model_df <- model_df[complete.cases(model_df[, needed, drop = FALSE]) & is.finite(model_df$outcome_log1p), , drop = FALSE]
  if (nrow(model_df) < 50) {
    return(data.frame(
      outcome = outcome,
      exposure = exposure_label,
      n = nrow(model_df),
      geometric_mean_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      estimable = FALSE,
      note = "Adjusted continuous-outcome model not estimated because fewer than 50 observations were available.",
      stringsAsFactors = FALSE
    ))
  }
  fit <- lm(formula_spec, data = model_df)
  coef_table <- summary(fit)$coefficients
  beta <- coef_table[exposure, "Estimate"]
  se <- coef_table[exposure, "Std. Error"]
  data.frame(
    outcome = outcome,
    exposure = exposure_label,
    n = nrow(model_df),
    geometric_mean_ratio = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    p_value = coef_table[exposure, "Pr(>|t|)"],
    estimable = TRUE,
    note = "Linear model of log1p outcome; estimate is a ratio on the transformed scale.",
    stringsAsFactors = FALSE
  )
}

ohca <- read.csv(ohca_path, stringsAsFactors = FALSE)
ohca$admission_date <- as.Date(ohca$admission_date)
ohca$admission_dttm <- as.POSIXct(ohca$admission_dttm, tz = "UTC")
ohca$discharge_dttm <- as.POSIXct(ohca$discharge_dttm, tz = "UTC")
ohca$age_at_admission <- as.numeric(ohca$age_at_admission)
ohca$hospital_los_days <- as.numeric(difftime(ohca$discharge_dttm, ohca$admission_dttm, units = "days"))
if (!"hospital_death" %in% names(ohca)) {
  ohca$hospital_death <- ifelse(ohca$discharge_category == "Expired", 1L, 0L)
}
if (!"death_or_hospice" %in% names(ohca)) {
  ohca$death_or_hospice <- ifelse(ohca$discharge_category == "Expired" | grepl("hospice", ohca$discharge_category, ignore.case = TRUE), 1L, 0L)
}
if (!"imv_any" %in% names(ohca)) ohca$imv_any <- NA_integer_
if (!"imv_duration_hours" %in% names(ohca)) ohca$imv_duration_hours <- NA_real_
if (!"vasopressor_any" %in% names(ohca)) ohca$vasopressor_any <- NA_integer_
ohca$county_fips <- normalize_county_fips(ohca$county_fips)
ohca$year <- as.integer(format(ohca$admission_date, "%Y"))
ohca$month <- as.integer(format(ohca$admission_date, "%m"))
ohca$dow <- factor(weekdays(ohca$admission_date), levels = c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"))
ohca$time_index <- as.integer(ohca$admission_date - min(ohca$admission_date, na.rm = TRUE)) + 1L

tmax <- arrow::read_parquet(tmax_path)
tmax$county_fips <- normalize_county_fips(tmax$geoid)
tmax$date <- as.Date(tmax$date)
tmax <- tmax[, c("county_fips", "date", "tmax_mean_c")]

rmax <- arrow::read_parquet(rmax_path)
rmax$county_fips <- normalize_county_fips(rmax$geoid)
rmax$date <- as.Date(rmax$date)
rmax <- rmax[, c("county_fips", "date", "rmax_mean_pct")]

no2 <- read.csv(no2_path, stringsAsFactors = FALSE)
no2$county_fips <- normalize_county_fips(no2$GEOID)
no2 <- no2[, c("county_fips", "year", "no2_mean")]

pm25 <- read.csv(pm25_path, stringsAsFactors = FALSE)
pm25$county_fips <- normalize_county_fips(pm25$GEOID)
pm25 <- pm25[, c("county_fips", "year", "pm25_mean")]

analysis_df <- merge(ohca, tmax, by.x = c("county_fips", "admission_date"), by.y = c("county_fips", "date"), all.x = TRUE)
analysis_df <- merge(analysis_df, rmax, by.x = c("county_fips", "admission_date"), by.y = c("county_fips", "date"), all.x = TRUE)
analysis_df <- merge(analysis_df, no2, by = c("county_fips", "year"), all.x = TRUE)
analysis_df <- merge(analysis_df, pm25, by = c("county_fips", "year"), all.x = TRUE)
analysis_df <- analysis_df[analysis_df$month %in% WARM_MONTHS & !is.na(analysis_df$tmax_mean_c), , drop = FALSE]
analysis_df$icu_los_days <- as.numeric(analysis_df$icu_los_hours) / 24
analysis_df$imv_duration_days <- as.numeric(analysis_df$imv_duration_hours) / 24

heat_threshold <- as.numeric(quantile(analysis_df$tmax_mean_c, 0.95, na.rm = TRUE, names = FALSE))
analysis_df$heat_95 <- ifelse(analysis_df$tmax_mean_c >= heat_threshold, 1L, 0L)

outcomes <- c("death_or_hospice", "hospital_death", "imv_any", "vasopressor_any")
model_rows <- list()
for (outcome in outcomes) {
  model_rows[[paste0(outcome, "_heat95")]] <- fit_outcome_model(analysis_df, outcome, "heat_95", "Tmax >= warm-season 95th percentile")
  model_rows[[paste0(outcome, "_linear")]] <- fit_outcome_model(analysis_df, outcome, "tmax_per_5c", "Tmax per 5 C")
}
model_results <- do.call(rbind, model_rows)

continuous_rows <- list(
  icu_los_heat95 = fit_continuous_outcome_model(analysis_df, "icu_los_days", "heat_95", "Tmax >= warm-season 95th percentile"),
  icu_los_linear = fit_continuous_outcome_model(analysis_df, "icu_los_days", "tmax_per_5c", "Tmax per 5 C"),
  imv_duration_heat95 = fit_continuous_outcome_model(analysis_df, "imv_duration_days", "heat_95", "Tmax >= warm-season 95th percentile", subset_var = "imv_any", subset_value = 1),
  imv_duration_linear = fit_continuous_outcome_model(analysis_df, "imv_duration_days", "tmax_per_5c", "Tmax per 5 C", subset_var = "imv_any", subset_value = 1)
)
continuous_results <- do.call(rbind, continuous_rows)

rate_rows <- do.call(rbind, lapply(outcomes, function(outcome) {
  by_heat <- aggregate(analysis_df[[outcome]], list(heat_95 = analysis_df$heat_95), function(x) c(events = sum(x == 1, na.rm = TRUE), n = sum(!is.na(x))))
  out <- do.call(data.frame, by_heat)
  names(out) <- c("heat_95", "events", "n")
  out$outcome <- outcome
  out$rate_pct <- safe_pct(out$events, out$n)
  out[, c("outcome", "heat_95", "events", "n", "rate_pct")]
}))

write.csv(analysis_df, file.path(intermediate_dir, "ohca_heat_adverse_outcome_analysis_dataset.csv"), row.names = FALSE)
write.csv(model_results, file.path(output_dir, "ohca_heat_adverse_outcome_models.csv"), row.names = FALSE)
write.csv(continuous_results, file.path(output_dir, "ohca_heat_continuous_outcome_models.csv"), row.names = FALSE)
write.csv(rate_rows, file.path(output_dir, "ohca_heat_adverse_outcome_rates.csv"), row.names = FALSE)
write.csv(data.frame(heat_threshold_tmax_c = heat_threshold), file.path(output_dir, "ohca_heat_adverse_outcome_threshold.csv"), row.names = FALSE)

message("Wrote OHCA heat adverse-outcome models to ", output_dir)

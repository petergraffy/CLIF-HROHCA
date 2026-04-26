#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) stop("Could not determine script path from commandArgs().")
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)

ensure_user_library <- function() {
  version_parts <- strsplit(as.character(getRversion()), "\\.")[[1]]
  version_stub <- paste(version_parts[1], version_parts[2], sep = ".")
  user_lib <- file.path(repo_root, ".r-user-lib", version_stub)
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}
ensure_user_library()

suppressPackageStartupMessages({
  library(dlnm)
  library(splines)
})

ohca_counts_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_daily_counts_2018_2024.csv")
ohca_cohort_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
daily_exposure_path <- file.path(repo_root, "output", "intermediate", "cohorts", "all_icu", "all_icu_daily_patient_address_tmax_2018_2024.csv")
output_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

START_DATE <- as.Date("2018-01-01")
END_DATE <- as.Date("2024-12-31")
WARM_MONTHS <- c(5L, 6L, 7L, 8L, 9L)
MAX_LAG <- 5L
TIME_DF_PER_YEAR <- 4L
VAR_DF <- 4L
LAG_DF <- 3L
RMAX_DF <- 3L

safe_quantile <- function(x, prob) stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7)

make_prediction_grid <- function(x) {
  unique(as.numeric(seq(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)), by = 0.5)))
}

normalize_county_fips <- function(x) {
  clean <- gsub("\\.0$", "", as.character(x))
  clean <- gsub("[^0-9]", "", clean)
  clean <- trimws(clean)
  clean <- ifelse(nchar(clean) == 0, NA_character_, clean)
  clean <- ifelse(is.na(clean), NA_character_, sprintf("%05s", clean))
  gsub(" ", "0", clean, fixed = TRUE)
}

build_formula <- function(extra_terms = character(), include_dow = TRUE, include_year = TRUE) {
  rhs <- c("cb_temp", "ns(time_index, df = time_df)")
  if (include_dow) rhs <- c(rhs, "dow")
  if (include_year) rhs <- c(rhs, "factor(year)")
  rhs <- c(rhs, extra_terms)
  as.formula(paste("ohca_admissions ~", paste(rhs, collapse = " + ")))
}

run_dlnm_spec <- function(
  df,
  label,
  model = "primary_humidity_adjusted",
  extra_terms = c("ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)"),
  reference = c("median", "mrt"),
  time_df_per_year = TIME_DF_PER_YEAR,
  include_dow = TRUE,
  include_year = TRUE,
  return_curve = FALSE
) {
  reference <- match.arg(reference)
  needed <- c("icu_patient_address_mean_tmax_c", "icu_patient_address_mean_rmax_pct")
  if (any(grepl("no2", extra_terms, fixed = TRUE))) needed <- c(needed, "icu_patient_address_mean_no2")
  if (any(grepl("pm25", extra_terms, fixed = TRUE))) needed <- c(needed, "icu_patient_address_mean_pm25")
  df <- df[complete.cases(df[, needed, drop = FALSE]), , drop = FALSE]
  df$time_index <- seq_len(nrow(df))
  cb_temp <- crossbasis(df$icu_patient_address_mean_tmax_c, lag = MAX_LAG, argvar = list(fun = "ns", df = VAR_DF), arglag = list(fun = "ns", df = LAG_DF))
  formula_spec <- build_formula(extra_terms, include_dow = include_dow, include_year = include_year)
  environment(formula_spec) <- environment()
  environment(formula_spec)$time_df <- length(unique(df$year)) * time_df_per_year
  fit <- glm(formula_spec, data = df, family = quasipoisson(link = "log"))
  dispersion <- summary(fit)$dispersion
  grid <- make_prediction_grid(df$icu_patient_address_mean_tmax_c)
  initial_center <- median(df$icu_patient_address_mean_tmax_c, na.rm = TRUE)
  pred_initial <- crosspred(cb_temp, fit, cen = initial_center, at = grid)
  center <- if (reference == "median") initial_center else grid[which.min(pred_initial$allRRfit)]
  pred <- crosspred(cb_temp, fit, cen = center, at = grid)
  reduced <- crossreduce(cb_temp, fit, cen = center)
  hot_temp <- grid[which.min(abs(grid - safe_quantile(df$icu_patient_address_mean_tmax_c, 0.95)))]
  hot_index <- which.min(abs(grid - hot_temp))
  result <- data.frame(
    stratum = label,
    model = model,
    n_days = nrow(df),
    n_ohca = sum(df$ohca_admissions),
    reference_type = reference,
    reference_temp_c = center,
    hot_temp_c = hot_temp,
    cumulative_rr = as.numeric(pred$allRRfit[hot_index]),
    cumulative_rr_low = as.numeric(pred$allRRlow[hot_index]),
    cumulative_rr_high = as.numeric(pred$allRRhigh[hot_index]),
    log_rr = log(as.numeric(pred$allRRfit[hot_index])),
    log_rr_se = (log(as.numeric(pred$allRRhigh[hot_index])) - log(as.numeric(pred$allRRlow[hot_index]))) / (2 * 1.96),
    aic = NA_real_,
    dispersion = dispersion,
    model_family = "quasipoisson",
    time_df_per_year = time_df_per_year,
    includes_day_of_week = include_dow,
    includes_year_fixed_effect = include_year,
    stringsAsFactors = FALSE
  )

  curve <- data.frame(
    stratum = label,
    model = model,
    reference_type = reference,
    reference_temp_c = center,
    tmax_mean_c = grid,
    cumulative_rr = as.numeric(pred$allRRfit),
    cumulative_rr_low = as.numeric(pred$allRRlow),
    cumulative_rr_high = as.numeric(pred$allRRhigh),
    log_rr = log(as.numeric(pred$allRRfit)),
    log_rr_se = (log(as.numeric(pred$allRRhigh)) - log(as.numeric(pred$allRRlow))) / (2 * 1.96),
    stringsAsFactors = FALSE
  )

  reduced_coef <- data.frame(
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
  reduced_vcov$stratum <- label
  reduced_vcov$model <- model
  reduced_vcov$reference_type <- reference
  reduced_vcov$reference_temp_c <- center
  reduced_vcov <- reduced_vcov[, c("stratum", "model", "reference_type", "reference_temp_c", "coefficient_row", "coefficient_col", "covariance")]

  if (return_curve) {
    return(list(result = result, curve = curve, reduced_coef = reduced_coef, reduced_vcov = reduced_vcov))
  }

  result
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
model_df$dow <- factor(weekdays(model_df$admission_date), levels = c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"))
model_df <- model_df[model_df$month %in% WARM_MONTHS & !is.na(model_df$icu_patient_address_mean_tmax_c), , drop = FALSE]

cohort <- read.csv(ohca_cohort_path, stringsAsFactors = FALSE)
cohort$admission_date <- as.Date(cohort$admission_date)
cohort$age_group <- ifelse(as.numeric(cohort$age_at_admission) >= 65, ">=65", "<65")
cohort$race_group <- ifelse(cohort$race_category == "Black or African American", "Black", "Non-Black")

make_stratum_daily_counts <- function(subset_df, label) {
  out <- aggregate(hospitalization_id ~ admission_date, data = subset_df, FUN = function(x) length(unique(x)))
  names(out)[2] <- "ohca_admissions"
  out$label <- label
  out
}

strata_counts <- list(
  overall = daily_counts,
  male = make_stratum_daily_counts(cohort[cohort$sex_category == "Male", ], "Male"),
  female = make_stratum_daily_counts(cohort[cohort$sex_category == "Female", ], "Female"),
  age_lt65 = make_stratum_daily_counts(cohort[cohort$age_group == "<65", ], "<65"),
  age_ge65 = make_stratum_daily_counts(cohort[cohort$age_group == ">=65", ], ">=65"),
  race_black = make_stratum_daily_counts(cohort[cohort$race_group == "Black", ], "Black"),
  race_nonblack = make_stratum_daily_counts(cohort[cohort$race_group == "Non-Black", ], "Non-Black")
)

results <- list()
curve_rows <- list()
reduced_coef_rows <- list()
reduced_vcov_rows <- list()

overall_primary <- run_dlnm_spec(model_df, "Overall", model = "primary_humidity_adjusted", reference = "median", return_curve = TRUE)
results[["overall_primary"]] <- overall_primary$result
curve_rows[["overall_primary"]] <- overall_primary$curve
reduced_coef_rows[["overall_primary"]] <- overall_primary$reduced_coef
reduced_vcov_rows[["overall_primary"]] <- overall_primary$reduced_vcov

overall_pollution <- run_dlnm_spec(
  model_df,
  "Overall",
  model = "sensitivity_humidity_pollution_adjusted",
  extra_terms = c("ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)", "icu_patient_address_mean_no2", "icu_patient_address_mean_pm25"),
  reference = "median",
  return_curve = TRUE
)
results[["overall_pollution"]] <- overall_pollution$result
curve_rows[["overall_pollution"]] <- overall_pollution$curve
reduced_coef_rows[["overall_pollution"]] <- overall_pollution$reduced_coef
reduced_vcov_rows[["overall_pollution"]] <- overall_pollution$reduced_vcov

overall_mrt <- run_dlnm_spec(model_df, "Overall", model = "sensitivity_mrt_reference", reference = "mrt", return_curve = TRUE)
results[["overall_mrt"]] <- overall_mrt$result
curve_rows[["overall_mrt"]] <- overall_mrt$curve
reduced_coef_rows[["overall_mrt"]] <- overall_mrt$reduced_coef
reduced_vcov_rows[["overall_mrt"]] <- overall_mrt$reduced_vcov

time_sensitivity_results <- list()
for (time_df_candidate in c(3L, 4L, 6L)) {
  time_sensitivity_results[[paste0("time_df_", time_df_candidate)]] <- run_dlnm_spec(
    model_df,
    "Overall",
    model = paste0("sensitivity_time_df_", time_df_candidate, "_per_year"),
    reference = "median",
    time_df_per_year = time_df_candidate
  )
}
time_sensitivity_results[["without_day_of_week"]] <- run_dlnm_spec(
  model_df,
  "Overall",
  model = "sensitivity_without_day_of_week",
  reference = "median",
  include_dow = FALSE
)
time_sensitivity_df <- do.call(rbind, time_sensitivity_results)

for (nm in c("male","female","age_lt65","age_ge65","race_black","race_nonblack")) {
  counts_df <- strata_counts[[nm]]
  counts_df$admission_date <- as.Date(counts_df$admission_date)
  merged <- merge(all_dates, counts_df[, c("admission_date","ohca_admissions")], by = "admission_date", all.x = TRUE, sort = TRUE)
  merged$ohca_admissions[is.na(merged$ohca_admissions)] <- 0L
  merged <- merge(merged, daily_exposure, by = "admission_date", all.x = TRUE, sort = TRUE)
  merged$year <- as.integer(format(merged$admission_date, "%Y"))
  merged$month <- as.integer(format(merged$admission_date, "%m"))
  merged$dow <- factor(weekdays(merged$admission_date), levels = c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"))
  merged <- merged[merged$month %in% WARM_MONTHS & !is.na(merged$icu_patient_address_mean_tmax_c), , drop = FALSE]
  label <- unique(strata_counts[[nm]]$label)[1]
  stratified <- run_dlnm_spec(
    merged,
    label,
    model = "stratified_humidity_adjusted",
    reference = "median",
    return_curve = TRUE
  )
  results[[paste0(nm, "_primary")]] <- stratified$result
  curve_rows[[paste0(nm, "_primary")]] <- stratified$curve
  reduced_coef_rows[[paste0(nm, "_primary")]] <- stratified$reduced_coef
  reduced_vcov_rows[[paste0(nm, "_primary")]] <- stratified$reduced_vcov
}

results_df <- do.call(rbind, results)
curves_df <- do.call(rbind, curve_rows)
reduced_coef_df <- do.call(rbind, reduced_coef_rows)
reduced_vcov_df <- do.call(rbind, reduced_vcov_rows)
write.csv(results_df, file.path(output_dir, "manuscript_dlnm_results.csv"), row.names = FALSE)
write.csv(curves_df, file.path(output_dir, "manuscript_dlnm_curves.csv"), row.names = FALSE)
write.csv(reduced_coef_df, file.path(output_dir, "manuscript_dlnm_reduced_coefficients.csv"), row.names = FALSE)
write.csv(reduced_vcov_df, file.path(output_dir, "manuscript_dlnm_reduced_vcov.csv"), row.names = FALSE)
write.csv(time_sensitivity_df, file.path(output_dir, "manuscript_dlnm_time_adjustment_sensitivity.csv"), row.names = FALSE)
message("Wrote manuscript-style DLNM results to ", output_dir)

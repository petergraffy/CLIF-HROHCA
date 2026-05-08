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
MIN_STRATUM_OHCA <- 30L
MIN_STRATUM_EVENT_DAYS <- 10L
MIN_FALLBACK_OHCA <- 10L
FALLBACK_TIME_DF_PER_YEAR <- 1L

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

finite_unique_n <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  length(unique(x[is.finite(x)]))
}

choose_spline_term <- function(df, variable, spline_term, requested_df, linear_term = variable) {
  n_unique <- finite_unique_n(df[[variable]])
  if (n_unique >= requested_df + 1L) {
    return(list(term = spline_term, note = NA_character_))
  }
  if (n_unique >= 2L) {
    return(list(term = linear_term, note = paste0(variable, "_spline_reduced_to_linear")))
  }
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
  if (n_days < 14L) {
    return(list(term = NA_character_, df = NA_integer_, note = "time_adjustment_omitted_too_few_days"))
  }

  capped_df <- min(as.integer(requested_df), max(1L, n_days - 2L))
  if (capped_df >= 2L) {
    note <- if (capped_df < requested_df) paste0("time_df_capped_from_", requested_df, "_to_", capped_df) else NA_character_
    return(list(term = "ns(time_index, df = time_df)", df = capped_df, note = note))
  }

  list(term = "time_index", df = NA_integer_, note = "time_spline_reduced_to_linear")
}

build_formula <- function(extra_terms = character(), time_term = "ns(time_index, df = time_df)", include_dow = TRUE, include_year = TRUE) {
  rhs <- c("cb_temp")
  if (!is.na(time_term)) rhs <- c(rhs, time_term)
  if (include_dow) rhs <- c(rhs, "dow")
  if (include_year) rhs <- c(rhs, "factor(year)")
  rhs <- c(rhs, extra_terms)
  as.formula(paste("ohca_admissions ~", paste(rhs, collapse = " + ")))
}

make_nonestimable_dlnm_result <- function(df, label, model, reference, time_df_per_year, include_dow, include_year, reason) {
  n_days <- nrow(df)
  n_ohca <- if ("ohca_admissions" %in% names(df)) sum(df$ohca_admissions, na.rm = TRUE) else 0L
  data.frame(
    stratum = label,
    model = model,
    n_days = n_days,
    n_ohca = n_ohca,
    reference_type = reference,
    reference_temp_c = NA_real_,
    hot_temp_c = NA_real_,
    cumulative_rr = NA_real_,
    cumulative_rr_low = NA_real_,
    cumulative_rr_high = NA_real_,
    log_rr = NA_real_,
    log_rr_se = NA_real_,
    aic = NA_real_,
    dispersion = NA_real_,
    model_family = "quasipoisson",
    time_df_per_year = time_df_per_year,
    includes_day_of_week = include_dow,
    includes_year_fixed_effect = include_year,
    estimable = FALSE,
    converged = FALSE,
    model_type = "not_estimable",
    fallback_reason = NA_character_,
    skip_reason = reason,
    stringsAsFactors = FALSE
  )
}

run_linear_fallback_spec <- function(
  df,
  label,
  model,
  reference,
  time_df_per_year,
  include_dow,
  reason
) {
  if (sum(df$ohca_admissions, na.rm = TRUE) < MIN_FALLBACK_OHCA) {
    return(make_nonestimable_dlnm_result(
      df, label, model, reference, time_df_per_year, include_dow, FALSE,
      paste0(reason, "; linear fallback not fit because fewer than ", MIN_FALLBACK_OHCA, " OHCA admissions")
    ))
  }

  fallback_model <- paste0(model, "_linear_fallback")
  df$time_index <- seq_len(nrow(df))
  df$tmax_per_5c <- df$icu_patient_address_mean_tmax_c / 5
  fallback_time_df <- max(1L, length(unique(df$year)) * FALLBACK_TIME_DF_PER_YEAR)
  formula_spec <- as.formula(
    "ohca_admissions ~ tmax_per_5c + ns(time_index, df = fallback_time_df) + dow + ns(icu_patient_address_mean_rmax_pct, df = RMAX_DF)"
  )
  environment(formula_spec) <- environment()

  fit <- tryCatch(
    glm(formula_spec, data = df, family = quasipoisson(link = "log")),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(make_nonestimable_dlnm_result(
      df, label, fallback_model, reference, FALLBACK_TIME_DF_PER_YEAR, include_dow, FALSE,
      paste0(reason, "; linear fallback failed: ", conditionMessage(fit))
    ))
  }

  beta <- unname(coef(fit)[["tmax_per_5c"]])
  beta_var <- vcov(fit)["tmax_per_5c", "tmax_per_5c"]
  beta_se <- sqrt(beta_var)
  stable <- isTRUE(fit$converged) && is.finite(beta) && is.finite(beta_se) && abs(beta) < 20 && beta_se < 10
  if (!stable) {
    return(make_nonestimable_dlnm_result(
      df, label, fallback_model, reference, FALLBACK_TIME_DF_PER_YEAR, include_dow, FALSE,
      paste0(reason, "; linear fallback did not converge or produced an unstable coefficient")
    ))
  }

  center <- median(df$icu_patient_address_mean_tmax_c, na.rm = TRUE)
  hot_temp <- safe_quantile(df$icu_patient_address_mean_tmax_c, 0.95)
  contrast <- (hot_temp - center) / 5
  log_rr <- beta * contrast
  log_rr_se <- sqrt(beta_var) * abs(contrast)

  data.frame(
    stratum = label,
    model = fallback_model,
    n_days = nrow(df),
    n_ohca = sum(df$ohca_admissions, na.rm = TRUE),
    reference_type = reference,
    reference_temp_c = center,
    hot_temp_c = hot_temp,
    cumulative_rr = exp(log_rr),
    cumulative_rr_low = exp(log_rr - 1.96 * log_rr_se),
    cumulative_rr_high = exp(log_rr + 1.96 * log_rr_se),
    log_rr = log_rr,
    log_rr_se = log_rr_se,
    aic = NA_real_,
    dispersion = summary(fit)$dispersion,
    model_family = "quasipoisson",
    time_df_per_year = FALLBACK_TIME_DF_PER_YEAR,
    includes_day_of_week = include_dow,
    includes_year_fixed_effect = FALSE,
    estimable = TRUE,
    converged = fit$converged,
    model_type = "linear_fallback",
    fallback_reason = reason,
    skip_reason = NA_character_,
    stringsAsFactors = FALSE
  )
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
  return_curve = FALSE,
  allow_skip = FALSE,
  min_ohca = 0L,
  min_event_days = 0L,
  fallback_linear = FALSE
) {
  reference <- match.arg(reference)
  needed <- c("icu_patient_address_mean_tmax_c", "icu_patient_address_mean_rmax_pct")
  if (any(grepl("no2", extra_terms, fixed = TRUE))) needed <- c(needed, "icu_patient_address_mean_no2")
  if (any(grepl("pm25", extra_terms, fixed = TRUE))) needed <- c(needed, "icu_patient_address_mean_pm25")
  df <- df[complete.cases(df[, needed, drop = FALSE]), , drop = FALSE]
  for (needed_col in needed) {
    df <- df[is.finite(suppressWarnings(as.numeric(df[[needed_col]]))), , drop = FALSE]
  }

  n_ohca <- sum(df$ohca_admissions, na.rm = TRUE)
  n_event_days <- sum(df$ohca_admissions > 0, na.rm = TRUE)
  skip_reason <- NULL
  if (nrow(df) <= MAX_LAG + 14L) {
    skip_reason <- paste0("Only ", nrow(df), " complete exposure days after filtering")
  } else if (n_ohca == 0L) {
    skip_reason <- "No OHCA admissions available"
  } else if (n_ohca < min_ohca) {
    skip_reason <- paste0("Only ", n_ohca, " OHCA admissions; minimum for this model is ", min_ohca)
  } else if (n_event_days < min_event_days) {
    skip_reason <- paste0("Only ", n_event_days, " event days; minimum for this model is ", min_event_days)
  } else if (finite_unique_n(df$icu_patient_address_mean_tmax_c) < 2L) {
    skip_reason <- "Temperature exposure is constant or unavailable"
  }
  if (!is.null(skip_reason)) {
    if (!allow_skip) stop(skip_reason, call. = FALSE)
    if (fallback_linear) {
      message("Falling back to linear model for DLNM model '", model, "' in stratum '", label, "': ", skip_reason)
      fallback <- run_linear_fallback_spec(df, label, model, reference, time_df_per_year, include_dow, skip_reason)
      if (return_curve) return(list(result = fallback, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
      return(fallback)
    } else {
      message("Skipping DLNM model '", model, "' in stratum '", label, "': ", skip_reason)
      skipped <- make_nonestimable_dlnm_result(df, label, model, reference, time_df_per_year, include_dow, include_year, skip_reason)
      if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
      return(skipped)
    }
  }

  df$time_index <- seq_len(nrow(df))
  temp_unique <- finite_unique_n(df$icu_patient_address_mean_tmax_c)
  temp_var_df <- min(VAR_DF, max(2L, temp_unique - 1L))
  if (temp_var_df < VAR_DF) {
    message("Reducing tmax spline df from ", VAR_DF, " to ", temp_var_df, " for ", label, " / ", model, ".")
  }

  cb_temp <- tryCatch(
    crossbasis(
      df$icu_patient_address_mean_tmax_c,
      lag = MAX_LAG,
      argvar = list(fun = "ns", df = temp_var_df),
      arglag = list(fun = "ns", df = LAG_DF)
    ),
    error = function(e) {
      if (!allow_skip) stop(e)
      e
    }
  )
  if (inherits(cb_temp, "error")) {
    reason <- paste0("DLNM crossbasis failed: ", conditionMessage(cb_temp))
    if (fallback_linear) {
      message("Falling back to linear model for DLNM model '", model, "' in stratum '", label, "': ", reason)
      fallback <- run_linear_fallback_spec(df, label, model, reference, time_df_per_year, include_dow, reason)
      if (return_curve) return(list(result = fallback, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
      return(fallback)
    }
    skipped <- make_nonestimable_dlnm_result(df, label, model, reference, time_df_per_year, include_dow, include_year, reason)
    message("Skipping DLNM model '", model, "' in stratum '", label, "': ", skipped$skip_reason)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
    return(skipped)
  }

  extra_info <- sanitize_extra_terms(df, extra_terms)
  if (length(extra_info$notes) > 0L) {
    message("Adjusted covariate terms for ", label, " / ", model, ": ", paste(extra_info$notes, collapse = "; "), ".")
  }

  requested_time_df <- length(unique(df$year)) * time_df_per_year
  time_info <- choose_time_adjustment(df, requested_time_df)
  if (!is.na(time_info$note)) {
    message("Adjusted time term for ", label, " / ", model, ": ", time_info$note, ".")
  }

  model_include_dow <- include_dow && length(unique(stats::na.omit(df$dow))) > 1L
  model_include_year <- include_year && length(unique(stats::na.omit(df$year))) > 1L
  if (include_dow && !model_include_dow) message("Omitting day-of-week term for ", label, " / ", model, " because it has fewer than 2 levels.")
  if (include_year && !model_include_year) message("Omitting year fixed effects for ", label, " / ", model, " because it has fewer than 2 levels.")

  formula_spec <- build_formula(
    extra_info$terms,
    time_term = time_info$term,
    include_dow = model_include_dow,
    include_year = model_include_year
  )
  environment(formula_spec) <- environment()
  environment(formula_spec)$time_df <- time_info$df
  design <- model.matrix(formula_spec, data = df)
  if (any(!is.finite(design))) {
    bad_cols <- names(which(colSums(!is.finite(design)) > 0L))
    reason <- paste0("Non-finite values found in DLNM design matrix. Columns: ", paste(bad_cols, collapse = ", "))
    if (!allow_skip) stop(reason, call. = FALSE)
    if (fallback_linear) {
      message("Falling back to linear model for DLNM model '", model, "' in stratum '", label, "': ", reason)
      fallback <- run_linear_fallback_spec(df, label, model, reference, time_df_per_year, include_dow, reason)
      if (return_curve) return(list(result = fallback, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
      return(fallback)
    }
    skipped <- make_nonestimable_dlnm_result(df, label, model, reference, time_df_per_year, include_dow, include_year, reason)
    message("Skipping DLNM model '", model, "' in stratum '", label, "': ", skipped$skip_reason)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
    return(skipped)
  }

  fit <- tryCatch(
    glm(formula_spec, data = df, family = quasipoisson(link = "log")),
    error = function(e) {
      if (!allow_skip) stop(e)
      e
    }
  )
  if (inherits(fit, "error")) {
    reason <- paste0("DLNM fit failed: ", conditionMessage(fit))
    if (fallback_linear) {
      message("Falling back to linear model for DLNM model '", model, "' in stratum '", label, "': ", reason)
      fallback <- run_linear_fallback_spec(df, label, model, reference, time_df_per_year, include_dow, reason)
      if (return_curve) return(list(result = fallback, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
      return(fallback)
    }
    skipped <- make_nonestimable_dlnm_result(df, label, model, reference, time_df_per_year, include_dow, include_year, reason)
    message("Skipping DLNM model '", model, "' in stratum '", label, "': ", skipped$skip_reason)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
    return(skipped)
  }
  dispersion <- summary(fit)$dispersion
  grid <- make_prediction_grid(df$icu_patient_address_mean_tmax_c)
  initial_center <- median(df$icu_patient_address_mean_tmax_c, na.rm = TRUE)
  prediction_objects <- tryCatch({
    pred_initial <- crosspred(cb_temp, fit, cen = initial_center, at = grid)
    center <- if (reference == "median") initial_center else grid[which.min(pred_initial$allRRfit)]
    list(
      center = center,
      pred = crosspred(cb_temp, fit, cen = center, at = grid),
      reduced = crossreduce(cb_temp, fit, cen = center)
    )
  }, error = function(e) {
    if (!allow_skip) stop(e)
    e
  })
  if (inherits(prediction_objects, "error")) {
    reason <- paste0("DLNM prediction failed: ", conditionMessage(prediction_objects))
    if (fallback_linear) {
      message("Falling back to linear model for DLNM model '", model, "' in stratum '", label, "': ", reason)
      fallback <- run_linear_fallback_spec(df, label, model, reference, time_df_per_year, include_dow, reason)
      if (return_curve) return(list(result = fallback, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
      return(fallback)
    }
    skipped <- make_nonestimable_dlnm_result(df, label, model, reference, time_df_per_year, include_dow, include_year, reason)
    message("Skipping DLNM model '", model, "' in stratum '", label, "': ", skipped$skip_reason)
    if (return_curve) return(list(result = skipped, curve = NULL, reduced_coef = NULL, reduced_vcov = NULL))
    return(skipped)
  }
  center <- prediction_objects$center
  pred <- prediction_objects$pred
  reduced <- prediction_objects$reduced
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
    includes_day_of_week = model_include_dow,
    includes_year_fixed_effect = model_include_year,
    estimable = TRUE,
    converged = fit$converged,
    model_type = "dlnm",
    fallback_reason = NA_character_,
    skip_reason = NA_character_,
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
cohort$race_group <- ifelse(is_black_race(cohort$race_category), "Black", "Non-Black")

make_stratum_daily_counts <- function(subset_df, label) {
  out <- aggregate(hospitalization_id ~ admission_date, data = subset_df, FUN = function(x) length(unique(x)))
  names(out)[2] <- "ohca_admissions"
  out$label <- label
  out
}

strata_counts <- list(
  overall = daily_counts,
  male = make_stratum_daily_counts(cohort[is_male(cohort$sex_category), ], "Male"),
  female = make_stratum_daily_counts(cohort[is_female(cohort$sex_category), ], "Female"),
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
    return_curve = TRUE,
    allow_skip = TRUE,
    min_ohca = MIN_STRATUM_OHCA,
    min_event_days = MIN_STRATUM_EVENT_DAYS,
    fallback_linear = TRUE
  )
  results[[paste0(nm, "_primary")]] <- stratified$result
  if (!is.null(stratified$curve)) curve_rows[[paste0(nm, "_primary")]] <- stratified$curve
  if (!is.null(stratified$reduced_coef)) reduced_coef_rows[[paste0(nm, "_primary")]] <- stratified$reduced_coef
  if (!is.null(stratified$reduced_vcov)) reduced_vcov_rows[[paste0(nm, "_primary")]] <- stratified$reduced_vcov
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

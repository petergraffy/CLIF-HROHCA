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
  library(cmprsk)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(readr)
  library(splines)
  library(stringr)
  library(survival)
  library(tidyr)
})

OHCA_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
TMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
RMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_rmax_2018_2024.parquet")
OUTPUT_DIR <- file.path(repo_root, "output", "final", "ohca_icu_phenotypes")
FIGURE_DIR <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(OUTPUT_DIR, c(
  "ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv",
  "ohca_icu_competing_risk_awake_extubated_cif_curves.csv",
  "ohca_icu_competing_risk_awake_extubated_fine_gray_models.csv",
  "ohca_icu_competing_risk_awake_extubated_summary.csv",
  "ohca_icu_competing_risk_awake_extubated_patient_audit.csv",
  "ohca_icu_competing_risk_extubation_cif_curves.csv",
  "ohca_icu_competing_risk_extubation_fine_gray_models.csv",
  "ohca_icu_competing_risk_extubation_summary.csv",
  "ohca_icu_competing_risk_extubation_patient_audit.csv"
)))
FOLLOWUP_HOURS <- 72

config <- load_project_config(repo_root)
validate_site_name(config, repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)
site_name <- config$site_name %||% Sys.getenv("CLIF_SITE_NAME", unset = "unknown_site")

safe_last <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  tail(x, 1)
}

clean_assessment_text <- function(x) {
  stringr::str_to_lower(stringr::str_squish(tidyr::replace_na(as.character(x), "")))
}

ensure_columns <- function(df, columns) {
  if (is.null(df)) df <- tibble::tibble()
  for (column in columns) {
    if (!column %in% names(df)) df[[column]] <- NA
  }
  df
}

choose_terms <- function(df, terms) {
  keep <- character()
  omitted <- character()
  for (term in terms) {
    if (!term %in% names(df)) {
      omitted <- c(omitted, paste0(term, "_missing"))
      next
    }
    x <- df[[term]]
    if (is.numeric(x)) {
      if (sum(is.finite(x)) > 0 && length(unique(x[is.finite(x)])) >= 2) {
        keep <- c(keep, term)
      } else {
        omitted <- c(omitted, paste0(term, "_constant"))
      }
    } else {
      if (length(unique(stats::na.omit(x))) >= 2) {
        keep <- c(keep, term)
      } else {
        omitted <- c(omitted, paste0(term, "_one_level"))
      }
    }
  }
  list(keep = keep, omitted = omitted)
}

make_model_matrix <- function(df, terms) {
  if (length(terms) == 0) return(NULL)
  stats::model.matrix(stats::as.formula(paste("~", paste(terms, collapse = " + "))), data = df)[, -1, drop = FALSE]
}

sanitize_covariate_matrix <- function(mat) {
  if (is.null(mat) || ncol(mat) == 0) return(list(matrix = mat, dropped = character()))
  dropped <- character()
  finite_cols <- vapply(seq_len(ncol(mat)), function(i) all(is.finite(mat[, i])), logical(1))
  if (any(!finite_cols)) dropped <- c(dropped, paste0(colnames(mat)[!finite_cols], "_nonfinite"))
  mat <- mat[, finite_cols, drop = FALSE]
  if (ncol(mat) == 0) return(list(matrix = mat, dropped = dropped))

  varying_cols <- vapply(seq_len(ncol(mat)), function(i) length(unique(mat[, i])) >= 2, logical(1))
  if (any(!varying_cols)) dropped <- c(dropped, paste0(colnames(mat)[!varying_cols], "_constant"))
  mat <- mat[, varying_cols, drop = FALSE]
  if (ncol(mat) == 0) return(list(matrix = mat, dropped = dropped))

  qr_mat <- qr(mat)
  if (qr_mat$rank < ncol(mat)) {
    independent <- sort(qr_mat$pivot[seq_len(qr_mat$rank)])
    dropped <- c(dropped, paste0(colnames(mat)[setdiff(seq_len(ncol(mat)), independent)], "_collinear"))
    mat <- mat[, independent, drop = FALSE]
  }
  list(matrix = mat, dropped = dropped)
}

wald_p <- function(beta, vcov) {
  if (length(beta) == 0 || any(!is.finite(beta)) || any(!is.finite(vcov))) return(NA_real_)
  inv <- tryCatch(solve(vcov), error = function(e) NULL)
  if (is.null(inv)) return(NA_real_)
  stat <- as.numeric(t(beta) %*% inv %*% beta)
  stats::pchisq(stat, df = length(beta), lower.tail = FALSE)
}

add_lag_window_means <- function(df, value_col, output_prefix) {
  value <- rlang::ensym(value_col)
  daily <- df |>
    dplyr::group_by(.data$county_fips, .data$admission_date) |>
    dplyr::summarise(!!rlang::as_name(value) := mean(!!value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(!!rlang::as_name(value) := ifelse(is.nan(!!value), NA_real_, !!value))
  daily |>
    dplyr::group_by(.data$county_fips) |>
    dplyr::arrange(.data$admission_date, .by_group = TRUE) |>
    dplyr::mutate(
      "{output_prefix}_lag0_1" := rowMeans(cbind(!!value, dplyr::lag(!!value, 1L)), na.rm = FALSE),
      "{output_prefix}_lag0_3" := rowMeans(cbind(!!value, dplyr::lag(!!value, 1L), dplyr::lag(!!value, 2L), dplyr::lag(!!value, 3L)), na.rm = FALSE),
      "{output_prefix}_lag0_5" := rowMeans(cbind(!!value, dplyr::lag(!!value, 1L), dplyr::lag(!!value, 2L), dplyr::lag(!!value, 3L), dplyr::lag(!!value, 4L), dplyr::lag(!!value, 5L)), na.rm = FALSE)
    ) |>
    dplyr::ungroup()
}

extract_crr_terms <- function(fit, model, exposure_label, exposure_regex, n, events, competing, covariates, omitted_terms, estimable = TRUE, skip_reason = NA_character_) {
  if (!estimable || inherits(fit, "error")) {
    return(tibble::tibble(
      model = model,
      exposure = exposure_label,
      term = exposure_label,
      n = n,
      events_awake_extubated = events,
      competing_deaths = competing,
      subdistribution_hr = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      overall_exposure_p_value = NA_real_,
      covariates = paste(covariates, collapse = " + "),
      omitted_terms = paste(omitted_terms, collapse = "; "),
      estimable = FALSE,
      skip_reason = if (inherits(fit, "error")) conditionMessage(fit) else skip_reason
    ))
  }

  coefs <- fit$coef
  vc <- fit$var
  se <- sqrt(diag(vc))
  keep <- stringr::str_detect(names(coefs), exposure_regex)
  if (!any(keep)) {
    return(tibble::tibble(
      model = model,
      exposure = exposure_label,
      term = exposure_label,
      n = n,
      events_awake_extubated = events,
      competing_deaths = competing,
      subdistribution_hr = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      overall_exposure_p_value = NA_real_,
      covariates = paste(covariates, collapse = " + "),
      omitted_terms = paste(omitted_terms, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Exposure term not present in fitted Fine-Gray model"
    ))
  }

  overall <- wald_p(coefs[keep], vc[keep, keep, drop = FALSE])
  tibble::tibble(
    model = model,
    exposure = exposure_label,
    term = names(coefs)[keep],
    n = n,
    events_awake_extubated = events,
    competing_deaths = competing,
    subdistribution_hr = exp(coefs[keep]),
    ci_low = exp(coefs[keep] - 1.96 * se[keep]),
    ci_high = exp(coefs[keep] + 1.96 * se[keep]),
    p_value = 2 * stats::pnorm(abs(coefs[keep] / se[keep]), lower.tail = FALSE),
    overall_exposure_p_value = overall,
    covariates = paste(covariates, collapse = " + "),
    omitted_terms = paste(omitted_terms, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
}

extract_crr_coefficients <- function(fit, model, n, events, competing, covariates, omitted_terms) {
  if (inherits(fit, "error") || is.null(fit$coef) || is.null(fit$var)) return(tibble::tibble())
  coefs <- fit$coef
  vc <- fit$var
  se <- sqrt(diag(vc))
  tibble::tibble(
    model = model,
    n = n,
    events_awake_extubated = events,
    competing_deaths = competing,
    coefficient = names(coefs),
    estimate = as.numeric(coefs),
    standard_error = as.numeric(se),
    subdistribution_hr = exp(as.numeric(coefs)),
    ci_low = exp(as.numeric(coefs) - 1.96 * as.numeric(se)),
    ci_high = exp(as.numeric(coefs) + 1.96 * as.numeric(se)),
    p_value = 2 * stats::pnorm(abs(as.numeric(coefs) / as.numeric(se)), lower.tail = FALSE),
    covariates = paste(covariates, collapse = " + "),
    omitted_terms = paste(omitted_terms, collapse = "; ")
  )
}

extract_crr_vcov <- function(fit, model, n, events, competing, covariates, omitted_terms) {
  if (inherits(fit, "error") || is.null(fit$coef) || is.null(fit$var)) return(tibble::tibble())
  vc <- fit$var
  if (is.null(rownames(vc))) rownames(vc) <- names(fit$coef)
  if (is.null(colnames(vc))) colnames(vc) <- names(fit$coef)
  as.data.frame(as.table(vc), stringsAsFactors = FALSE) |>
    transmute(
      model = model,
      n = n,
      events_awake_extubated = events,
      competing_deaths = competing,
      coefficient_row = as.character(.data$Var1),
      coefficient_col = as.character(.data$Var2),
      covariance = as.numeric(.data$Freq),
      covariates = paste(covariates, collapse = " + "),
      omitted_terms = paste(omitted_terms, collapse = "; ")
    )
}

extract_cif_curve <- function(df, model = "awake_extubated_72h_cif", time_grid_hours = 0:FOLLOWUP_HOURS,
                              stratification = "Overall", stratum = "Overall", stratum_order = 0,
                              stratum_min_temp_c = NA_real_, stratum_max_temp_c = NA_real_) {
  base <- df |>
    filter(!is.na(.data$fg_time_hours), is.finite(.data$fg_time_hours), .data$fg_time_hours >= 0, !is.na(.data$fg_status))
  if (nrow(base) == 0 || !any(base$fg_status %in% c(1L, 2L))) {
    return(tibble::tibble())
  }

  cif <- cmprsk::cuminc(
    ftime = base$fg_time_hours,
    fstatus = base$fg_status,
    cencode = 0
  )

  event_labels <- c(
    "1" = "Awake/extubated",
    "2" = "Death before awake/extubated"
  )

  purrr::map_dfr(names(event_labels), function(event_code) {
    object_name <- names(cif)[stringr::str_detect(names(cif), paste0("(^| )", event_code, "$"))][1]
    if (is.na(object_name)) return(tibble::tibble())
    curve <- cif[[object_name]]
    curve_times <- c(0, curve$time)
    curve_est <- c(0, curve$est)
    curve_var <- c(0, curve$var)
    unique_times <- sort(unique(curve_times))
    curve_est <- vapply(unique_times, function(x) tail(curve_est[curve_times == x], 1), numeric(1))
    curve_var <- vapply(unique_times, function(x) tail(curve_var[curve_times == x], 1), numeric(1))
    curve_times <- unique_times
    out <- tibble::tibble(
      model = model,
      event_code = as.integer(event_code),
      event_type = unname(event_labels[[event_code]]),
      time_hours = time_grid_hours,
      cif = stats::approx(curve_times, curve_est, xout = time_grid_hours, method = "constant", f = 0, rule = 2)$y,
      variance = stats::approx(curve_times, curve_var, xout = time_grid_hours, method = "constant", f = 0, rule = 2)$y,
      n = nrow(base),
      events_awake_extubated = sum(base$fg_status == 1L, na.rm = TRUE),
      competing_deaths = sum(base$fg_status == 2L, na.rm = TRUE)
    )
    out$standard_error <- sqrt(pmax(out$variance, 0))
    out$stratification <- stratification
    out$stratum <- stratum
    out$stratum_order <- stratum_order
    out$stratum_min_temp_c <- stratum_min_temp_c
    out$stratum_max_temp_c <- stratum_max_temp_c
    out[, c(
      "model", "stratification", "stratum", "stratum_order", "stratum_min_temp_c", "stratum_max_temp_c",
      "event_code", "event_type", "time_hours", "cif", "variance", "standard_error",
      "n", "events_awake_extubated", "competing_deaths"
    )]
  })
}

run_fine_gray <- function(df, exposure_terms, model, exposure_label, exposure_regex, adjust_terms = character()) {
  base <- df |>
    filter(!is.na(.data$fg_time_days), is.finite(.data$fg_time_days), .data$fg_time_days > 0, !is.na(.data$fg_status))
  events <- sum(base$fg_status == 1, na.rm = TRUE)
  competing <- sum(base$fg_status == 2, na.rm = TRUE)
  term_info <- choose_terms(base, c(exposure_terms, adjust_terms))
  exposure_keep <- intersect(exposure_terms, term_info$keep)
  adjust_keep <- setdiff(term_info$keep, exposure_keep)

  if (nrow(base) < 50 || events < 5 || competing < 5 || length(exposure_keep) == 0) {
    return(list(
      summary = extract_crr_terms(
        NULL, model, exposure_label, exposure_regex, nrow(base), events, competing,
        c(exposure_keep, adjust_keep), term_info$omitted, FALSE,
        "Insufficient sample size, awake/extubated events, competing deaths, or exposure variation"
      ),
      coefficients = tibble::tibble(),
      vcov = tibble::tibble()
    ))
  }

  model_df <- base |>
    select("fg_time_days", "fg_status", dplyr::all_of(c(exposure_keep, adjust_keep))) |>
    tidyr::drop_na()
  events <- sum(model_df$fg_status == 1, na.rm = TRUE)
  competing <- sum(model_df$fg_status == 2, na.rm = TRUE)
  matrix_info <- sanitize_covariate_matrix(make_model_matrix(model_df, c(exposure_keep, adjust_keep)))
  cov_matrix <- matrix_info$matrix
  matrix_dropped <- matrix_info$dropped
  if (is.null(cov_matrix) || ncol(cov_matrix) == 0 || !any(stringr::str_detect(colnames(cov_matrix), exposure_regex))) {
    return(list(
      summary = extract_crr_terms(
        NULL, model, exposure_label, exposure_regex, nrow(model_df), events, competing,
        c(exposure_keep, adjust_keep), c(term_info$omitted, matrix_dropped), FALSE,
        "Exposure terms were removed from the Fine-Gray design matrix"
      ),
      coefficients = tibble::tibble(),
      vcov = tibble::tibble()
    ))
  }
  fit <- tryCatch(
    cmprsk::crr(
      ftime = model_df$fg_time_days,
      fstatus = model_df$fg_status,
      cov1 = cov_matrix,
      failcode = 1,
      cencode = 0
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    fg_df <- model_df |>
      mutate(
        event_factor = factor(
          dplyr::case_when(
            .data$fg_status == 1L ~ "awake_extubated",
            .data$fg_status == 2L ~ "death",
            TRUE ~ "censor"
          ),
          levels = c("censor", "awake_extubated", "death")
        )
      )
    fg_formula <- stats::as.formula(paste("survival::Surv(fg_time_days, event_factor) ~", paste(c(exposure_keep, adjust_keep), collapse = " + ")))
    fg_expanded <- tryCatch(survival::finegray(fg_formula, data = fg_df, etype = "awake_extubated"), error = function(e) e)
    if (!inherits(fg_expanded, "error")) {
      cox_formula <- stats::as.formula(paste("survival::Surv(fgstart, fgstop, fgstatus) ~", paste(c(exposure_keep, adjust_keep), collapse = " + ")))
      fit <- tryCatch(
        survival::coxph(cox_formula, data = fg_expanded, weights = fgwt, ties = "breslow"),
        error = function(e) e
      )
      if (!inherits(fit, "error")) {
        fit <- list(coef = stats::coef(fit), var = stats::vcov(fit))
      }
    } else {
      fit <- fg_expanded
    }
  }

  all_covariates <- c(exposure_keep, adjust_keep)
  all_omitted <- c(term_info$omitted, matrix_dropped)
  list(
    summary = extract_crr_terms(
      fit, model, exposure_label, exposure_regex, nrow(model_df), events, competing,
      all_covariates, all_omitted
    ),
    coefficients = extract_crr_coefficients(
      fit, model, nrow(model_df), events, competing, all_covariates, all_omitted
    ),
    vcov = extract_crr_vcov(
      fit, model, nrow(model_df), events, competing, all_covariates, all_omitted
    )
  )
}

ohca <- readr::read_csv(OHCA_PATH, show_col_types = FALSE) |>
  mutate(
    hospitalization_id = as.character(.data$hospitalization_id),
    patient_id = as.character(.data$patient_id),
    county_fips = normalize_county_fips(.data$county_fips),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    first_icu_in = as_utc_datetime(.data$first_icu_in),
    admission_date = as.Date(.data$admission_date),
    admit_year = lubridate::year(.data$admission_date),
    race_group = ifelse(is_black_race(.data$race_category), "Black", "Non-Black"),
    sex_group = dplyr::case_when(
      is_male(.data$sex_category) ~ "Male",
      is_female(.data$sex_category) ~ "Female",
      TRUE ~ "Other/Unknown"
    ),
    hospital_death = ifelse(.data$hospital_death == 1L | is_expired_discharge(.data$discharge_category), 1L, 0L)
  )

tmax_raw <- arrow::read_parquet(TMAX_PATH)
tmax_county_col <- if ("county_fips" %in% names(tmax_raw)) "county_fips" else if ("geoid" %in% names(tmax_raw)) "geoid" else NA_character_
if (is.na(tmax_county_col)) stop("Tmax file must contain county_fips or geoid.", call. = FALSE)
tmax <- tmax_raw |>
  transmute(
    county_fips = normalize_county_fips(.data[[tmax_county_col]]),
    admission_date = as.Date(.data$date),
    tmax_mean_c = suppressWarnings(as.numeric(.data$tmax_mean_c))
  ) |>
  add_lag_window_means("tmax_mean_c", "tmax_mean_c")

rmax_raw <- arrow::read_parquet(RMAX_PATH)
rmax_county_col <- if ("county_fips" %in% names(rmax_raw)) "county_fips" else if ("geoid" %in% names(rmax_raw)) "geoid" else NA_character_
if (is.na(rmax_county_col)) stop("Rmax file must contain county_fips or geoid.", call. = FALSE)
rmax <- rmax_raw |>
  transmute(
    county_fips = normalize_county_fips(.data[[rmax_county_col]]),
    admission_date = as.Date(.data$date),
    rmax_mean_pct = suppressWarnings(as.numeric(.data$rmax_mean_pct))
  ) |>
  add_lag_window_means("rmax_mean_pct", "rmax_mean_pct")

cohort <- ohca |>
  left_join(tmax, by = c("county_fips", "admission_date")) |>
  left_join(rmax, by = c("county_fips", "admission_date")) |>
  filter(!is.na(.data$tmax_mean_c), !is.na(.data$rmax_mean_pct), !is.na(.data$first_icu_in)) |>
  mutate(
    tmax_per5c = .data$tmax_mean_c / 5,
    rmax_per10pct = .data$rmax_mean_pct / 10
  )

cohort_ids <- unique(cohort$hospitalization_id)

diagnosis_source <- NULL
for (candidate in c("hospital_diagnosis", "admission_diagnosis")) {
  if (file.exists(file.path(tables_path, sprintf("clif_%s.%s", candidate, file_type)))) {
    diagnosis_source <- candidate
    break
  }
}
if (is.null(diagnosis_source)) stop("Could not find clif_hospital_diagnosis or clif_admission_diagnosis.")

diagnosis_raw <- read_clif_table(tables_path, file_type, diagnosis_source, required = TRUE)
ohca_mechanism <- derive_ohca_mechanism(diagnosis_raw, cohort_ids, diagnosis_source)

patient <- read_clif_table(tables_path, file_type, "patient", columns = c("patient_id", "death_dttm"), required = FALSE)
if (is.null(patient) || nrow(patient) == 0) {
  patient_min <- tibble::tibble(patient_id = character(), patient_death_dttm = as.POSIXct(character(), tz = "UTC"))
} else {
  patient_min <- patient |>
    transmute(
      patient_id = as.character(.data$patient_id),
      patient_death_dttm = as_utc_datetime(.data$death_dttm)
    )
}

vitals <- read_clif_table(tables_path, file_type, "vitals", columns = c("hospitalization_id", "recorded_dttm"), required = FALSE)
if (is.null(vitals) || nrow(vitals) == 0) {
  vitals_last <- tibble::tibble(hospitalization_id = character(), last_vital_dttm = as.POSIXct(character(), tz = "UTC"))
} else {
  vitals_last <- vitals |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm)
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    group_by(.data$hospitalization_id) |>
    summarise(last_vital_dttm = max(.data$recorded_dttm, na.rm = TRUE), .groups = "drop")
}

cohort <- cohort |>
  left_join(ohca_mechanism, by = "hospitalization_id") |>
  left_join(patient_min, by = "patient_id") |>
  left_join(vitals_last, by = "hospitalization_id") |>
  mutate(
    ohca_mechanism = tidyr::replace_na(.data$ohca_mechanism, "unclear_other"),
    death_dttm_final = dplyr::case_when(
      .data$hospital_death == 1L & !is.na(.data$patient_death_dttm) ~ .data$patient_death_dttm,
      .data$hospital_death == 1L & is.na(.data$patient_death_dttm) & !is.na(.data$last_vital_dttm) ~ .data$last_vital_dttm,
      TRUE ~ as.POSIXct(NA, tz = "UTC")
    ),
    death_source = dplyr::case_when(
      .data$hospital_death == 1L & !is.na(.data$patient_death_dttm) ~ "patient_death_dttm",
      .data$hospital_death == 1L & is.na(.data$patient_death_dttm) & !is.na(.data$last_vital_dttm) ~ "expired_discharge_fallback_last_vital",
      .data$hospital_death == 1L ~ "expired_discharge_missing_death_time",
      TRUE ~ "no_death"
    )
  )

respiratory <- read_clif_table(
  tables_path,
  file_type,
  "respiratory_support",
  columns = c("hospitalization_id", "recorded_dttm", "device_category", "artificial_airway", "tracheostomy"),
  required = FALSE
)
if (is.null(respiratory) || nrow(respiratory) == 0) {
  resp_events <- tibble::tibble(hospitalization_id = cohort_ids)
} else if (!all(c("hospitalization_id", "recorded_dttm") %in% names(respiratory))) {
  resp_events <- tibble::tibble(hospitalization_id = cohort_ids)
} else {
  respiratory <- ensure_columns(respiratory, c("device_category", "artificial_airway", "tracheostomy"))
  resp_long <- respiratory |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      device_category = stringr::str_to_upper(tidyr::replace_na(as.character(.data$device_category), "")),
      artificial_airway = stringr::str_to_upper(tidyr::replace_na(as.character(.data$artificial_airway), "")),
      tracheostomy = stringr::str_to_upper(tidyr::replace_na(as.character(.data$tracheostomy), ""))
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    inner_join(cohort |> select("hospitalization_id", "first_icu_in"), by = "hospitalization_id") |>
    mutate(
      icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")),
      imv = .data$device_category %in% c("IMV", "VENT") |
        .data$artificial_airway %in% c("ETT", "NASAL ETT", "TRACH") |
        .data$tracheostomy %in% c("1", "TRUE", "YES")
    ) |>
    filter(.data$icu_hour >= 0) |>
    arrange(.data$hospitalization_id, .data$recorded_dttm)

  resp_events <- resp_long |>
    group_by(.data$hospitalization_id) |>
    summarise(
      cr_first_imv_dttm = suppressWarnings(min(.data$recorded_dttm[.data$imv], na.rm = TRUE)),
      cr_last_imv_dttm = {
        imv_times <- .data$recorded_dttm[.data$imv]
        if (length(imv_times) == 0L) as.POSIXct(NA, tz = "UTC") else max(imv_times, na.rm = TRUE)
      },
      cr_first_extubation_dttm = {
        imv_times <- .data$recorded_dttm[.data$imv]
        if (length(imv_times) == 0L) {
          as.POSIXct(NA, tz = "UTC")
        } else {
          last_imv <- max(imv_times, na.rm = TRUE)
          ext_times <- .data$recorded_dttm[!.data$imv & .data$recorded_dttm > last_imv]
          if (length(ext_times) == 0L) as.POSIXct(NA, tz = "UTC") else min(ext_times, na.rm = TRUE)
        }
      },
      cr_last_resp_dttm = max(.data$recorded_dttm, na.rm = TRUE),
      cr_n_resp_records = dplyr::n(),
      cr_n_imv_records = sum(.data$imv, na.rm = TRUE),
      cr_last_resp_device = safe_last(.data$device_category),
      cr_last_resp_imv = safe_last(.data$imv),
      .groups = "drop"
    ) |>
    mutate(
      cr_first_imv_dttm = ifelse(is.infinite(.data$cr_first_imv_dttm), as.POSIXct(NA, tz = "UTC"), .data$cr_first_imv_dttm),
      cr_first_imv_dttm = as_utc_datetime(.data$cr_first_imv_dttm),
      cr_last_imv_dttm = as_utc_datetime(.data$cr_last_imv_dttm)
    )
}

assessments <- read_clif_table(
  tables_path,
  file_type,
  "patient_assessments",
  columns = c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value", "categorical_value", "text_value"),
  required = FALSE
)
if (is.null(assessments) || nrow(assessments) == 0) {
  awake_events <- tibble::tibble(hospitalization_id = cohort_ids)
} else if (!all(c("hospitalization_id", "recorded_dttm", "assessment_category") %in% names(assessments))) {
  awake_events <- tibble::tibble(hospitalization_id = cohort_ids)
} else {
  assessments <- ensure_columns(assessments, c("numerical_value", "categorical_value", "text_value"))
  awake_events <- assessments |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      assessment_category_clean = clean_assessment_text(.data$assessment_category),
      numerical_value = suppressWarnings(as.numeric(.data$numerical_value)),
      categorical_value = clean_assessment_text(.data$categorical_value),
      text_value = clean_assessment_text(.data$text_value)
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    inner_join(cohort |> select("hospitalization_id", "first_icu_in"), by = "hospitalization_id") |>
    mutate(
      icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")),
      awake_signal = (
        .data$icu_hour >= 24 &
          .data$assessment_category_clean == "gcs_total" & !is.na(.data$numerical_value) & .data$numerical_value >= 13
      ) | (
        .data$icu_hour >= 24 &
          .data$assessment_category_clean == "gcs_motor" & !is.na(.data$numerical_value) & .data$numerical_value >= 6
      ) | (
        .data$icu_hour >= 24 &
          .data$assessment_category_clean == "rass" & !is.na(.data$numerical_value) & .data$numerical_value >= -1 & .data$numerical_value <= 4
      ) | (
        .data$assessment_category_clean == "avpu" & stringr::str_detect(paste(.data$categorical_value, .data$text_value), "alert")
      ) | (
        .data$assessment_category_clean %in% c("sat_delivery_pass_fail", "sbt_delivery_pass_fail") &
          .data$categorical_value == "pass"
      )
    ) |>
    filter(.data$icu_hour >= 0, .data$awake_signal) |>
    group_by(.data$hospitalization_id) |>
    summarise(
      first_awake_dttm = min(.data$recorded_dttm, na.rm = TRUE),
      n_awake_signal_records = dplyr::n(),
      .groups = "drop"
    )
}

analysis <- cohort |>
  left_join(resp_events, by = "hospitalization_id") |>
  left_join(awake_events, by = "hospitalization_id") |>
  filter(!is.na(.data$cr_first_imv_dttm)) |>
  mutate(
    cr_first_extubation_dttm = as_utc_datetime(.data$cr_first_extubation_dttm),
    first_awake_dttm = as_utc_datetime(.data$first_awake_dttm),
    awake_extubated_dttm = dplyr::case_when(
      !is.na(.data$cr_first_extubation_dttm) & !is.na(.data$first_awake_dttm) ~ pmax(.data$cr_first_extubation_dttm, .data$first_awake_dttm),
      TRUE ~ as.POSIXct(NA, tz = "UTC")
    ),
    administrative_censor_dttm = .data$first_icu_in + lubridate::hours(FOLLOWUP_HOURS),
    observed_censor_dttm = dplyr::coalesce(.data$discharge_dttm, .data$cr_last_resp_dttm, .data$last_vital_dttm, .data$administrative_censor_dttm),
    censor_dttm = pmin(.data$observed_censor_dttm, .data$administrative_censor_dttm, na.rm = TRUE),
    awake_extubated_before_death = !is.na(.data$awake_extubated_dttm) &
      .data$awake_extubated_dttm <= .data$administrative_censor_dttm &
      (is.na(.data$death_dttm_final) | .data$awake_extubated_dttm <= .data$death_dttm_final),
    death_before_awake_extubated = !is.na(.data$death_dttm_final) &
      .data$death_dttm_final <= .data$administrative_censor_dttm &
      (is.na(.data$awake_extubated_dttm) | .data$death_dttm_final < .data$awake_extubated_dttm),
    event_dttm = dplyr::case_when(
      .data$awake_extubated_before_death ~ .data$awake_extubated_dttm,
      .data$death_before_awake_extubated ~ .data$death_dttm_final,
      TRUE ~ .data$censor_dttm
    ),
    fg_status = dplyr::case_when(
      .data$awake_extubated_before_death ~ 1L,
      .data$death_before_awake_extubated ~ 2L,
      TRUE ~ 0L
    ),
    fg_time_days = as.numeric(difftime(.data$event_dttm, .data$first_icu_in, units = "days")),
    fg_time_hours = as.numeric(difftime(.data$event_dttm, .data$first_icu_in, units = "hours"))
  ) |>
  filter(!is.na(.data$fg_time_days), .data$fg_time_days > 0) |>
  mutate(
    admission_tmax_quartile_id = dplyr::ntile(.data$tmax_mean_c, 4L),
    admission_tmax_quartile = factor(
      dplyr::case_when(
        .data$admission_tmax_quartile_id == 1L ~ "Q1 coolest",
        .data$admission_tmax_quartile_id == 2L ~ "Q2",
        .data$admission_tmax_quartile_id == 3L ~ "Q3",
        .data$admission_tmax_quartile_id == 4L ~ "Q4 hottest",
        TRUE ~ NA_character_
      ),
      levels = c("Q1 coolest", "Q2", "Q3", "Q4 hottest")
    )
  )

summary_tbl <- analysis |>
  mutate(status_label = dplyr::case_when(
    .data$fg_status == 1L ~ "Awake and extubated by 72h",
    .data$fg_status == 2L ~ "Death before awake/extubated by 72h",
    TRUE ~ "Neither by 72h"
  )) |>
  count(.data$status_label, name = "n") |>
  mutate(
    site_name = site_name,
    pct = 100 * .data$n / sum(.data$n),
    definition = "Fine-Gray at-risk cohort is OHCA ICU admissions with invasive ventilation after ICU entry. Follow-up is administratively censored at 72 ICU hours. Event of interest is being both extubated and awake by 72 hours, defined as the later of first post-final-IMV non-IMV respiratory support record and first awake signal. Awake signal mirrors the regained-consciousness phenotype: GCS total >=13, GCS motor 6, or RASS >=-1 after ICU hour 24; AVPU alert, SAT pass, or SBT pass after ICU entry. Competing event is death before awake/extubated by 72 hours. Patients with neither event by 72 hours are censored as neither by 72h. Death time uses patient death_dttm when present, otherwise expired/death discharge with last recorded vital as fallback.",
    .before = 1
  )

cif_overall <- extract_cif_curve(analysis)
cif_by_tmax_quartile <- purrr::map_dfr(seq_len(4), function(quartile_id) {
  quartile_df <- analysis |>
    filter(.data$admission_tmax_quartile_id == quartile_id)
  if (nrow(quartile_df) == 0) return(tibble::tibble())
  quartile_label <- as.character(stats::na.omit(unique(quartile_df$admission_tmax_quartile))[1])
  extract_cif_curve(
    quartile_df,
    stratification = "Admission Tmax quartile",
    stratum = quartile_label,
    stratum_order = quartile_id,
    stratum_min_temp_c = min(quartile_df$tmax_mean_c, na.rm = TRUE),
    stratum_max_temp_c = max(quartile_df$tmax_mean_c, na.rm = TRUE)
  )
})
cif_curves <- bind_rows(cif_overall, cif_by_tmax_quartile) |>
  mutate(site_name = site_name, .before = 1)

death_source_summary <- analysis |>
  count(.data$death_source, name = "n") |>
  mutate(site_name = site_name, .before = 1)

mechanism_summary <- analysis |>
  mutate(status_label = dplyr::case_when(
    .data$fg_status == 1L ~ "Awake and extubated by 72h",
    .data$fg_status == 2L ~ "Death before awake/extubated by 72h",
    TRUE ~ "Neither by 72h"
  )) |>
  count(.data$ohca_mechanism, .data$status_label, name = "n") |>
  group_by(.data$ohca_mechanism) |>
  mutate(pct_within_mechanism = 100 * .data$n / sum(.data$n)) |>
  ungroup() |>
  group_by(.data$status_label) |>
  mutate(pct_within_status = 100 * .data$n / sum(.data$n)) |>
  ungroup() |>
  mutate(site_name = site_name, .before = 1)

base_adjust <- c("age_at_admission", "sex_group", "race_group", "admit_year")
mechanism_adjust <- c(base_adjust, "ohca_mechanism")

add_exposure_splines <- function(df, temp_var, humidity_var, prefix = "") {
  temp_spline <- splines::ns(df[[temp_var]], df = 3)
  humidity_spline <- splines::ns(df[[humidity_var]], df = 3)
  if (nzchar(prefix)) prefix <- paste0(prefix, "_")
  df[[paste0(prefix, "tmax_spline_1")]] <- temp_spline[, 1]
  df[[paste0(prefix, "tmax_spline_2")]] <- temp_spline[, 2]
  df[[paste0(prefix, "tmax_spline_3")]] <- temp_spline[, 3]
  df[[paste0(prefix, "rmax_spline_1")]] <- humidity_spline[, 1]
  df[[paste0(prefix, "rmax_spline_2")]] <- humidity_spline[, 2]
  df[[paste0(prefix, "rmax_spline_3")]] <- humidity_spline[, 3]
  df
}

analysis <- add_exposure_splines(analysis, "tmax_mean_c", "rmax_mean_pct")

fine_gray_results <- list(
  run_fine_gray(
    analysis,
    c("tmax_spline_1", "tmax_spline_2", "tmax_spline_3", "rmax_spline_1", "rmax_spline_2", "rmax_spline_3"),
    "awake_extubated_72h_vs_death_72h_temperature_humidity_demographics",
    "Natural spline Tmax and humidity, df=3",
    "^tmax_spline_|^rmax_spline_",
    base_adjust
  ),
  run_fine_gray(
    analysis,
    c("tmax_spline_1", "tmax_spline_2", "tmax_spline_3", "rmax_spline_1", "rmax_spline_2", "rmax_spline_3"),
    "awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_mechanism",
    "Natural spline Tmax and humidity, df=3",
    "^tmax_spline_|^rmax_spline_",
    mechanism_adjust
  )
)

lag_window_specs <- tibble::tribble(
  ~exposure_window, ~temp_var, ~humidity_var,
  "lag0_1", "tmax_mean_c_lag0_1", "rmax_mean_pct_lag0_1",
  "lag0_3", "tmax_mean_c_lag0_3", "rmax_mean_pct_lag0_3",
  "lag0_5", "tmax_mean_c_lag0_5", "rmax_mean_pct_lag0_5"
)

fine_gray_lag_results <- purrr::map(seq_len(nrow(lag_window_specs)), function(i) {
  spec <- lag_window_specs[i, ]
  lag_df <- analysis |>
    filter(!is.na(.data[[spec$temp_var]]), is.finite(.data[[spec$temp_var]]), !is.na(.data[[spec$humidity_var]]), is.finite(.data[[spec$humidity_var]]))
  if (nrow(lag_df) > 0) {
    lag_df <- add_exposure_splines(lag_df, spec$temp_var, spec$humidity_var, prefix = spec$exposure_window)
  }
  exposure_terms <- paste0(spec$exposure_window, "_", c(
    "tmax_spline_1", "tmax_spline_2", "tmax_spline_3",
    "rmax_spline_1", "rmax_spline_2", "rmax_spline_3"
  ))
  exposure_regex <- paste0("^", spec$exposure_window, "_tmax_spline_|^", spec$exposure_window, "_rmax_spline_")
  primary <- run_fine_gray(
    lag_df,
    exposure_terms,
    paste0("awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_", spec$exposure_window),
    paste0("Natural spline mean Tmax and humidity ", spec$exposure_window, ", df=3"),
    exposure_regex,
    base_adjust
  )
  mechanism <- run_fine_gray(
    lag_df,
    exposure_terms,
    paste0("awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_mechanism_", spec$exposure_window),
    paste0("Natural spline mean Tmax and humidity ", spec$exposure_window, ", df=3"),
    exposure_regex,
    mechanism_adjust
  )
  list(primary = primary, mechanism = mechanism)
})

models <- purrr::map_dfr(fine_gray_results, "summary") |>
  mutate(site_name = site_name, .before = 1)
model_coefficients <- purrr::map_dfr(fine_gray_results, "coefficients") |>
  mutate(site_name = site_name, .before = 1)
model_vcov <- purrr::map_dfr(fine_gray_results, "vcov") |>
  mutate(site_name = site_name, .before = 1)

fine_gray_lag_models <- purrr::map_dfr(fine_gray_lag_results, function(x) {
  bind_rows(x$primary$summary, x$mechanism$summary)
}) |>
  mutate(site_name = site_name, .before = 1)
fine_gray_lag_coefficients <- purrr::map_dfr(fine_gray_lag_results, function(x) {
  bind_rows(x$primary$coefficients, x$mechanism$coefficients)
}) |>
  mutate(site_name = site_name, .before = 1)
fine_gray_lag_vcov <- purrr::map_dfr(fine_gray_lag_results, function(x) {
  bind_rows(x$primary$vcov, x$mechanism$vcov)
}) |>
  mutate(site_name = site_name, .before = 1)

patient_audit <- analysis |>
  select(
    "hospitalization_id", "patient_id", "admission_date", "first_icu_in", "cr_first_imv_dttm", "cr_last_imv_dttm",
    "ohca_mechanism", "ohca_mechanism_detail",
    "cr_first_extubation_dttm", "first_awake_dttm", "awake_extubated_dttm", "administrative_censor_dttm",
    "death_dttm_final", "death_source", "discharge_dttm",
    "last_vital_dttm", "cr_last_resp_dttm", "fg_status", "fg_time_hours", "fg_time_days",
    "tmax_mean_c", "rmax_mean_pct", "cr_n_resp_records", "cr_n_imv_records",
    "cr_last_resp_device", "cr_last_resp_imv", "n_awake_signal_records"
  ) |>
  mutate(site_name = site_name, .before = 1)

readr::write_csv(summary_tbl, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_summary.csv"))
readr::write_csv(cif_curves, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv"))
readr::write_csv(death_source_summary, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_death_source_summary.csv"))
readr::write_csv(mechanism_summary, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_ohca_mechanism_summary.csv"))
readr::write_csv(models, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_fine_gray_models.csv"))
readr::write_csv(model_coefficients, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv"))
readr::write_csv(model_vcov, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_vcov.csv"))
readr::write_csv(fine_gray_lag_models, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_fine_gray_models.csv"))
readr::write_csv(fine_gray_lag_coefficients, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_coefficients.csv"))
readr::write_csv(fine_gray_lag_vcov, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_vcov.csv"))
readr::write_csv(patient_audit, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_patient_audit.csv"))

if (nrow(cif_curves) > 0) {
  cif_plot <- ggplot(
    filter(cif_curves, .data$stratification == "Overall"),
    aes(x = .data$time_hours, y = .data$cif, color = .data$event_type)
  ) +
    geom_step(linewidth = 1.1) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, NA)) +
    scale_color_manual(values = c("Awake/extubated" = "#0f6b78", "Death before awake/extubated" = "#8c3b3b")) +
    labs(
      title = "72-hour competing-risk cumulative incidence after OHCA ICU admission",
      subtitle = "Event of interest is awake/extubated; death before awake/extubated is the competing event.",
      x = "Hours since ICU admission",
      y = "Cumulative incidence",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey35"),
      legend.position = "bottom"
    )
  ggsave(file.path(FIGURE_DIR, "figure_ohca_icu_competing_risk_awake_extubated_72h_cif.png"), cif_plot, width = 8, height = 5, dpi = 300)

  quartile_cif <- cif_curves |>
    filter(.data$stratification == "Admission Tmax quartile") |>
    mutate(stratum = factor(.data$stratum, levels = c("Q1 coolest", "Q2", "Q3", "Q4 hottest")))
  if (nrow(quartile_cif) > 0) {
    quartile_plot <- ggplot(quartile_cif, aes(x = .data$time_hours, y = .data$cif, color = .data$stratum)) +
      geom_step(linewidth = 1.05) +
      facet_wrap(vars(.data$event_type), nrow = 1) +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, NA)) +
      scale_x_continuous(breaks = seq(0, FOLLOWUP_HOURS, by = 12), limits = c(0, FOLLOWUP_HOURS)) +
      scale_color_manual(values = c("Q1 coolest" = "#335c67", "Q2" = "#8f7a2f", "Q3" = "#c06c3e", "Q4 hottest" = "#9b2f37")) +
      labs(
        title = "72-hour competing-risk cumulative incidence by admission temperature quartile",
        subtitle = "Quartiles are defined within site from admission daily maximum temperature.",
        x = "Hours since ICU admission",
        y = "Cumulative incidence",
        color = "Tmax quartile"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey35"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      )
    ggsave(file.path(FIGURE_DIR, "figure_ohca_icu_competing_risk_awake_extubated_72h_cif_by_tmax_quartile.png"), quartile_plot, width = 10, height = 5, dpi = 300)
  }
}

print(summary_tbl)
print(cif_curves |> filter(.data$time_hours %in% c(0, 24, 48, 72)))
print(death_source_summary)
print(models)

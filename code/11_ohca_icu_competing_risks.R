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
  library(jsonlite)
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

run_fine_gray <- function(df, exposure_terms, model, exposure_label, exposure_regex, adjust_terms = character()) {
  base <- df |>
    filter(!is.na(.data$fg_time_days), is.finite(.data$fg_time_days), .data$fg_time_days > 0, !is.na(.data$fg_status))
  events <- sum(base$fg_status == 1, na.rm = TRUE)
  competing <- sum(base$fg_status == 2, na.rm = TRUE)
  term_info <- choose_terms(base, c(exposure_terms, adjust_terms))
  exposure_keep <- intersect(exposure_terms, term_info$keep)
  adjust_keep <- setdiff(term_info$keep, exposure_keep)

  if (nrow(base) < 50 || events < 5 || competing < 5 || length(exposure_keep) == 0) {
    return(extract_crr_terms(
      NULL, model, exposure_label, exposure_regex, nrow(base), events, competing,
      c(exposure_keep, adjust_keep), term_info$omitted, FALSE,
      "Insufficient sample size, awake/extubated events, competing deaths, or exposure variation"
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
    return(extract_crr_terms(
      NULL, model, exposure_label, exposure_regex, nrow(model_df), events, competing,
      c(exposure_keep, adjust_keep), c(term_info$omitted, matrix_dropped), FALSE,
      "Exposure terms were removed from the Fine-Gray design matrix"
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

  extract_crr_terms(
    fit, model, exposure_label, exposure_regex, nrow(model_df), events, competing,
    c(exposure_keep, adjust_keep), c(term_info$omitted, matrix_dropped)
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
  )

rmax_raw <- arrow::read_parquet(RMAX_PATH)
rmax_county_col <- if ("county_fips" %in% names(rmax_raw)) "county_fips" else if ("geoid" %in% names(rmax_raw)) "geoid" else NA_character_
if (is.na(rmax_county_col)) stop("Rmax file must contain county_fips or geoid.", call. = FALSE)
rmax <- rmax_raw |>
  transmute(
    county_fips = normalize_county_fips(.data[[rmax_county_col]]),
    admission_date = as.Date(.data$date),
    rmax_mean_pct = suppressWarnings(as.numeric(.data$rmax_mean_pct))
  )

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
} else {
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
} else {
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
  filter(!is.na(.data$fg_time_days), .data$fg_time_days > 0)

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
analysis <- analysis |>
  mutate(
    tmax_spline_1 = splines::ns(.data$tmax_mean_c, df = 3)[, 1],
    tmax_spline_2 = splines::ns(.data$tmax_mean_c, df = 3)[, 2],
    tmax_spline_3 = splines::ns(.data$tmax_mean_c, df = 3)[, 3],
    rmax_spline_1 = splines::ns(.data$rmax_mean_pct, df = 3)[, 1],
    rmax_spline_2 = splines::ns(.data$rmax_mean_pct, df = 3)[, 2],
    rmax_spline_3 = splines::ns(.data$rmax_mean_pct, df = 3)[, 3]
  )

models <- dplyr::bind_rows(
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
) |>
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
readr::write_csv(death_source_summary, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_death_source_summary.csv"))
readr::write_csv(mechanism_summary, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_ohca_mechanism_summary.csv"))
readr::write_csv(models, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_fine_gray_models.csv"))
readr::write_csv(patient_audit, file.path(OUTPUT_DIR, "ohca_icu_competing_risk_awake_extubated_72h_patient_audit.csv"))

print(summary_tbl)
print(death_source_summary)
print(models)

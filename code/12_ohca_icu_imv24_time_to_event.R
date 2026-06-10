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
  library(lubridate)
  library(readr)
  library(stringr)
  library(tidyr)
})

OHCA_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
TMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
OUTPUT_DIR <- file.path(repo_root, "output", "final", "ohca_icu_phenotypes")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

unlink(file.path(OUTPUT_DIR, c(
  "ohca_icu_imv24_time_to_event_summary.csv",
  "ohca_icu_imv24_time_to_event_landmark_flow.csv",
  "ohca_icu_imv24_time_to_event_cif_curves.csv",
  "ohca_icu_imv24_time_to_event_evidence_summary.csv",
  "ohca_icu_imv24_time_to_event_table_by_outcome.csv",
  "ohca_icu_imv24_time_to_event_fine_gray_models.csv",
  "ohca_icu_imv24_time_to_event_coefficients.csv",
  "ohca_icu_imv24_time_to_event_vcov.csv",
  "ohca_icu_imv24_time_to_event_patient_audit.csv"
)))

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

empty_time <- function() as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")

min_time <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(empty_time())
  min(x)
}

max_time <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(empty_time())
  max(x)
}

ensure_columns <- function(df, columns) {
  if (is.null(df)) df <- tibble::tibble()
  for (column in columns) {
    if (!column %in% names(df)) df[[column]] <- NA
  }
  df
}

is_hospice_discharge <- function(x) {
  stringr::str_detect(normalize_category(as.character(x)), "hospice|comfort")
}

yesno <- function(x) {
  tidyr::replace_na(as.logical(x), FALSE)
}

norm_code <- function(x) {
  x |>
    tidyr::replace_na("") |>
    as.character() |>
    stringr::str_to_upper() |>
    stringr::str_replace_all("[^A-Z0-9]", "") |>
    trimws()
}

infer_diagnosis_format <- function(dx_clean, diagnosis_code_format) {
  explicit <- diagnosis_code_format |>
    tidyr::replace_na("") |>
    as.character() |>
    stringr::str_to_upper()
  inferred <- ifelse(
    stringr::str_detect(dx_clean, "^[A-Z]"),
    "ICD10",
    ifelse(nchar(dx_clean) >= 3, "ICD9", NA_character_)
  )
  ifelse(nzchar(explicit), explicit, inferred)
}

format_n_pct <- function(x, denom) {
  fmt_n_pct(sum(x, na.rm = TRUE), denom)
}

format_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  fmt_median_iqr(x)
}

add_outcome_continuous <- function(df, section, characteristic, variable, outcome_levels) {
  row <- tibble::tibble(section = section, characteristic = characteristic, level = "Median [IQR]")
  for (outcome in outcome_levels) {
    x <- df[[variable]][df$outcome_status == outcome]
    row[[outcome]] <- format_median_iqr(x)
  }
  row
}

add_outcome_binary <- function(df, section, characteristic, variable, outcome_levels) {
  row <- tibble::tibble(section = section, characteristic = characteristic, level = "Yes")
  for (outcome in outcome_levels) {
    keep <- df$outcome_status == outcome
    row[[outcome]] <- format_n_pct(df[[variable]][keep] %in% TRUE | df[[variable]][keep] == 1L, sum(keep, na.rm = TRUE))
  }
  row
}

add_outcome_categorical <- function(df, section, characteristic, variable, outcome_levels) {
  levels <- sort(unique(stats::na.omit(as.character(df[[variable]]))))
  if (length(levels) == 0L) return(tibble::tibble())
  rows <- lapply(levels, function(level_value) {
    row <- tibble::tibble(section = section, characteristic = characteristic, level = level_value)
    for (outcome in outcome_levels) {
      keep <- df$outcome_status == outcome
      row[[outcome]] <- fmt_n_pct(sum(as.character(df[[variable]][keep]) == level_value, na.rm = TRUE), sum(keep, na.rm = TRUE))
    }
    row
  })
  dplyr::bind_rows(rows)
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
  if (length(terms) == 0L) return(NULL)
  stats::model.matrix(stats::as.formula(paste("~", paste(terms, collapse = " + "))), data = df)[, -1, drop = FALSE]
}

sanitize_covariate_matrix <- function(mat) {
  if (is.null(mat) || ncol(mat) == 0L) return(list(matrix = mat, dropped = character()))
  dropped <- character()
  finite_cols <- vapply(seq_len(ncol(mat)), function(i) all(is.finite(mat[, i])), logical(1))
  if (any(!finite_cols)) dropped <- c(dropped, paste0(colnames(mat)[!finite_cols], "_nonfinite"))
  mat <- mat[, finite_cols, drop = FALSE]
  if (ncol(mat) == 0L) return(list(matrix = mat, dropped = dropped))

  varying_cols <- vapply(seq_len(ncol(mat)), function(i) length(unique(mat[, i])) >= 2, logical(1))
  if (any(!varying_cols)) dropped <- c(dropped, paste0(colnames(mat)[!varying_cols], "_constant"))
  mat <- mat[, varying_cols, drop = FALSE]
  if (ncol(mat) == 0L) return(list(matrix = mat, dropped = dropped))

  qr_mat <- qr(mat)
  if (qr_mat$rank < ncol(mat)) {
    independent <- sort(qr_mat$pivot[seq_len(qr_mat$rank)])
    dropped <- c(dropped, paste0(colnames(mat)[setdiff(seq_len(ncol(mat)), independent)], "_collinear"))
    mat <- mat[, independent, drop = FALSE]
  }
  list(matrix = mat, dropped = dropped)
}

extract_crr_summary <- function(fit, model, event_code, event_type, n, events, competing, covariates, omitted_terms, estimable = TRUE, skip_reason = NA_character_) {
  if (!estimable || inherits(fit, "error") || is.null(fit$coef) || is.null(fit$var) || !"tmax_per5c" %in% names(fit$coef)) {
    return(tibble::tibble(
      model = model,
      event_code = event_code,
      event_type = event_type,
      exposure = "Admission Tmax per 5 C",
      n = n,
      events = events,
      competing_events = competing,
      subdistribution_hr = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(covariates, collapse = " + "),
      omitted_terms = paste(omitted_terms, collapse = "; "),
      estimable = FALSE,
      skip_reason = if (inherits(fit, "error")) conditionMessage(fit) else skip_reason
    ))
  }

  beta <- fit$coef[["tmax_per5c"]]
  se <- sqrt(diag(fit$var))
  if (is.null(names(se))) names(se) <- names(fit$coef)
  if (!"tmax_per5c" %in% names(se)) {
    return(tibble::tibble(
      model = model,
      event_code = event_code,
      event_type = event_type,
      exposure = "Admission Tmax per 5 C",
      n = n,
      events = events,
      competing_events = competing,
      subdistribution_hr = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(covariates, collapse = " + "),
      omitted_terms = paste(omitted_terms, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Fine-Gray variance matrix did not preserve the temperature coefficient name"
    ))
  }
  se <- se[["tmax_per5c"]]
  tibble::tibble(
    model = model,
    event_code = event_code,
    event_type = event_type,
    exposure = "Admission Tmax per 5 C",
    n = n,
    events = events,
    competing_events = competing,
    subdistribution_hr = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    p_value = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE),
    covariates = paste(covariates, collapse = " + "),
    omitted_terms = paste(omitted_terms, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
}

extract_crr_coefficients <- function(fit, model, event_code, event_type, n, events, competing, covariates, omitted_terms) {
  if (inherits(fit, "error") || is.null(fit$coef) || is.null(fit$var)) return(tibble::tibble())
  coefs <- fit$coef
  vc <- fit$var
  se <- sqrt(diag(vc))
  if (is.null(names(se))) names(se) <- names(coefs)
  tibble::tibble(
    model = model,
    event_code = event_code,
    event_type = event_type,
    n = n,
    events = events,
    competing_events = competing,
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

extract_crr_vcov <- function(fit, model, event_code, event_type, n, events, competing, covariates, omitted_terms) {
  if (inherits(fit, "error") || is.null(fit$coef) || is.null(fit$var)) return(tibble::tibble())
  vc <- fit$var
  if (is.null(rownames(vc))) rownames(vc) <- names(fit$coef)
  if (is.null(colnames(vc))) colnames(vc) <- names(fit$coef)
  as.data.frame(as.table(vc), stringsAsFactors = FALSE) |>
    transmute(
      model = model,
      event_code = event_code,
      event_type = event_type,
      n = n,
      events = events,
      competing_events = competing,
      coefficient_row = as.character(.data$Var1),
      coefficient_col = as.character(.data$Var2),
      covariance = as.numeric(.data$Freq),
      covariates = paste(covariates, collapse = " + "),
      omitted_terms = paste(omitted_terms, collapse = "; ")
    )
}

run_fine_gray <- function(df, event_code, event_type, adjust_terms) {
  model <- paste0("imv24_", event_type, "_temperature_demographics")
  model <- stringr::str_replace_all(model, "[^A-Za-z0-9]+", "_")
  base <- df |>
    filter(!is.na(.data$fg_time_days), is.finite(.data$fg_time_days), .data$fg_time_days > 0, !is.na(.data$fg_status))
  events <- sum(base$fg_status == event_code, na.rm = TRUE)
  competing <- sum(base$fg_status != 0L & base$fg_status != event_code, na.rm = TRUE)
  term_info <- choose_terms(base, c("tmax_per5c", adjust_terms))

  if (!"tmax_per5c" %in% term_info$keep || nrow(base) < 50 || events < 5 || competing < 5) {
    return(list(
      summary = extract_crr_summary(
        NULL, model, event_code, event_type, nrow(base), events, competing,
        term_info$keep, term_info$omitted, FALSE,
        "Insufficient sample size, event counts, competing events, or exposure variation"
      ),
      coefficients = tibble::tibble(),
      vcov = tibble::tibble()
    ))
  }

  model_df <- base |>
    select("fg_time_days", "fg_status", dplyr::all_of(term_info$keep)) |>
    tidyr::drop_na()
  events <- sum(model_df$fg_status == event_code, na.rm = TRUE)
  competing <- sum(model_df$fg_status != 0L & model_df$fg_status != event_code, na.rm = TRUE)
  matrix_info <- sanitize_covariate_matrix(make_model_matrix(model_df, term_info$keep))
  cov_matrix <- matrix_info$matrix
  all_omitted <- c(term_info$omitted, matrix_info$dropped)
  if (is.null(cov_matrix) || ncol(cov_matrix) == 0L || !"tmax_per5c" %in% colnames(cov_matrix)) {
    return(list(
      summary = extract_crr_summary(
        NULL, model, event_code, event_type, nrow(model_df), events, competing,
        term_info$keep, all_omitted, FALSE,
        "Temperature term was removed from the Fine-Gray design matrix"
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
      failcode = event_code,
      cencode = 0
    ),
    error = function(e) e
  )

  list(
    summary = extract_crr_summary(fit, model, event_code, event_type, nrow(model_df), events, competing, term_info$keep, all_omitted),
    coefficients = extract_crr_coefficients(fit, model, event_code, event_type, nrow(model_df), events, competing, term_info$keep, all_omitted),
    vcov = extract_crr_vcov(fit, model, event_code, event_type, nrow(model_df), events, competing, term_info$keep, all_omitted)
  )
}

extract_cif_curve <- function(df, time_grid_days = seq(0, 60, by = 1),
                              stratification = "Overall", stratum = "Overall", stratum_order = 0,
                              stratum_min_temp_c = NA_real_, stratum_max_temp_c = NA_real_) {
  base <- df |>
    filter(!is.na(.data$fg_time_days), is.finite(.data$fg_time_days), .data$fg_time_days >= 0, !is.na(.data$fg_status))
  if (nrow(base) == 0 || !any(base$fg_status %in% c(1L, 2L, 3L))) return(tibble::tibble())

  cif <- cmprsk::cuminc(ftime = base$fg_time_days, fstatus = base$fg_status, cencode = 0)
  event_labels <- c(
    "1" = "Successful extubation",
    "2" = "Death/hospice discharge",
    "3" = "Tracheostomy"
  )

  purrr::map_dfr(names(event_labels), function(event_code) {
    code_int <- as.integer(event_code)
    object_name <- names(cif)[stringr::str_detect(names(cif), paste0("(^| )", event_code, "$"))][1]
    if (is.na(object_name)) return(tibble::tibble())
    curve <- cif[[object_name]]
    curve_times <- c(0, curve$time)
    curve_est <- c(0, curve$est)
    curve_var <- c(0, curve$var)
    unique_times <- sort(unique(curve_times))
    curve_est <- vapply(unique_times, function(x) tail(curve_est[curve_times == x], 1), numeric(1))
    curve_var <- vapply(unique_times, function(x) tail(curve_var[curve_times == x], 1), numeric(1))
    tibble::tibble(
      model = "imv24_competing_risks_cif",
      stratification = stratification,
      stratum = stratum,
      stratum_order = stratum_order,
      stratum_min_temp_c = stratum_min_temp_c,
      stratum_max_temp_c = stratum_max_temp_c,
      event_code = code_int,
      event_type = unname(event_labels[[event_code]]),
      time_days = time_grid_days,
      time_hours = time_grid_days * 24,
      cif = stats::approx(unique_times, curve_est, xout = time_grid_days, method = "constant", f = 0, rule = 2)$y,
      variance = stats::approx(unique_times, curve_var, xout = time_grid_days, method = "constant", f = 0, rule = 2)$y,
      n = nrow(base),
      events = sum(base$fg_status == code_int, na.rm = TRUE)
    ) |>
      mutate(standard_error = sqrt(pmax(.data$variance, 0)))
  })
}

ohca <- readr::read_csv(OHCA_PATH, show_col_types = FALSE) |>
  ensure_columns(c("discharge_name")) |>
  mutate(
    hospitalization_id = as.character(.data$hospitalization_id),
    patient_id = as.character(.data$patient_id),
    county_fips = normalize_county_fips(.data$county_fips),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    first_icu_in = as_utc_datetime(.data$first_icu_in),
    admission_date = as.Date(.data$admission_date),
    race_group = ifelse(is_black_race(.data$race_category), "Black", "Non-Black"),
    sex_group = dplyr::case_when(
      is_male(.data$sex_category) ~ "Male",
      is_female(.data$sex_category) ~ "Female",
      TRUE ~ "Other/Unknown"
    ),
    hospital_death = ifelse(.data$hospital_death == 1L | is_expired_discharge(.data$discharge_category), 1L, 0L),
    hospice_discharge = is_hospice_discharge(.data$discharge_category) | is_hospice_discharge(.data$discharge_name)
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

cohort <- ohca |>
  left_join(tmax, by = c("county_fips", "admission_date")) |>
  filter(!is.na(.data$first_icu_in), !is.na(.data$tmax_mean_c)) |>
  mutate(
    tmax_per5c = .data$tmax_mean_c / 5,
    landmark_dttm = .data$first_icu_in + lubridate::hours(24)
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
diagnosis <- diagnosis_raw |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    diagnosis_code = if ("diagnosis_code" %in% names(diagnosis_raw)) as.character(.data$diagnosis_code) else NA_character_,
    diagnosis_code_format = if ("diagnosis_code_format" %in% names(diagnosis_raw)) as.character(.data$diagnosis_code_format) else NA_character_
  ) |>
  filter(.data$hospitalization_id %in% cohort_ids) |>
  mutate(
    diagnosis_code_clean = norm_code(.data$diagnosis_code),
    diagnosis_code_format = infer_diagnosis_format(.data$diagnosis_code_clean, .data$diagnosis_code_format),
    anoxic_brain_injury_dx = (
      stringr::str_detect(.data$diagnosis_code_format, "10") & stringr::str_starts(.data$diagnosis_code_clean, "G931")
    ) | (
      stringr::str_detect(.data$diagnosis_code_format, "9") & stringr::str_starts(.data$diagnosis_code_clean, "3481")
    ),
    brain_death_dx = stringr::str_detect(.data$diagnosis_code_format, "10") & stringr::str_starts(.data$diagnosis_code_clean, "G9382")
  )

neuro_dx <- diagnosis |>
  group_by(.data$hospitalization_id) |>
  summarise(
    anoxic_brain_injury_dx = any(.data$anoxic_brain_injury_dx, na.rm = TRUE),
    brain_death_dx = any(.data$brain_death_dx, na.rm = TRUE),
    anoxic_brain_injury_codes = paste(sort(unique(.data$diagnosis_code_clean[.data$anoxic_brain_injury_dx | .data$brain_death_dx])), collapse = " | "),
    .groups = "drop"
  )

patient <- read_clif_table(tables_path, file_type, "patient", columns = c("patient_id", "death_dttm"), required = FALSE)
patient_min <- if (is.null(patient) || nrow(patient) == 0) {
  tibble::tibble(patient_id = character(), patient_death_dttm = as.POSIXct(character(), tz = "UTC"))
} else {
  patient |>
    transmute(
      patient_id = as.character(.data$patient_id),
      patient_death_dttm = as_utc_datetime(.data$death_dttm)
    )
}

vitals <- read_clif_table(tables_path, file_type, "vitals", columns = c("hospitalization_id", "recorded_dttm"), required = FALSE)
vitals_last <- if (is.null(vitals) || nrow(vitals) == 0) {
  tibble::tibble(hospitalization_id = character(), last_vital_dttm = as.POSIXct(character(), tz = "UTC"))
} else {
  vitals |>
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
  left_join(neuro_dx, by = "hospitalization_id") |>
  left_join(patient_min, by = "patient_id") |>
  left_join(vitals_last, by = "hospitalization_id") |>
  mutate(
    ohca_mechanism = tidyr::replace_na(.data$ohca_mechanism, "unclear_other"),
    anoxic_brain_injury_dx = yesno(.data$anoxic_brain_injury_dx),
    brain_death_dx = yesno(.data$brain_death_dx),
    death_dttm_final = dplyr::case_when(
      .data$hospital_death == 1L & !is.na(.data$patient_death_dttm) ~ .data$patient_death_dttm,
      .data$hospital_death == 1L & is.na(.data$patient_death_dttm) & !is.na(.data$last_vital_dttm) ~ .data$last_vital_dttm,
      TRUE ~ as.POSIXct(NA, tz = "UTC")
    ),
    hospice_dttm = dplyr::case_when(
      .data$hospice_discharge & !is.na(.data$discharge_dttm) ~ .data$discharge_dttm,
      TRUE ~ as.POSIXct(NA, tz = "UTC")
    ),
    death_hospice_dttm = pmin(.data$death_dttm_final, .data$hospice_dttm, na.rm = TRUE),
    death_hospice_dttm = ifelse(is.infinite(.data$death_hospice_dttm), as.POSIXct(NA, tz = "UTC"), .data$death_hospice_dttm),
    death_hospice_dttm = as_utc_datetime(.data$death_hospice_dttm)
  )

respiratory <- read_clif_table(
  tables_path,
  file_type,
  "respiratory_support",
  columns = c("hospitalization_id", "recorded_dttm", "device_category", "tracheostomy"),
  required = FALSE
)

if (is.null(respiratory) || nrow(respiratory) == 0 || !all(c("hospitalization_id", "recorded_dttm") %in% names(respiratory))) {
  resp_events <- tibble::tibble(hospitalization_id = cohort_ids)
} else {
  respiratory <- ensure_columns(respiratory, c("device_category", "tracheostomy"))
  resp_long <- respiratory |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      device_category = stringr::str_to_upper(tidyr::replace_na(as.character(.data$device_category), "")),
      tracheostomy = stringr::str_to_upper(tidyr::replace_na(as.character(.data$tracheostomy), ""))
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    inner_join(cohort |> select("hospitalization_id", "first_icu_in", "landmark_dttm"), by = "hospitalization_id") |>
    mutate(
      icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")),
      trach = .data$tracheostomy %in% c("1", "TRUE", "YES"),
      imv = .data$device_category %in% c("IMV", "VENT")
    ) |>
    filter(.data$icu_hour >= 0) |>
    arrange(.data$hospitalization_id, .data$recorded_dttm)

  resp_events <- resp_long |>
    group_by(.data$hospitalization_id) |>
    summarise(
      last_resp_before_landmark_dttm = max_time(.data$recorded_dttm[.data$recorded_dttm <= safe_last(.data$landmark_dttm)]),
      imv_at_landmark = {
        landmark <- safe_last(.data$landmark_dttm)
        before <- .data$recorded_dttm <= landmark
        if (!any(before, na.rm = TRUE)) FALSE else safe_last(.data$imv[before])
      },
      first_imv_dttm = min_time(.data$recorded_dttm[.data$imv %in% TRUE]),
      last_resp_dttm = max_time(.data$recorded_dttm),
      trach_before_or_at_landmark = {
        landmark <- safe_last(.data$landmark_dttm)
        any(.data$trach & .data$recorded_dttm <= landmark, na.rm = TRUE)
      },
      first_trach_after_landmark_dttm = {
        landmark <- safe_last(.data$landmark_dttm)
        min_time(.data$recorded_dttm[.data$trach & .data$recorded_dttm > landmark])
      },
      final_imv_after_landmark_dttm = {
        landmark <- safe_last(.data$landmark_dttm)
        max_time(c(landmark, .data$recorded_dttm[.data$imv %in% TRUE & .data$recorded_dttm > landmark]))
      },
      first_non_imv_after_final_imv_dttm = {
        landmark <- safe_last(.data$landmark_dttm)
        final_imv <- max_time(c(landmark, .data$recorded_dttm[.data$imv %in% TRUE & .data$recorded_dttm > landmark]))
        min_time(.data$recorded_dttm[!(.data$imv %in% TRUE) & .data$recorded_dttm > final_imv])
      },
      n_resp_records = dplyr::n(),
      n_imv_records = sum(.data$imv, na.rm = TRUE),
      n_trach_records = sum(.data$trach, na.rm = TRUE),
      .groups = "drop"
    )
}

analysis_base <- cohort |>
  left_join(resp_events, by = "hospitalization_id") |>
  mutate(
    pre_landmark_death_hospice = !is.na(.data$death_hospice_dttm) & .data$death_hospice_dttm <= .data$landmark_dttm,
    landmark_eligible = .data$imv_at_landmark %in% TRUE &
      !(.data$trach_before_or_at_landmark %in% TRUE) &
      !(.data$pre_landmark_death_hospice %in% TRUE),
    successful_extubation_dttm = dplyr::case_when(
      .data$landmark_eligible & !is.na(.data$first_non_imv_after_final_imv_dttm) & is.na(.data$death_hospice_dttm) ~ .data$first_non_imv_after_final_imv_dttm,
      TRUE ~ as.POSIXct(NA, tz = "UTC")
    ),
    censor_dttm = dplyr::coalesce(.data$discharge_dttm, .data$last_resp_dttm, .data$last_vital_dttm),
    censor_dttm = dplyr::case_when(
      is.na(.data$censor_dttm) | .data$censor_dttm < .data$landmark_dttm ~ .data$landmark_dttm,
      TRUE ~ .data$censor_dttm
    )
  )

analysis_base <- ensure_columns(analysis_base, c(
  "first_imv_dttm",
  "final_imv_after_landmark_dttm",
  "n_resp_records",
  "n_imv_records",
  "n_trach_records"
)) |>
  mutate(
    first_imv_dttm = as_utc_datetime(.data$first_imv_dttm),
    final_imv_after_landmark_dttm = as_utc_datetime(.data$final_imv_after_landmark_dttm)
  )

analysis <- analysis_base |>
  filter(.data$landmark_eligible) |>
  rowwise() |>
  mutate(
    event_dttm = min_time(c(.data$successful_extubation_dttm, .data$death_hospice_dttm, .data$first_trach_after_landmark_dttm)),
    event_dttm = ifelse(!is.na(.data$event_dttm) & .data$event_dttm <= .data$censor_dttm, .data$event_dttm, as.POSIXct(NA, tz = "UTC")),
    event_dttm = as_utc_datetime(.data$event_dttm),
    fg_status = dplyr::case_when(
      !is.na(.data$event_dttm) & !is.na(.data$successful_extubation_dttm) & .data$event_dttm == .data$successful_extubation_dttm ~ 1L,
      !is.na(.data$event_dttm) & !is.na(.data$death_hospice_dttm) & .data$event_dttm == .data$death_hospice_dttm ~ 2L,
      !is.na(.data$event_dttm) & !is.na(.data$first_trach_after_landmark_dttm) & .data$event_dttm == .data$first_trach_after_landmark_dttm ~ 3L,
      TRUE ~ 0L
    ),
    observed_dttm = dplyr::coalesce(.data$event_dttm, .data$censor_dttm),
    fg_time_days = as.numeric(difftime(.data$observed_dttm, .data$landmark_dttm, units = "days")),
    fg_time_hours = as.numeric(difftime(.data$observed_dttm, .data$landmark_dttm, units = "hours"))
  ) |>
  ungroup() |>
  filter(!is.na(.data$fg_time_days), .data$fg_time_days >= 0) |>
  mutate(
    fg_time_days = pmax(.data$fg_time_days, 1 / (24 * 60)),
    fg_time_hours = pmax(.data$fg_time_hours, 1 / 60),
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
    ),
    outcome_status = factor(
      dplyr::case_when(
        .data$fg_status == 1L ~ "successful_extubation",
        .data$fg_status == 2L ~ "death_hospice_discharge",
        .data$fg_status == 3L ~ "tracheostomy",
        TRUE ~ "censored_no_event_before_discharge"
      ),
      levels = c(
        "successful_extubation",
        "death_hospice_discharge",
        "tracheostomy",
        "censored_no_event_before_discharge"
      )
    ),
    hours_first_imv_to_landmark = as.numeric(difftime(.data$landmark_dttm, .data$first_imv_dttm, units = "hours")),
    hours_landmark_to_final_imv = as.numeric(difftime(.data$final_imv_after_landmark_dttm, .data$landmark_dttm, units = "hours")),
    hours_landmark_to_event_or_censor = .data$fg_time_hours
  )

landmark_flow <- tibble::tibble(
  site_name = site_name,
  step = c(
    "OHCA ICU admissions with admission Tmax",
    "Any respiratory support record",
    "Latest respiratory support at or before ICU hour 24 is IMV",
    "Exclude tracheostomy before or at ICU hour 24",
    "Exclude death/hospice before or at ICU hour 24",
    "Final IMV-24 landmark cohort"
  ),
  n = c(
    nrow(cohort),
    sum(!is.na(analysis_base$n_resp_records)),
    sum(analysis_base$imv_at_landmark %in% TRUE, na.rm = TRUE),
    sum(analysis_base$imv_at_landmark %in% TRUE & !(analysis_base$trach_before_or_at_landmark %in% TRUE), na.rm = TRUE),
    sum(analysis_base$imv_at_landmark %in% TRUE & !(analysis_base$trach_before_or_at_landmark %in% TRUE) & !(analysis_base$pre_landmark_death_hospice %in% TRUE), na.rm = TRUE),
    nrow(analysis)
  ),
  definition = "Time zero is ICU hour 24 among OHCA ICU admissions whose latest respiratory-support record at or before ICU hour 24 is IMV. Prevalent tracheostomy and death/hospice before or at the landmark are excluded."
)

summary_tbl <- analysis |>
  mutate(status_label = dplyr::case_when(
    .data$fg_status == 1L ~ "Successful extubation",
    .data$fg_status == 2L ~ "Death/hospice discharge",
    .data$fg_status == 3L ~ "Tracheostomy",
    TRUE ~ "Censored/no event before discharge"
  )) |>
  count(.data$status_label, name = "n") |>
  mutate(
    site_name = site_name,
    pct = 100 * .data$n / sum(.data$n),
    definition = "Successful extubation is final IMV stop followed by a non-IMV respiratory-support record, with no subsequent IMV and no death/hospice during hospitalization. Death/hospice uses death_dttm when available, expired/death discharge with last vital fallback, or hospice/comfort discharge at discharge time. Tracheostomy uses the CLIF respiratory-support tracheostomy flag only. Follow-up starts at ICU hour 24 and ends at first event or hospital discharge/last observed support or vital.",
    .before = 1
  )

evidence_summary <- analysis |>
  group_by(.data$outcome_status) |>
  summarise(
    n = dplyr::n(),
    anoxic_brain_injury_dx_n = sum(.data$anoxic_brain_injury_dx, na.rm = TRUE),
    brain_death_dx_n = sum(.data$brain_death_dx, na.rm = TRUE),
    hospital_death_n = sum(.data$hospital_death == 1L, na.rm = TRUE),
    hospice_discharge_n = sum(.data$hospice_discharge, na.rm = TRUE),
    median_tmax_c = median(.data$tmax_mean_c, na.rm = TRUE),
    median_hours_landmark_to_event_or_censor = median(.data$hours_landmark_to_event_or_censor, na.rm = TRUE),
    median_hours_landmark_to_final_imv = median(.data$hours_landmark_to_final_imv, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    anoxic_brain_injury_dx_pct = safe_pct(.data$anoxic_brain_injury_dx_n, .data$n),
    brain_death_dx_pct = safe_pct(.data$brain_death_dx_n, .data$n),
    hospital_death_pct = safe_pct(.data$hospital_death_n, .data$n),
    hospice_discharge_pct = safe_pct(.data$hospice_discharge_n, .data$n),
    site_name = site_name,
    .before = 1
  )

outcome_levels <- levels(analysis$outcome_status)
outcome_table <- bind_rows(
  add_outcome_continuous(analysis, "Demographics", "Age, years", "age_at_admission", outcome_levels),
  add_outcome_categorical(analysis, "Demographics", "Sex", "sex_group", outcome_levels),
  add_outcome_categorical(analysis, "Demographics", "Race", "race_group", outcome_levels),
  add_outcome_continuous(analysis, "Admission exposure", "Admission-day Tmax, C", "tmax_mean_c", outcome_levels),
  add_outcome_categorical(analysis, "OHCA mechanism", "POA diagnosis mechanism", "ohca_mechanism", outcome_levels),
  add_outcome_continuous(analysis, "IMV course", "Hours from landmark to final IMV record", "hours_landmark_to_final_imv", outcome_levels),
  add_outcome_continuous(analysis, "Follow-up", "Hours from landmark to event/censor", "hours_landmark_to_event_or_censor", outcome_levels),
  add_outcome_continuous(analysis, "Respiratory support data", "Respiratory support records, n", "n_resp_records", outcome_levels),
  add_outcome_continuous(analysis, "Respiratory support data", "IMV records, n", "n_imv_records", outcome_levels),
  add_outcome_continuous(analysis, "Respiratory support data", "Tracheostomy records, n", "n_trach_records", outcome_levels),
  add_outcome_binary(analysis, "ICD evidence", "Anoxic brain injury ICD", "anoxic_brain_injury_dx", outcome_levels),
  add_outcome_binary(analysis, "ICD evidence", "Brain death ICD", "brain_death_dx", outcome_levels),
  add_outcome_binary(analysis, "Outcome", "In-hospital mortality", "hospital_death", outcome_levels),
  add_outcome_binary(analysis, "Outcome", "Hospice/comfort discharge", "hospice_discharge", outcome_levels)
) |>
  mutate(site_name = site_name, .before = 1)

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

adjust_terms <- c("age_at_admission", "sex_group", "race_group")
fg_results <- list(
  run_fine_gray(analysis, 1L, "Successful extubation", adjust_terms),
  run_fine_gray(analysis, 2L, "Death/hospice discharge", adjust_terms),
  run_fine_gray(analysis, 3L, "Tracheostomy", adjust_terms)
)

models <- purrr::map_dfr(fg_results, "summary") |>
  mutate(site_name = site_name, .before = 1)
model_coefficients <- purrr::map_dfr(fg_results, "coefficients") |>
  mutate(site_name = site_name, .before = 1)
model_vcov <- purrr::map_dfr(fg_results, "vcov") |>
  mutate(site_name = site_name, .before = 1)

patient_audit <- analysis |>
  select(dplyr::any_of(c(
    "hospitalization_id", "patient_id", "admission_date", "first_icu_in", "landmark_dttm",
    "tmax_mean_c", "age_at_admission", "sex_group", "race_group", "ohca_mechanism",
    "anoxic_brain_injury_dx", "brain_death_dx", "anoxic_brain_injury_codes",
    "first_imv_dttm", "last_resp_before_landmark_dttm", "imv_at_landmark",
    "final_imv_after_landmark_dttm", "first_non_imv_after_final_imv_dttm",
    "successful_extubation_dttm", "death_dttm_final", "hospice_dttm", "death_hospice_dttm",
    "first_trach_after_landmark_dttm", "censor_dttm", "event_dttm", "fg_status",
    "fg_time_hours", "fg_time_days", "n_resp_records", "n_imv_records", "n_trach_records"
  ))) |>
  mutate(site_name = site_name, .before = 1)

readr::write_csv(summary_tbl, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_summary.csv"))
readr::write_csv(landmark_flow, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_landmark_flow.csv"))
readr::write_csv(cif_curves, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_cif_curves.csv"))
readr::write_csv(evidence_summary, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_evidence_summary.csv"))
readr::write_csv(outcome_table, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_table_by_outcome.csv"))
readr::write_csv(models, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_fine_gray_models.csv"))
readr::write_csv(model_coefficients, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_coefficients.csv"))
readr::write_csv(model_vcov, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_vcov.csv"))
readr::write_csv(patient_audit, file.path(OUTPUT_DIR, "ohca_icu_imv24_time_to_event_patient_audit.csv"))

print(landmark_flow)
print(summary_tbl)
print(models)

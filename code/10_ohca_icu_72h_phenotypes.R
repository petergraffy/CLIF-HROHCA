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
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(lubridate)
  library(readr)
  library(splines)
  library(stringr)
  library(tidyr)
})

OHCA_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
OHCA_SUMMARY_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024_summary.csv")
OHCA_RESTRICTION_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_pathway_timing_restriction_summary.csv")
TMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
RMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_rmax_2018_2024.parquet")
OUTPUT_DIR <- file.path(repo_root, "output", "final", "ohca_icu_phenotypes")
FIGURE_DIR <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(OUTPUT_DIR, c(
  "ohca_icu_72h_phenotype_heat_thresholds.csv",
  "ohca_icu_72h_heat_by_phenotype_summary.csv",
  "ohca_icu_72h_phenotype_heat_models.csv",
  "ohca_icu_72h_phenotype_multinomial_heat_models.csv",
  "ohca_icu_72h_phenotype_spline_heat_models.csv",
  "ohca_icu_72h_phenotype_spline_heat_curves.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_model.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_curves.csv",
  "ohca_icu_72h_phenotype_heat_models_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_multinomial_heat_models_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_spline_heat_models_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_spline_heat_curves_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_model_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_curves_mechanism_adjusted.csv",
  "ohca_icu_72h_middle_imv_mental_state_summary.csv",
  "ohca_icu_72h_middle_imv_mental_state_heat_summary.csv",
  "ohca_icu_72h_middle_imv_mental_state_binary_models.csv",
  "ohca_icu_72h_middle_imv_mental_state_continuous_models.csv",
  "ohca_icu_72h_middle_imv_mental_state_spline_models.csv",
  "ohca_icu_72h_middle_imv_mental_state_spline_curves.csv",
  "ohca_icu_72h_middle_imv_mental_state_patient_audit.csv",
  "ohca_icu_72h_gcs_hourly_by_phenotype.csv"
)))

WINDOW_HOURS <- 72
LATE_START_HOUR <- 24
EXTUBATION_ASSESSMENT_START_HOUR <- 48

config <- load_project_config(repo_root)
validate_site_name(config, repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)
site_name <- config$site_name %||% Sys.getenv("CLIF_SITE_NAME", unset = "unknown_site")

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
    ifelse(stringr::str_detect(dx_clean, "^[0-9]"), "ICD9", "")
  )
  ifelse(nchar(explicit) > 0, explicit, inferred)
}

safe_last <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  tail(x, 1)
}

yesno <- function(x) {
  ifelse(is.na(x), FALSE, x)
}

clean_binary_text <- function(x) {
  stringr::str_to_lower(stringr::str_squish(tidyr::replace_na(as.character(x), "")))
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(stringr::str_squish(x))] <- "Missing"
  stringr::str_squish(x)
}

ensure_columns <- function(df, columns) {
  if (is.null(df)) df <- tibble::tibble()
  for (column in columns) {
    if (!column %in% names(df)) df[[column]] <- NA
  }
  df
}

fmt_n_pct <- function(n, denom) {
  if (is.na(denom) || denom <= 0) return(NA_character_)
  sprintf("%s (%.1f%%)", format(as.integer(n), big.mark = ","), 100 * n / denom)
}

fmt_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  sprintf(
    "%.1f [%.1f, %.1f]",
    stats::median(x, na.rm = TRUE),
    stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  )
}

fmt_mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
}

fmt_median_iqr_n <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  sprintf("%s; n=%s", fmt_median_iqr(x), format(length(x), big.mark = ","))
}

fmt_yes <- function(x) {
  x <- ifelse(is.na(x), NA, as.integer(x %in% c(TRUE, 1L, "1", "TRUE", "true", "Yes", "yes")))
  fmt_n_pct(sum(x == 1L, na.rm = TRUE), sum(!is.na(x)))
}

add_table1_continuous <- function(df, section, characteristic, variable) {
  if (!variable %in% names(df)) return(tibble::tibble())
  tibble::tibble(section = section, characteristic = characteristic, level = "", value = fmt_median_iqr(df[[variable]]))
}

add_table1_binary <- function(df, section, characteristic, variable) {
  if (!variable %in% names(df)) return(tibble::tibble())
  tibble::tibble(section = section, characteristic = characteristic, level = "Yes", value = fmt_yes(df[[variable]]))
}

add_table1_categorical <- function(df, section, characteristic, variable) {
  if (!variable %in% names(df)) return(tibble::tibble())
  counts <- sort(table(safe_chr(df[[variable]])), decreasing = TRUE)
  tibble::tibble(
    section = section,
    characteristic = characteristic,
    level = names(counts),
    value = vapply(as.integer(counts), fmt_n_pct, character(1), denom = nrow(df))
  )
}

add_table2_continuous <- function(df, section, characteristic, variable, phenotype_levels) {
  if (!variable %in% names(df)) return(tibble::tibble())
  vals <- vapply(phenotype_levels, function(ph) fmt_median_iqr(df[[variable]][df$phenotype == ph]), character(1))
  tibble::tibble(section = section, characteristic = characteristic, level = "") |>
    bind_cols(as_tibble(as.list(vals)))
}

add_table2_mean_sd <- function(df, section, characteristic, variable, phenotype_levels) {
  if (!variable %in% names(df)) return(tibble::tibble())
  vals <- vapply(phenotype_levels, function(ph) fmt_mean_sd(df[[variable]][df$phenotype == ph]), character(1))
  tibble::tibble(section = section, characteristic = characteristic, level = "") |>
    bind_cols(as_tibble(as.list(vals)))
}

add_table2_binary <- function(df, section, characteristic, variable, phenotype_levels) {
  if (!variable %in% names(df)) return(tibble::tibble())
  vals <- vapply(phenotype_levels, function(ph) fmt_yes(df[[variable]][df$phenotype == ph]), character(1))
  tibble::tibble(section = section, characteristic = characteristic, level = "Yes") |>
    bind_cols(as_tibble(as.list(vals)))
}

add_table2_categorical <- function(df, section, characteristic, variable, phenotype_levels) {
  if (!variable %in% names(df)) return(tibble::tibble())
  levels_all <- names(sort(table(safe_chr(df[[variable]])), decreasing = TRUE))
  rows <- lapply(levels_all, function(lvl) {
    vals <- vapply(phenotype_levels, function(ph) {
      sub <- safe_chr(df[[variable]][df$phenotype == ph])
      fmt_n_pct(sum(sub == lvl, na.rm = TRUE), length(sub))
    }, character(1))
    tibble::tibble(section = section, characteristic = characteristic, level = lvl) |>
      bind_cols(as_tibble(as.list(vals)))
  })
  bind_rows(rows)
}

calc_pao2_from_spo2 <- function(spo2) {
  s <- suppressWarnings(as.numeric(spo2)) / 100
  out <- rep(NA_real_, length(s))
  ok <- is.finite(s) & s > 0 & s < 1
  a <- 11700 / ((1 / s[ok]) - 1)
  b <- sqrt(50^3 + a^2)
  out[ok] <- ((b + a)^(1 / 3)) - ((b - a)^(1 / 3))
  out
}

ensure_numeric_cols <- function(df, cols) {
  for (col in cols) {
    if (!col %in% names(df)) df[[col]] <- NA_real_
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  df
}

standardize_vaso_dose <- function(med_category, med_dose, med_dose_unit, weight_kg) {
  med_category <- stringr::str_to_lower(as.character(med_category))
  med_dose_unit <- stringr::str_to_lower(as.character(med_dose_unit))
  med_dose <- suppressWarnings(as.numeric(med_dose))
  weight_kg <- dplyr::coalesce(suppressWarnings(as.numeric(weight_kg)), 80)
  factor <- dplyr::case_when(
    med_category %in% c("norepinephrine", "epinephrine", "phenylephrine", "dopamine", "dobutamine", "milrinone") &
      med_dose_unit == "mcg/kg/min" ~ 1,
    med_category %in% c("norepinephrine", "epinephrine", "phenylephrine", "dopamine", "dobutamine", "milrinone") &
      med_dose_unit %in% c("mcg/kg/hr", "mcg/kg/hour") ~ 1 / 60,
    med_category %in% c("norepinephrine", "epinephrine", "phenylephrine", "dopamine", "dobutamine", "milrinone") &
      med_dose_unit %in% c("mg/kg/hr", "mg/kg/hour") ~ 1000 / 60,
    med_category %in% c("norepinephrine", "epinephrine", "phenylephrine", "dopamine", "dobutamine", "milrinone") &
      med_dose_unit == "mcg/min" ~ 1 / weight_kg,
    med_category %in% c("norepinephrine", "epinephrine", "phenylephrine", "dopamine", "dobutamine", "milrinone") &
      med_dose_unit %in% c("mg/hr", "mg/hour") ~ (1000 / 60) / weight_kg,
    med_category == "vasopressin" & med_dose_unit == "units/min" ~ 1,
    med_category == "vasopressin" & med_dose_unit %in% c("units/hr", "units/hour") ~ 1 / 60,
    med_category == "vasopressin" & med_dose_unit == "milliunits/min" ~ 1 / 1000,
    med_category == "vasopressin" & med_dose_unit %in% c("milliunits/hr", "milliunits/hour") ~ 1 / 1000 / 60,
    TRUE ~ NA_real_
  )
  med_dose * factor
}

calculate_sofa_windows <- function(cohort_df, vitals_df, labs_df, support_df, med_admin_df, scores_df, windows = c(24, 48, 72)) {
  ids <- unique(as.character(cohort_df$hospitalization_id))
  base <- cohort_df |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      first_icu_in = as_utc_datetime(.data$first_icu_in)
    )

  if (is.null(vitals_df) || nrow(vitals_df) == 0) vitals_df <- tibble::tibble()
  if (is.null(labs_df) || nrow(labs_df) == 0) labs_df <- tibble::tibble()
  if (is.null(support_df) || nrow(support_df) == 0) support_df <- tibble::tibble()
  if (is.null(med_admin_df) || nrow(med_admin_df) == 0) med_admin_df <- tibble::tibble()
  if (is.null(scores_df) || nrow(scores_df) == 0) scores_df <- tibble::tibble()
  support_df <- ensure_columns(support_df, c("device_category", "fio2_set"))
  med_admin_df <- ensure_columns(med_admin_df, c("med_dose", "med_dose_unit"))

  vitals_min <- if (all(c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value") %in% names(vitals_df))) {
    vitals_df |>
      transmute(
        hospitalization_id = as.character(.data$hospitalization_id),
        recorded_dttm = as_utc_datetime(.data$recorded_dttm),
        vital_category = stringr::str_to_lower(as.character(.data$vital_category)),
        vital_value = suppressWarnings(as.numeric(.data$vital_value))
      ) |>
      filter(.data$hospitalization_id %in% ids, !is.na(.data$recorded_dttm))
  } else tibble::tibble(hospitalization_id = character(), recorded_dttm = as.POSIXct(character(), tz = "UTC"), vital_category = character(), vital_value = numeric())

  weights <- vitals_min |>
    filter(.data$vital_category == "weight_kg", .data$vital_value > 10, .data$vital_value < 500) |>
    group_by(.data$hospitalization_id) |>
    summarise(weight_kg = median(.data$vital_value, na.rm = TRUE), .groups = "drop")

  support_min <- if (all(c("hospitalization_id", "recorded_dttm") %in% names(support_df))) {
    support_df |>
      transmute(
        hospitalization_id = as.character(.data$hospitalization_id),
        recorded_dttm = as_utc_datetime(.data$recorded_dttm),
        device_category = stringr::str_to_lower(tidyr::replace_na(as.character(.data$device_category), "")),
        fio2_set = suppressWarnings(as.numeric(.data$fio2_set))
      ) |>
      filter(.data$hospitalization_id %in% ids, !is.na(.data$recorded_dttm))
  } else tibble::tibble(hospitalization_id = character(), recorded_dttm = as.POSIXct(character(), tz = "UTC"), device_category = character(), fio2_set = numeric())

  labs_min <- if (all(c("hospitalization_id", "lab_result_dttm", "lab_category", "lab_value_numeric") %in% names(labs_df))) {
    labs_df |>
      transmute(
        hospitalization_id = as.character(.data$hospitalization_id),
        lab_result_dttm = as_utc_datetime(.data$lab_result_dttm),
        lab_category = stringr::str_to_lower(as.character(.data$lab_category)),
        lab_value_numeric = suppressWarnings(as.numeric(.data$lab_value_numeric))
      ) |>
      filter(.data$hospitalization_id %in% ids, !is.na(.data$lab_result_dttm))
  } else tibble::tibble(hospitalization_id = character(), lab_result_dttm = as.POSIXct(character(), tz = "UTC"), lab_category = character(), lab_value_numeric = numeric())

  med_min <- if (all(c("hospitalization_id", "admin_dttm", "med_category", "med_group", "mar_action_group") %in% names(med_admin_df))) {
    med_admin_df |>
      transmute(
        hospitalization_id = as.character(.data$hospitalization_id),
        admin_dttm = as_utc_datetime(.data$admin_dttm),
        med_category = stringr::str_to_lower(as.character(.data$med_category)),
        med_group = stringr::str_to_lower(as.character(.data$med_group)),
        mar_action_group = stringr::str_to_lower(as.character(.data$mar_action_group)),
        med_dose = suppressWarnings(as.numeric(.data$med_dose)),
        med_dose_unit = stringr::str_to_lower(as.character(.data$med_dose_unit))
      ) |>
      filter(.data$hospitalization_id %in% ids, !is.na(.data$admin_dttm))
  } else tibble::tibble(hospitalization_id = character(), admin_dttm = as.POSIXct(character(), tz = "UTC"), med_category = character(), med_group = character(), mar_action_group = character(), med_dose = numeric(), med_dose_unit = character())

  scores_min <- if (all(c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value") %in% names(scores_df))) {
    scores_df |>
      transmute(
        hospitalization_id = as.character(.data$hospitalization_id),
        recorded_dttm = as_utc_datetime(.data$recorded_dttm),
        assessment_category = stringr::str_to_lower(as.character(.data$assessment_category)),
        numerical_value = suppressWarnings(as.numeric(.data$numerical_value))
      ) |>
      filter(.data$hospitalization_id %in% ids, !is.na(.data$recorded_dttm))
  } else tibble::tibble(hospitalization_id = character(), recorded_dttm = as.POSIXct(character(), tz = "UTC"), assessment_category = character(), numerical_value = numeric())

  dplyr::bind_rows(lapply(windows, function(window_hour) {
    window_start <- max(0, window_hour - 24)
    window_end <- window_hour
    vitals_window <- vitals_min |>
      inner_join(base, by = "hospitalization_id") |>
      mutate(icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours"))) |>
      filter(.data$icu_hour >= window_start, .data$icu_hour <= window_end) |>
      filter(
        (.data$vital_category == "map" & .data$vital_value >= 20 & .data$vital_value <= 250) |
          (.data$vital_category == "spo2" & .data$vital_value >= 50 & .data$vital_value <= 100)
      ) |>
      group_by(.data$hospitalization_id, .data$vital_category) |>
      summarise(worst_val = min(.data$vital_value, na.rm = TRUE), .groups = "drop") |>
      tidyr::pivot_wider(names_from = "vital_category", values_from = "worst_val") |>
      mutate(pao2_imputed = calc_pao2_from_spo2(.data$spo2))

    resp_window <- support_min |>
      inner_join(base, by = "hospitalization_id") |>
      mutate(
        icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")),
        fio2_std = dplyr::case_when(
          is.na(.data$fio2_set) ~ NA_real_,
          .data$fio2_set > 1 & .data$fio2_set <= 100 ~ .data$fio2_set / 100,
          .data$fio2_set >= 0.21 & .data$fio2_set <= 1 ~ .data$fio2_set,
          TRUE ~ NA_real_
        )
      ) |>
      filter(.data$icu_hour >= window_start, .data$icu_hour <= window_end) |>
      group_by(.data$hospitalization_id) |>
      summarise(
        fio2_max = suppressWarnings(max(.data$fio2_std, na.rm = TRUE)),
        has_imv = any(stringr::str_detect(.data$device_category, "imv|vent"), na.rm = TRUE),
        has_nippv = any(stringr::str_detect(.data$device_category, "nippv"), na.rm = TRUE),
        has_cpap = any(stringr::str_detect(.data$device_category, "cpap"), na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        fio2_max = ifelse(is.finite(.data$fio2_max), .data$fio2_max, NA_real_),
        resp_support_max = dplyr::case_when(
          .data$has_imv ~ "Vent",
          .data$has_nippv ~ "NIPPV",
          .data$has_cpap ~ "CPAP",
          TRUE ~ "Other"
        )
      )

    labs_window <- labs_min |>
      inner_join(base, by = "hospitalization_id") |>
      mutate(icu_hour = as.numeric(difftime(.data$lab_result_dttm, .data$first_icu_in, units = "hours"))) |>
      filter(.data$icu_hour >= window_start, .data$icu_hour <= window_end) |>
      filter(
        (.data$lab_category == "creatinine" & .data$lab_value_numeric >= 0 & .data$lab_value_numeric <= 25) |
          (.data$lab_category == "bilirubin_total" & .data$lab_value_numeric >= 0 & .data$lab_value_numeric <= 100) |
          (.data$lab_category == "po2_arterial" & .data$lab_value_numeric >= 20 & .data$lab_value_numeric <= 800) |
          (.data$lab_category == "platelet_count" & .data$lab_value_numeric >= 0 & .data$lab_value_numeric <= 3000)
      )

    labs_combined <- base |>
      select("hospitalization_id") |>
      left_join(labs_window |> filter(.data$lab_category == "creatinine") |> group_by(.data$hospitalization_id) |> summarise(creatinine = max(.data$lab_value_numeric, na.rm = TRUE), .groups = "drop"), by = "hospitalization_id") |>
      left_join(labs_window |> filter(.data$lab_category == "bilirubin_total") |> group_by(.data$hospitalization_id) |> summarise(bilirubin_total = max(.data$lab_value_numeric, na.rm = TRUE), .groups = "drop"), by = "hospitalization_id") |>
      left_join(labs_window |> filter(.data$lab_category == "po2_arterial") |> group_by(.data$hospitalization_id) |> summarise(po2_arterial = min(.data$lab_value_numeric, na.rm = TRUE), .groups = "drop"), by = "hospitalization_id") |>
      left_join(labs_window |> filter(.data$lab_category == "platelet_count") |> group_by(.data$hospitalization_id) |> summarise(platelet_count = min(.data$lab_value_numeric, na.rm = TRUE), .groups = "drop"), by = "hospitalization_id") |>
      mutate(across(c("creatinine", "bilirubin_total", "po2_arterial", "platelet_count"), ~ ifelse(is.infinite(.x), NA_real_, .x)))

    gcs_window <- scores_min |>
      inner_join(base, by = "hospitalization_id") |>
      mutate(icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours"))) |>
      filter(
        .data$icu_hour >= window_start, .data$icu_hour <= window_end,
        .data$assessment_category == "gcs_total",
        .data$numerical_value >= 3, .data$numerical_value <= 15
      ) |>
      group_by(.data$hospitalization_id) |>
      summarise(min_gcs_score = min(.data$numerical_value, na.rm = TRUE), .groups = "drop") |>
      mutate(min_gcs_score = ifelse(is.infinite(.data$min_gcs_score), NA_real_, .data$min_gcs_score))

    meds_window <- med_min |>
      inner_join(base, by = "hospitalization_id") |>
      mutate(icu_hour = as.numeric(difftime(.data$admin_dttm, .data$first_icu_in, units = "hours"))) |>
      filter(
        .data$icu_hour >= window_start, .data$icu_hour <= window_end,
        .data$mar_action_group == "administered",
        .data$med_category %in% c("norepinephrine", "epinephrine", "phenylephrine", "vasopressin", "dopamine", "dobutamine", "milrinone"),
        .data$med_dose > 0
      ) |>
      left_join(weights, by = "hospitalization_id") |>
      mutate(dose_converted = standardize_vaso_dose(.data$med_category, .data$med_dose, .data$med_dose_unit, .data$weight_kg)) |>
      filter(!is.na(.data$dose_converted)) |>
      group_by(.data$hospitalization_id, .data$med_category) |>
      summarise(max_dose = max(.data$dose_converted, na.rm = TRUE), .groups = "drop") |>
      tidyr::pivot_wider(names_from = "med_category", values_from = "max_dose")

    sofa_data <- base |>
      select("hospitalization_id") |>
      left_join(vitals_window, by = "hospitalization_id") |>
      left_join(resp_window, by = "hospitalization_id") |>
      left_join(labs_combined, by = "hospitalization_id") |>
      left_join(gcs_window, by = "hospitalization_id") |>
      left_join(meds_window, by = "hospitalization_id") |>
      ensure_numeric_cols(c("map", "spo2", "pao2_imputed", "fio2_max", "po2_arterial", "platelet_count", "bilirubin_total", "creatinine", "min_gcs_score", "dopamine", "epinephrine", "norepinephrine", "dobutamine"))

    sofa_data |>
      mutate(
        p_f = ifelse(!is.na(.data$po2_arterial) & !is.na(.data$fio2_max) & .data$fio2_max > 0, .data$po2_arterial / .data$fio2_max, NA_real_),
        p_f_imputed = ifelse(!is.na(.data$pao2_imputed) & !is.na(.data$fio2_max) & .data$fio2_max > 0, .data$pao2_imputed / .data$fio2_max, NA_real_),
        sofa_cv = dplyr::case_when(
          .data$dopamine > 15 | .data$epinephrine > 0.1 | .data$norepinephrine > 0.1 ~ 4,
          .data$dopamine > 5 | (.data$epinephrine > 0 & .data$epinephrine <= 0.1) | (.data$norepinephrine > 0 & .data$norepinephrine <= 0.1) ~ 3,
          (.data$dopamine > 0 & .data$dopamine <= 5) | .data$dobutamine > 0 ~ 2,
          .data$map < 70 ~ 1,
          TRUE ~ 0
        ),
        sofa_coag = dplyr::case_when(.data$platelet_count < 20 ~ 4, .data$platelet_count < 50 ~ 3, .data$platelet_count < 100 ~ 2, .data$platelet_count < 150 ~ 1, TRUE ~ 0),
        sofa_liver = dplyr::case_when(.data$bilirubin_total >= 12 ~ 4, .data$bilirubin_total >= 6 ~ 3, .data$bilirubin_total >= 2 ~ 2, .data$bilirubin_total >= 1.2 ~ 1, TRUE ~ 0),
        sofa_renal = dplyr::case_when(.data$creatinine >= 5 ~ 4, .data$creatinine >= 3.5 ~ 3, .data$creatinine >= 2 ~ 2, .data$creatinine >= 1.2 ~ 1, TRUE ~ 0),
        sofa_resp = dplyr::case_when(
          !is.na(.data$p_f) & .data$p_f < 100 & .data$resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 4,
          !is.na(.data$p_f) & .data$p_f < 200 & .data$resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 3,
          !is.na(.data$p_f) & .data$p_f < 300 ~ 2,
          !is.na(.data$p_f) & .data$p_f < 400 ~ 1,
          !is.na(.data$p_f_imputed) & .data$p_f_imputed < 100 & .data$resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 4,
          !is.na(.data$p_f_imputed) & .data$p_f_imputed < 200 & .data$resp_support_max %in% c("Vent", "NIPPV", "CPAP") ~ 3,
          !is.na(.data$p_f_imputed) & .data$p_f_imputed < 300 ~ 2,
          !is.na(.data$p_f_imputed) & .data$p_f_imputed < 400 ~ 1,
          TRUE ~ 0
        ),
        sofa_cns = dplyr::case_when(.data$min_gcs_score < 6 ~ 4, .data$min_gcs_score <= 9 ~ 3, .data$min_gcs_score <= 12 ~ 2, .data$min_gcs_score <= 14 ~ 1, TRUE ~ 0),
        sofa_any_data = !is.na(.data$map) |
          !is.na(.data$spo2) |
          !is.na(.data$fio2_max) |
          !is.na(.data$po2_arterial) |
          !is.na(.data$platelet_count) |
          !is.na(.data$bilirubin_total) |
          !is.na(.data$creatinine) |
          !is.na(.data$min_gcs_score) |
          !is.na(.data$dopamine) |
          !is.na(.data$epinephrine) |
          !is.na(.data$norepinephrine) |
          !is.na(.data$dobutamine),
        sofa_total = ifelse(
          .data$sofa_any_data,
          .data$sofa_cv + .data$sofa_coag + .data$sofa_liver + .data$sofa_renal + .data$sofa_resp + .data$sofa_cns,
          NA_real_
        ),
        sofa_window_hours = window_hour
      ) |>
      select("hospitalization_id", "sofa_window_hours", "sofa_total", "sofa_cv", "sofa_coag", "sofa_liver", "sofa_renal", "sofa_resp", "sofa_cns")
  }))
}

fmt_or <- function(beta, se) {
  if (!is.finite(beta) || !is.finite(se)) return(c(NA_real_, NA_real_, NA_real_))
  c(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
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
      if (sum(is.finite(x)) > 0 && length(unique(x[is.finite(x)])) >= 2) keep <- c(keep, term) else omitted <- c(omitted, paste0(term, "_constant"))
    } else {
      if (length(unique(stats::na.omit(x))) >= 2) keep <- c(keep, term) else omitted <- c(omitted, paste0(term, "_one_level"))
    }
  }
  list(keep = keep, omitted = omitted)
}

run_logistic <- function(df, outcome, exposure, label, adjust_terms = character()) {
  base <- df |>
    filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]]))
  events <- sum(base[[outcome]] == 1, na.rm = TRUE)
  nonevents <- sum(base[[outcome]] == 0, na.rm = TRUE)
  if (nrow(base) < 25 || events < 5 || nonevents < 5 || length(unique(base[[exposure]])) < 2) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      events = events,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Insufficient events, non-events, sample size, or exposure variation"
    ))
  }
  term_info <- choose_terms(base, adjust_terms)
  formula_txt <- paste(outcome, "~", paste(c(exposure, term_info$keep), collapse = " + "))
  fit <- tryCatch(
    stats::glm(stats::as.formula(formula_txt), data = base, family = stats::binomial()),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      events = events,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = conditionMessage(fit)
    ))
  }
  coef_tbl <- summary(fit)$coefficients
  if (!exposure %in% rownames(coef_tbl)) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(stats::model.frame(fit)),
      events = sum(stats::model.frame(fit)[[outcome]] == 1, na.rm = TRUE),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Exposure term unavailable after model fit"
    ))
  }
  if (
    !is.finite(coef_tbl[exposure, "Std. Error"]) ||
      coef_tbl[exposure, "Std. Error"] <= 1e-6 ||
      coef_tbl[exposure, "Std. Error"] > 10 ||
      !is.finite(coef_tbl[exposure, "Estimate"])
  ) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(stats::model.frame(fit)),
      events = sum(stats::model.frame(fit)[[outcome]] == 1, na.rm = TRUE),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Potential separation or numerical instability in logistic exposure estimate"
    ))
  }
  ors <- fmt_or(coef_tbl[exposure, "Estimate"], coef_tbl[exposure, "Std. Error"])
  tibble::tibble(
    outcome = outcome,
    model = label,
    exposure = exposure,
    n = nrow(stats::model.frame(fit)),
    events = sum(stats::model.frame(fit)[[outcome]] == 1, na.rm = TRUE),
    odds_ratio = ors[[1]],
    ci_low = ors[[2]],
    ci_high = ors[[3]],
    p_value = coef_tbl[exposure, "Pr(>|z|)"],
    covariates = paste(term_info$keep, collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
}

run_linear <- function(df, outcome, exposure, label, adjust_terms = character()) {
  base <- df |>
    filter(!is.na(.data[[outcome]]), is.finite(.data[[outcome]]), !is.na(.data[[exposure]]), is.finite(.data[[exposure]]))
  if (nrow(base) < 25 || length(unique(base[[outcome]])) < 2 || length(unique(base[[exposure]])) < 2) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      outcome_median = median(base[[outcome]], na.rm = TRUE),
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Insufficient sample size or outcome/exposure variation"
    ))
  }
  term_info <- choose_terms(base, adjust_terms)
  formula_txt <- paste(outcome, "~", paste(c(exposure, term_info$keep), collapse = " + "))
  fit <- tryCatch(
    stats::lm(stats::as.formula(formula_txt), data = base),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      outcome_median = median(base[[outcome]], na.rm = TRUE),
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = conditionMessage(fit)
    ))
  }
  coef_tbl <- summary(fit)$coefficients
  if (!exposure %in% rownames(coef_tbl)) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(stats::model.frame(fit)),
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      outcome_median = median(stats::model.frame(fit)[[outcome]], na.rm = TRUE),
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Exposure term unavailable after model fit"
    ))
  }
  beta <- coef_tbl[exposure, "Estimate"]
  se <- coef_tbl[exposure, "Std. Error"]
  tibble::tibble(
    outcome = outcome,
    model = label,
    exposure = exposure,
    n = nrow(stats::model.frame(fit)),
    estimate = beta,
    ci_low = beta - 1.96 * se,
    ci_high = beta + 1.96 * se,
    p_value = coef_tbl[exposure, "Pr(>|t|)"],
    outcome_median = median(stats::model.frame(fit)[[outcome]], na.rm = TRUE),
    covariates = paste(term_info$keep, collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
}

run_multinomial <- function(df, exposure, label, adjust_terms = character()) {
  model_levels <- c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury")
  reference_level <- "regained_consciousness_extubated"
  if (!requireNamespace("nnet", quietly = TRUE)) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = NA_character_,
      reference_level = reference_level,
      n = 0L,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Package nnet unavailable"
    ))
  }
  base <- df |>
    filter(.data$phenotype %in% model_levels, !is.na(.data[[exposure]])) |>
    mutate(phenotype = stats::relevel(factor(.data$phenotype, levels = model_levels), ref = reference_level))
  if (nrow(base) < 50 || length(unique(base$phenotype)) < 3 || any(table(base$phenotype) < 5) || length(unique(base[[exposure]])) < 2) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = setdiff(model_levels, reference_level),
      reference_level = reference_level,
      n = nrow(base),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Insufficient phenotype counts, sample size, or exposure variation"
    ))
  }
  term_info <- choose_terms(base, adjust_terms)
  formula_txt <- paste("phenotype ~", paste(c(exposure, term_info$keep), collapse = " + "))
  fit <- tryCatch(
    nnet::multinom(stats::as.formula(formula_txt), data = base, trace = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = setdiff(model_levels, reference_level),
      reference_level = reference_level,
      n = nrow(base),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = conditionMessage(fit)
    ))
  }
  s <- summary(fit)
  coef_mat <- s$coefficients
  se_mat <- s$standard.errors
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, nrow = 1, dimnames = list(names(coef_mat)[[1]], names(coef_mat)))
    se_mat <- matrix(se_mat, nrow = 1, dimnames = dimnames(coef_mat))
  }
  if (
    !exposure %in% colnames(coef_mat) ||
      any(!is.finite(se_mat[, exposure])) ||
      any(se_mat[, exposure] <= 1e-6) ||
      any(abs(coef_mat[, exposure] / se_mat[, exposure]) > 8)
  ) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = rownames(coef_mat),
      reference_level = reference_level,
      n = nrow(stats::model.frame(fit)),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Potential separation or numerical instability in multinomial exposure estimate"
    ))
  }
  rows <- rownames(coef_mat)
  out <- lapply(rows, function(level) {
    if (!exposure %in% colnames(coef_mat)) return(NULL)
    beta <- coef_mat[level, exposure]
    se <- se_mat[level, exposure]
    ors <- fmt_or(beta, se)
    z <- beta / se
    tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = level,
      reference_level = reference_level,
      n = nrow(stats::model.frame(fit)),
      odds_ratio = ors[[1]],
      ci_low = ors[[2]],
      ci_high = ors[[3]],
      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = TRUE,
      skip_reason = NA_character_
    )
  })
  dplyr::bind_rows(out)
}

mode_value <- function(x) {
  x <- stats::na.omit(x)
  if (length(x) == 0L) return(NA)
  names(sort(table(x), decreasing = TRUE))[[1]]
}

make_reference_newdata <- function(df, temp_grid, adjust_terms = character()) {
  out <- tibble::tibble(tmax_mean_c = temp_grid, tmax_per5c = temp_grid / 5)
  for (term in adjust_terms) {
    if (!term %in% names(df)) next
    if (is.numeric(df[[term]])) {
      out[[term]] <- mean(df[[term]], na.rm = TRUE)
    } else {
      out[[term]] <- mode_value(df[[term]])
    }
  }
  out
}

run_spline_logistic <- function(df, outcome, label, adjust_terms = character(), spline_df = 3L) {
  base <- df |>
    filter(!is.na(.data[[outcome]]), !is.na(.data$tmax_mean_c), is.finite(.data$tmax_mean_c))
  events <- sum(base[[outcome]] == 1, na.rm = TRUE)
  nonevents <- sum(base[[outcome]] == 0, na.rm = TRUE)
  empty_summary <- tibble::tibble(
    outcome = outcome,
    model = label,
    n = nrow(base),
    events = events,
    spline_df = spline_df,
    reference_temp_c = NA_real_,
    min_temp_c = NA_real_,
    max_temp_c = NA_real_,
    linear_lrt_p_value = NA_real_,
    spline_lrt_p_value = NA_real_,
    nonlinear_lrt_p_value = NA_real_,
    aic = NA_real_,
    covariates = paste(adjust_terms, collapse = " + "),
    omitted_terms = NA_character_,
    estimable = FALSE,
    skip_reason = "Insufficient events, non-events, sample size, or temperature variation"
  )
  if (nrow(base) < 50 || events < 5 || nonevents < 5 || length(unique(base$tmax_mean_c)) < 5) {
    return(list(summary = empty_summary, curve = tibble::tibble()))
  }

  term_info <- choose_terms(base, adjust_terms)
  adjust_txt <- if (length(term_info$keep) > 0) paste("+", paste(term_info$keep, collapse = " + ")) else ""
  null_formula <- stats::as.formula(paste(outcome, "~", if (length(term_info$keep) > 0) paste(term_info$keep, collapse = " + ") else "1"))
  linear_formula <- stats::as.formula(paste(outcome, "~ tmax_per5c", adjust_txt))
  spline_formula <- stats::as.formula(paste0(outcome, " ~ splines::ns(tmax_mean_c, df = ", spline_df, ") ", adjust_txt))

  null_fit <- tryCatch(stats::glm(null_formula, data = base, family = stats::binomial()), error = function(e) e)
  linear_fit <- tryCatch(stats::glm(linear_formula, data = base, family = stats::binomial()), error = function(e) e)
  spline_fit <- tryCatch(stats::glm(spline_formula, data = base, family = stats::binomial()), error = function(e) e)
  if (inherits(null_fit, "error") || inherits(linear_fit, "error") || inherits(spline_fit, "error")) {
    reason <- paste(
      stats::na.omit(c(
        if (inherits(null_fit, "error")) paste0("null: ", conditionMessage(null_fit)),
        if (inherits(linear_fit, "error")) paste0("linear: ", conditionMessage(linear_fit)),
        if (inherits(spline_fit, "error")) paste0("spline: ", conditionMessage(spline_fit))
      )),
      collapse = "; "
    )
    empty_summary$skip_reason <- reason
    empty_summary$omitted_terms <- paste(term_info$omitted, collapse = "; ")
    return(list(summary = empty_summary, curve = tibble::tibble()))
  }

  temp_grid <- seq(
    stats::quantile(base$tmax_mean_c, 0.01, na.rm = TRUE, names = FALSE),
    stats::quantile(base$tmax_mean_c, 0.99, na.rm = TRUE, names = FALSE),
    length.out = 100
  )
  ref_temp <- stats::median(base$tmax_mean_c, na.rm = TRUE)
  newdata <- make_reference_newdata(base, temp_grid, term_info$keep)
  pred <- predict(spline_fit, newdata = newdata, type = "link", se.fit = TRUE)
  curve <- tibble::tibble(
    model = label,
    outcome = outcome,
    tmax_mean_c = temp_grid,
    predicted_probability = stats::plogis(pred$fit),
    predicted_probability_low = stats::plogis(pred$fit - 1.96 * pred$se.fit),
    predicted_probability_high = stats::plogis(pred$fit + 1.96 * pred$se.fit),
    reference_temp_c = ref_temp,
    spline_df = spline_df
  )

  linear_lrt <- tryCatch(stats::anova(null_fit, linear_fit, test = "LRT"), error = function(e) NULL)
  spline_lrt <- tryCatch(stats::anova(null_fit, spline_fit, test = "LRT"), error = function(e) NULL)
  nonlinear_lrt <- tryCatch(stats::anova(linear_fit, spline_fit, test = "LRT"), error = function(e) NULL)

  summary <- tibble::tibble(
    outcome = outcome,
    model = label,
    n = nrow(stats::model.frame(spline_fit)),
    events = sum(stats::model.frame(spline_fit)[[outcome]] == 1, na.rm = TRUE),
    spline_df = spline_df,
    reference_temp_c = ref_temp,
    min_temp_c = min(base$tmax_mean_c, na.rm = TRUE),
    max_temp_c = max(base$tmax_mean_c, na.rm = TRUE),
    linear_lrt_p_value = if (!is.null(linear_lrt)) linear_lrt$`Pr(>Chi)`[[2]] else NA_real_,
    spline_lrt_p_value = if (!is.null(spline_lrt)) spline_lrt$`Pr(>Chi)`[[2]] else NA_real_,
    nonlinear_lrt_p_value = if (!is.null(nonlinear_lrt)) nonlinear_lrt$`Pr(>Chi)`[[2]] else NA_real_,
    aic = stats::AIC(spline_fit),
    covariates = paste(term_info$keep, collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
  list(summary = summary, curve = curve)
}

run_phenotype_assignment_model <- function(df, label, adjust_terms = character(), spline_df = 3L) {
  keep_levels <- c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury")
  reference_level <- "regained_consciousness_extubated"
  base <- df |>
    filter(
      .data$phenotype %in% keep_levels,
      !is.na(.data$tmax_mean_c), is.finite(.data$tmax_mean_c),
      !is.na(.data$rmax_mean_pct), is.finite(.data$rmax_mean_pct)
    ) |>
    mutate(phenotype = stats::relevel(factor(.data$phenotype, levels = keep_levels), ref = reference_level))

  empty_summary <- tibble::tibble(
    model = label,
    n = nrow(base),
    spline_df = spline_df,
    reference_level = reference_level,
    reference_temp_c = NA_real_,
    reference_humidity_pct = NA_real_,
    min_temp_c = NA_real_,
    max_temp_c = NA_real_,
    min_humidity_pct = NA_real_,
    max_humidity_pct = NA_real_,
    temperature_lrt_p_value = NA_real_,
    humidity_lrt_p_value = NA_real_,
    temperature_nonlinear_lrt_p_value = NA_real_,
    aic = NA_real_,
    covariates = paste(c("splines::ns(tmax_mean_c, df = 3)", "splines::ns(rmax_mean_pct, df = 3)", adjust_terms), collapse = " + "),
    modeled_phenotypes = paste(keep_levels, collapse = " | "),
    excluded_phenotypes = paste(setdiff(levels(factor(df$phenotype)), keep_levels), collapse = " | "),
    omitted_terms = NA_character_,
    estimable = FALSE,
    skip_reason = "Insufficient phenotype counts, sample size, temperature variation, or humidity variation"
  )
  if (!requireNamespace("nnet", quietly = TRUE)) {
    empty_summary$skip_reason <- "Package nnet unavailable"
    return(list(summary = empty_summary, curve = tibble::tibble(), coefficients = tibble::tibble(), vcov = tibble::tibble()))
  }
  if (
    nrow(base) < 50 ||
      length(unique(base$phenotype)) < 3 ||
      any(table(base$phenotype) < 5) ||
      length(unique(base$tmax_mean_c)) < 5 ||
      length(unique(base$rmax_mean_pct)) < 5
  ) {
    return(list(summary = empty_summary, curve = tibble::tibble(), coefficients = tibble::tibble(), vcov = tibble::tibble()))
  }

  term_info <- choose_terms(base, adjust_terms)
  adjust_txt <- if (length(term_info$keep) > 0) paste("+", paste(term_info$keep, collapse = " + ")) else ""
  covar_rhs <- if (length(term_info$keep) > 0) paste(term_info$keep, collapse = " + ") else "1"
  humidity_term <- paste0("splines::ns(rmax_mean_pct, df = ", spline_df, ")")
  temp_term <- paste0("splines::ns(tmax_mean_c, df = ", spline_df, ")")
  full_formula <- stats::as.formula(paste("phenotype ~", paste(c(temp_term, humidity_term, term_info$keep), collapse = " + ")))
  no_temp_formula <- stats::as.formula(paste("phenotype ~", paste(c(humidity_term, term_info$keep), collapse = " + ")))
  no_humidity_formula <- stats::as.formula(paste("phenotype ~", paste(c(temp_term, term_info$keep), collapse = " + ")))
  linear_temp_formula <- stats::as.formula(paste("phenotype ~", paste(c("tmax_per5c", humidity_term, term_info$keep), collapse = " + ")))

  full_fit <- tryCatch(nnet::multinom(full_formula, data = base, trace = FALSE), error = function(e) e)
  no_temp_fit <- tryCatch(nnet::multinom(no_temp_formula, data = base, trace = FALSE), error = function(e) e)
  no_humidity_fit <- tryCatch(nnet::multinom(no_humidity_formula, data = base, trace = FALSE), error = function(e) e)
  linear_temp_fit <- tryCatch(nnet::multinom(linear_temp_formula, data = base, trace = FALSE), error = function(e) e)
  if (inherits(full_fit, "error") || inherits(no_temp_fit, "error") || inherits(no_humidity_fit, "error") || inherits(linear_temp_fit, "error")) {
    reason <- paste(
      stats::na.omit(c(
        if (inherits(full_fit, "error")) paste0("full: ", conditionMessage(full_fit)),
        if (inherits(no_temp_fit, "error")) paste0("no_temp: ", conditionMessage(no_temp_fit)),
        if (inherits(no_humidity_fit, "error")) paste0("no_humidity: ", conditionMessage(no_humidity_fit)),
        if (inherits(linear_temp_fit, "error")) paste0("linear_temp: ", conditionMessage(linear_temp_fit))
      )),
      collapse = "; "
    )
    empty_summary$skip_reason <- reason
    empty_summary$omitted_terms <- paste(term_info$omitted, collapse = "; ")
    return(list(summary = empty_summary, curve = tibble::tibble(), coefficients = tibble::tibble(), vcov = tibble::tibble()))
  }

  lrt_p <- function(reduced, full) {
    ll_reduced <- stats::logLik(reduced)
    ll_full <- stats::logLik(full)
    chi <- 2 * (as.numeric(ll_full) - as.numeric(ll_reduced))
    df_diff <- attr(ll_full, "df") - attr(ll_reduced, "df")
    if (!is.finite(chi) || !is.finite(df_diff) || df_diff <= 0) return(NA_real_)
    stats::pchisq(chi, df = df_diff, lower.tail = FALSE)
  }

  temp_grid <- seq(
    stats::quantile(base$tmax_mean_c, 0.01, na.rm = TRUE, names = FALSE),
    stats::quantile(base$tmax_mean_c, 0.99, na.rm = TRUE, names = FALSE),
    length.out = 100
  )
  ref_temp <- stats::median(base$tmax_mean_c, na.rm = TRUE)
  ref_humidity <- stats::median(base$rmax_mean_pct, na.rm = TRUE)
  newdata <- make_reference_newdata(base, temp_grid, term_info$keep) |>
    mutate(rmax_mean_pct = ref_humidity)
  pred <- as.data.frame(predict(full_fit, newdata = newdata, type = "probs"), stringsAsFactors = FALSE)
  if (is.null(dim(pred))) pred <- as.data.frame(t(pred), stringsAsFactors = FALSE)
  pred$tmax_mean_c <- temp_grid
  pred$rmax_mean_pct <- ref_humidity
  curve <- pred |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(keep_levels),
      names_to = "phenotype",
      values_to = "predicted_probability"
    ) |>
    mutate(
      model = label,
      reference_temp_c = ref_temp,
      reference_humidity_pct = ref_humidity,
      spline_df = spline_df,
      .before = 1
    )

  summary <- tibble::tibble(
    model = label,
    n = nrow(stats::model.frame(full_fit)),
    spline_df = spline_df,
    reference_level = reference_level,
    reference_temp_c = ref_temp,
    reference_humidity_pct = ref_humidity,
    min_temp_c = min(base$tmax_mean_c, na.rm = TRUE),
    max_temp_c = max(base$tmax_mean_c, na.rm = TRUE),
    min_humidity_pct = min(base$rmax_mean_pct, na.rm = TRUE),
    max_humidity_pct = max(base$rmax_mean_pct, na.rm = TRUE),
    temperature_lrt_p_value = lrt_p(no_temp_fit, full_fit),
    humidity_lrt_p_value = lrt_p(no_humidity_fit, full_fit),
    temperature_nonlinear_lrt_p_value = lrt_p(linear_temp_fit, full_fit),
    aic = stats::AIC(full_fit),
    covariates = paste(c(temp_term, humidity_term, term_info$keep), collapse = " + "),
    modeled_phenotypes = paste(keep_levels, collapse = " | "),
    excluded_phenotypes = paste(setdiff(levels(factor(df$phenotype)), keep_levels), collapse = " | "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )

  coef_mat <- stats::coef(full_fit)
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, nrow = 1, dimnames = list(names(coef_mat)[[1]], names(coef_mat)))
  }
  coefficients <- as.data.frame(coef_mat, stringsAsFactors = FALSE) |>
    tibble::rownames_to_column("outcome_level") |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of("outcome_level"),
      names_to = "coefficient",
      values_to = "estimate"
    ) |>
    mutate(
      model = label,
      reference_level = reference_level,
      n = nrow(stats::model.frame(full_fit)),
      spline_df = spline_df,
      covariates = paste(c(temp_term, humidity_term, term_info$keep), collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      .before = 1
    )

  vc <- tryCatch(stats::vcov(full_fit), error = function(e) NULL)
  vcov_tbl <- tibble::tibble()
  if (!is.null(vc) && length(vc) > 0) {
    parse_multinom_vcov_name <- function(x) {
      pieces <- strsplit(x, ":", fixed = TRUE)[[1]]
      if (length(pieces) < 2) {
        tibble::tibble(outcome_level = NA_character_, coefficient = x)
      } else {
        tibble::tibble(
          outcome_level = pieces[[1]],
          coefficient = paste(pieces[-1], collapse = ":")
        )
      }
    }
    vc_names <- colnames(vc)
    if (is.null(vc_names)) {
      vc_names <- as.vector(outer(rownames(coef_mat), colnames(coef_mat), paste, sep = ":"))
      colnames(vc) <- vc_names
      rownames(vc) <- vc_names
    }
    row_info <- purrr::map_dfr(rownames(vc), parse_multinom_vcov_name) |>
      rename(outcome_level_row = "outcome_level", coefficient_row = "coefficient") |>
      mutate(row_name = rownames(vc), .before = 1)
    col_info <- purrr::map_dfr(colnames(vc), parse_multinom_vcov_name) |>
      rename(outcome_level_col = "outcome_level", coefficient_col = "coefficient") |>
      mutate(col_name = colnames(vc), .before = 1)
    vcov_tbl <- as.data.frame(as.table(vc), stringsAsFactors = FALSE) |>
      transmute(
        row_name = as.character(.data$Var1),
        col_name = as.character(.data$Var2),
        covariance = as.numeric(.data$Freq)
      ) |>
      left_join(row_info, by = "row_name") |>
      left_join(col_info, by = "col_name") |>
      transmute(
        model = label,
        reference_level = reference_level,
        n = nrow(stats::model.frame(full_fit)),
        spline_df = spline_df,
        outcome_level_row,
        coefficient_row,
        outcome_level_col,
        coefficient_col,
        covariance,
        covariates = paste(c(temp_term, humidity_term, term_info$keep), collapse = " + "),
        omitted_terms = paste(term_info$omitted, collapse = "; ")
      )
  }

  list(summary = summary, curve = curve, coefficients = coefficients, vcov = vcov_tbl)
}

ohca <- readr::read_csv(OHCA_PATH, show_col_types = FALSE) |>
  mutate(
    hospitalization_id = as.character(.data$hospitalization_id),
    patient_id = as.character(.data$patient_id),
    county_fips = normalize_county_fips(.data$county_fips),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    first_icu_in = as_utc_datetime(.data$first_icu_in),
    last_icu_out = as_utc_datetime(.data$last_icu_out),
    admission_date = as.Date(.data$admission_date),
    admit_year = lubridate::year(.data$admission_date),
    age_group = ifelse(as.numeric(.data$age_at_admission) >= 65, ">=65", "<65"),
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
    rmax_per10pct = .data$rmax_mean_pct / 10,
    death_within_72h = !is.na(.data$discharge_dttm) &
      .data$hospital_death == 1L &
      as.numeric(difftime(.data$discharge_dttm, .data$first_icu_in, units = "hours")) <= WINDOW_HOURS
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

respiratory <- read_clif_table(tables_path, file_type, "respiratory_support", required = FALSE)
if (is.null(respiratory) || nrow(respiratory) == 0) {
  resp_summary <- tibble::tibble(hospitalization_id = cohort_ids)
} else if (!all(c("hospitalization_id", "recorded_dttm") %in% names(respiratory))) {
  resp_summary <- tibble::tibble(hospitalization_id = cohort_ids)
} else {
  respiratory <- ensure_columns(respiratory, c("device_category", "artificial_airway", "tracheostomy"))
  resp_summary <- respiratory |>
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
    filter(.data$icu_hour >= 0, .data$icu_hour <= WINDOW_HOURS) |>
    arrange(.data$hospitalization_id, .data$recorded_dttm) |>
    group_by(.data$hospitalization_id) |>
    summarise(
      any_imv_0_72h = any(.data$imv, na.rm = TRUE),
      any_imv_48_72h = any(.data$imv & .data$icu_hour >= EXTUBATION_ASSESSMENT_START_HOUR, na.rm = TRUE),
      any_non_imv_after_imv = {
        imv_times <- .data$icu_hour[.data$imv]
        if (length(imv_times) == 0L) FALSE else any(!.data$imv & .data$icu_hour > max(imv_times, na.rm = TRUE), na.rm = TRUE)
      },
      last_resp_hour_0_72h = max(.data$icu_hour, na.rm = TRUE),
      last_resp_device_0_72h = safe_last(.data$device_category),
      last_resp_imv_0_72h = safe_last(.data$imv),
      .groups = "drop"
    )
}

assessments <- read_clif_table(tables_path, file_type, "patient_assessments", required = FALSE)
if (is.null(assessments) || nrow(assessments) == 0) {
  neuro_summary <- tibble::tibble(
    hospitalization_id = cohort_ids,
    n_neuro_assessments_0_72h = 0L,
    awake_gcs_total_ge13_24_72h = FALSE,
    awake_gcs_motor_6_24_72h = FALSE,
    awake_rass_ge_minus1_24_72h = FALSE,
    any_avpu_alert_0_72h = FALSE,
    any_avpu_unresponsive_0_72h = FALSE,
    sat_pass_0_72h = FALSE,
    sbt_pass_0_72h = FALSE
  )
  gcs_landmark_patient <- tibble::tibble(hospitalization_id = character(), gcs_landmark_hour = integer(), gcs_total_landmark = numeric())
} else if (!all(c("hospitalization_id", "recorded_dttm", "assessment_category") %in% names(assessments))) {
  neuro_summary <- tibble::tibble(
    hospitalization_id = cohort_ids,
    n_neuro_assessments_0_72h = 0L,
    awake_gcs_total_ge13_24_72h = FALSE,
    awake_gcs_motor_6_24_72h = FALSE,
    awake_rass_ge_minus1_24_72h = FALSE,
    any_avpu_alert_0_72h = FALSE,
    any_avpu_unresponsive_0_72h = FALSE,
    sat_pass_0_72h = FALSE,
    sbt_pass_0_72h = FALSE
  )
  gcs_landmark_patient <- tibble::tibble(hospitalization_id = character(), gcs_landmark_hour = integer(), gcs_total_landmark = numeric())
} else {
  assessments <- ensure_columns(assessments, c("numerical_value", "categorical_value", "text_value"))
  assessment_long <- assessments |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      assessment_category = as.character(.data$assessment_category),
      numerical_value = suppressWarnings(as.numeric(.data$numerical_value)),
      categorical_value = clean_binary_text(.data$categorical_value),
      text_value = clean_binary_text(.data$text_value)
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    inner_join(cohort |> select("hospitalization_id", "first_icu_in"), by = "hospitalization_id") |>
    mutate(
      assessment_category_clean = stringr::str_to_lower(.data$assessment_category),
      icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours"))
    ) |>
    filter(.data$icu_hour >= 0, .data$icu_hour <= WINDOW_HOURS)

  neuro_summary <- assessment_long |>
    group_by(.data$hospitalization_id) |>
    summarise(
      n_neuro_assessments_0_72h = sum(.data$assessment_category_clean %in% c("gcs_total", "gcs_motor", "rass", "avpu"), na.rm = TRUE),
      min_gcs_0_72h = suppressWarnings(min(.data$numerical_value[.data$assessment_category_clean == "gcs_total" & .data$numerical_value >= 3 & .data$numerical_value <= 15], na.rm = TRUE)),
      best_gcs_24_72h = suppressWarnings(max(.data$numerical_value[.data$assessment_category_clean == "gcs_total" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 3 & .data$numerical_value <= 15], na.rm = TRUE)),
      last_gcs_24_72h = safe_last(.data$numerical_value[.data$assessment_category_clean == "gcs_total" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 3 & .data$numerical_value <= 15]),
      min_gcs_motor_0_72h = suppressWarnings(min(.data$numerical_value[.data$assessment_category_clean == "gcs_motor" & .data$numerical_value >= 1 & .data$numerical_value <= 6], na.rm = TRUE)),
      best_gcs_motor_24_72h = suppressWarnings(max(.data$numerical_value[.data$assessment_category_clean == "gcs_motor" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 1 & .data$numerical_value <= 6], na.rm = TRUE)),
      last_gcs_motor_24_72h = safe_last(.data$numerical_value[.data$assessment_category_clean == "gcs_motor" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 1 & .data$numerical_value <= 6]),
      min_rass_0_72h = suppressWarnings(min(.data$numerical_value[.data$assessment_category_clean == "rass" & .data$numerical_value >= -5 & .data$numerical_value <= 4], na.rm = TRUE)),
      best_rass_24_72h = suppressWarnings(max(.data$numerical_value[.data$assessment_category_clean == "rass" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= -5 & .data$numerical_value <= 4], na.rm = TRUE)),
      last_rass_24_72h = safe_last(.data$numerical_value[.data$assessment_category_clean == "rass" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= -5 & .data$numerical_value <= 4]),
      awake_gcs_total_ge13_24_72h = any(.data$assessment_category_clean == "gcs_total" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 13 & .data$numerical_value <= 15, na.rm = TRUE),
      awake_gcs_motor_6_24_72h = any(.data$assessment_category_clean == "gcs_motor" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 6 & .data$numerical_value <= 6, na.rm = TRUE),
      awake_rass_ge_minus1_24_72h = any(.data$assessment_category_clean == "rass" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= -1 & .data$numerical_value <= 4, na.rm = TRUE),
      any_avpu_alert_0_72h = any(.data$assessment_category_clean == "avpu" & stringr::str_detect(paste(.data$categorical_value, .data$text_value), "alert"), na.rm = TRUE),
      any_avpu_unresponsive_0_72h = any(.data$assessment_category_clean == "avpu" & stringr::str_detect(paste(.data$categorical_value, .data$text_value), "unresponsive|pain"), na.rm = TRUE),
      sat_pass_0_72h = any(.data$assessment_category_clean == "sat_delivery_pass_fail" & .data$categorical_value == "pass", na.rm = TRUE),
      sbt_pass_0_72h = any(.data$assessment_category_clean == "sbt_delivery_pass_fail" & .data$categorical_value == "pass", na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      across(c(min_gcs_0_72h, best_gcs_24_72h, min_gcs_motor_0_72h, best_gcs_motor_24_72h, min_rass_0_72h, best_rass_24_72h), ~ ifelse(is.infinite(.x), NA_real_, .x))
    )

  gcs_landmark_patient <- assessment_long |>
    filter(
      .data$assessment_category_clean == "gcs_total",
      .data$numerical_value >= 3,
      .data$numerical_value <= 15
    ) |>
    mutate(gcs_landmark_hour = dplyr::case_when(
      .data$icu_hour > 12 & .data$icu_hour <= 24 ~ 24L,
      .data$icu_hour > 24 & .data$icu_hour <= 48 ~ 48L,
      .data$icu_hour > 48 & .data$icu_hour <= 72 ~ 72L,
      TRUE ~ NA_integer_
    )) |>
    filter(!is.na(.data$gcs_landmark_hour)) |>
    group_by(.data$hospitalization_id, .data$gcs_landmark_hour) |>
    summarise(gcs_total_landmark = median(.data$numerical_value, na.rm = TRUE), .groups = "drop")
}

medication_continuous <- read_clif_table(
  tables_path,
  file_type,
  "medication_admin_continuous",
  columns = c("hospitalization_id", "admin_dttm", "med_category", "med_group", "mar_action_group", "med_dose", "med_dose_unit"),
  required = FALSE
)
if (is.null(medication_continuous) || nrow(medication_continuous) == 0) {
  vaso_summary <- tibble::tibble(hospitalization_id = cohort_ids, vaso_any_0_72h = 0L)
} else {
  vaso_summary <- medication_continuous |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      admin_dttm = as_utc_datetime(.data$admin_dttm),
      med_category = stringr::str_to_lower(as.character(.data$med_category)),
      med_group = stringr::str_to_lower(as.character(.data$med_group)),
      mar_action_group = stringr::str_to_lower(as.character(.data$mar_action_group))
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$admin_dttm)) |>
    inner_join(cohort |> select("hospitalization_id", "first_icu_in"), by = "hospitalization_id") |>
    mutate(
      icu_hour = as.numeric(difftime(.data$admin_dttm, .data$first_icu_in, units = "hours")),
      vasoactive = .data$med_group == "vasoactives" |
        .data$med_category %in% c("norepinephrine", "epinephrine", "phenylephrine", "vasopressin", "dopamine", "angiotensin", "dobutamine", "milrinone")
    ) |>
    filter(.data$icu_hour >= 0, .data$icu_hour <= WINDOW_HOURS, .data$mar_action_group == "administered", .data$vasoactive) |>
    distinct(.data$hospitalization_id) |>
    mutate(vaso_any_0_72h = 1L) |>
    right_join(tibble::tibble(hospitalization_id = cohort_ids), by = "hospitalization_id") |>
    mutate(vaso_any_0_72h = tidyr::replace_na(.data$vaso_any_0_72h, 0L))
}

vitals_for_sofa <- read_clif_table(tables_path, file_type, "vitals", columns = c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value"), required = FALSE)
labs_for_sofa <- read_clif_table(tables_path, file_type, "labs", columns = c("hospitalization_id", "lab_result_dttm", "lab_category", "lab_value_numeric"), required = FALSE)
sofa_long <- calculate_sofa_windows(
  cohort_df = cohort,
  vitals_df = vitals_for_sofa,
  labs_df = labs_for_sofa,
  support_df = respiratory,
  med_admin_df = medication_continuous,
  scores_df = assessments,
  windows = c(24, 48, 72)
)
sofa_wide <- sofa_long |>
  select("hospitalization_id", "sofa_window_hours", "sofa_total") |>
  tidyr::pivot_wider(names_from = "sofa_window_hours", values_from = "sofa_total", names_prefix = "sofa_total_") |>
  rename(
    sofa_total_24h = "sofa_total_24",
    sofa_total_48h = "sofa_total_48",
    sofa_total_72h = "sofa_total_72"
  )

phenotype_cohort <- cohort |>
  left_join(ohca_mechanism, by = "hospitalization_id") |>
  left_join(neuro_dx, by = "hospitalization_id") |>
  left_join(resp_summary, by = "hospitalization_id") |>
  left_join(neuro_summary, by = "hospitalization_id") |>
  left_join(vaso_summary, by = "hospitalization_id") |>
  left_join(sofa_wide, by = "hospitalization_id") |>
  mutate(
    ohca_mechanism = tidyr::replace_na(.data$ohca_mechanism, "unclear_other"),
    anoxic_brain_injury_dx = yesno(.data$anoxic_brain_injury_dx),
    brain_death_dx = yesno(.data$brain_death_dx),
    any_imv_0_72h = yesno(.data$any_imv_0_72h),
    any_imv_48_72h = yesno(.data$any_imv_48_72h),
    any_non_imv_after_imv = yesno(.data$any_non_imv_after_imv),
    last_resp_imv_0_72h = yesno(.data$last_resp_imv_0_72h),
    vaso_any_0_72h = tidyr::replace_na(.data$vaso_any_0_72h, 0L),
    n_neuro_assessments_0_72h = tidyr::replace_na(.data$n_neuro_assessments_0_72h, 0L),
    awake_gcs_total_ge13_24_72h = yesno(.data$awake_gcs_total_ge13_24_72h),
    awake_gcs_motor_6_24_72h = yesno(.data$awake_gcs_motor_6_24_72h),
    awake_rass_ge_minus1_24_72h = yesno(.data$awake_rass_ge_minus1_24_72h),
    any_avpu_alert_0_72h = yesno(.data$any_avpu_alert_0_72h),
    sat_pass_0_72h = yesno(.data$sat_pass_0_72h),
    sbt_pass_0_72h = yesno(.data$sbt_pass_0_72h),
    severe_neuro_signal = (
      (!is.na(.data$min_gcs_0_72h) & .data$min_gcs_0_72h <= 5) |
        (!is.na(.data$min_gcs_motor_0_72h) & .data$min_gcs_motor_0_72h <= 2) |
        (!is.na(.data$min_rass_0_72h) & .data$min_rass_0_72h <= -4)
    ),
    awake_signal = (
      (!is.na(.data$best_gcs_24_72h) & .data$best_gcs_24_72h >= 13) |
        (!is.na(.data$last_gcs_24_72h) & .data$last_gcs_24_72h >= 13) |
        (!is.na(.data$best_gcs_motor_24_72h) & .data$best_gcs_motor_24_72h >= 6) |
        (!is.na(.data$best_rass_24_72h) & .data$best_rass_24_72h >= -1) |
        yesno(.data$any_avpu_alert_0_72h) |
        yesno(.data$sat_pass_0_72h) |
        yesno(.data$sbt_pass_0_72h)
    ),
    impaired_neuro_signal = (
      (!is.na(.data$last_gcs_24_72h) & .data$last_gcs_24_72h < 13) |
        (!is.na(.data$min_gcs_0_72h) & .data$min_gcs_0_72h <= 8) |
        (!is.na(.data$last_gcs_motor_24_72h) & .data$last_gcs_motor_24_72h < 6) |
        (!is.na(.data$last_rass_24_72h) & .data$last_rass_24_72h <= -2) |
        yesno(.data$any_avpu_unresponsive_0_72h)
    ),
    extubated_by_72h = .data$any_imv_0_72h & (
      .data$any_non_imv_after_imv |
        (!.data$any_imv_48_72h & !.data$last_resp_imv_0_72h)
    ),
    phenotype = dplyr::case_when(
      .data$death_within_72h ~ "anoxic_brain_injury",
      !.data$any_imv_0_72h ~ "alive_no_imv",
      .data$extubated_by_72h & .data$awake_signal ~ "regained_consciousness_extubated",
      .data$severe_neuro_signal | .data$impaired_neuro_signal | .data$any_imv_48_72h | (.data$any_imv_0_72h & !.data$extubated_by_72h) ~ "limited_brain_function",
      TRUE ~ "unclassified"
    ),
    phenotype = factor(.data$phenotype, levels = c(
      "alive_no_imv",
      "regained_consciousness_extubated",
      "limited_brain_function",
      "anoxic_brain_injury",
      "unclassified"
    )),
    unclassified_reason = dplyr::case_when(
      .data$phenotype != "unclassified" ~ NA_character_,
      !.data$any_imv_0_72h & .data$n_neuro_assessments_0_72h == 0 ~ "Alive at 72h; no IMV in first 72h and no neurologic assessment evidence",
      !.data$any_imv_0_72h ~ "Alive at 72h; no IMV in first 72h and no qualifying impaired/awake neurologic phenotype signal",
      .data$extubated_by_72h & !.data$awake_signal ~ "Extubated/no ongoing IMV by 72h but no awake signal",
      .data$awake_signal & !.data$extubated_by_72h & !.data$impaired_neuro_signal & !.data$any_imv_48_72h ~ "Awake signal without extubation and without limited-brain-function criteria",
      TRUE ~ "Alive at 72h; structured criteria did not uniquely identify a primary phenotype"
    ),
    anoxic_outcome = as.integer(.data$phenotype == "anoxic_brain_injury"),
    limited_vs_regained = dplyr::case_when(
      .data$phenotype == "limited_brain_function" ~ 1L,
      .data$phenotype == "regained_consciousness_extubated" ~ 0L,
      TRUE ~ NA_integer_
    ),
    extubated_awake_outcome = as.integer(.data$phenotype == "regained_consciousness_extubated")
  )

phenotype_summary <- phenotype_cohort |>
  count(.data$phenotype, name = "n") |>
  mutate(
    site_name = site_name,
    pct = 100 * .data$n / sum(.data$n),
    .before = 1
  )

mechanism_summary <- phenotype_cohort |>
  count(.data$ohca_mechanism, .data$phenotype, name = "n") |>
  group_by(.data$ohca_mechanism) |>
  mutate(pct_within_mechanism = 100 * .data$n / sum(.data$n)) |>
  ungroup() |>
  group_by(.data$phenotype) |>
  mutate(pct_within_phenotype = 100 * .data$n / sum(.data$n)) |>
  ungroup() |>
  mutate(site_name = site_name, .before = 1)

evidence_summary <- phenotype_cohort |>
  group_by(.data$phenotype) |>
  summarise(
    n = n(),
    anoxic_dx_pct = safe_pct(sum(.data$anoxic_brain_injury_dx, na.rm = TRUE), n()),
    brain_death_dx_pct = safe_pct(sum(.data$brain_death_dx, na.rm = TRUE), n()),
    severe_neuro_signal_pct = safe_pct(sum(.data$severe_neuro_signal, na.rm = TRUE), n()),
    awake_signal_pct = safe_pct(sum(.data$awake_signal, na.rm = TRUE), n()),
    impaired_neuro_signal_pct = safe_pct(sum(.data$impaired_neuro_signal, na.rm = TRUE), n()),
    extubated_by_72h_pct = safe_pct(sum(.data$extubated_by_72h, na.rm = TRUE), n()),
    imv_48_72h_pct = safe_pct(sum(.data$any_imv_48_72h, na.rm = TRUE), n()),
    hospital_death_pct = safe_pct(sum(.data$hospital_death == 1, na.rm = TRUE), n()),
    median_tmax_c = median(.data$tmax_mean_c, na.rm = TRUE),
    median_rmax_pct = median(.data$rmax_mean_pct, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(site_name = site_name, .before = 1)

base_adjust <- c("age_at_admission", "sex_group", "race_group", "admit_year")
mechanism_adjust <- c(base_adjust, "ohca_mechanism")

phenotype_assignment_primary <- run_phenotype_assignment_model(
  phenotype_cohort,
  "phenotype_assignment_temperature_humidity_demographics",
  base_adjust
)
phenotype_assignment_model <- phenotype_assignment_primary$summary |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_curve <- phenotype_assignment_primary$curve |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_coefficients <- phenotype_assignment_primary$coefficients |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_vcov <- phenotype_assignment_primary$vcov |>
  mutate(site_name = site_name, .before = 1)

phenotype_assignment_mechanism <- run_phenotype_assignment_model(
  phenotype_cohort,
  "phenotype_assignment_temperature_humidity_demographics_mechanism",
  mechanism_adjust
)
phenotype_assignment_model_mechanism <- phenotype_assignment_mechanism$summary |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_curve_mechanism <- phenotype_assignment_mechanism$curve |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_coefficients_mechanism <- phenotype_assignment_mechanism$coefficients |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_vcov_mechanism <- phenotype_assignment_mechanism$vcov |>
  mutate(site_name = site_name, .before = 1)

phenotype_definitions <- tibble::tibble(
  phenotype = c(
    "alive_no_imv",
    "anoxic_brain_injury",
    "regained_consciousness_extubated",
    "limited_brain_function",
    "unclassified"
  ),
  hierarchy_order = c(2L, 1L, 3L, 4L, 5L),
  definition = c(
    "Alive at 72 hours and no invasive mechanical ventilation evidence during the first 72 ICU hours.",
    "Death within 72 hours of ICU entry. No diagnosis-code criterion is applied.",
    "Not anoxic phenotype; invasive ventilation occurred in first 72 ICU hours; evidence of extubation/no ongoing IMV by 72h; and awake signal from GCS >=13, GCS motor 6, RASS >=-1, AVPU alert, SAT pass, or SBT pass.",
    "Not anoxic or extubated-awake phenotype; severe or impaired neurologic signal in first 72 ICU hours, ongoing IMV at 48-72h, or IMV in first 72h without evidence of extubation.",
    "Does not meet structured-data criteria for the four primary phenotypes."
  )
) |>
  mutate(site_name = site_name, .before = 1)

primary_phenotype_levels <- c("alive_no_imv", "regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury")
table2_phenotype_levels <- c(primary_phenotype_levels, "unclassified")
phenotype_analysis_cohort <- phenotype_cohort |>
  filter(.data$phenotype %in% primary_phenotype_levels) |>
  mutate(phenotype = factor(as.character(.data$phenotype), levels = primary_phenotype_levels))

phenotype_table_cohort <- phenotype_cohort |>
  mutate(phenotype = factor(as.character(.data$phenotype), levels = table2_phenotype_levels))

gcs_landmark_by_phenotype <- gcs_landmark_patient |>
  left_join(phenotype_table_cohort |> select("hospitalization_id", "phenotype"), by = "hospitalization_id") |>
  filter(!is.na(.data$phenotype)) |>
  group_by(.data$phenotype, .data$gcs_landmark_hour) |>
  summarise(
    n_with_gcs = sum(!is.na(.data$gcs_total_landmark)),
    mean_gcs = mean(.data$gcs_total_landmark, na.rm = TRUE),
    median_gcs = median(.data$gcs_total_landmark, na.rm = TRUE),
    q1_gcs = stats::quantile(.data$gcs_total_landmark, 0.25, na.rm = TRUE, names = FALSE),
    q3_gcs = stats::quantile(.data$gcs_total_landmark, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) |>
  mutate(site_name = site_name, .before = 1)

gcs_landmark_table_rows <- dplyr::bind_rows(lapply(c(24L, 48L, 72L), function(landmark_hour) {
  vals <- vapply(table2_phenotype_levels, function(ph) {
    fmt_median_iqr_n(gcs_landmark_patient$gcs_total_landmark[
      gcs_landmark_patient$hospitalization_id %in% phenotype_table_cohort$hospitalization_id[phenotype_table_cohort$phenotype == ph] &
        gcs_landmark_patient$gcs_landmark_hour == landmark_hour
    ])
  }, character(1))
  tibble::tibble(
    section = "Neurologic status",
    characteristic = "GCS total by ICU landmark, median [IQR]",
    level = sprintf("Hour %d", landmark_hour)
  ) |>
    bind_cols(as_tibble(as.list(vals)))
}))

ohca_summary_raw <- if (file.exists(OHCA_SUMMARY_PATH)) readr::read_csv(OHCA_SUMMARY_PATH, show_col_types = FALSE) else tibble::tibble()
summary_value <- function(col) {
  if (col %in% names(ohca_summary_raw) && nrow(ohca_summary_raw) > 0) suppressWarnings(as.numeric(ohca_summary_raw[[col]][[1]])) else NA_real_
}

consort_counts <- tibble::tibble(
  step_order = seq_len(9),
  step = c(
    "All ICU admissions in CLIF, 2018-2024",
    "OHCA present-on-admission ICU admissions before pathway/timing restriction",
    "OHCA ICU admissions with allowed ED/procedure/direct ICU pathway and ICU entry <24h",
    "OHCA ICU admissions with admission-day county Tmax and humidity",
    "Classified into one of four primary 72h phenotypes",
    "Alive at 72h with no IMV in first 72h",
    "Regained consciousness and extubated",
    "Limited brain function",
    "Anoxic/severe group: death within 72h"
  ),
  n = c(
    summary_value("n_all_icu_admissions"),
    summary_value("n_ohca_admissions_before_pathway_timing_restriction"),
    summary_value("n_ohca_admissions"),
    nrow(phenotype_cohort),
    nrow(phenotype_analysis_cohort),
    sum(phenotype_analysis_cohort$phenotype == "alive_no_imv", na.rm = TRUE),
    sum(phenotype_analysis_cohort$phenotype == "regained_consciousness_extubated", na.rm = TRUE),
    sum(phenotype_analysis_cohort$phenotype == "limited_brain_function", na.rm = TRUE),
    sum(phenotype_analysis_cohort$phenotype == "anoxic_brain_injury", na.rm = TRUE)
  ),
  notes = c(
    "Source: clif_hospitalization/ADT ICU cohort summary",
    "POA cardiac arrest/OHCA-proxy diagnosis plus ICU admission",
    "Cohort definition used for all OHCA ICU heat analyses",
    "Required exposure fields for 72h phenotype and competing-risk models",
    "Excludes unclassified structured-data phenotype",
    "Alive at 72h and no invasive mechanical ventilation evidence in first 72 ICU hours",
    "Not death within 72h; IMV in first 72h; extubated/no ongoing IMV by 72h; awake signal",
    "Not death within 72h or regained/extubated; impaired neurologic signal and/or ongoing IMV evidence",
    "Death within 72h of ICU entry"
  )
) |>
  mutate(
    site_name = site_name,
    n_excluded_since_previous = dplyr::lag(.data$n) - .data$n,
    n_excluded_since_previous = ifelse(.data$step_order == 1L | .data$step_order >= 6L, NA_real_, .data$n_excluded_since_previous),
    .before = 1
  )

table1 <- bind_rows(
  tibble::tibble(section = "Cohort", characteristic = "Four primary phenotype cohort", level = "", value = format(nrow(phenotype_analysis_cohort), big.mark = ",")),
  add_table1_continuous(phenotype_analysis_cohort, "Demographics", "Age, years, median [IQR]", "age_at_admission"),
  add_table1_categorical(phenotype_analysis_cohort, "Demographics", "Sex", "sex_category"),
  add_table1_categorical(phenotype_analysis_cohort, "Demographics", "Race", "race_category"),
  add_table1_categorical(phenotype_analysis_cohort, "Demographics", "Ethnicity", "ethnicity_category"),
  add_table1_continuous(phenotype_analysis_cohort, "Admission exposure", "Admission-day Tmax, C, median [IQR]", "tmax_mean_c"),
  add_table1_continuous(phenotype_analysis_cohort, "Admission exposure", "Admission-day relative humidity, %, median [IQR]", "rmax_mean_pct"),
  add_table1_categorical(phenotype_analysis_cohort, "OHCA mechanism", "POA diagnosis mechanism", "ohca_mechanism"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Any IMV in first 72h", "any_imv_0_72h"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "IMV evidence at ICU hours 48-72", "any_imv_48_72h"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Extubated/no ongoing IMV by 72h", "extubated_by_72h"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Awake signal", "awake_signal"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Impaired neurologic signal", "impaired_neuro_signal"),
  add_table1_binary(phenotype_analysis_cohort, "Outcome", "Hospital death", "hospital_death"),
  add_table1_binary(phenotype_analysis_cohort, "Outcome", "Death within 72h", "death_within_72h")
) |>
  mutate(site_name = site_name, .before = 1)

table2 <- bind_rows(
  tibble::tibble(section = "Cohort", characteristic = "Phenotype cohort size", level = "") |>
    bind_cols(as_tibble(as.list(setNames(vapply(table2_phenotype_levels, function(ph) format(sum(phenotype_table_cohort$phenotype == ph), big.mark = ","), character(1)), table2_phenotype_levels)))),
  add_table2_continuous(phenotype_table_cohort, "Demographics", "Age, years, median [IQR]", "age_at_admission", table2_phenotype_levels),
  add_table2_categorical(phenotype_table_cohort, "Demographics", "Sex", "sex_category", table2_phenotype_levels),
  add_table2_categorical(phenotype_table_cohort, "Demographics", "Race", "race_category", table2_phenotype_levels),
  add_table2_categorical(phenotype_table_cohort, "Demographics", "Ethnicity", "ethnicity_category", table2_phenotype_levels),
  add_table2_continuous(phenotype_table_cohort, "Admission exposure", "Admission-day Tmax, C, median [IQR]", "tmax_mean_c", table2_phenotype_levels),
  add_table2_continuous(phenotype_table_cohort, "Admission exposure", "Admission-day relative humidity, %, median [IQR]", "rmax_mean_pct", table2_phenotype_levels),
  add_table2_categorical(phenotype_table_cohort, "OHCA mechanism", "POA diagnosis mechanism", "ohca_mechanism", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Organ support", "Any IMV in first 72h", "any_imv_0_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Organ support", "Any vasopressor in first 72h", "vaso_any_0_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Organ support", "IMV evidence at ICU hours 48-72", "any_imv_48_72h", table2_phenotype_levels),
  add_table2_continuous(phenotype_table_cohort, "Organ dysfunction", "SOFA total at 24h, median [IQR]", "sofa_total_24h", table2_phenotype_levels),
  add_table2_continuous(phenotype_table_cohort, "Organ dysfunction", "SOFA total at 48h, median [IQR]", "sofa_total_48h", table2_phenotype_levels),
  add_table2_continuous(phenotype_table_cohort, "Organ dysfunction", "SOFA total at 72h, median [IQR]", "sofa_total_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "72h phenotype evidence", "Extubated/no ongoing IMV by 72h", "extubated_by_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "72h phenotype evidence", "Composite awake signal", "awake_signal", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Awake signal components", "GCS total >=13 after ICU hour 24", "awake_gcs_total_ge13_24_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Awake signal components", "GCS motor = 6 after ICU hour 24", "awake_gcs_motor_6_24_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Awake signal components", "RASS >= -1 after ICU hour 24", "awake_rass_ge_minus1_24_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Awake signal components", "AVPU alert in first 72h", "any_avpu_alert_0_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Awake signal components", "SAT pass in first 72h", "sat_pass_0_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Awake signal components", "SBT pass in first 72h", "sbt_pass_0_72h", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "72h phenotype evidence", "Impaired neurologic signal", "impaired_neuro_signal", table2_phenotype_levels),
  gcs_landmark_table_rows,
  add_table2_binary(phenotype_table_cohort, "Outcome", "In-hospital mortality", "hospital_death", table2_phenotype_levels),
  add_table2_binary(phenotype_table_cohort, "Outcome", "Death within 72h", "death_within_72h", table2_phenotype_levels),
  add_table2_categorical(phenotype_table_cohort |> filter(.data$phenotype == "unclassified"), "Unclassified audit", "Reason unclassified", "unclassified_reason", table2_phenotype_levels)
) |>
  mutate(site_name = site_name, .before = 1)

readr::write_csv(
  phenotype_cohort |> select(
    "hospitalization_id", "patient_id", "admission_date", "tmax_mean_c", "rmax_mean_pct", "phenotype",
    "unclassified_reason",
    "ohca_mechanism", "ohca_mechanism_detail",
    "anoxic_brain_injury_dx", "brain_death_dx", "severe_neuro_signal", "awake_signal",
    "awake_gcs_total_ge13_24_72h", "awake_gcs_motor_6_24_72h", "awake_rass_ge_minus1_24_72h",
    "any_avpu_alert_0_72h", "sat_pass_0_72h", "sbt_pass_0_72h",
    "impaired_neuro_signal", "extubated_by_72h", "any_imv_0_72h", "any_imv_48_72h",
    "vaso_any_0_72h", "sofa_total_24h", "sofa_total_48h", "sofa_total_72h",
    "min_gcs_0_72h", "best_gcs_24_72h", "last_gcs_24_72h", "best_rass_24_72h",
    "last_rass_24_72h", "hospital_death", "death_within_72h"
  ),
  file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_patient_audit.csv")
)
readr::write_csv(phenotype_summary, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_summary.csv"))
readr::write_csv(consort_counts, file.path(OUTPUT_DIR, "ohca_icu_72h_consort_flow.csv"))
readr::write_csv(table1, file.path(OUTPUT_DIR, "ohca_icu_72h_table1.csv"))
readr::write_csv(table2, file.path(OUTPUT_DIR, "ohca_icu_72h_table2_by_phenotype.csv"))
readr::write_csv(gcs_landmark_by_phenotype, file.path(OUTPUT_DIR, "ohca_icu_72h_gcs_landmark_by_phenotype.csv"))
readr::write_csv(mechanism_summary, file.path(OUTPUT_DIR, "ohca_icu_72h_ohca_mechanism_summary.csv"))
readr::write_csv(evidence_summary, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_evidence_summary.csv"))
readr::write_csv(phenotype_assignment_model, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_model.csv"))
readr::write_csv(phenotype_assignment_curve, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_temperature_curve.csv"))
readr::write_csv(phenotype_assignment_coefficients, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_coefficients.csv"))
readr::write_csv(phenotype_assignment_vcov, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_vcov.csv"))
readr::write_csv(phenotype_assignment_model_mechanism, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted.csv"))
readr::write_csv(phenotype_assignment_curve_mechanism, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted.csv"))
readr::write_csv(phenotype_assignment_coefficients_mechanism, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted.csv"))
readr::write_csv(phenotype_assignment_vcov_mechanism, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_vcov_mechanism_adjusted.csv"))
readr::write_csv(phenotype_definitions, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_definitions.csv"))

plot_df <- phenotype_summary |>
  filter(.data$phenotype != "unclassified") |>
  mutate(phenotype = factor(.data$phenotype, levels = rev(c("alive_no_imv", "regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury"))))

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = .data$n, y = .data$phenotype, fill = .data$phenotype)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%s (%.1f%%)", .data$n, .data$pct)), hjust = -0.05, size = 3.5) +
    scale_fill_manual(values = c(
      "alive_no_imv" = "#457B9D",
      "regained_consciousness_extubated" = "#2A9D8F",
      "limited_brain_function" = "#E9C46A",
      "anoxic_brain_injury" = "#C44536"
    ), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "OHCA ICU 72-hour Structured Phenotypes",
      x = "Patients",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIGURE_DIR, "figure_ohca_icu_72h_phenotype_counts.png"), p, width = 8.5, height = 4.8, dpi = 300)
}

if (nrow(phenotype_assignment_curve) > 0 && isTRUE(phenotype_assignment_model$estimable[[1]])) {
  curve_plot <- phenotype_assignment_curve |>
    mutate(
      phenotype = factor(
        .data$phenotype,
        levels = c("alive_no_imv", "regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury"),
        labels = c("Alive/no IMV", "Regained/extubated", "Limited brain function", "Death within 72h")
      )
    )

  p_curve <- ggplot(curve_plot, aes(x = .data$tmax_mean_c, y = .data$predicted_probability, color = .data$phenotype)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(
      "Alive/no IMV" = "#457B9D",
      "Regained/extubated" = "#2A9D8F",
      "Limited brain function" = "#E9C46A",
      "Death within 72h" = "#C44536"
    )) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(
      title = "Adjusted OHCA ICU Phenotype Probability by Admission Temperature",
      x = "Admission-day assigned-county Tmax, C",
      y = "Predicted phenotype probability",
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  ggsave(file.path(FIGURE_DIR, "figure_ohca_icu_72h_phenotype_temperature_curves.png"), p_curve, width = 8.5, height = 5.2, dpi = 300)
}

print(phenotype_summary)
print(phenotype_assignment_model)
print(phenotype_assignment_model_mechanism)
message("Wrote OHCA ICU 72-hour phenotype outputs to ", OUTPUT_DIR)

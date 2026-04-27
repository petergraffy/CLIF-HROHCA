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
  library(dplyr)
  library(jsonlite)
  library(lubridate)
  library(readr)
  library(stringr)
  library(tidyr)
})

START_DATE <- as_utc_datetime("2018-01-01")
END_DATE <- as_utc_datetime("2025-01-01")
ADULT_AGE_YEARS <- 18
OHCA_ICD10_PREFIXES <- c("I46")
OHCA_ICD9_PREFIXES <- c("4275")

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

derive_poa_flag <- function(df, diagnosis_source) {
  if ("poa_present" %in% names(df)) {
    raw <- df$poa_present
    numeric_flag <- suppressWarnings(as.integer(as.character(raw)))
    char_flag <- stringr::str_to_upper(as.character(raw))
    return(ifelse(!is.na(numeric_flag), numeric_flag == 1L, char_flag %in% c("Y", "YES", "TRUE", "T", "1")))
  }
  if (identical(diagnosis_source, "admission_diagnosis")) return(rep(TRUE, nrow(df)))
  rep(FALSE, nrow(df))
}

config <- load_project_config(repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)
site_name <- config$site_name %||% Sys.getenv("CLIF_SITE_NAME", unset = "unknown_site")

output_dir <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

patient <- read_clif_table(tables_path, file_type, "patient")
hospitalization <- read_clif_table(tables_path, file_type, "hospitalization")
adt <- read_clif_table(tables_path, file_type, "adt")

diagnosis_source <- NULL
for (candidate in c("hospital_diagnosis", "admission_diagnosis")) {
  if (file.exists(file.path(tables_path, sprintf("clif_%s.%s", candidate, file_type)))) {
    diagnosis_source <- candidate
    break
  }
}
if (is.null(diagnosis_source)) stop("Could not find clif_hospital_diagnosis or clif_admission_diagnosis.")
diagnosis <- read_clif_table(tables_path, file_type, diagnosis_source)

patient_min <- patient |>
  mutate(patient_id = as.character(.data$patient_id))
for (col_name in c("sex_category", "race_category", "ethnicity_category")) {
  if (!col_name %in% names(patient_min)) patient_min[[col_name]] <- NA_character_
}
patient_min <- patient_min |>
  transmute(
    patient_id = .data$patient_id,
    sex_category = as.character(.data$sex_category),
    race_category = as.character(.data$race_category),
    ethnicity_category = as.character(.data$ethnicity_category)
  )

hosp <- hospitalization |>
  transmute(
    patient_id = as.character(.data$patient_id),
    hospitalization_id = as.character(.data$hospitalization_id),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    age_at_admission = suppressWarnings(as.numeric(.data$age_at_admission)),
    discharge_category = if ("discharge_category" %in% names(hospitalization)) as.character(.data$discharge_category) else NA_character_,
    county_code = if ("county_code" %in% names(hospitalization)) as.character(.data$county_code) else NA_character_,
    state_code = if ("state_code" %in% names(hospitalization)) as.character(.data$state_code) else NA_character_
  ) |>
  mutate(
    county_fips = normalize_county_fips(.data$county_code),
    admission_date = as.Date(.data$admission_dttm)
  ) |>
  filter(
    !is.na(.data$admission_dttm),
    .data$admission_dttm >= START_DATE,
    .data$admission_dttm < END_DATE,
    .data$age_at_admission >= ADULT_AGE_YEARS
  )

adt_min <- adt |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    in_dttm = as_utc_datetime(.data$in_dttm),
    out_dttm = as_utc_datetime(.data$out_dttm),
    location_category = if ("location_category" %in% names(adt)) as.character(.data$location_category) else NA_character_,
    hospital_id = if ("hospital_id" %in% names(adt)) as.character(.data$hospital_id) else NA_character_
  )

icu_bounds <- adt_min |>
  filter(stringr::str_detect(stringr::str_to_lower(.data$location_category), "icu"), !is.na(.data$in_dttm)) |>
  group_by(.data$hospitalization_id) |>
  summarise(
    first_icu_in = min(.data$in_dttm, na.rm = TRUE),
    last_icu_out = if (all(is.na(.data$out_dttm))) max(.data$in_dttm, na.rm = TRUE) else max(.data$out_dttm, na.rm = TRUE),
    n_icu_segments = n(),
    .groups = "drop"
  ) |>
  mutate(icu_los_hours = as.numeric(difftime(.data$last_icu_out, .data$first_icu_in, units = "hours")))

hospital_lookup <- build_hospital_id_lookup(adt_min)

all_icu <- hosp |>
  inner_join(icu_bounds, by = "hospitalization_id") |>
  left_join(hospital_lookup, by = "hospitalization_id") |>
  left_join(patient_min, by = "patient_id") |>
  arrange(.data$admission_dttm, .data$hospitalization_id)

dx <- diagnosis |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    diagnosis_code = if ("diagnosis_code" %in% names(diagnosis)) as.character(.data$diagnosis_code) else NA_character_,
    diagnosis_code_format = if ("diagnosis_code_format" %in% names(diagnosis)) as.character(.data$diagnosis_code_format) else NA_character_
  ) |>
  mutate(
    diagnosis_code_clean = norm_code(.data$diagnosis_code),
    diagnosis_code_format = infer_diagnosis_format(.data$diagnosis_code_clean, .data$diagnosis_code_format),
    poa_flag = derive_poa_flag(diagnosis, diagnosis_source)
  )

dx_ohca <- dx |>
  filter(
    .data$poa_flag,
    (
      stringr::str_detect(.data$diagnosis_code_format, "10") &
        stringr::str_starts(.data$diagnosis_code_clean, OHCA_ICD10_PREFIXES)
    ) |
      (
        stringr::str_detect(.data$diagnosis_code_format, "9") &
          stringr::str_starts(.data$diagnosis_code_clean, OHCA_ICD9_PREFIXES)
      )
  ) |>
  group_by(.data$hospitalization_id) |>
  summarise(
    ohca_poa = 1L,
    ohca_codes = paste(sort(unique(.data$diagnosis_code_clean)), collapse = " | "),
    .groups = "drop"
  )

ohca <- all_icu |>
  inner_join(dx_ohca, by = "hospitalization_id") |>
  arrange(.data$admission_dttm, .data$hospitalization_id)

ohca <- apply_site_county_assignment(ohca, repo_root, config)

respiratory <- read_clif_table(tables_path, file_type, "respiratory_support", required = FALSE)
ohca$imv_any <- 0L
ohca$imv_duration_hours <- NA_real_
ohca$imv_record_count <- 0L
if (!is.null(respiratory) && nrow(respiratory) > 0) {
  resp <- respiratory |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      device_category = if ("device_category" %in% names(respiratory)) as.character(.data$device_category) else NA_character_
    ) |>
    semi_join(ohca |> select("hospitalization_id"), by = "hospitalization_id") |>
    mutate(
      imv_row = stringr::str_to_upper(tidyr::replace_na(.data$device_category, "")) %in% c("IMV", "VENT")
    )
  imv_summary <- resp |>
    filter(.data$imv_row, !is.na(.data$recorded_dttm)) |>
    group_by(.data$hospitalization_id) |>
    summarise(
      first_imv_dttm = min(.data$recorded_dttm),
      last_imv_dttm = max(.data$recorded_dttm),
      imv_record_count = n(),
      imv_duration_hours = as.numeric(difftime(max(.data$recorded_dttm), min(.data$recorded_dttm), units = "hours")),
      .groups = "drop"
    )
  ohca <- ohca |>
    left_join(imv_summary, by = "hospitalization_id", suffix = c("", "_calc")) |>
    mutate(
      imv_any = ifelse(!is.na(.data$first_imv_dttm), 1L, 0L),
      imv_record_count = tidyr::replace_na(.data$imv_record_count_calc, 0L),
      imv_duration_hours = .data$imv_duration_hours_calc
    ) |>
    select(-dplyr::any_of(c("imv_record_count_calc", "imv_duration_hours_calc")))
}

medication <- read_clif_table(tables_path, file_type, "medication_admin_continuous", required = FALSE)
ohca$vasopressor_any <- 0L
if (!is.null(medication) && nrow(medication) > 0) {
  vasoactive_categories <- c("angiotensin", "dobutamine", "dopamine", "epinephrine", "milrinone", "norepinephrine", "phenylephrine", "vasopressin")
  vasopressor_ids <- medication |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      med_category = stringr::str_to_lower(if ("med_category" %in% names(medication)) as.character(.data$med_category) else NA_character_),
      med_group = stringr::str_to_lower(if ("med_group" %in% names(medication)) as.character(.data$med_group) else NA_character_),
      mar_action_group = stringr::str_to_lower(if ("mar_action_group" %in% names(medication)) as.character(.data$mar_action_group) else NA_character_)
    ) |>
    semi_join(ohca |> select("hospitalization_id"), by = "hospitalization_id") |>
    filter(
      .data$med_group == "vasoactives",
      .data$med_category %in% vasoactive_categories,
      is.na(.data$mar_action_group) | .data$mar_action_group != "not_administered"
    ) |>
    distinct(.data$hospitalization_id)
  ohca <- ohca |>
    mutate(vasopressor_any = ifelse(.data$hospitalization_id %in% vasopressor_ids$hospitalization_id, 1L, 0L))
}

ohca <- ohca |>
  mutate(
    discharge_category_clean = stringr::str_to_lower(tidyr::replace_na(.data$discharge_category, "")),
    hospital_death = ifelse(stringr::str_detect(.data$discharge_category_clean, "expired|death|dead"), 1L, 0L),
    death_or_hospice = ifelse(.data$hospital_death == 1L | stringr::str_detect(.data$discharge_category_clean, "hospice"), 1L, 0L)
  )

summary_tbl <- tibble::tibble(
  site_name = site_name,
  diagnosis_source = diagnosis_source,
  n_all_icu_admissions = dplyr::n_distinct(all_icu$hospitalization_id),
  n_ohca_admissions = dplyr::n_distinct(ohca$hospitalization_id),
  n_ohca_patients = dplyr::n_distinct(ohca$patient_id),
  n_with_home_county_fips = sum(!is.na(ohca$home_county_fips)),
  n_with_assigned_county_fips = sum(!is.na(ohca$county_fips)),
  n_county_fips_overridden = sum(ohca$county_fips_was_overridden == 1, na.rm = TRUE),
  n_assigned_hospitals = dplyr::n_distinct(ohca$assigned_hospital_id),
  assigned_hospital_county_fips = paste(sort(unique(stats::na.omit(ohca$assigned_hospital_county_fips))), collapse = " | "),
  admission_start = as.character(min(ohca$admission_date)),
  admission_end = as.character(max(ohca$admission_date))
)

daily_counts <- ohca |>
  count(.data$admission_date, name = "ohca_admissions") |>
  arrange(.data$admission_date)

arrow::write_parquet(ohca, sink = file.path(output_dir, "ohca_poa_icu_2018_2024.parquet"))
readr::write_csv(ohca, file.path(output_dir, "ohca_poa_icu_2018_2024.csv"))
readr::write_csv(summary_tbl, file.path(output_dir, "ohca_poa_icu_2018_2024_summary.csv"))
readr::write_csv(daily_counts, file.path(output_dir, "ohca_daily_counts_2018_2024.csv"))

print(summary_tbl)
message("Wrote OHCA cohort outputs to ", output_dir)

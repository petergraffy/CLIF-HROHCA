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

config <- load_project_config(repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)
output_dir <- file.path(repo_root, "output", "intermediate", "cohorts", "all_icu")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hospitalization <- read_clif_table(tables_path, file_type, "hospitalization")
adt <- read_clif_table(tables_path, file_type, "adt")

tmax <- read_exposome_daymet(repo_root, "daymet_county_tmax_2018_2024_conus.parquet", "tmax_mean_c")
rmax <- read_exposome_daymet(repo_root, "daymet_county_rmax_2018_2024.parquet", "rmax_mean_pct")
no2 <- read_exposome_pollution(repo_root, "no2_county_year.csv", "no2_mean")
pm25 <- read_exposome_pollution(repo_root, "pm25_county_year.csv", "pm25_mean")

hosp <- hospitalization |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    age_at_admission = suppressWarnings(as.numeric(.data$age_at_admission)),
    county_fips = normalize_county_fips(if ("county_code" %in% names(hospitalization)) .data$county_code else NA_character_),
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
    location_category = if ("location_category" %in% names(adt)) as.character(.data$location_category) else NA_character_,
    hospital_id = if ("hospital_id" %in% names(adt)) as.character(.data$hospital_id) else NA_character_
  )

icu_ids <- adt_min |>
  filter(stringr::str_detect(stringr::str_to_lower(.data$location_category), "icu")) |>
  distinct(.data$hospitalization_id)

hospital_lookup <- build_hospital_id_lookup(adt_min)

all_icu <- hosp |>
  inner_join(icu_ids, by = "hospitalization_id") |>
  left_join(hospital_lookup, by = "hospitalization_id") |>
  distinct(.data$hospitalization_id, .keep_all = TRUE) |>
  arrange(.data$admission_dttm, .data$hospitalization_id)

all_icu <- apply_site_county_assignment(all_icu, repo_root, config) |>
  mutate(year = as.integer(format(.data$admission_date, "%Y")))

merged <- all_icu |>
  filter(!is.na(.data$county_fips)) |>
  select("hospitalization_id", "admission_date", "county_fips", "year") |>
  left_join(tmax, by = c("admission_date", "county_fips")) |>
  left_join(rmax, by = c("admission_date", "county_fips")) |>
  left_join(no2, by = c("county_fips", "year")) |>
  left_join(pm25, by = c("county_fips", "year"))

daily <- merged |>
  group_by(.data$admission_date) |>
  summarise(
    n_icu_admissions = n_distinct(.data$hospitalization_id),
    n_icu_admissions_with_tmax = sum(!is.na(.data$tmax_mean_c)),
    n_icu_admissions_with_rmax = sum(!is.na(.data$rmax_mean_pct)),
    n_icu_admissions_with_no2 = sum(!is.na(.data$no2_mean)),
    n_icu_admissions_with_pm25 = sum(!is.na(.data$pm25_mean)),
    icu_patient_address_mean_tmax_c = mean(.data$tmax_mean_c, na.rm = TRUE),
    icu_patient_address_p95_tmax_c = quantile(.data$tmax_mean_c, 0.95, na.rm = TRUE, names = FALSE),
    icu_patient_address_mean_rmax_pct = mean(.data$rmax_mean_pct, na.rm = TRUE),
    icu_patient_address_p95_rmax_pct = quantile(.data$rmax_mean_pct, 0.95, na.rm = TRUE, names = FALSE),
    icu_patient_address_mean_no2 = mean(.data$no2_mean, na.rm = TRUE),
    icu_patient_address_mean_pm25 = mean(.data$pm25_mean, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x))) |>
  arrange(.data$admission_date)

summary_tbl <- tibble::tibble(
  n_all_icu_admissions = n_distinct(all_icu$hospitalization_id),
  n_with_home_county_fips = sum(!is.na(all_icu$home_county_fips)),
  n_with_assigned_county_fips = sum(!is.na(all_icu$county_fips)),
  n_county_fips_overridden = sum(all_icu$county_fips_was_overridden == 1, na.rm = TRUE),
  n_assigned_hospitals = n_distinct(all_icu$assigned_hospital_id),
  assigned_hospital_county_fips = paste(sort(unique(stats::na.omit(all_icu$assigned_hospital_county_fips))), collapse = " | "),
  n_days_with_icu_admissions = n_distinct(daily$admission_date),
  n_days_with_linked_tmax = sum(!is.na(daily$icu_patient_address_mean_tmax_c)),
  n_days_with_linked_rmax = sum(!is.na(daily$icu_patient_address_mean_rmax_pct)),
  n_days_with_linked_no2 = sum(!is.na(daily$icu_patient_address_mean_no2)),
  n_days_with_linked_pm25 = sum(!is.na(daily$icu_patient_address_mean_pm25)),
  date_start = as.character(min(daily$admission_date)),
  date_end = as.character(max(daily$admission_date))
)

readr::write_csv(daily, file.path(output_dir, "all_icu_daily_patient_address_tmax_2018_2024.csv"))
readr::write_csv(summary_tbl, file.path(output_dir, "all_icu_daily_patient_address_tmax_2018_2024_summary.csv"))

print(summary_tbl)
message("Wrote all-ICU daily exposure outputs to ", output_dir)

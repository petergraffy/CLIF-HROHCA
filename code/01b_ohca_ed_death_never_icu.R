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

normalize_location_category <- function(x) {
  clean <- x |>
    tidyr::replace_na("Unknown") |>
    as.character() |>
    stringr::str_squish() |>
    stringr::str_to_lower()
  dplyr::na_if(clean, "")
}

is_ed_location <- function(x) {
  clean <- normalize_location_category(x)
  stringr::str_detect(clean, "(^|[^a-z])ed([^a-z]|$)|emerg")
}

collapse_pathway <- function(x) {
  x <- tidyr::replace_na(x, "unknown")
  if (length(x) == 0L) return("unknown")
  x <- x[c(TRUE, x[-1] != x[-length(x)])]
  paste(x, collapse = " -> ")
}

config <- load_project_config(repo_root)
validate_site_name(config, repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)
site_name <- config$site_name %||% Sys.getenv("CLIF_SITE_NAME", unset = "unknown_site")

output_dir <- file.path(repo_root, "output", "final", "descriptive")
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
  transmute(patient_id = as.character(.data$patient_id))

hosp <- hospitalization |>
  transmute(
    patient_id = as.character(.data$patient_id),
    hospitalization_id = as.character(.data$hospitalization_id),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    age_at_admission = suppressWarnings(as.numeric(.data$age_at_admission)),
    discharge_category = if ("discharge_category" %in% names(hospitalization)) as.character(.data$discharge_category) else NA_character_
  ) |>
  filter(
    !is.na(.data$admission_dttm),
    .data$admission_dttm >= START_DATE,
    .data$admission_dttm < END_DATE,
    .data$age_at_admission >= ADULT_AGE_YEARS
  ) |>
  left_join(patient_min, by = "patient_id")

dx <- diagnosis |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    diagnosis_code = if ("diagnosis_code" %in% names(diagnosis)) as.character(.data$diagnosis_code) else NA_character_,
    diagnosis_code_format = if ("diagnosis_code_format" %in% names(diagnosis)) as.character(.data$diagnosis_code_format) else NA_character_
  ) |>
  mutate(
    diagnosis_code_clean = norm_code(.data$diagnosis_code),
    diagnosis_code_format = infer_diagnosis_format(.data$diagnosis_code_clean, .data$diagnosis_code_format)
  )

dx_ohca_any <- dx |>
  filter(
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
    ohca_dx = 1L,
    ohca_codes = paste(sort(unique(.data$diagnosis_code_clean)), collapse = " | "),
    .groups = "drop"
  )

adt_location_summary <- adt |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    in_dttm = as_utc_datetime(.data$in_dttm),
    out_dttm = as_utc_datetime(.data$out_dttm),
    location_category = if ("location_category" %in% names(adt)) as.character(.data$location_category) else NA_character_
  ) |>
  filter(!is.na(.data$in_dttm)) |>
  mutate(
    location_category_clean = normalize_location_category(.data$location_category),
    location_is_ed = is_ed_location(.data$location_category),
    location_is_icu = stringr::str_detect(stringr::str_to_lower(tidyr::replace_na(.data$location_category, "")), "icu")
  ) |>
  arrange(.data$hospitalization_id, .data$in_dttm, .data$out_dttm) |>
  group_by(.data$hospitalization_id) |>
  summarise(
    n_adt_segments = n(),
    first_location_category = first(.data$location_category_clean),
    final_location_category = dplyr::last(.data$location_category_clean),
    care_pathway = collapse_pathway(.data$location_category_clean),
    first_location_is_ed = first(.data$location_is_ed),
    final_location_is_ed = dplyr::last(.data$location_is_ed),
    all_locations_are_ed = all(.data$location_is_ed, na.rm = TRUE),
    any_icu_location = any(.data$location_is_icu, na.rm = TRUE),
    .groups = "drop"
  )

ohca_all <- hosp |>
  inner_join(dx_ohca_any, by = "hospitalization_id") |>
  left_join(adt_location_summary, by = "hospitalization_id") |>
  mutate(
    first_location_is_ed = tidyr::replace_na(.data$first_location_is_ed, FALSE),
    final_location_is_ed = tidyr::replace_na(.data$final_location_is_ed, FALSE),
    all_locations_are_ed = tidyr::replace_na(.data$all_locations_are_ed, FALSE),
    any_icu_location = tidyr::replace_na(.data$any_icu_location, FALSE),
    hospital_death = ifelse(is_expired_discharge(.data$discharge_category), 1L, 0L),
    ed_first_only_never_icu = .data$first_location_is_ed &
      .data$all_locations_are_ed &
      !.data$any_icu_location,
    ed_first_only_death_never_icu = .data$ed_first_only_never_icu &
      .data$final_location_is_ed &
      .data$hospital_death == 1L
  )

summary_tbl <- tibble::tibble(
  site_name = site_name,
  diagnosis_source = diagnosis_source,
  study_start = as.character(as.Date(START_DATE)),
  study_end = as.character(as.Date(END_DATE) - 1),
  n_adult_hospitalizations = dplyr::n_distinct(hosp$hospitalization_id),
  n_ohca_dx_hospitalizations = dplyr::n_distinct(ohca_all$hospitalization_id),
  n_ohca_dx_patients = dplyr::n_distinct(ohca_all$patient_id),
  n_ohca_dx_first_location_ed = dplyr::n_distinct(ohca_all$hospitalization_id[ohca_all$first_location_is_ed]),
  n_ohca_dx_ed_first_only_never_icu = dplyr::n_distinct(ohca_all$hospitalization_id[ohca_all$ed_first_only_never_icu]),
  n_ohca_dx_ed_first_only_death_never_icu_hospitalizations = dplyr::n_distinct(ohca_all$hospitalization_id[ohca_all$ed_first_only_death_never_icu]),
  n_ohca_dx_ed_first_only_death_never_icu_patients = dplyr::n_distinct(ohca_all$patient_id[ohca_all$ed_first_only_death_never_icu]),
  definition = "Adult hospitalization with any OHCA diagnosis code, no POA restriction, first ADT location ED/emergency, all recorded ADT locations ED/emergency, no ICU ADT location, final ADT location ED/emergency, and expired/death discharge category."
)

pathway_summary <- ohca_all |>
  count(.data$care_pathway, .data$hospital_death, name = "n", sort = TRUE)

readr::write_csv(summary_tbl, file.path(output_dir, "ohca_ed_only_death_never_icu_summary.csv"))
readr::write_csv(pathway_summary, file.path(output_dir, "ohca_ed_only_death_never_icu_pathway_audit.csv"))

print(summary_tbl)
message("Wrote OHCA ED-only death never-ICU outputs to ", output_dir)

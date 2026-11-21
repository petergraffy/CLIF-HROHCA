# Run through 01_build_ohca_cohort.R --aim1-export.
export_aim1_counts <- function(path, hosp_source, hosp, all_icu, probable, ohca, adt,
                               hospitalization, diagnosis_source, repo_root, config, start, end) {
  geography <- get_site_geography(repo_root, config)$geography
  if (anyNA(geography$hospital_county_fips)) stop("Missing county in the hospital crosswalk.")
  key_number <- function(x) sub("[.]0+$", "", x)
  stopifnot(!anyDuplicated(tolower(geography$hospital_id)),
            !anyDuplicated(key_number(geography$hospital_number)))
  first_icu <- adt |>
    mutate(hospital_id = normalize_hospital_id(.data$hospital_id)) |>
    filter(stringr::str_detect(tolower(.data$location_category), "icu")) |>
    inner_join(select(all_icu, hospitalization_id, first_icu_in),
               by = c("hospitalization_id", "in_dttm" = "first_icu_in")) |>
    group_by(hospitalization_id) |>
    summarise(
      ambiguous = n_distinct(hospital_id, na.rm = TRUE) > 1L,
      hospital_id = if (!ambiguous) first(hospital_id[!is.na(hospital_id)]) else NA_character_,
      .groups = "drop"
    )
  all <- all_icu |>
    select(hospitalization_id, admission_date, age_at_admission, sex_category, race_category,
           ethnicity_category, discharge_category) |>
    left_join(first_icu, by = "hospitalization_id")
  stopifnot(!anyDuplicated(all$hospitalization_id), !anyDuplicated(ohca$hospitalization_id))
  key <- match(tolower(all$hospital_id), tolower(geography$hospital_id))
  key <- coalesce(key, match(key_number(all$hospital_id), key_number(geography$hospital_number)))
  county <- geography$hospital_county_fips[key]
  unique_counties <- unique(geography$hospital_county_fips)
  fallback <- is.na(county) & length(unique_counties) == 1L
  county[fallback] <- unique_counties[[1L]]
  ed_path <- ohca$care_pathway_to_icu %in% c("ed -> icu", "ed -> procedural -> icu")
  all <- all |>
    mutate(site_name = config$site_name, county_fips = county, fallback = fallback,
           ohca = hospitalization_id %in% ohca$hospitalization_id,
           ohca_ed = hospitalization_id %in% .env$ohca$hospitalization_id[ed_path])
  audit <- all |>
    group_by(site_name) |>
    summarise(n_icu_admissions = n(), ohca_admissions = sum(ohca),
              unknown_county_icu = sum(is.na(county_fips)),
              unknown_county_ohca = sum(ohca & is.na(county_fips)),
              ambiguous_first_icu_hospital = sum(ambiguous, na.rm = TRUE),
              single_county_fallback_icu = sum(fallback),
              single_county_fallback_ohca = sum(fallback & ohca))
  daily <- all |>
    filter(!is.na(county_fips)) |>
    group_by(site_name, county_fips, admission_date) |>
    summarise(n_icu_admissions = n(), ohca_admissions = sum(ohca),
              ohca_ed_admissions = sum(ohca_ed), .groups = "drop")

  flow_step <- function(step, x) {
    ids <- distinct(x, hospitalization_id, .keep_all = TRUE)
    known <- !is.na(ids$patient_id) & nzchar(ids$patient_id)
    unique_patients <- n_distinct(ids$patient_id[known])
    has_icu <- "n_icu_segments" %in% names(ids)
    tibble::tibble(step, n_hospitalizations = nrow(ids), n_unique_patients = unique_patients,
                   n_repeat_hospitalizations = sum(known) - unique_patients,
                   n_missing_patient_id = sum(!known),
                   n_adt_icu_segments = if (has_icu) sum(ids$n_icu_segments) else NA_integer_,
                   n_with_multiple_adt_icu_segments = if (has_icu) sum(ids$n_icu_segments > 1L) else NA_integer_)
  }
  assigned_ids <- all$hospitalization_id[!is.na(county)]
  window <- filter(hosp_source, !is.na(admission_dttm), as.Date(admission_dttm) >= start,
                   as.Date(admission_dttm) <= end)
  flow <- bind_rows(
    flow_step("hospitalizations_in_extract", hosp_source),
    flow_step("hospitalizations_missing_admission_timestamp", filter(hosp_source, is.na(admission_dttm))),
    flow_step("hospitalizations_in_study_period", window),
    flow_step("hospitalizations_in_study_period_missing_age", filter(window, is.na(age_at_admission))),
    flow_step("hospitalizations_in_study_period_age_under_18", filter(window, age_at_admission < 18)),
    flow_step("adult_hospitalizations_in_study_period", hosp),
    flow_step("adult_icu_hospitalizations", all_icu),
    flow_step("adult_icu_hospitalizations_with_county", filter(all_icu, hospitalization_id %in% assigned_ids)),
    flow_step("present_on_admission_cardiac_arrest", probable),
    flow_step("eligible_care_pathway", filter(probable, allowed_ohca_care_pathway)),
    flow_step("eligible_ohca_icu_entry_lt_24h", ohca),
    flow_step("eligible_ohca_with_county", filter(ohca, hospitalization_id %in% assigned_ids)),
    flow_step("ed_pathway_ohca_with_county", filter(ohca, hospitalization_id %in% assigned_ids, ed_path)),
    flow_step("poa_arrest_disallowed_pathway", filter(probable, !allowed_ohca_care_pathway)),
    flow_step("allowed_pathway_icu_entry_ge_24h_or_missing",
              filter(probable, allowed_ohca_care_pathway, !icu_entry_within_24h)),
    flow_step("eligible_ohca_icu_entry_before_hospital_admission", filter(ohca, hours_admission_to_icu < 0))
  ) |>
    mutate(site_name = config$site_name, diagnosis_source = diagnosis_source, .before = 1)

  coverage <- tidyr::expand_grid(site_name = config$site_name, county_fips = sort(unique_counties),
                                admission_date = seq(start, end, by = "day")) |>
    left_join(select(daily, site_name, county_fips, admission_date, n_icu_admissions),
              by = c("site_name", "county_fips", "admission_date")) |>
    mutate(observed_status = if_else(is.na(n_icu_admissions), "no_admissions_observed", "admissions_observed"),
           site_status = NA_character_)

  admission_source <- if ("admission_type_category" %in% names(hospitalization))
    normalize_category(hospitalization$admission_type_category) else rep("", nrow(hospitalization))
  admission_type <- tibble::tibble(
    hospitalization_id = as.character(hospitalization$hospitalization_id),
    admission_source = case_when(
      admission_source == "osh" ~ "outside_hospital",
      admission_source %in% c("ed", "facility", "direct", "elective") ~ "other_known",
      TRUE ~ "unknown"
    )
  )
  stopifnot(!anyDuplicated(admission_type$hospitalization_id))
  hospital_moves <- adt |>
    semi_join(all, by = "hospitalization_id") |>
    group_by(hospitalization_id) |>
    summarise(n_hospitals = n_distinct(tolower(normalize_hospital_id(hospital_id)), na.rm = TRUE),
              .groups = "drop")
  transfer <- all |>
    left_join(admission_type, by = "hospitalization_id") |>
    left_join(hospital_moves, by = "hospitalization_id") |>
    mutate(discharge = normalize_category(discharge_category),
           discharge_acute_transfer = discharge == "acute care hospital",
           discharge_unknown = discharge %in% c("", "unknown", "missing"),
           multiple_hospitals = n_hospitals > 1L,
           adt_hospital_unknown = n_hospitals == 0L) |>
    group_by(site_name, county_fips, admission_date, admission_source) |>
    summarise(n_icu_admissions = n(), ohca_admissions = sum(ohca),
              ohca_ed_admissions = sum(ohca_ed),
              n_icu_acute_transfer_out = sum(discharge_acute_transfer),
              n_ohca_acute_transfer_out = sum(ohca & discharge_acute_transfer),
              n_icu_discharge_unknown = sum(discharge_unknown),
              n_ohca_discharge_unknown = sum(ohca & discharge_unknown),
              n_icu_multiple_hospitals_in_adt = sum(multiple_hospitals),
              n_ohca_multiple_hospitals_in_adt = sum(ohca & multiple_hospitals),
              n_icu_adt_hospital_unknown = sum(adt_hospital_unknown),
              n_ohca_adt_hospital_unknown = sum(ohca & adt_hospital_unknown), .groups = "drop")

  subgroup <- all |>
    mutate(age_group = if_else(age_at_admission < 65, "18-64", "65+"),
           sex = normalize_category(sex_category), race = normalize_category(race_category),
           ethnicity = normalize_category(ethnicity_category),
           sex = case_when(sex %in% c("male", "m") ~ "Male", sex %in% c("female", "f") ~ "Female",
                           sex %in% c("", "unknown") ~ "Unknown", TRUE ~ "Other"),
           race = case_when(is_black_race(race) ~ "Black or African American",
                            race == "white" ~ "White", race == "asian" ~ "Asian",
                            race == "american indian or alaska native" ~ "American Indian or Alaska Native",
                            race == "native hawaiian or other pacific islander" ~ "Native Hawaiian or Other Pacific Islander",
                            race %in% c("", "unknown") ~ "Unknown", TRUE ~ "Other"),
           ethnicity = case_when(ethnicity == "hispanic" ~ "Hispanic",
                                 ethnicity == "non-hispanic" ~ "Non-Hispanic", TRUE ~ "Unknown")) |>
    select(site_name, county_fips, admission_date, ohca, ohca_ed, age_group, sex, race, ethnicity) |>
    tidyr::pivot_longer(c(age_group, sex, race, ethnicity), names_to = "subgroup_variable",
                        values_to = "subgroup_level") |>
    group_by(site_name, county_fips, admission_date, subgroup_variable, subgroup_level) |>
    summarise(n_icu_admissions = n(), ohca_admissions = sum(ohca),
              ohca_ed_admissions = sum(ohca_ed), .groups = "drop")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(daily, path)
  readr::write_csv(audit, paste0(path, ".audit.csv"))
  readr::write_csv(flow, paste0(path, ".flow.csv"))
  readr::write_csv(coverage, paste0(path, ".coverage.csv"))
  readr::write_csv(transfer, paste0(path, ".transfers.csv"))
  readr::write_csv(subgroup, paste0(path, ".subgroups.csv"))
}

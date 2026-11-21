# Run through 01_build_ohca_cohort.R --aim1-export.
export_aim1_counts <- function(path, all_icu, ohca, adt, repo_root, config) {
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
    select(hospitalization_id, admission_date) |>
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
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(daily, path)
  readr::write_csv(audit, paste0(path, ".audit.csv"))
}

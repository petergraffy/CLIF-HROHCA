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
  library(dplyr)
  library(grid)
  library(jsonlite)
  library(lubridate)
  library(readr)
  library(stringr)
  library(tidyr)
})

config <- load_project_config(repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)

START_DATE <- as_utc_datetime("2018-01-01")
END_DATE <- as_utc_datetime("2025-01-01")
ADULT_AGE_YEARS <- 18
WARM_MONTHS <- c(5L, 6L, 7L, 8L, 9L)

ohca_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
outcome_dataset_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca_outcomes", "ohca_heat_adverse_outcome_analysis_dataset.csv")
output_dir <- file.path(repo_root, "output", "final", "quality_checks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

format_count <- function(x) {
  if (is.na(x)) return("Not available")
  format(as.integer(x), big.mark = ",", scientific = FALSE)
}

plot_cohort_flow <- function(cohort_flow, path) {
  png(path, width = 1800, height = 2400, res = 220)
  on.exit(dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::grid.text(
    "OHCA ICU Cohort Flow",
    x = 0.5,
    y = 0.965,
    gp = grid::gpar(fontface = "bold", fontsize = 18, col = "#1b2838")
  )
  grid::grid.text(
    "Aggregate site-level inclusion counts; no patient-level rows exported",
    x = 0.5,
    y = 0.94,
    gp = grid::gpar(fontsize = 9, col = "#576574")
  )

  box_x <- 0.5
  box_width <- 0.74
  box_height <- 0.082
  y_values <- seq(0.86, 0.15, length.out = nrow(cohort_flow))

  for (i in seq_len(nrow(cohort_flow))) {
    row <- cohort_flow[i, ]
    y <- y_values[[i]]
    excluded_text <- if (is.na(row$n_excluded_since_previous)) {
      ""
    } else {
      paste0("\nExcluded since prior step: ", format_count(row$n_excluded_since_previous))
    }
    note_text <- if (!is.na(row$notes) && nzchar(row$notes)) paste0("\n", row$notes) else ""
    label <- paste0(row$step, "\nN = ", format_count(row$n), excluded_text, note_text)

    grid::grid.roundrect(
      x = box_x,
      y = y,
      width = box_width,
      height = box_height,
      r = grid::unit(0.02, "npc"),
      gp = grid::gpar(fill = if (i %% 2 == 0) "#eef6f8" else "#f8fbf5", col = "#2d3436", lwd = 1.2)
    )
    grid::grid.text(
      label,
      x = box_x,
      y = y,
      gp = grid::gpar(fontsize = 8.5, lineheight = 0.95, col = "#1b2838")
    )
    if (i < nrow(cohort_flow)) {
      grid::grid.lines(
        x = c(box_x, box_x),
        y = c(y - box_height / 2 - 0.008, y_values[[i + 1]] + box_height / 2 + 0.008),
        arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed"),
        gp = grid::gpar(col = "#576574", lwd = 1)
      )
    }
  }
}

ohca <- readr::read_csv(ohca_path, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_guess(), hospitalization_id = readr::col_character()))
ohca <- ohca |>
  mutate(
    admission_date = as.Date(.data$admission_date),
    month = as.integer(format(.data$admission_date, "%m")),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    first_icu_in = as_utc_datetime(.data$first_icu_in),
    hours_admission_to_icu = as.numeric(difftime(.data$first_icu_in, .data$admission_dttm, units = "hours"))
  )

hospitalization <- read_clif_table(
  tables_path,
  file_type,
  "hospitalization",
  columns = c("patient_id", "hospitalization_id", "admission_dttm", "age_at_admission")
)
adt <- read_clif_table(
  tables_path,
  file_type,
  "adt",
  columns = c("hospitalization_id", "in_dttm", "out_dttm", "location_category")
)

if (file.exists(outcome_dataset_path)) {
  outcome_dataset <- readr::read_csv(outcome_dataset_path, show_col_types = FALSE)
  complete_model_cases <- sum(stats::complete.cases(outcome_dataset[, c("tmax_mean_c", "rmax_mean_pct", "no2_mean", "pm25_mean", "age_at_admission", "sex_category", "race_category", "dow", "year", "time_index")]))
  ventilated_with_duration <- sum(outcome_dataset$imv_any == 1 & !is.na(outcome_dataset$imv_duration_days))
} else {
  complete_model_cases <- NA_integer_
  ventilated_with_duration <- NA_integer_
}

hospitalization_flow <- hospitalization |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    age_at_admission = suppressWarnings(as.numeric(.data$age_at_admission))
  ) |>
  distinct(.data$hospitalization_id, .keep_all = TRUE)

adult_hospitalizations <- hospitalization_flow |>
  filter(.data$age_at_admission >= ADULT_AGE_YEARS)

study_period_adult_hospitalizations <- adult_hospitalizations |>
  filter(
    !is.na(.data$admission_dttm),
    .data$admission_dttm >= START_DATE,
    .data$admission_dttm < END_DATE
  )

icu_hospitalization_ids <- adt |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    location_category = as.character(.data$location_category)
  ) |>
  filter(stringr::str_detect(stringr::str_to_lower(tidyr::replace_na(.data$location_category, "")), "icu")) |>
  distinct(.data$hospitalization_id)

study_period_adult_icu_hospitalizations <- study_period_adult_hospitalizations |>
  semi_join(icu_hospitalization_ids, by = "hospitalization_id")

warm_season_ohca_n <- sum(ohca$month %in% WARM_MONTHS, na.rm = TRUE)
county_assigned_ohca_n <- if ("county_fips" %in% names(ohca)) {
  sum(!is.na(ohca$county_fips) & nzchar(as.character(ohca$county_fips)))
} else {
  NA_integer_
}

cohort_flow <- tibble::tibble(
  step_order = seq_len(8),
  step = c(
    "All CLIF hospitalizations",
    "Adult hospitalizations (age >= 18 years)",
    "Adult hospitalizations admitted 2018-2024",
    "Adult 2018-2024 hospitalizations with ICU stay",
    "Present-on-admission cardiac arrest / OHCA-proxy ICU admissions",
    "OHCA ICU admissions with assigned county exposure geography",
    "Warm-season OHCA ICU admissions (May-September)",
    "Complete cases for primary heat/outcome models"
  ),
  n = c(
    nrow(hospitalization_flow),
    nrow(adult_hospitalizations),
    nrow(study_period_adult_hospitalizations),
    nrow(study_period_adult_icu_hospitalizations),
    nrow(ohca),
    county_assigned_ohca_n,
    warm_season_ohca_n,
    complete_model_cases
  ),
  notes = c(
    "Source: clif_hospitalization",
    "",
    "",
    "ICU identified from clif_adt location_category",
    "Present-on-admission cardiac arrest diagnosis plus ICU stay",
    "Patient county retained if local/adjacent; otherwise hospital county",
    "Primary heat-phenotype season",
    "Nonmissing primary exposure/covariate fields"
  )
) |>
  mutate(
    n_excluded_since_previous = dplyr::lag(.data$n) - .data$n,
    n_excluded_since_previous = ifelse(.data$step_order == 1L, NA_integer_, .data$n_excluded_since_previous)
  )

denominator_audit <- tibble::tibble(
  step = c(
    "All OHCA ICU admissions, 2018-2024",
    "Warm-season OHCA admissions, May-September",
    "Complete cases for heat/outcome models",
    "Ventilated warm-season OHCA with IMV duration"
  ),
  n = c(
    nrow(ohca),
    sum(ohca$month %in% c(5L, 6L, 7L, 8L, 9L), na.rm = TRUE),
    complete_model_cases,
    ventilated_with_duration
  )
)

timing_summary <- tibble::tibble(
  n = sum(!is.na(ohca$hours_admission_to_icu)),
  median_hours = median(ohca$hours_admission_to_icu, na.rm = TRUE),
  iqr_low_hours = quantile(ohca$hours_admission_to_icu, 0.25, na.rm = TRUE, names = FALSE),
  iqr_high_hours = quantile(ohca$hours_admission_to_icu, 0.75, na.rm = TRUE, names = FALSE),
  p90_hours = quantile(ohca$hours_admission_to_icu, 0.90, na.rm = TRUE, names = FALSE),
  p95_hours = quantile(ohca$hours_admission_to_icu, 0.95, na.rm = TRUE, names = FALSE),
  within_6h_pct = mean(ohca$hours_admission_to_icu <= 6, na.rm = TRUE) * 100,
  within_12h_pct = mean(ohca$hours_admission_to_icu <= 12, na.rm = TRUE) * 100,
  within_24h_pct = mean(ohca$hours_admission_to_icu <= 24, na.rm = TRUE) * 100,
  within_48h_pct = mean(ohca$hours_admission_to_icu <= 48, na.rm = TRUE) * 100
)

timing_bins <- ohca |>
  mutate(
    icu_entry_timing = cut(
      .data$hours_admission_to_icu,
      breaks = c(-Inf, 0, 6, 12, 24, 48, 72, Inf),
      labels = c("<=0 hours", "0-6 hours", "6-12 hours", "12-24 hours", "24-48 hours", "48-72 hours", ">72 hours"),
      right = TRUE
    )
  ) |>
  count(.data$icu_entry_timing, name = "n") |>
  mutate(pct = 100 * .data$n / nrow(ohca))

adt_pre_icu <- adt |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    in_dttm = as_utc_datetime(.data$in_dttm),
    out_dttm = as_utc_datetime(.data$out_dttm),
    location_category = if ("location_category" %in% names(adt)) as.character(.data$location_category) else NA_character_
  ) |>
  inner_join(ohca |> select("hospitalization_id", "first_icu_in"), by = "hospitalization_id") |>
  filter(!is.na(.data$in_dttm), !is.na(.data$first_icu_in), .data$in_dttm <= .data$first_icu_in) |>
  arrange(.data$hospitalization_id, .data$in_dttm, .data$out_dttm) |>
  mutate(location_category_clean = dplyr::na_if(stringr::str_trim(tidyr::replace_na(.data$location_category, "Unknown")), ""))

collapse_pathway <- function(x) {
  x <- tidyr::replace_na(x, "Unknown")
  x <- x[c(TRUE, x[-1] != x[-length(x)])]
  paste(x, collapse = " -> ")
}

pathway_summary <- adt_pre_icu |>
  group_by(.data$hospitalization_id) |>
  summarise(
    care_pathway = collapse_pathway(.data$location_category_clean),
    first_location_category = first(.data$location_category_clean),
    .groups = "drop"
  ) |>
  count(.data$care_pathway, name = "n", sort = TRUE) |>
  mutate(pct = 100 * .data$n / n_distinct(ohca$hospitalization_id))

first_location_summary <- adt_pre_icu |>
  group_by(.data$hospitalization_id) |>
  summarise(first_location_category = first(.data$location_category_clean), .groups = "drop") |>
  count(.data$first_location_category, name = "n", sort = TRUE) |>
  mutate(pct = 100 * .data$n / n_distinct(ohca$hospitalization_id))

readr::write_csv(cohort_flow, file.path(output_dir, "cohort_flow.csv"))
readr::write_csv(denominator_audit, file.path(output_dir, "model_denominator_audit.csv"))
readr::write_csv(timing_summary, file.path(output_dir, "ohca_admission_to_icu_timing_summary.csv"))
readr::write_csv(timing_bins, file.path(output_dir, "ohca_admission_to_icu_timing_bins.csv"))
readr::write_csv(pathway_summary, file.path(output_dir, "ohca_pre_icu_care_pathway_summary.csv"))
readr::write_csv(first_location_summary, file.path(output_dir, "ohca_first_location_summary.csv"))
plot_cohort_flow(cohort_flow, file.path(output_dir, "cohort_flow_diagram.png"))

print(cohort_flow)
print(timing_summary)
print(head(pathway_summary, 15))
message("Wrote aggregate quality-check outputs to ", output_dir)

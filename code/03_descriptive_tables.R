#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    stop("Could not determine script path from commandArgs().")
  }
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(
  file.path(dirname(get_script_path()), ".."),
  winslash = "/",
  mustWork = TRUE
)

ensure_user_library <- function() {
  version_parts <- strsplit(as.character(getRversion()), "\\.")[[1]]
  version_stub <- paste(version_parts[1], version_parts[2], sep = ".")
  user_lib <- file.path(repo_root, ".r-user-lib", version_stub)
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}

ensure_user_library()

ohca_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
output_dir <- file.path(repo_root, "output", "final", "descriptive")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

safe_pct <- function(x, denom) {
  ifelse(denom > 0, 100 * x / denom, NA_real_)
}

fmt_n_pct <- function(n, denom) {
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), safe_pct(n, denom))
}

add_group_rows <- function(df, var, label) {
  counts <- sort(table(df[[var]], useNA = "ifany"), decreasing = TRUE)
  out <- data.frame(
    section = label,
    characteristic = names(counts),
    value = vapply(as.numeric(counts), fmt_n_pct, character(1), denom = nrow(df)),
    stringsAsFactors = FALSE
  )
  out$characteristic[is.na(out$characteristic) | out$characteristic == ""] <- "Missing"
  out
}

ohca <- read.csv(ohca_path, stringsAsFactors = FALSE)
ohca$admission_dttm <- as.POSIXct(ohca$admission_dttm, tz = "UTC")
ohca$discharge_dttm <- as.POSIXct(ohca$discharge_dttm, tz = "UTC")
ohca$first_icu_in <- as.POSIXct(ohca$first_icu_in, tz = "UTC")
ohca$last_icu_out <- as.POSIXct(ohca$last_icu_out, tz = "UTC")
ohca$age_at_admission <- as.numeric(ohca$age_at_admission)
ohca$county_fips_was_overridden <- as.integer(ohca$county_fips_was_overridden)
ohca$hospital_los_days <- as.numeric(difftime(ohca$discharge_dttm, ohca$admission_dttm, units = "days"))
ohca$time_to_icu_days <- as.numeric(difftime(ohca$first_icu_in, ohca$admission_dttm, units = "days"))
if (!"hospital_death" %in% names(ohca)) {
  ohca$hospital_death <- ifelse(ohca$discharge_category == "Expired", 1L, 0L)
}
if (!"death_or_hospice" %in% names(ohca)) {
  ohca$death_or_hospice <- ifelse(ohca$discharge_category == "Expired" | grepl("hospice", ohca$discharge_category, ignore.case = TRUE), 1L, 0L)
}
if (!"imv_any" %in% names(ohca)) ohca$imv_any <- NA_integer_
if (!"imv_duration_hours" %in% names(ohca)) ohca$imv_duration_hours <- NA_real_
if (!"vasopressor_any" %in% names(ohca)) ohca$vasopressor_any <- NA_integer_
ohca$age_group <- ifelse(ohca$age_at_admission >= 65, ">=65", "<65")

table1_rows <- list(
  data.frame(section = "Cohort", characteristic = "OHCA ICU admissions", value = format(nrow(ohca), big.mark = ","), stringsAsFactors = FALSE),
  data.frame(section = "Age", characteristic = "Age, median [IQR]", value = sprintf("%.1f [%.1f, %.1f]", median(ohca$age_at_admission, na.rm = TRUE), quantile(ohca$age_at_admission, 0.25, na.rm = TRUE), quantile(ohca$age_at_admission, 0.75, na.rm = TRUE)), stringsAsFactors = FALSE),
  add_group_rows(ohca, "age_group", "Age group"),
  add_group_rows(ohca, "sex_category", "Sex"),
  add_group_rows(ohca, "race_category", "Race"),
  add_group_rows(ohca, "ethnicity_category", "Ethnicity"),
  add_group_rows(ohca, "discharge_category", "Discharge category"),
  data.frame(section = "Geography", characteristic = "County reassigned to hospital county", value = fmt_n_pct(sum(ohca$county_fips_was_overridden == 1, na.rm = TRUE), nrow(ohca)), stringsAsFactors = FALSE),
  data.frame(section = "ICU utilization", characteristic = "ICU LOS days, median [IQR]", value = sprintf("%.1f [%.1f, %.1f]", median(ohca$icu_los_hours / 24, na.rm = TRUE), quantile(ohca$icu_los_hours / 24, 0.25, na.rm = TRUE), quantile(ohca$icu_los_hours / 24, 0.75, na.rm = TRUE)), stringsAsFactors = FALSE),
  data.frame(section = "ICU therapies", characteristic = "Invasive mechanical ventilation", value = fmt_n_pct(sum(ohca$imv_any == 1, na.rm = TRUE), sum(!is.na(ohca$imv_any))), stringsAsFactors = FALSE),
  data.frame(section = "ICU therapies", characteristic = "IMV duration days among ventilated, median [IQR]", value = sprintf("%.1f [%.1f, %.1f]", median(ohca$imv_duration_hours[ohca$imv_any == 1] / 24, na.rm = TRUE), quantile(ohca$imv_duration_hours[ohca$imv_any == 1] / 24, 0.25, na.rm = TRUE), quantile(ohca$imv_duration_hours[ohca$imv_any == 1] / 24, 0.75, na.rm = TRUE)), stringsAsFactors = FALSE),
  data.frame(section = "ICU therapies", characteristic = "Vasopressor use", value = fmt_n_pct(sum(ohca$vasopressor_any == 1, na.rm = TRUE), sum(!is.na(ohca$vasopressor_any))), stringsAsFactors = FALSE),
  data.frame(section = "Hospital utilization", characteristic = "Hospital LOS days, median [IQR]", value = sprintf("%.1f [%.1f, %.1f]", median(ohca$hospital_los_days, na.rm = TRUE), quantile(ohca$hospital_los_days, 0.25, na.rm = TRUE), quantile(ohca$hospital_los_days, 0.75, na.rm = TRUE)), stringsAsFactors = FALSE),
  data.frame(section = "Outcomes", characteristic = "Hospital mortality", value = fmt_n_pct(sum(ohca$hospital_death == 1, na.rm = TRUE), nrow(ohca)), stringsAsFactors = FALSE),
  data.frame(section = "Outcomes", characteristic = "Death or discharge to hospice", value = fmt_n_pct(sum(ohca$death_or_hospice == 1, na.rm = TRUE), nrow(ohca)), stringsAsFactors = FALSE),
  data.frame(section = "Outcomes", characteristic = "Time to death among decedents, days median [IQR]", value = sprintf("%.1f [%.1f, %.1f]", median(ohca$hospital_los_days[ohca$hospital_death == 1], na.rm = TRUE), quantile(ohca$hospital_los_days[ohca$hospital_death == 1], 0.25, na.rm = TRUE), quantile(ohca$hospital_los_days[ohca$hospital_death == 1], 0.75, na.rm = TRUE)), stringsAsFactors = FALSE)
)

table1 <- do.call(rbind, table1_rows)

outcomes <- data.frame(
  metric = c(
    "Hospital mortality",
    "Death or hospice",
    "Invasive mechanical ventilation",
    "IMV duration days median among ventilated",
    "IMV duration days IQR low among ventilated",
    "IMV duration days IQR high among ventilated",
    "Vasopressor use",
    "ICU LOS days median",
    "ICU LOS days IQR low",
    "ICU LOS days IQR high",
    "Hospital LOS days median",
    "Hospital LOS days IQR low",
    "Hospital LOS days IQR high",
    "Time to death days median",
    "Time to death days IQR low",
    "Time to death days IQR high"
  ),
  value = c(
    safe_pct(sum(ohca$hospital_death == 1, na.rm = TRUE), nrow(ohca)),
    safe_pct(sum(ohca$death_or_hospice == 1, na.rm = TRUE), nrow(ohca)),
    safe_pct(sum(ohca$imv_any == 1, na.rm = TRUE), sum(!is.na(ohca$imv_any))),
    median(ohca$imv_duration_hours[ohca$imv_any == 1] / 24, na.rm = TRUE),
    quantile(ohca$imv_duration_hours[ohca$imv_any == 1] / 24, 0.25, na.rm = TRUE),
    quantile(ohca$imv_duration_hours[ohca$imv_any == 1] / 24, 0.75, na.rm = TRUE),
    safe_pct(sum(ohca$vasopressor_any == 1, na.rm = TRUE), sum(!is.na(ohca$vasopressor_any))),
    median(ohca$icu_los_hours / 24, na.rm = TRUE),
    quantile(ohca$icu_los_hours / 24, 0.25, na.rm = TRUE),
    quantile(ohca$icu_los_hours / 24, 0.75, na.rm = TRUE),
    median(ohca$hospital_los_days, na.rm = TRUE),
    quantile(ohca$hospital_los_days, 0.25, na.rm = TRUE),
    quantile(ohca$hospital_los_days, 0.75, na.rm = TRUE),
    median(ohca$hospital_los_days[ohca$hospital_death == 1], na.rm = TRUE),
    quantile(ohca$hospital_los_days[ohca$hospital_death == 1], 0.25, na.rm = TRUE),
    quantile(ohca$hospital_los_days[ohca$hospital_death == 1], 0.75, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(table1, file.path(output_dir, "table1_ohca_cohort_characteristics.csv"), row.names = FALSE)
write.csv(outcomes, file.path(output_dir, "ohca_outcomes_summary.csv"), row.names = FALSE)

message("Wrote descriptive OHCA outputs to ", output_dir)

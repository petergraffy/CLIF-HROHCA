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

ohca_path <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
tmax_path <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
output_dir <- file.path(repo_root, "output", "final", "descriptive")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_county_fips <- function(x) {
  clean <- gsub("\\.0$", "", as.character(x))
  clean <- gsub("[^0-9]", "", clean)
  clean <- trimws(clean)
  clean <- ifelse(nchar(clean) == 0, NA_character_, clean)
  clean <- ifelse(is.na(clean), NA_character_, sprintf("%05s", clean))
  gsub(" ", "0", clean, fixed = TRUE)
}

safe_pct <- function(x, denom) ifelse(denom > 0, 100 * x / denom, NA_real_)

fmt_n_pct <- function(n, denom) {
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), safe_pct(n, denom))
}

fmt_median_iqr <- function(x) {
  sprintf(
    "%.1f [%.1f, %.1f]",
    median(x, na.rm = TRUE),
    stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  )
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

safe_fisher_p <- function(group, value) {
  ok <- !is.na(group) & !is.na(value)
  if (length(unique(group[ok])) < 2 || length(unique(value[ok])) < 2) return(NA_real_)
  stats::fisher.test(table(group[ok], value[ok]), simulate.p.value = TRUE, B = 10000)$p.value
}

safe_wilcox_p <- function(group, value) {
  ok <- !is.na(group) & !is.na(value)
  if (length(unique(group[ok])) < 2 || any(table(group[ok]) == 0)) return(NA_real_)
  suppressWarnings(stats::wilcox.test(value[ok] ~ group[ok])$p.value)
}

ohca <- read.csv(ohca_path, stringsAsFactors = FALSE)
ohca$admission_date <- as.Date(ohca$admission_date)
ohca$admission_dttm <- as.POSIXct(ohca$admission_dttm, tz = "UTC")
ohca$discharge_dttm <- as.POSIXct(ohca$discharge_dttm, tz = "UTC")
ohca$hospital_los_days <- as.numeric(difftime(ohca$discharge_dttm, ohca$admission_dttm, units = "days"))
if (!"hospital_death" %in% names(ohca)) {
  ohca$hospital_death <- ifelse(ohca$discharge_category == "Expired", 1L, 0L)
}
if (!"death_or_hospice" %in% names(ohca)) {
  ohca$death_or_hospice <- ifelse(ohca$discharge_category == "Expired" | grepl("hospice", ohca$discharge_category, ignore.case = TRUE), 1L, 0L)
}
if (!"imv_any" %in% names(ohca)) ohca$imv_any <- NA_integer_
if (!"imv_duration_hours" %in% names(ohca)) ohca$imv_duration_hours <- NA_real_
if (!"vasopressor_any" %in% names(ohca)) ohca$vasopressor_any <- NA_integer_

suppressPackageStartupMessages(library(arrow))
tmax <- arrow::read_parquet(tmax_path)
tmax$county_fips <- normalize_county_fips(tmax$geoid)
tmax$date <- as.Date(tmax$date)

warm_tmax <- tmax[as.integer(format(tmax$date, "%m")) %in% c(5L,6L,7L,8L,9L), c("county_fips","date","tmax_mean_c")]
ohca2 <- merge(
  ohca,
  warm_tmax,
  by.x = c("county_fips", "admission_date"),
  by.y = c("county_fips", "date"),
  all.x = TRUE
)
threshold <- as.numeric(stats::quantile(ohca2$tmax_mean_c, 0.95, na.rm = TRUE, names = FALSE))
ohca2$heat_related_ohca <- ifelse(!is.na(ohca2$tmax_mean_c) & ohca2$tmax_mean_c >= threshold, "Heat-related OHCA", "Non-heat-related OHCA")
ohca2$heat_related_ohca <- factor(ohca2$heat_related_ohca, levels = c("Heat-related OHCA", "Non-heat-related OHCA"))
ohca2$age_at_admission <- as.numeric(ohca2$age_at_admission)
ohca2$age_group <- ifelse(ohca2$age_at_admission >= 65, ">=65", "<65")
ohca2$race_group <- ifelse(ohca2$race_category == "Black or African American", "Black", "Non-Black")
ohca2$icu_los_days <- as.numeric(ohca2$icu_los_hours) / 24
ohca2$imv_duration_days <- as.numeric(ohca2$imv_duration_hours) / 24

summarize_group <- function(df, label) {
  data.frame(
    group = label,
    n = nrow(df),
    median_age = median(as.numeric(df$age_at_admission), na.rm = TRUE),
    hospital_mortality_pct = safe_pct(sum(df$hospital_death == 1, na.rm = TRUE), nrow(df)),
    death_or_hospice_pct = safe_pct(sum(df$death_or_hospice == 1, na.rm = TRUE), nrow(df)),
    imv_pct = safe_pct(sum(df$imv_any == 1, na.rm = TRUE), sum(!is.na(df$imv_any))),
    median_imv_duration_days = median(df$imv_duration_hours[df$imv_any == 1] / 24, na.rm = TRUE),
    vasopressor_pct = safe_pct(sum(df$vasopressor_any == 1, na.rm = TRUE), sum(!is.na(df$vasopressor_any))),
    median_icu_los_days = median(as.numeric(df$icu_los_hours) / 24, na.rm = TRUE),
    median_hospital_los_days = median(df$hospital_los_days, na.rm = TRUE),
    median_time_to_death_days = median(df$hospital_los_days[df$hospital_death == 1], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(split(ohca2, ohca2$heat_related_ohca), function(x) summarize_group(x, unique(x$heat_related_ohca)[1])))

by_discharge <- as.data.frame.matrix(table(ohca2$heat_related_ohca, ohca2$discharge_category, useNA = "ifany"))
by_discharge <- cbind(group = rownames(by_discharge), by_discharge, row.names = NULL)

make_binary_row <- function(label, var) {
  heat <- ohca2[ohca2$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- ohca2[ohca2$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  data.frame(
    characteristic = label,
    heat_related_ohca = fmt_n_pct(sum(heat[[var]] == 1, na.rm = TRUE), sum(!is.na(heat[[var]]))),
    non_heat_related_ohca = fmt_n_pct(sum(non_heat[[var]] == 1, na.rm = TRUE), sum(!is.na(non_heat[[var]]))),
    p_value = fmt_p(safe_fisher_p(ohca2$heat_related_ohca, ohca2[[var]])),
    stringsAsFactors = FALSE
  )
}

make_continuous_row <- function(label, var, subset_var = NULL) {
  df <- ohca2
  if (!is.null(subset_var)) df <- df[df[[subset_var]] == 1, , drop = FALSE]
  heat <- df[df$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- df[df$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  data.frame(
    characteristic = label,
    heat_related_ohca = fmt_median_iqr(heat[[var]]),
    non_heat_related_ohca = fmt_median_iqr(non_heat[[var]]),
    p_value = fmt_p(safe_wilcox_p(df$heat_related_ohca, df[[var]])),
    stringsAsFactors = FALSE
  )
}

add_tmp_category_rows <- function(label, var) {
  levels_found <- names(sort(table(ohca2[[var]], useNA = "no"), decreasing = TRUE))
  p <- safe_fisher_p(ohca2$heat_related_ohca, ohca2[[var]])
  do.call(rbind, lapply(seq_along(levels_found), function(i) {
    level_value <- levels_found[[i]]
    heat <- ohca2[ohca2$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
    non_heat <- ohca2[ohca2$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
    data.frame(
      characteristic = paste0(label, ": ", level_value),
      heat_related_ohca = fmt_n_pct(sum(heat[[var]] == level_value, na.rm = TRUE), sum(!is.na(heat[[var]]))),
      non_heat_related_ohca = fmt_n_pct(sum(non_heat[[var]] == level_value, na.rm = TRUE), sum(!is.na(non_heat[[var]]))),
      p_value = ifelse(i == 1, fmt_p(p), ""),
      stringsAsFactors = FALSE
    )
  }))
}

heat <- ohca2[ohca2$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
non_heat <- ohca2[ohca2$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
table2_rows <- list(
  data.frame(
    characteristic = "N",
    heat_related_ohca = format(nrow(heat), big.mark = ","),
    non_heat_related_ohca = format(nrow(non_heat), big.mark = ","),
    p_value = "",
    stringsAsFactors = FALSE
  ),
  make_continuous_row("Age, years", "age_at_admission"),
  add_tmp_category_rows("Age group", "age_group"),
  add_tmp_category_rows("Sex", "sex_category"),
  add_tmp_category_rows("Race", "race_group"),
  add_tmp_category_rows("Ethnicity", "ethnicity_category"),
  make_binary_row("County reassigned to hospital county", "county_fips_was_overridden"),
  make_continuous_row("Assigned-county Tmax, C", "tmax_mean_c"),
  make_binary_row("Invasive mechanical ventilation", "imv_any"),
  make_continuous_row("IMV duration among ventilated, days", "imv_duration_days", subset_var = "imv_any"),
  make_binary_row("Vasopressor use", "vasopressor_any"),
  make_continuous_row("ICU length of stay, days", "icu_los_days"),
  make_continuous_row("Hospital length of stay, days", "hospital_los_days"),
  make_binary_row("Hospital death", "hospital_death"),
  make_binary_row("Death or hospice", "death_or_hospice"),
  make_continuous_row("Time to death among decedents, days", "hospital_los_days", subset_var = "hospital_death")
)
table2 <- do.call(rbind, table2_rows)

write.csv(summary_df, file.path(output_dir, "heat_related_vs_non_heat_related_ohca_outcomes.csv"), row.names = FALSE)
write.csv(by_discharge, file.path(output_dir, "heat_related_vs_non_heat_related_discharge_categories.csv"), row.names = FALSE)
write.csv(table2, file.path(output_dir, "table2_heat_related_vs_non_heat_related_ohca.csv"), row.names = FALSE)
write.csv(data.frame(heat_threshold_tmax_c = threshold), file.path(output_dir, "heat_related_ohca_threshold.csv"), row.names = FALSE)

message("Wrote heat-related OHCA outcome comparison outputs to ", output_dir)

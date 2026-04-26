#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) stop("Could not determine script path from commandArgs().")
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)

ensure_user_library <- function() {
  version_parts <- strsplit(as.character(getRversion()), "\\.")[[1]]
  version_stub <- paste(version_parts[1], version_parts[2], sep = ".")
  user_lib <- file.path(repo_root, ".r-user-lib", version_stub)
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}
ensure_user_library()

config_path <- file.path(repo_root, "config", "config.json")
output_dir <- file.path(repo_root, "output", "final", "federated_exports")
figure_source_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
figure_export_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_export_dir, recursive = TRUE, showWarnings = FALSE)

config <- jsonlite::fromJSON(config_path)
site_name <- config$site_name

results_path <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript", "manuscript_dlnm_results.csv")
curves_path <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript", "manuscript_dlnm_curves.csv")
reduced_coef_path <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript", "manuscript_dlnm_reduced_coefficients.csv")
reduced_vcov_path <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript", "manuscript_dlnm_reduced_vcov.csv")
table1_path <- file.path(repo_root, "output", "final", "descriptive", "table1_ohca_cohort_characteristics.csv")
outcomes_path <- file.path(repo_root, "output", "final", "descriptive", "ohca_outcomes_summary.csv")
heat_table2_path <- file.path(repo_root, "output", "final", "descriptive", "table2_heat_related_vs_non_heat_related_ohca.csv")
heat_table2_all_path <- file.path(repo_root, "output", "final", "descriptive", "table2_heat_related_vs_non_heat_related_ohca_all_definitions.csv")
heat90_table2_path <- file.path(repo_root, "output", "final", "descriptive", "table2_heat90_vs_non_heat90_ohca.csv")
heat_thresholds_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_thresholds.csv")
heat_summary_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_vs_non_heat_related_ohca_outcomes.csv")
heat_discharge_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_vs_non_heat_related_discharge_categories.csv")
heat_hourly_vitals_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_hourly_vital_trajectories.csv")
heat_hourly_labs_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_hourly_lab_trajectories.csv")
heat_hourly_support_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_hourly_support_trajectories.csv")
heat_hourly_vitals_smoothed_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_hourly_vital_trajectories_smoothed.csv")
heat_hourly_labs_smoothed_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_hourly_lab_trajectories_smoothed.csv")
heat_hourly_support_smoothed_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_hourly_support_trajectories_smoothed.csv")
heat_hourly_cumulative_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_hourly_cumulative_incidence.csv")
heat_renal_marker_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_renal_metabolic_marker_summary.csv")
heat_crrt_window_path <- file.path(repo_root, "output", "final", "descriptive", "heat_related_ohca_crrt_window_summary.csv")
time_sensitivity_path <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript", "manuscript_dlnm_time_adjustment_sensitivity.csv")
adverse_models_path <- file.path(repo_root, "output", "final", "ohca_outcomes", "ohca_heat_adverse_outcome_models.csv")
continuous_models_path <- file.path(repo_root, "output", "final", "ohca_outcomes", "ohca_heat_continuous_outcome_models.csv")
pollution_binary_models_path <- file.path(repo_root, "output", "final", "ohca_outcomes", "ohca_pollution_12m_binary_outcome_models.csv")
pollution_continuous_models_path <- file.path(repo_root, "output", "final", "ohca_outcomes", "ohca_pollution_12m_continuous_outcome_models.csv")
adverse_rates_path <- file.path(repo_root, "output", "final", "ohca_outcomes", "ohca_heat_adverse_outcome_rates.csv")
denominator_audit_path <- file.path(repo_root, "output", "final", "quality_checks", "model_denominator_audit.csv")
icu_timing_path <- file.path(repo_root, "output", "final", "quality_checks", "ohca_admission_to_icu_timing_summary.csv")
icu_timing_bins_path <- file.path(repo_root, "output", "final", "quality_checks", "ohca_admission_to_icu_timing_bins.csv")
care_pathway_path <- file.path(repo_root, "output", "final", "quality_checks", "ohca_pre_icu_care_pathway_summary.csv")

results <- read.csv(results_path, stringsAsFactors = FALSE)
results$site_name <- site_name
results <- results[, c("site_name", setdiff(names(results), "site_name"))]

write.csv(results, file.path(output_dir, paste0(site_name, "_dlnm_site_estimates.csv")), row.names = FALSE)

if (file.exists(curves_path)) {
  curves <- read.csv(curves_path, stringsAsFactors = FALSE)
  curves$site_name <- site_name
  curves <- curves[, c("site_name", setdiff(names(curves), "site_name"))]
  write.csv(curves, file.path(output_dir, paste0(site_name, "_dlnm_curves.csv")), row.names = FALSE)
}

if (file.exists(reduced_coef_path)) {
  reduced_coef <- read.csv(reduced_coef_path, stringsAsFactors = FALSE)
  reduced_coef$site_name <- site_name
  reduced_coef <- reduced_coef[, c("site_name", setdiff(names(reduced_coef), "site_name"))]
  write.csv(reduced_coef, file.path(output_dir, paste0(site_name, "_dlnm_reduced_coefficients.csv")), row.names = FALSE)
}

if (file.exists(reduced_vcov_path)) {
  reduced_vcov <- read.csv(reduced_vcov_path, stringsAsFactors = FALSE)
  reduced_vcov$site_name <- site_name
  reduced_vcov <- reduced_vcov[, c("site_name", setdiff(names(reduced_vcov), "site_name"))]
  write.csv(reduced_vcov, file.path(output_dir, paste0(site_name, "_dlnm_reduced_vcov.csv")), row.names = FALSE)
}

if (file.exists(table1_path)) {
  table1 <- read.csv(table1_path, stringsAsFactors = FALSE)
  table1$site_name <- site_name
  write.csv(table1, file.path(output_dir, paste0(site_name, "_table1.csv")), row.names = FALSE)
}

if (file.exists(outcomes_path)) {
  outcomes <- read.csv(outcomes_path, stringsAsFactors = FALSE)
  outcomes$site_name <- site_name
  write.csv(outcomes, file.path(output_dir, paste0(site_name, "_outcomes.csv")), row.names = FALSE)
}

if (file.exists(heat_table2_path)) {
  heat_table2 <- read.csv(heat_table2_path, stringsAsFactors = FALSE)
  heat_table2$site_name <- site_name
  write.csv(heat_table2, file.path(output_dir, paste0(site_name, "_heat_related_vs_non_heat_related_table.csv")), row.names = FALSE)
}

for (item in list(
  list(path = heat_table2_all_path, suffix = "heat_related_vs_non_heat_related_table_all_definitions"),
  list(path = heat90_table2_path, suffix = "heat90_vs_non_heat90_table"),
  list(path = heat_thresholds_path, suffix = "heat_related_ohca_thresholds"),
  list(path = heat_summary_path, suffix = "heat_related_vs_non_heat_related_outcomes"),
  list(path = heat_discharge_path, suffix = "heat_related_vs_non_heat_related_discharge_categories"),
  list(path = heat_hourly_vitals_path, suffix = "heat_related_hourly_vital_trajectories"),
  list(path = heat_hourly_labs_path, suffix = "heat_related_hourly_lab_trajectories"),
  list(path = heat_hourly_support_path, suffix = "heat_related_hourly_support_trajectories"),
  list(path = heat_hourly_vitals_smoothed_path, suffix = "heat_related_hourly_vital_trajectories_smoothed"),
  list(path = heat_hourly_labs_smoothed_path, suffix = "heat_related_hourly_lab_trajectories_smoothed"),
  list(path = heat_hourly_support_smoothed_path, suffix = "heat_related_hourly_support_trajectories_smoothed"),
  list(path = heat_hourly_cumulative_path, suffix = "heat_related_hourly_cumulative_incidence"),
  list(path = heat_renal_marker_path, suffix = "heat_related_renal_metabolic_marker_summary"),
  list(path = heat_crrt_window_path, suffix = "heat_related_crrt_window_summary")
)) {
  if (file.exists(item$path)) {
    dat <- read.csv(item$path, stringsAsFactors = FALSE)
    dat$site_name <- site_name
    write.csv(dat, file.path(output_dir, paste0(site_name, "_", item$suffix, ".csv")), row.names = FALSE)
  }
}

if (file.exists(time_sensitivity_path)) {
  time_sensitivity <- read.csv(time_sensitivity_path, stringsAsFactors = FALSE)
  time_sensitivity$site_name <- site_name
  write.csv(time_sensitivity, file.path(output_dir, paste0(site_name, "_dlnm_time_sensitivity.csv")), row.names = FALSE)
}

if (file.exists(adverse_models_path)) {
  adverse_models <- read.csv(adverse_models_path, stringsAsFactors = FALSE)
  adverse_models$site_name <- site_name
  write.csv(adverse_models, file.path(output_dir, paste0(site_name, "_adverse_outcome_models.csv")), row.names = FALSE)
}

if (file.exists(continuous_models_path)) {
  continuous_models <- read.csv(continuous_models_path, stringsAsFactors = FALSE)
  continuous_models$site_name <- site_name
  write.csv(continuous_models, file.path(output_dir, paste0(site_name, "_continuous_outcome_models.csv")), row.names = FALSE)
}

if (file.exists(pollution_binary_models_path)) {
  pollution_binary_models <- read.csv(pollution_binary_models_path, stringsAsFactors = FALSE)
  pollution_binary_models$site_name <- site_name
  write.csv(pollution_binary_models, file.path(output_dir, paste0(site_name, "_pollution_12m_binary_outcome_models.csv")), row.names = FALSE)
}

if (file.exists(pollution_continuous_models_path)) {
  pollution_continuous_models <- read.csv(pollution_continuous_models_path, stringsAsFactors = FALSE)
  pollution_continuous_models$site_name <- site_name
  write.csv(pollution_continuous_models, file.path(output_dir, paste0(site_name, "_pollution_12m_continuous_outcome_models.csv")), row.names = FALSE)
}

if (file.exists(adverse_rates_path)) {
  adverse_rates <- read.csv(adverse_rates_path, stringsAsFactors = FALSE)
  adverse_rates$site_name <- site_name
  write.csv(adverse_rates, file.path(output_dir, paste0(site_name, "_adverse_outcome_rates.csv")), row.names = FALSE)
}

for (item in list(
  list(path = denominator_audit_path, suffix = "denominator_audit"),
  list(path = icu_timing_path, suffix = "icu_timing_summary"),
  list(path = icu_timing_bins_path, suffix = "icu_timing_bins"),
  list(path = care_pathway_path, suffix = "care_pathway_summary")
)) {
  if (file.exists(item$path)) {
    dat <- read.csv(item$path, stringsAsFactors = FALSE)
    dat$site_name <- site_name
    write.csv(dat, file.path(output_dir, paste0(site_name, "_", item$suffix, ".csv")), row.names = FALSE)
  }
}

if (dir.exists(figure_source_dir)) {
  old_site_figures <- list.files(
    figure_export_dir,
    pattern = paste0("^", site_name, "_.*\\.png$"),
    full.names = TRUE
  )
  if (length(old_site_figures) > 0) unlink(old_site_figures)

  figure_files <- list.files(
    figure_source_dir,
    pattern = "\\.png$",
    full.names = TRUE
  )
  for (figure_path in figure_files) {
    dest <- file.path(figure_export_dir, paste0(site_name, "_", basename(figure_path)))
    file.copy(figure_path, dest, overwrite = TRUE)
  }
}

message("Wrote federated site export files to ", output_dir)

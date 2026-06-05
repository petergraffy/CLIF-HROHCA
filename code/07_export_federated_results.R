#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    ofiles <- vapply(sys.frames(), function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile, character(1))
    ofiles <- stats::na.omit(ofiles)
    if (length(ofiles) > 0) return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = TRUE))
    stop("Could not determine script path. Run with Rscript.")
  }
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))

config <- load_project_config(repo_root)
site_name <- validate_site_name(config, repo_root)

output_dir <- file.path(repo_root, "output", "final", "federated_exports")
figure_export_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_export_dir, recursive = TRUE, showWarnings = FALSE)

unlink(list.files(output_dir, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE))
unlink(file.path(output_dir, ".DS_Store"))
unlink(list.files(figure_export_dir, full.names = TRUE, recursive = FALSE))
dir.create(figure_export_dir, recursive = TRUE, showWarnings = FALSE)

copy_csv <- function(path, suffix) {
  if (!file.exists(path)) return(FALSE)
  dat <- read.csv(path, stringsAsFactors = FALSE)
  dat$site_name <- site_name
  dat <- dat[, c("site_name", setdiff(names(dat), "site_name"))]
  write.csv(dat, file.path(output_dir, paste0(site_name, "_", suffix, ".csv")), row.names = FALSE)
  TRUE
}

copy_figure <- function(path) {
  if (!file.exists(path)) return(FALSE)
  file.copy(path, file.path(figure_export_dir, paste0(site_name, "_", basename(path))), overwrite = TRUE)
}

manuscript_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript")
lag_diag_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "lag_diagnostics")
rate_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "icu_admission_rate")
case_crossover_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "case_crossover")
phenotype_dir <- file.path(repo_root, "output", "final", "ohca_icu_phenotypes")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")

exports <- list(
  list(file.path(manuscript_dir, "manuscript_dlnm_results.csv"), "dlnm_site_estimates"),
  list(file.path(manuscript_dir, "manuscript_dlnm_curves.csv"), "dlnm_curves"),
  list(file.path(manuscript_dir, "manuscript_dlnm_lag_summaries.csv"), "dlnm_lag_summaries"),
  list(file.path(manuscript_dir, "manuscript_dlnm_lag_specific_summaries.csv"), "dlnm_lag_specific_summaries"),
  list(file.path(manuscript_dir, "manuscript_dlnm_reduced_coefficients.csv"), "dlnm_reduced_coefficients"),
  list(file.path(manuscript_dir, "manuscript_dlnm_reduced_vcov.csv"), "dlnm_reduced_vcov"),
  list(file.path(manuscript_dir, "manuscript_dlnm_time_adjustment_sensitivity.csv"), "dlnm_time_sensitivity"),
  list(file.path(manuscript_dir, "manuscript_dlnm_time_adjustment_lag_summaries.csv"), "dlnm_time_sensitivity_lag_summaries"),
  list(file.path(manuscript_dir, "manuscript_dlnm_time_adjustment_lag_specific_summaries.csv"), "dlnm_time_sensitivity_lag_specific_summaries"),
  list(file.path(lag_diag_dir, "lag30_diagnostic_summary.csv"), "lag30_diagnostic_summary"),
  list(file.path(lag_diag_dir, "lag30_temperature_lag_rr_surface.csv"), "lag30_temperature_lag_rr_surface"),
  list(file.path(lag_diag_dir, "lag30_hot_temperature_lag_specific_rr.csv"), "lag30_hot_temperature_lag_specific_rr"),
  list(file.path(lag_diag_dir, "lag30_hot_temperature_cumulative_rr_by_lag.csv"), "lag30_hot_temperature_cumulative_rr_by_lag"),
  list(file.path(lag_diag_dir, "lag30_temperature_distribution_source.csv"), "lag30_temperature_distribution"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_dlnm_results.csv"), "ohca_icu_admission_rate_dlnm_site_estimates"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_dlnm_curves.csv"), "ohca_icu_admission_rate_dlnm_curves"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_dlnm_lag_summaries.csv"), "ohca_icu_admission_rate_dlnm_lag_summaries"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_dlnm_lag_specific_summaries.csv"), "ohca_icu_admission_rate_dlnm_lag_specific_summaries"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_dlnm_reduced_coefficients.csv"), "ohca_icu_admission_rate_dlnm_reduced_coefficients"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_dlnm_reduced_vcov.csv"), "ohca_icu_admission_rate_dlnm_reduced_vcov"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_denominator_summary.csv"), "ohca_icu_admission_rate_denominator_summary"),
  list(file.path(rate_dir, "ohca_icu_admission_rate_daily_timeseries.csv"), "ohca_icu_admission_rate_daily_timeseries"),
  list(file.path(case_crossover_dir, "case_crossover_dlnm_results.csv"), "case_crossover_dlnm_site_estimates"),
  list(file.path(case_crossover_dir, "case_crossover_dlnm_curves.csv"), "case_crossover_dlnm_curves"),
  list(file.path(case_crossover_dir, "case_crossover_dlnm_lag_summaries.csv"), "case_crossover_dlnm_lag_summaries"),
  list(file.path(case_crossover_dir, "case_crossover_dlnm_lag_specific_summaries.csv"), "case_crossover_dlnm_lag_specific_summaries"),
  list(file.path(case_crossover_dir, "case_crossover_dlnm_reduced_coefficients.csv"), "case_crossover_dlnm_reduced_coefficients"),
  list(file.path(case_crossover_dir, "case_crossover_dlnm_reduced_vcov.csv"), "case_crossover_dlnm_reduced_vcov"),
  list(file.path(case_crossover_dir, "case_crossover_referent_set_summary.csv"), "case_crossover_referent_set_summary"),
  list(file.path(phenotype_dir, "ohca_icu_72h_consort_flow.csv"), "ohca_icu_72h_consort_flow"),
  list(file.path(phenotype_dir, "ohca_icu_72h_table1.csv"), "ohca_icu_72h_table1"),
  list(file.path(phenotype_dir, "ohca_icu_72h_table2_by_phenotype.csv"), "ohca_icu_72h_table2_by_phenotype"),
  list(file.path(phenotype_dir, "ohca_icu_72h_gcs_landmark_by_phenotype.csv"), "ohca_icu_72h_gcs_landmark_by_phenotype"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_summary.csv"), "ohca_icu_72h_phenotype_summary"),
  list(file.path(phenotype_dir, "ohca_icu_72h_ohca_mechanism_summary.csv"), "ohca_icu_72h_ohca_mechanism_summary"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_evidence_summary.csv"), "ohca_icu_72h_phenotype_evidence_summary"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_definitions.csv"), "ohca_icu_72h_phenotype_definitions"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_model.csv"), "ohca_icu_72h_phenotype_assignment_model"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_temperature_curve.csv"), "ohca_icu_72h_phenotype_assignment_temperature_curve"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_coefficients.csv"), "ohca_icu_72h_phenotype_assignment_coefficients"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_vcov.csv"), "ohca_icu_72h_phenotype_assignment_vcov"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted.csv"), "ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted.csv"), "ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted.csv"), "ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted"),
  list(file.path(phenotype_dir, "ohca_icu_72h_phenotype_assignment_vcov_mechanism_adjusted.csv"), "ohca_icu_72h_phenotype_assignment_vcov_mechanism_adjusted"),
  list(file.path(phenotype_dir, "ohca_icu_competing_risk_awake_extubated_72h_summary.csv"), "ohca_icu_competing_risk_awake_extubated_72h_summary"),
  list(file.path(phenotype_dir, "ohca_icu_competing_risk_death_source_summary.csv"), "ohca_icu_competing_risk_death_source_summary"),
  list(file.path(phenotype_dir, "ohca_icu_competing_risk_awake_extubated_72h_ohca_mechanism_summary.csv"), "ohca_icu_competing_risk_awake_extubated_72h_ohca_mechanism_summary"),
  list(file.path(phenotype_dir, "ohca_icu_competing_risk_awake_extubated_72h_fine_gray_models.csv"), "ohca_icu_competing_risk_awake_extubated_72h_fine_gray_models"),
  list(file.path(phenotype_dir, "ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv"), "ohca_icu_competing_risk_awake_extubated_72h_coefficients"),
  list(file.path(phenotype_dir, "ohca_icu_competing_risk_awake_extubated_72h_vcov.csv"), "ohca_icu_competing_risk_awake_extubated_72h_vcov")
)

invisible(lapply(exports, function(item) copy_csv(item[[1]], item[[2]])))

figure_paths <- file.path(figure_dir, c(
  "figure_primary_dlnm_curve.png",
  "figure_regional_dlnm_mrt_curves.png",
  "figure_all_year_dlnm_mrt_median_comparison.png",
  "figure_lag30_temperature_distribution.png",
  "figure_lag30_temperature_lag_rr_surface.png",
  "figure_lag30_hot_temperature_lag_specific_rr.png",
  "figure_lag30_hot_temperature_cumulative_rr_by_lag.png",
  "figure_ohca_icu_admission_rate_temperature_timeseries.png",
  "figure_ohca_icu_72h_phenotype_counts.png",
  "figure_ohca_icu_72h_phenotype_temperature_curves.png"
))
invisible(lapply(figure_paths, copy_figure))

message("Wrote DLNM and 72-hour phenotype federated exports to ", output_dir)

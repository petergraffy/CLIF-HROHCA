#!/usr/bin/env Rscript

# Michigan-specific recovery helper.
#
# Purpose:
# 1. Scan output/final recursively for the aggregate CSV files needed for
#    federated pooling.
# 2. Write a manifest of found, missing, ambiguous, and unsafe files.
# 3. Copy only safe aggregate CSVs into output/final/federated_exports with a
#    Michigan_ prefix and a site_name column.
#
# This script intentionally does not copy any file with patient/row-level
# language in the filename, and it skips any candidate with identifier-like
# column names.

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) > 0) {
    return(normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE))
  }
  ofiles <- vapply(
    sys.frames(),
    function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile,
    character(1)
  )
  ofiles <- stats::na.omit(ofiles)
  if (length(ofiles) > 0) {
    return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = TRUE))
  }
  stop("Could not determine script path. Run with Rscript.", call. = FALSE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
site_name <- Sys.getenv("MICHIGAN_SITE_NAME", unset = "Michigan")

final_dir <- file.path(repo_root, "output", "final")
export_dir <- file.path(final_dir, "federated_exports")
diagnostic_dir <- file.path(final_dir, "federated_export_diagnostics")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

expected_exports <- c(
  "case_crossover_dlnm_contrast_summaries.csv",
  "case_crossover_dlnm_curves.csv",
  "case_crossover_dlnm_lag_specific_summaries.csv",
  "case_crossover_dlnm_lag_summaries.csv",
  "case_crossover_dlnm_lag_temperature_surface.csv",
  "case_crossover_dlnm_reduced_coefficients.csv",
  "case_crossover_dlnm_reduced_vcov.csv",
  "case_crossover_dlnm_site_estimates.csv",
  "case_crossover_lag0_coefficients.csv",
  "case_crossover_lag0_contrast_summaries.csv",
  "case_crossover_lag0_curves.csv",
  "case_crossover_lag0_site_estimates.csv",
  "case_crossover_lag0_vcov.csv",
  "case_crossover_lag30_diagnostic_summary.csv",
  "case_crossover_lag30_hot_temperature_cumulative_rr_by_lag.csv",
  "case_crossover_lag30_hot_temperature_lag_specific_rr.csv",
  "case_crossover_referent_set_summary.csv",
  "dlnm_contrast_summaries.csv",
  "dlnm_curves.csv",
  "dlnm_lag_specific_summaries.csv",
  "dlnm_lag_summaries.csv",
  "dlnm_lag_temperature_surface.csv",
  "dlnm_reduced_coefficients.csv",
  "dlnm_reduced_vcov.csv",
  "dlnm_site_estimates.csv",
  "dlnm_time_sensitivity_contrast_summaries.csv",
  "dlnm_time_sensitivity_lag_specific_summaries.csv",
  "dlnm_time_sensitivity_lag_summaries.csv",
  "dlnm_time_sensitivity_lag_temperature_surface.csv",
  "dlnm_time_sensitivity.csv",
  "lag30_diagnostic_summary.csv",
  "lag30_hot_temperature_cumulative_rr_by_lag.csv",
  "lag30_hot_temperature_lag_specific_rr.csv",
  "lag30_temperature_distribution.csv",
  "lag30_temperature_lag_rr_surface.csv",
  "ohca_diagnosis_code_audit.csv",
  "ohca_ed_only_death_never_icu_los_distribution.csv",
  "ohca_ed_only_death_never_icu_pathway_audit.csv",
  "ohca_ed_only_death_never_icu_poa_audit.csv",
  "ohca_ed_only_death_never_icu_summary.csv",
  "ohca_ed_only_death_never_icu_timing_summary.csv",
  "ohca_ed_only_never_icu_discharge_category_audit.csv",
  "ohca_icu_72h_admission_temperature_density.csv",
  "ohca_icu_72h_admission_temperature_distribution_summary.csv",
  "ohca_icu_72h_consort_flow.csv",
  "ohca_icu_72h_gcs_landmark_by_phenotype.csv",
  "ohca_icu_72h_ohca_mechanism_summary.csv",
  "ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_assignment_coefficients.csv",
  "ohca_icu_72h_phenotype_assignment_lag_sensitivity_coefficients.csv",
  "ohca_icu_72h_phenotype_assignment_lag_sensitivity_model.csv",
  "ohca_icu_72h_phenotype_assignment_lag_sensitivity_temperature_curve.csv",
  "ohca_icu_72h_phenotype_assignment_lag_sensitivity_vcov.csv",
  "ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_assignment_model.csv",
  "ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_assignment_temperature_curve.csv",
  "ohca_icu_72h_phenotype_assignment_vcov_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_assignment_vcov.csv",
  "ohca_icu_72h_phenotype_definitions.csv",
  "ohca_icu_72h_phenotype_evidence_summary.csv",
  "ohca_icu_72h_phenotype_summary.csv",
  "ohca_icu_72h_table1.csv",
  "ohca_icu_72h_table2_by_phenotype.csv",
  "ohca_icu_admission_rate_daily_timeseries.csv",
  "ohca_icu_admission_rate_denominator_summary.csv",
  "ohca_icu_admission_rate_dlnm_contrast_summaries.csv",
  "ohca_icu_admission_rate_dlnm_curves.csv",
  "ohca_icu_admission_rate_dlnm_lag_specific_summaries.csv",
  "ohca_icu_admission_rate_dlnm_lag_summaries.csv",
  "ohca_icu_admission_rate_dlnm_lag_temperature_surface.csv",
  "ohca_icu_admission_rate_dlnm_reduced_coefficients.csv",
  "ohca_icu_admission_rate_dlnm_reduced_vcov.csv",
  "ohca_icu_admission_rate_dlnm_site_estimates.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_fine_gray_models.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_coefficients.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_fine_gray_models.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_vcov.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_ohca_mechanism_summary.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_summary.csv",
  "ohca_icu_competing_risk_awake_extubated_72h_vcov.csv",
  "ohca_icu_competing_risk_death_source_summary.csv",
  "ohca_icu_imv12_discharge_outcome_coefficients.csv",
  "ohca_icu_imv12_discharge_outcome_flow.csv",
  "ohca_icu_imv12_discharge_outcome_model.csv",
  "ohca_icu_imv12_discharge_outcome_summary.csv",
  "ohca_icu_imv12_discharge_outcome_table_by_outcome.csv",
  "ohca_icu_imv12_discharge_outcome_vcov.csv",
  "ohca_icu_imv24_time_to_event_cif_curves.csv",
  "ohca_icu_imv24_time_to_event_coefficients.csv",
  "ohca_icu_imv24_time_to_event_evidence_summary.csv",
  "ohca_icu_imv24_time_to_event_fine_gray_models.csv",
  "ohca_icu_imv24_time_to_event_landmark_flow.csv",
  "ohca_icu_imv24_time_to_event_summary.csv",
  "ohca_icu_imv24_time_to_event_table_by_outcome.csv",
  "ohca_icu_imv24_time_to_event_vcov.csv",
  "rate_lag30_diagnostic_summary.csv",
  "rate_lag30_hot_temperature_cumulative_rr_by_lag.csv",
  "rate_lag30_hot_temperature_lag_specific_rr.csv",
  "rate_lag30_temperature_distribution.csv",
  "rate_lag30_temperature_lag_rr_surface.csv"
)

source_aliases <- list(
  dlnm_site_estimates.csv = "manuscript_dlnm_results.csv",
  dlnm_curves.csv = "manuscript_dlnm_curves.csv",
  dlnm_lag_summaries.csv = "manuscript_dlnm_lag_summaries.csv",
  dlnm_lag_specific_summaries.csv = "manuscript_dlnm_lag_specific_summaries.csv",
  dlnm_lag_temperature_surface.csv = "manuscript_dlnm_lag_temperature_surface.csv",
  dlnm_contrast_summaries.csv = "manuscript_dlnm_contrast_summaries.csv",
  dlnm_reduced_coefficients.csv = "manuscript_dlnm_reduced_coefficients.csv",
  dlnm_reduced_vcov.csv = "manuscript_dlnm_reduced_vcov.csv",
  dlnm_time_sensitivity.csv = "manuscript_dlnm_time_adjustment_sensitivity.csv",
  dlnm_time_sensitivity_lag_summaries.csv = "manuscript_dlnm_time_adjustment_lag_summaries.csv",
  dlnm_time_sensitivity_lag_specific_summaries.csv = "manuscript_dlnm_time_adjustment_lag_specific_summaries.csv",
  dlnm_time_sensitivity_lag_temperature_surface.csv = "manuscript_dlnm_time_adjustment_lag_temperature_surface.csv",
  dlnm_time_sensitivity_contrast_summaries.csv = "manuscript_dlnm_time_adjustment_contrast_summaries.csv",
  lag30_temperature_distribution.csv = "lag30_temperature_distribution_source.csv",
  rate_lag30_temperature_distribution.csv = "rate_lag30_temperature_distribution_source.csv",
  case_crossover_dlnm_site_estimates.csv = "case_crossover_dlnm_results.csv",
  case_crossover_lag0_site_estimates.csv = "case_crossover_lag0_results.csv",
  ohca_icu_admission_rate_dlnm_site_estimates.csv = "ohca_icu_admission_rate_dlnm_results.csv"
)

make_aliases <- function(expected_file) {
  aliases <- c(expected_file, sub("[.]csv$", "_source.csv", expected_file))
  extra <- source_aliases[[expected_file]]
  if (!is.null(extra)) aliases <- c(aliases, extra)
  unique(aliases)
}

unsafe_file_pattern <- paste(
  c(
    "patient_audit",
    "patient_level",
    "row_level",
    "row_audit",
    "encounter_level",
    "hospitalization_level",
    "mrn"
  ),
  collapse = "|"
)

unsafe_column_names <- c(
  "patient_id",
  "person_id",
  "hospitalization_id",
  "encounter_id",
  "mrn",
  "medical_record_number",
  "birth_date",
  "date_of_birth",
  "dob",
  "first_name",
  "last_name"
)

is_unsafe_filename <- function(path) {
  grepl(unsafe_file_pattern, tolower(basename(path)))
}

is_excluded_search_path <- function(path) {
  lower_path <- tolower(normalizePath(path, winslash = "/", mustWork = FALSE))
  grepl("/federated_exports", lower_path, fixed = TRUE) ||
    grepl("/federated_export_diagnostics", lower_path, fixed = TRUE) ||
    grepl("/federated_pooled", lower_path, fixed = TRUE)
}

read_header <- function(path) {
  tryCatch(
    names(utils::read.csv(path, nrows = 0, check.names = FALSE)),
    error = function(e) character()
  )
}

identifier_columns <- function(path) {
  lower_header <- tolower(read_header(path))
  lower_header[lower_header %in% unsafe_column_names]
}

all_csv <- list.files(final_dir, pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
all_csv <- all_csv[!vapply(all_csv, is_excluded_search_path, logical(1))]

find_candidates <- function(expected_file) {
  aliases <- make_aliases(expected_file)
  candidate_paths <- all_csv[basename(all_csv) %in% aliases]
  candidate_paths <- candidate_paths[!vapply(candidate_paths, is_unsafe_filename, logical(1))]
  candidate_paths
}

choose_candidate <- function(paths, expected_file) {
  if (length(paths) == 0) return(NA_character_)
  expected_basename <- expected_file
  exact <- paths[basename(paths) == expected_basename]
  if (length(exact) > 0) paths <- exact
  path_order <- order(file.info(paths)$mtime, decreasing = TRUE, na.last = TRUE)
  paths[[path_order[[1]]]]
}

manifest_rows <- lapply(expected_exports, function(expected_file) {
  candidates <- find_candidates(expected_file)
  chosen <- choose_candidate(candidates, expected_file)
  id_cols <- if (is.na(chosen)) character() else identifier_columns(chosen)
  status <- if (length(candidates) == 0) {
    "missing"
  } else if (length(id_cols) > 0) {
    "unsafe_identifier_columns"
  } else if (length(candidates) > 1) {
    "found_multiple_using_newest"
  } else {
    "found"
  }
  data.frame(
    expected_export = expected_file,
    status = status,
    chosen_source = ifelse(is.na(chosen), "", chosen),
    n_candidates = length(candidates),
    candidate_basenames = paste(unique(basename(candidates)), collapse = "; "),
    identifier_columns = paste(id_cols, collapse = "; "),
    stringsAsFactors = FALSE
  )
})
manifest <- do.call(rbind, manifest_rows)

write_with_site_name <- function(source_path, expected_file) {
  dat <- utils::read.csv(source_path, stringsAsFactors = FALSE, check.names = FALSE)
  id_cols <- tolower(names(dat))[tolower(names(dat)) %in% unsafe_column_names]
  if (length(id_cols) > 0) {
    stop("Refusing to copy ", source_path, "; identifier-like columns found: ", paste(id_cols, collapse = ", "))
  }
  dat$site_name <- site_name
  dat <- dat[, c("site_name", setdiff(names(dat), "site_name")), drop = FALSE]
  target <- file.path(export_dir, paste0(site_name, "_", expected_file))
  utils::write.csv(dat, target, row.names = FALSE)
  target
}

copyable <- manifest$status %in% c("found", "found_multiple_using_newest")
targets <- character()
if (!dry_run) {
  targets <- mapply(
    write_with_site_name,
    source_path = manifest$chosen_source[copyable],
    expected_file = manifest$expected_export[copyable],
    USE.NAMES = FALSE
  )
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
manifest_path <- file.path(diagnostic_dir, paste0(site_name, "_federated_export_recovery_manifest_", timestamp, ".csv"))
missing_path <- file.path(diagnostic_dir, paste0(site_name, "_missing_federated_exports_", timestamp, ".csv"))
unsafe_path <- file.path(diagnostic_dir, paste0(site_name, "_unsafe_federated_export_candidates_", timestamp, ".csv"))

utils::write.csv(manifest, manifest_path, row.names = FALSE)
utils::write.csv(manifest[manifest$status == "missing", , drop = FALSE], missing_path, row.names = FALSE)
utils::write.csv(manifest[grepl("^unsafe", manifest$status), , drop = FALSE], unsafe_path, row.names = FALSE)

message("Michigan federated export recovery complete.")
message("Mode: ", ifelse(dry_run, "dry run; no files copied", "copied safe files"))
message("Expected CSV exports: ", length(expected_exports))
message("Found safe exports: ", sum(copyable))
message("Missing exports: ", sum(manifest$status == "missing"))
message("Unsafe candidates skipped: ", sum(grepl("^unsafe", manifest$status)))
message("Federated export folder: ", export_dir)
message("Manifest: ", manifest_path)
message("Missing manifest: ", missing_path)
message("Unsafe manifest: ", unsafe_path)

if (!dry_run && length(targets) > 0) {
  message("Wrote ", length(targets), " Michigan-prefixed aggregate CSV(s).")
}

if (sum(manifest$status == "missing") > 0) {
  message("Review the missing manifest before upload; this folder may still be incomplete.")
}

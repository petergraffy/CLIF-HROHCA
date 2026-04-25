#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) stop("Could not determine script path from commandArgs().")
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)

scripts <- c(
  "00_install_or_restore_packages.R",
  "01_build_ohca_cohort.R",
  "02_build_icu_exposure_series.R",
  "03_descriptive_tables.R",
  "04_dlnm_primary_and_sensitivity.R",
  "05_heat_related_vs_non_heat_related_table.R",
  "09_supplementary_ohca_outcome_models.R",
  "08_quality_checks.R",
  "06_manuscript_tables_figures.R",
  "07_export_federated_results.R"
)

for (script in scripts) {
  path <- file.path(repo_root, "code", script)
  message("\n=== Running ", script, " ===")
  status <- system2("Rscript", shQuote(path))
  if (is.null(status)) status <- 0
  if (status != 0) stop("Script failed: ", script, call. = FALSE)
}

message("\nSite analysis complete. Share only files in output/final/federated_exports with the coordinating center.")

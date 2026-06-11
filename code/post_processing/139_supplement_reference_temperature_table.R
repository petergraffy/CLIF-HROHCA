#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(file_name) {
  path <- file.path(pooled_dir, file_name)
  if (!file.exists(path)) stop("Missing required pooled file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

c_to_f <- function(x) x * 9 / 5 + 32
fmt_temp <- function(celsius) {
  ifelse(
    is.na(celsius),
    NA_character_,
    sprintf("%.1f C (%.1f F)", celsius, c_to_f(celsius))
  )
}

pretty_family <- function(x) {
  dplyr::recode(
    x,
    case_crossover_dlnm = "Time-stratified case-crossover DLNM",
    count_dlnm = "Daily count DLNM",
    rate_dlnm = "Daily OHCA/ICU-admission rate DLNM",
    .default = x
  )
}

pretty_model <- function(x) {
  x |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_replace("^primary humidity adjusted$", "primary humidity adjusted") |>
    stringr::str_replace("^sensitivity humidity pollution adjusted$", "humidity + pollution adjusted") |>
    stringr::str_replace("^sensitivity mrt reference$", "MRT reference") |>
    stringr::str_replace("^stratified humidity adjusted$", "stratified humidity adjusted") |>
    stringr::str_replace("^stratified mrt reference$", "stratified MRT reference") |>
    stringr::str_replace("^rate median humidity plus pollution$", "rate median humidity + pollution adjusted") |>
    stringr::str_replace("^rate mrt humidity plus pollution$", "rate MRT humidity + pollution adjusted") |>
    stringr::str_replace("^rate stratified humidity adjusted$", "rate stratified humidity adjusted") |>
    stringr::str_replace("^rate stratified mrt reference$", "rate stratified MRT reference")
}

dlnm_estimates <- read_required("all_site_dlnm_estimates.csv")

reference_table_full <- dlnm_estimates |>
  dplyr::filter(.data$estimable %in% c(TRUE, "TRUE", NA)) |>
  dplyr::distinct(
    .data$export_family,
    .data$site_name,
    .data$stratum,
    .data$model,
    .data$reference_type,
    .data$reference_temp_c,
    .data$hot_temp_c
  ) |>
  dplyr::mutate(
    analysis = pretty_family(.data$export_family),
    model_label = pretty_model(.data$model),
    reference_label = dplyr::if_else(.data$reference_type == "mrt", "MRT", "Median"),
    reference_temp_f = c_to_f(.data$reference_temp_c),
    hot_temp_f = c_to_f(.data$hot_temp_c),
    reference_temperature = fmt_temp(.data$reference_temp_c),
    percentile_95_temperature = fmt_temp(.data$hot_temp_c)
  ) |>
  dplyr::arrange(.data$analysis, .data$stratum, .data$model_label, .data$reference_label, .data$site_name) |>
  dplyr::select(
    analysis,
    site_name,
    stratum,
    model = "model_label",
    reference = "reference_label",
    reference_temp_c,
    reference_temp_f,
    hot_temp_c = "hot_temp_c",
    hot_temp_f,
    reference_temperature,
    percentile_95_temperature
  )

reference_table_primary <- reference_table_full |>
  dplyr::filter(
    .data$stratum == "Overall",
    .data$reference == "Median",
    .data$model %in% c(
      "primary humidity adjusted",
      "rate median humidity + pollution adjusted"
    )
  ) |>
  dplyr::arrange(.data$analysis, .data$site_name)

full_csv <- file.path(table_dir, "supplement_table_dlnm_reference_and_95th_temperatures_full.csv")
primary_csv <- file.path(table_dir, "supplement_table_dlnm_reference_and_95th_temperatures_primary.csv")
primary_md <- file.path(table_dir, "supplement_table_dlnm_reference_and_95th_temperatures_primary.md")

readr::write_csv(reference_table_full, full_csv)
readr::write_csv(reference_table_primary, primary_csv)

md_cols <- c("analysis", "site_name", "stratum", "model", "reference", "reference_temperature", "percentile_95_temperature")
md_lines <- c(
  paste0("| ", paste(c("Analysis", "Site", "Stratum", "Model", "Reference", "Reference temperature", "95th percentile temperature"), collapse = " | "), " |"),
  paste0("|", paste(rep("---", length(md_cols)), collapse = "|"), "|")
)
for (i in seq_len(nrow(reference_table_primary))) {
  vals <- reference_table_primary[i, md_cols]
  vals <- vapply(vals, function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    stringr::str_replace_all(x, "\\|", "/")
  }, character(1))
  md_lines <- c(md_lines, paste0("| ", paste(vals, collapse = " | "), " |"))
}
writeLines(md_lines, primary_md)

message("Wrote DLNM reference-temperature supplement tables to ", table_dir)

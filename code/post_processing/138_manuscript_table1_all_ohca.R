#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "stringr", "tidyr")
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

parse_count <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(stringr::str_replace_all(stringr::str_extract(x, "^[0-9,]+"), ",", "")))
  out[is.na(out)] <- 0
  out
}

parse_median_iqr <- function(x) {
  x <- as.character(x)
  matches <- stringr::str_match(x, "^\\s*(-?[0-9.]+)\\s*\\[\\s*(-?[0-9.]+)\\s*,\\s*(-?[0-9.]+)\\s*\\]")
  tibble::tibble(
    median = suppressWarnings(as.numeric(matches[, 2])),
    q1 = suppressWarnings(as.numeric(matches[, 3])),
    q3 = suppressWarnings(as.numeric(matches[, 4]))
  )
}

format_n <- function(n) {
  format(round(n), big.mark = ",", scientific = FALSE)
}

format_n_pct <- function(n, denom) {
  if (!is.finite(denom) || denom <= 0) return(NA_character_)
  sprintf("%s (%.1f%%)", format_n(n), 100 * n / denom)
}

lookup_count <- function(x, name) {
  if (!name %in% names(x)) return(0)
  value <- x[[name]]
  if (length(value) == 0 || is.na(value)) 0 else value
}

format_site_median_range <- function(values) {
  parsed <- parse_median_iqr(values)
  if (all(is.na(parsed$median))) return(NA_character_)
  sprintf(
    "site medians %.1f-%.1f; site IQRs %.1f-%.1f to %.1f-%.1f",
    min(parsed$median, na.rm = TRUE),
    max(parsed$median, na.rm = TRUE),
    min(parsed$q1, na.rm = TRUE),
    max(parsed$q1, na.rm = TRUE),
    min(parsed$q3, na.rm = TRUE),
    max(parsed$q3, na.rm = TRUE)
  )
}

pretty_phenotype <- function(x) {
  dplyr::recode(
    x,
    alive_no_imv = "Alive/no IMV",
    regained_consciousness_extubated = "Regained consciousness/extubated",
    limited_brain_function = "Limited brain function",
    anoxic_brain_injury = "Anoxic brain injury",
    unclassified = "Unclassified",
    .default = x
  )
}

table1_site <- read_required("all_site_ohca_icu_72h_table1.csv")
table2_site <- read_required("all_site_ohca_icu_72h_table2_by_phenotype.csv")
phenotype_summary <- read_required("all_site_ohca_icu_72h_phenotype_summary.csv")
ed_death_summary <- read_required("all_site_ohca_ed_only_death_never_icu_summary.csv")

site_order <- unique(phenotype_summary$site_name)

site_denoms <- phenotype_summary |>
  dplyr::group_by(.data$site_name) |>
  dplyr::summarise(n = sum(.data$n, na.rm = TRUE), .groups = "drop")
site_denoms_vec <- stats::setNames(site_denoms$n, site_denoms$site_name)
overall_n <- sum(site_denoms$n, na.rm = TRUE)
phenotype_cols <- c("alive_no_imv", "regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury", "unclassified")

site_value_lookup <- function(dat) {
  vals <- stats::setNames(dat$value, dat$site_name)
  out <- vals[site_order]
  out[is.na(out)] <- ""
  out
}

cohort_row <- tibble::tibble(
  section = "Cohort",
  characteristic = "OHCA ICU cohort",
  level = "N",
  overall = format_n(overall_n)
)
for (site in site_order) cohort_row[[site]] <- format_n(site_denoms_vec[[site]])

ed_death_counts <- stats::setNames(
  ed_death_summary$n_ohca_dx_ed_first_only_death_never_icu_patients,
  ed_death_summary$site_name
)
ed_death_row <- tibble::tibble(
  section = "Cohort",
  characteristic = "OHCA died in ED, never ICU",
  level = "N patients",
  overall = format_n(sum(ed_death_counts, na.rm = TRUE))
)
for (site in site_order) ed_death_row[[site]] <- format_n(lookup_count(ed_death_counts, site))

continuous_keys <- table1_site |>
  dplyr::filter((is.na(.data$level) | .data$level == "") & .data$section != "Cohort") |>
  dplyr::distinct(.data$section, .data$characteristic, .data$level)

continuous_rows <- lapply(seq_len(nrow(continuous_keys)), function(i) {
  key <- continuous_keys[i, ]
  dat <- table1_site |>
    dplyr::filter(
      .data$section == key$section,
      .data$characteristic == key$characteristic,
      (is.na(.data$level) | .data$level == "")
    )
  row <- tibble::tibble(
    section = key$section,
    characteristic = key$characteristic,
    level = key$level,
    overall = format_site_median_range(dat$value)
  )
  site_vals <- site_value_lookup(dat)
  for (site in site_order) row[[site]] <- site_vals[[site]]
  row
}) |>
  dplyr::bind_rows()

categorical_sections <- c("Demographics", "OHCA mechanism", "Organ support", "Outcome", "Discharge outcomes")
categorical_keys <- table2_site |>
  dplyr::filter(.data$section %in% categorical_sections, !is.na(.data$level) & .data$level != "") |>
  dplyr::distinct(.data$section, .data$characteristic, .data$level)

categorical_rows <- lapply(seq_len(nrow(categorical_keys)), function(i) {
  key <- categorical_keys[i, ]
  dat <- table2_site |>
    dplyr::filter(
      .data$section == key$section,
      .data$characteristic == key$characteristic,
      .data$level == key$level
    )
  dat$site_count <- rowSums(
    as.data.frame(lapply(phenotype_cols, function(col) parse_count(dat[[col]]))),
    na.rm = TRUE
  )
  site_counts <- stats::setNames(dat$site_count, dat$site_name)
  row <- tibble::tibble(
    section = key$section,
    characteristic = key$characteristic,
    level = key$level,
    overall = format_n_pct(sum(site_counts, na.rm = TRUE), overall_n)
  )
  for (site in site_order) row[[site]] <- format_n_pct(lookup_count(site_counts, site), site_denoms_vec[[site]])
  row
}) |>
  dplyr::bind_rows()

phenotype_rows <- phenotype_summary |>
  dplyr::mutate(
    section = "72h phenotype",
    characteristic = "Assigned phenotype",
    level = pretty_phenotype(.data$phenotype)
  )

phenotype_levels <- c(
  "Alive/no IMV",
  "Regained consciousness/extubated",
  "Limited brain function",
  "Anoxic brain injury",
  "Unclassified"
)

phenotype_out <- lapply(phenotype_levels, function(level_value) {
  dat <- phenotype_rows |>
    dplyr::filter(.data$level == level_value)
  site_counts <- stats::setNames(dat$n, dat$site_name)
  row <- tibble::tibble(
    section = "72h phenotype",
    characteristic = "Assigned phenotype",
    level = level_value,
    overall = format_n_pct(sum(site_counts, na.rm = TRUE), overall_n)
  )
  for (site in site_order) row[[site]] <- format_n_pct(lookup_count(site_counts, site), site_denoms_vec[[site]])
  row
}) |>
  dplyr::bind_rows()

section_order <- c(
  "Cohort",
  "Demographics",
  "Admission exposure",
  "OHCA mechanism",
  "Organ support",
  "72h phenotype",
  "Outcome",
  "Discharge outcomes"
)

table1 <- dplyr::bind_rows(cohort_row, ed_death_row, continuous_rows, categorical_rows, phenotype_out) |>
  dplyr::filter(.data$section %in% section_order) |>
  dplyr::mutate(section = factor(.data$section, levels = section_order)) |>
  dplyr::arrange(.data$section, .data$characteristic, .data$level) |>
  dplyr::mutate(
    section = as.character(.data$section),
    level = tidyr::replace_na(.data$level, "")
  )

csv_path <- file.path(table_dir, "table1_all_ohca_icu_cohort.csv")
md_path <- file.path(table_dir, "table1_all_ohca_icu_cohort.md")
readr::write_csv(table1, csv_path)

md_cols <- c("section", "characteristic", "level", "overall", site_order)
md_lines <- c(
  paste0("| ", paste(c("Section", "Characteristic", "Level", "Overall", site_order), collapse = " | "), " |"),
  paste0("|", paste(rep("---", length(md_cols)), collapse = "|"), "|")
)
for (i in seq_len(nrow(table1))) {
  vals <- table1[i, md_cols]
  vals <- vapply(vals, function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    stringr::str_replace_all(x, "\\|", "/")
  }, character(1))
  md_lines <- c(md_lines, paste0("| ", paste(vals, collapse = " | "), " |"))
}
writeLines(md_lines, md_path)

message("Wrote full OHCA ICU Table 1 to ", csv_path)

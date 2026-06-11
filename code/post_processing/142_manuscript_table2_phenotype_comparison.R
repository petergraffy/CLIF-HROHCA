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

table2_path <- file.path(pooled_dir, "all_site_ohca_icu_72h_table2_by_phenotype.csv")
if (!file.exists(table2_path)) stop("Missing pooled Table 2 source: ", table2_path, call. = FALSE)

table2_site <- readr::read_csv(table2_path, show_col_types = FALSE)

phenotype_map <- c(
  alive_no_imv = "No IMV in first 72h",
  regained_consciousness_extubated = "Extubated by 72h",
  limited_brain_function = "On IMV at 72h",
  anoxic_brain_injury = "Death within 72h"
)
phenotype_cols <- names(phenotype_map)
phenotype_labels <- unname(phenotype_map)

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

format_n <- function(n) format(round(n), big.mark = ",", scientific = FALSE)

format_n_pct <- function(n, denom) {
  if (!is.finite(denom) || denom <= 0) return(NA_character_)
  sprintf("%s (%.1f%%)", format_n(n), 100 * n / denom)
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

cohort_size_rows <- table2_site |>
  dplyr::filter(.data$section == "Cohort", .data$characteristic == "Phenotype cohort size")

denoms <- vapply(phenotype_cols, function(col) {
  if (!col %in% names(cohort_size_rows)) return(0)
  sum(parse_count(cohort_size_rows[[col]]), na.rm = TRUE)
}, numeric(1))

summarise_row <- function(dat) {
  row <- tibble::tibble(
    section = dat$section[[1]],
    characteristic = dat$characteristic[[1]],
    level = tidyr::replace_na(dat$level[[1]], "")
  )

  for (col in phenotype_cols) {
    values <- dat[[col]]
    if (row$section == "Cohort" || (!is.na(row$level) && nzchar(row$level))) {
      n <- sum(parse_count(values), na.rm = TRUE)
      row[[phenotype_map[[col]]]] <- if (row$section == "Cohort") format_n(n) else format_n_pct(n, denoms[[col]])
    } else {
      row[[phenotype_map[[col]]]] <- format_site_median_range(values)
    }
  }

  row
}

group_keys <- table2_site |>
  dplyr::select("section", "characteristic", "level") |>
  dplyr::distinct()

table2 <- lapply(seq_len(nrow(group_keys)), function(i) {
  key <- group_keys[i, ]
  dat <- table2_site |>
    dplyr::filter(
      .data$section == key$section,
      .data$characteristic == key$characteristic,
      tidyr::replace_na(.data$level, "") == tidyr::replace_na(key$level, "")
    )
  summarise_row(dat)
}) |>
  dplyr::bind_rows()

section_order <- c(
  "Cohort",
  "Demographics",
  "Admission exposure",
  "OHCA mechanism",
  "Organ support",
  "Organ dysfunction",
  "72h phenotype evidence",
  "Awake signal components",
  "Neurologic injury",
  "Outcome",
  "Unclassified audit"
)

table2 <- table2 |>
  dplyr::mutate(
    section = factor(.data$section, levels = section_order),
    level = tidyr::replace_na(.data$level, "")
  ) |>
  dplyr::arrange(.data$section, .data$characteristic, .data$level) |>
  dplyr::mutate(section = as.character(.data$section)) |>
  dplyr::select("section", "characteristic", "level", dplyr::all_of(phenotype_labels))

table2_site_labeled <- table2_site |>
  dplyr::rename(
    "No IMV in first 72h" = "alive_no_imv",
    "Extubated by 72h" = "regained_consciousness_extubated",
    "On IMV at 72h" = "limited_brain_function",
    "Death within 72h" = "anoxic_brain_injury"
  )

csv_path <- file.path(table_dir, "table2_ohca_icu_72h_by_phenotype.csv")
md_path <- file.path(table_dir, "table2_ohca_icu_72h_by_phenotype.md")
source_path <- file.path(table_dir, "table2_ohca_icu_72h_by_phenotype_site_source.csv")

readr::write_csv(table2, csv_path)
readr::write_csv(table2_site_labeled, source_path)

md_cols <- names(table2)
md_lines <- c(
  paste0("| ", paste(c("Section", "Characteristic", "Level", phenotype_labels), collapse = " | "), " |"),
  paste0("|", paste(rep("---", length(md_cols)), collapse = "|"), "|")
)
for (i in seq_len(nrow(table2))) {
  vals <- table2[i, md_cols]
  vals <- vapply(vals, function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    stringr::str_replace_all(x, "\\|", "/")
  }, character(1))
  md_lines <- c(md_lines, paste0("| ", paste(vals, collapse = " | "), " |"))
}
writeLines(md_lines, md_path)

message("Wrote OHCA ICU 72h phenotype Table 2 to ", csv_path, " and ", md_path)

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
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(openxlsx)
  library(readr)
  library(tidyr)
})

definition_label <- c("heat95" = "Heat95", "heat90" = "Heat90")

suffix_for <- function(definition) {
  if (identical(definition, "heat95")) "" else paste0("_", definition)
}

format_p <- function(p) {
  out <- ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  out
}

format_ci <- function(est, lo, hi, digits = 1, suffix = "") {
  ifelse(
    is.na(est),
    "",
    paste0(
      sprintf(paste0("%.", digits, "f"), est), suffix,
      " (",
      sprintf(paste0("%.", digits, "f"), lo), ", ",
      sprintf(paste0("%.", digits, "f"), hi), suffix,
      ")"
    )
  )
}

format_est <- function(est, digits = 1, suffix = "") {
  ifelse(is.na(est), "", paste0(sprintf(paste0("%.", digits, "f"), est), suffix))
}

read_table1_source <- function(definition) {
  read.csv(
    file.path(table_dir, paste0("table1_pooled_hrohca_vs_non_hrohca", suffix_for(definition), "_source.csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    mutate(definition = definition)
}

table1_rows <- bind_rows(read_table1_source("heat95"), read_table1_source("heat90")) |>
  group_by(.data$definition) |>
  mutate(row_index = row_number()) |>
  ungroup() |>
  mutate(
    Characteristic = .data$characteristic,
    Section = .data$section,
    `HROHCA vs non-HROHCA` = if_else(
      .data$characteristic == "Patients, n",
      paste0(.data$heat_related_ohca, " HROHCA; ", .data$non_heat_related_ohca, " non-HROHCA"),
      paste0(.data$heat_related_ohca, " vs ", .data$non_heat_related_ohca)
    ),
    `p value` = .data$p_value
  ) |>
  select("row_index", "Section", "Characteristic", "definition", "HROHCA vs non-HROHCA", "p value") |>
  pivot_wider(
    id_cols = c("row_index", "Section", "Characteristic"),
    names_from = "definition",
    values_from = c("HROHCA vs non-HROHCA", "p value"),
    names_glue = "{definition}_{.value}"
  ) |>
  transmute(
    Section = .data$Section,
    Characteristic = .data$Characteristic,
    `Heat95 HROHCA vs non-HROHCA` = .data$`heat95_HROHCA vs non-HROHCA`,
    `Heat95 p value` = .data$`heat95_p value`,
    `Heat90 HROHCA vs non-HROHCA` = .data$`heat90_HROHCA vs non-HROHCA`,
    `Heat90 p value` = .data$`heat90_p value`
  ) |>
  mutate(across(everything(), ~tidyr::replace_na(as.character(.x), "")))

read_lab_or_vital <- function(kind, definition) {
  file <- file.path(
    figure_dir,
    if (kind == "Laboratory") {
      paste0("figure_hrohca_lab_trajectories", suffix_for(definition), "_source.csv")
    } else {
      paste0("figure_hrohca_vital_trajectory_differences_source", suffix_for(definition), ".csv")
    }
  )
  dat <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  label_col <- if (kind == "Laboratory") "lab" else "vital"

  dat |>
    mutate(
      domain = kind,
      variable_label = .data[[label_col]],
      abs_diff = abs(.data$approx_difference),
      significant = .data$approx_p_value < 0.05
    ) |>
    group_by(.data$domain, .data$variable_label) |>
    arrange(desc(.data$abs_diff), .by_group = TRUE) |>
    summarise(
      definition = definition,
      `Largest HROHCA - non-HROHCA difference` = format_ci(
        first(.data$approx_difference),
        first(.data$approx_ci_low),
        first(.data$approx_ci_high),
        digits = if (kind == "Laboratory" && first(.data$variable) == "ph") 2 else 1
      ),
      `Hour of largest difference` = first(.data$icu_hour),
      `p value at largest difference` = format_p(first(.data$approx_p_value)),
      `Significant time points` = paste0(sum(.data$significant, na.rm = TRUE), "/", n()),
      .groups = "drop"
    )
}

read_support_prevalence <- function(definition) {
  read.csv(
    file.path(figure_dir, paste0("figure_hrohca_support_trajectories", suffix_for(definition), "_source.csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    mutate(
      domain = "Organ support prevalence",
      variable_label = .data$support,
      abs_diff = abs(.data$absolute_difference_pct),
      significant = .data$p_value < 0.05
    ) |>
    group_by(.data$domain, .data$variable_label) |>
    arrange(desc(.data$abs_diff), .by_group = TRUE) |>
    summarise(
      definition = definition,
      `Largest HROHCA - non-HROHCA difference` = format_est(
        first(.data$absolute_difference_pct),
        digits = 1,
        suffix = " pp"
      ),
      `Hour of largest difference` = first(.data$icu_hour),
      `p value at largest difference` = format_p(first(.data$p_value)),
      `Significant time points` = paste0(sum(.data$significant, na.rm = TRUE), "/", n()),
      .groups = "drop"
    )
}

read_outcome_incidence <- function(definition) {
  read.csv(
    file.path(figure_dir, paste0("figure_hrohca_support_outcome_combined_difference_source", suffix_for(definition), ".csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    mutate(
      domain = "Cumulative outcome incidence",
      variable_label = as.character(.data$event),
      significant = .data$p_value < 0.05
    ) |>
    group_by(.data$domain, .data$variable_label) |>
    summarise(
      definition = definition,
      `Largest HROHCA - non-HROHCA difference` = format_est(
        diff_pct[which.max(abs(diff_pct))],
        digits = 1,
        suffix = " pp"
      ),
      `Hour of largest difference` = icu_hour[which.max(abs(diff_pct))],
      `p value at largest difference` = format_p(p_value[which.max(abs(diff_pct))]),
      `72-hour HROHCA - non-HROHCA difference` = paste0(sprintf("%.1f", diff_pct[icu_hour == 72]), " pp"),
      `72-hour p value` = format_p(p_value[icu_hour == 72]),
      `Significant time points` = paste0(sum(.data$significant, na.rm = TRUE), "/", n()),
      .groups = "drop"
    )
}

trajectory_long <- bind_rows(
  read_lab_or_vital("Laboratory", "heat95"),
  read_lab_or_vital("Laboratory", "heat90"),
  read_lab_or_vital("Vital signs", "heat95"),
  read_lab_or_vital("Vital signs", "heat90"),
  read_support_prevalence("heat95"),
  read_support_prevalence("heat90"),
  read_outcome_incidence("heat95"),
  read_outcome_incidence("heat90")
)

trajectory_summary <- trajectory_long |>
  mutate(
    row_id = paste(.data$domain, .data$variable_label, sep = "||"),
    row_order = row_number()
  ) |>
  pivot_wider(
    id_cols = c("domain", "variable_label"),
    names_from = "definition",
    values_from = c(
      "Largest HROHCA - non-HROHCA difference",
      "Hour of largest difference",
      "p value at largest difference",
      "72-hour HROHCA - non-HROHCA difference",
      "72-hour p value",
      "Significant time points"
    ),
    names_glue = "{definition}_{.value}"
  ) |>
  transmute(
    Domain = .data$domain,
    Measure = .data$variable_label,
    `Heat95 largest difference` = .data$`heat95_Largest HROHCA - non-HROHCA difference`,
    `Heat95 hour` = .data$`heat95_Hour of largest difference`,
    `Heat95 p value` = .data$`heat95_p value at largest difference`,
    `Heat95 significant time points` = .data$`heat95_Significant time points`,
    `Heat90 largest difference` = .data$`heat90_Largest HROHCA - non-HROHCA difference`,
    `Heat90 hour` = .data$`heat90_Hour of largest difference`,
    `Heat90 p value` = .data$`heat90_p value at largest difference`,
    `Heat90 significant time points` = .data$`heat90_Significant time points`,
    `Heat95 72-hour difference` = .data$`heat95_72-hour HROHCA - non-HROHCA difference`,
    `Heat95 72-hour p value` = .data$`heat95_72-hour p value`,
    `Heat90 72-hour difference` = .data$`heat90_72-hour HROHCA - non-HROHCA difference`,
    `Heat90 72-hour p value` = .data$`heat90_72-hour p value`
  ) |>
  mutate(across(everything(), ~tidyr::replace_na(as.character(.x), "")))

readr::write_csv(table1_rows, file.path(table_dir, "supplement_table_heat90_vs_heat95_clinical_phenotype.csv"))
readr::write_csv(trajectory_summary, file.path(table_dir, "supplement_table_heat90_vs_heat95_trajectory_summary.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Clinical phenotype")
openxlsx::writeData(wb, "Clinical phenotype", table1_rows)
openxlsx::setColWidths(wb, "Clinical phenotype", cols = 1:ncol(table1_rows), widths = "auto")
openxlsx::addWorksheet(wb, "Trajectory summary")
openxlsx::writeData(wb, "Trajectory summary", trajectory_summary)
openxlsx::setColWidths(wb, "Trajectory summary", cols = 1:ncol(trajectory_summary), widths = "auto")
openxlsx::saveWorkbook(
  wb,
  file.path(table_dir, "supplement_table_heat90_vs_heat95_phenotype_comparison.xlsx"),
  overwrite = TRUE
)

ft_clinical <- flextable(table1_rows) |>
  autofit() |>
  align(align = "left", part = "all") |>
  theme_booktabs() |>
  bold(part = "header")

ft_trajectory <- flextable(trajectory_summary) |>
  autofit() |>
  align(align = "left", part = "all") |>
  theme_booktabs() |>
  bold(part = "header") |>
  add_footer_lines(
    values = paste(
      "Largest differences are HROHCA minus non-HROHCA over the 0-72 hour ICU window.",
      "Laboratory and vital sign differences use approximate random-effects pooled median differences.",
      "Support and cumulative outcome differences are percentage-point differences using pooled counts."
    )
  )

flextable::save_as_docx(
  `Supplemental Table. Clinical phenotype under heat95 and heat90 definitions` = ft_clinical,
  `Supplemental Table. Trajectory signal summary under heat95 and heat90 definitions` = ft_trajectory,
  path = file.path(table_dir, "supplement_table_heat90_vs_heat95_phenotype_comparison.docx")
)

message("Wrote heat90 vs heat95 phenotype comparison tables to ", table_dir)

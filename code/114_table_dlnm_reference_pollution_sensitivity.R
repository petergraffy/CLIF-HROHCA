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
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(openxlsx)
  library(readr)
})

fmt_n <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_rr <- function(rr, low, high) sprintf("%.2f (%.2f-%.2f)", rr, low, high)
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_pct <- function(x) sprintf("%.0f%%", x)
fmt_temp <- function(x) sprintf("%.1f", x)

model_labels <- c(
  "primary_humidity_adjusted" = "Median reference, humidity adjusted",
  "sensitivity_humidity_pollution_adjusted" = "Median reference, humidity and air-pollution adjusted"
)

adjustment_labels <- c(
  "primary_humidity_adjusted" = "Relative humidity; time trend; day of week; calendar year",
  "sensitivity_humidity_pollution_adjusted" = "Relative humidity; 12-month NO2 and PM2.5; time trend; day of week; calendar year"
)

pooled <- readr::read_csv(
  file.path(pooled_dir, "pooled_dlnm_random_effects_results.csv"),
  show_col_types = FALSE
) |>
  filter(
    .data$stratum == "Overall",
    .data$model %in% names(model_labels),
    .data$reference_type == "median"
  )

site <- readr::read_csv(
  file.path(pooled_dir, "all_sites_dlnm_site_estimates.csv"),
  show_col_types = FALSE
) |>
  filter(
    .data$stratum == "Overall",
    .data$model %in% names(model_labels),
    .data$reference_type == "median"
  ) |>
  mutate(
    n_days = as.numeric(.data$n_days),
    n_ohca = as.numeric(.data$n_ohca),
    reference_temp_c = as.numeric(.data$reference_temp_c),
    hot_temp_c = as.numeric(.data$hot_temp_c)
  )

site_summary <- site |>
  group_by(.data$model) |>
  summarise(
    total_site_days = sum(.data$n_days, na.rm = TRUE),
    total_ohca_admissions = sum(.data$n_ohca, na.rm = TRUE),
    reference_temp_c_weighted = weighted.mean(.data$reference_temp_c, .data$n_days, na.rm = TRUE),
    hot_temp_c_weighted = weighted.mean(.data$hot_temp_c, .data$n_days, na.rm = TRUE),
    site_specific_reference_range = paste0(fmt_temp(min(.data$reference_temp_c, na.rm = TRUE)), "-", fmt_temp(max(.data$reference_temp_c, na.rm = TRUE))),
    site_specific_hot_range = paste0(fmt_temp(min(.data$hot_temp_c, na.rm = TRUE)), "-", fmt_temp(max(.data$hot_temp_c, na.rm = TRUE))),
    .groups = "drop"
  )

source_table <- pooled |>
  left_join(site_summary, by = "model") |>
  mutate(
    sensitivity_specification = unname(model_labels[.data$model]),
    model_adjustment = unname(adjustment_labels[.data$model])
  ) |>
  arrange(match(.data$model, names(model_labels)))

display_table <- source_table |>
  transmute(
    `Sensitivity analysis` = .data$sensitivity_specification,
    `Reference` = "Site-specific median warm-season Tmax",
    `Model adjustment` = .data$model_adjustment,
    `Sites, n` = .data$k_sites,
    `OHCA admissions, n` = fmt_n(.data$total_ohca_admissions),
    `Weighted reference Tmax, °C` = fmt_temp(.data$reference_temp_c_weighted),
    `Reference range, °C` = .data$site_specific_reference_range,
    `Weighted hot Tmax, °C` = fmt_temp(.data$hot_temp_c_weighted),
    `Hot range, °C` = .data$site_specific_hot_range,
    `Cumulative RR (95% CI)` = fmt_rr(.data$ratio, .data$ci_low, .data$ci_high),
    `P value` = fmt_p(.data$p_value),
    `I²` = fmt_pct(.data$i2)
  )

readr::write_csv(
  source_table,
  file.path(table_dir, "supplement_table_dlnm_reference_pollution_sensitivity_source.csv")
)
readr::write_csv(
  display_table,
  file.path(table_dir, "supplement_table_dlnm_reference_pollution_sensitivity.csv")
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "DLNM sensitivity")
openxlsx::writeData(wb, "DLNM sensitivity", display_table)
openxlsx::setColWidths(wb, "DLNM sensitivity", cols = seq_along(display_table), widths = c(34, 34, 54, 9, 18, 20, 17, 17, 14, 22, 10, 8))
openxlsx::addStyle(
  wb,
  "DLNM sensitivity",
  style = openxlsx::createStyle(textDecoration = "bold", fgFill = "#E9EDF2", border = "Bottom", wrapText = TRUE),
  rows = 1,
  cols = seq_along(display_table),
  gridExpand = TRUE
)
openxlsx::saveWorkbook(
  wb,
  file.path(table_dir, "supplement_table_dlnm_reference_pollution_sensitivity.xlsx"),
  overwrite = TRUE
)

ft <- flextable(display_table) |>
  theme_booktabs() |>
  autofit() |>
  fontsize(size = 8.5, part = "all") |>
  bold(part = "header") |>
  align(align = "center", part = "all") |>
  align(j = c("Sensitivity analysis", "Reference", "Model adjustment"), align = "left", part = "all") |>
  set_caption(
    "Supplementary Table. Median-referenced and air-pollution-adjusted DLNM sensitivity analyses for daily maximum temperature and OHCA ICU admission."
  ) |>
  add_footer_lines(
    "Random-effects pooled cumulative relative risks compare site-specific hot temperature with the site-specific median warm-season reference temperature. Weighted temperatures are weighted by site days. The air-pollution-adjusted model additionally includes county-year 12-month NO2 and PM2.5. I² indicates between-site heterogeneity."
  )

flextable::save_as_docx(
  ft,
  path = file.path(table_dir, "supplement_table_dlnm_reference_pollution_sensitivity.docx")
)

message("Wrote DLNM reference/pollution sensitivity table to ", table_dir)

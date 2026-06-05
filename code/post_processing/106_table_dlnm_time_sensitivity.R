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

repo_root <- normalizePath(file.path(dirname(get_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
box_root <- "/Users/saborpete/Library/CloudStorage/Box-Box/CLIF/Projects/CLIF-Heat-Related-OHCA"
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(pooled_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(openxlsx)
  library(readr)
  library(stringr)
})

read_export <- function(path) {
  out <- readr::read_csv(path, show_col_types = FALSE)
  if (!"site_name" %in% names(out)) {
    parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
    idx <- match("CLIF-Heat-Related-OHCA", parts)
    out$site_name <- if (!is.na(idx) && length(parts) > idx) parts[[idx + 1L]] else sub("_.*$", "", basename(path))
  }
  out$source_file <- path
  out
}

pool_der_simonian_laird <- function(est, se) {
  keep <- is.finite(est) & is.finite(se) & se >= 0
  est <- est[keep]
  se <- se[keep]
  k <- length(est)
  if (k == 0) return(NULL)
  vi <- pmax(se^2, .Machine$double.eps)
  wi <- 1 / vi
  fixed <- sum(wi * est) / sum(wi)
  q <- sum(wi * (est - fixed)^2)
  c_val <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (k > 1) max(0, (q - (k - 1)) / c_val) else 0
  w_re <- 1 / (vi + tau2)
  pooled <- sum(w_re * est) / sum(w_re)
  pooled_se <- sqrt(1 / sum(w_re))
  tibble(
    k_sites = k,
    log_estimate = pooled,
    log_estimate_se = pooled_se,
    cumulative_rr = exp(pooled),
    ci_low = exp(pooled - 1.96 * pooled_se),
    ci_high = exp(pooled + 1.96 * pooled_se),
    p_value = 2 * stats::pnorm(abs(pooled / pooled_se), lower.tail = FALSE),
    tau2 = tau2,
    q = q,
    i2 = ifelse(k > 1 && q > (k - 1) && q > 0, 100 * (q - (k - 1)) / q, 0)
  )
}

fmt_n <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_rr <- function(rr, low, high) sprintf("%.2f (%.2f-%.2f)", rr, low, high)
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_pct <- function(x) sprintf("%.0f%%", x)
fmt_temp <- function(x) sprintf("%.1f", x)

model_labels <- c(
  "sensitivity_time_df_3_per_year" = "Time spline: 3 df/year",
  "sensitivity_time_df_4_per_year" = "Time spline: 4 df/year",
  "sensitivity_time_df_6_per_year" = "Time spline: 6 df/year",
  "sensitivity_without_day_of_week" = "Without day-of-week adjustment"
)

time_files <- list.files(
  box_root,
  pattern = "^[A-Za-z0-9]+_dlnm_time_sensitivity[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
time_files <- time_files[grepl("/federated_exports/", time_files)]
if (length(time_files) == 0) stop("No DLNM time-sensitivity exports found.")

site_time <- bind_rows(lapply(time_files, read_export)) |>
  mutate(
    n_days = suppressWarnings(as.numeric(.data$n_days)),
    n_ohca = suppressWarnings(as.numeric(.data$n_ohca)),
    reference_temp_c = suppressWarnings(as.numeric(.data$reference_temp_c)),
    hot_temp_c = suppressWarnings(as.numeric(.data$hot_temp_c)),
    cumulative_rr = suppressWarnings(as.numeric(.data$cumulative_rr)),
    cumulative_rr_low = suppressWarnings(as.numeric(.data$cumulative_rr_low)),
    cumulative_rr_high = suppressWarnings(as.numeric(.data$cumulative_rr_high)),
    log_rr = suppressWarnings(as.numeric(.data$log_rr)),
    log_rr_se = suppressWarnings(as.numeric(.data$log_rr_se)),
    time_df_per_year = suppressWarnings(as.numeric(.data$time_df_per_year)),
    includes_day_of_week = .data$includes_day_of_week %in% c(TRUE, "TRUE", "true", "1"),
    includes_year_fixed_effect = .data$includes_year_fixed_effect %in% c(TRUE, "TRUE", "true", "1")
  ) |>
  filter(.data$stratum == "Overall", .data$reference_type == "median")

pooled_rows <- lapply(split(site_time, site_time$model), function(dat) {
  pooled <- pool_der_simonian_laird(dat$log_rr, dat$log_rr_se)
  if (is.null(pooled)) return(NULL)
  tibble(
    model = dat$model[1],
    specification = unname(model_labels[dat$model[1]]),
    time_df_per_year = unique(dat$time_df_per_year)[1],
    day_of_week_adjusted = unique(dat$includes_day_of_week)[1],
    year_fixed_effect = unique(dat$includes_year_fixed_effect)[1],
    total_site_days = sum(dat$n_days, na.rm = TRUE),
    total_ohca_admissions = sum(dat$n_ohca, na.rm = TRUE),
    reference_temp_c_weighted = weighted.mean(dat$reference_temp_c, dat$n_days, na.rm = TRUE),
    hot_temp_c_weighted = weighted.mean(dat$hot_temp_c, dat$n_days, na.rm = TRUE)
  ) |>
    bind_cols(pooled)
})

pooled <- bind_rows(pooled_rows) |>
  mutate(
    specification = if_else(is.na(.data$specification), .data$model, .data$specification),
    day_of_week_adjusted = if_else(.data$day_of_week_adjusted, "Yes", "No"),
    year_fixed_effect = if_else(.data$year_fixed_effect, "Yes", "No")
  ) |>
  arrange(match(.data$model, names(model_labels)))

display_table <- pooled |>
  transmute(
    `Sensitivity specification` = .data$specification,
    `Time df/year` = if_else(is.na(.data$time_df_per_year), "", sprintf("%.0f", .data$time_df_per_year)),
    `Day-of-week adjusted` = .data$day_of_week_adjusted,
    `Year fixed effects` = .data$year_fixed_effect,
    `Sites, n` = .data$k_sites,
    `OHCA admissions, n` = fmt_n(.data$total_ohca_admissions),
    `Reference Tmax, °C` = fmt_temp(.data$reference_temp_c_weighted),
    `Hot Tmax, °C` = fmt_temp(.data$hot_temp_c_weighted),
    `Cumulative RR (95% CI)` = fmt_rr(.data$cumulative_rr, .data$ci_low, .data$ci_high),
    `P value` = fmt_p(.data$p_value),
    `I²` = fmt_pct(.data$i2)
  )

readr::write_csv(site_time, file.path(pooled_dir, "all_sites_dlnm_time_sensitivity.csv"))
readr::write_csv(pooled, file.path(pooled_dir, "pooled_dlnm_time_sensitivity.csv"))
readr::write_csv(display_table, file.path(table_dir, "supplement_table_dlnm_time_sensitivity.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "DLNM time sensitivity")
openxlsx::writeData(wb, "DLNM time sensitivity", display_table)
openxlsx::setColWidths(wb, "DLNM time sensitivity", cols = seq_along(display_table), widths = c(34, 13, 20, 18, 10, 18, 18, 15, 23, 11, 9))
openxlsx::addStyle(
  wb,
  "DLNM time sensitivity",
  style = openxlsx::createStyle(textDecoration = "bold", fgFill = "#E9EDF2", border = "Bottom"),
  rows = 1,
  cols = seq_along(display_table),
  gridExpand = TRUE
)
openxlsx::saveWorkbook(wb, file.path(table_dir, "supplement_table_dlnm_time_sensitivity.xlsx"), overwrite = TRUE)

ft <- flextable(display_table) |>
  theme_booktabs() |>
  autofit() |>
  fontsize(size = 9, part = "all") |>
  bold(part = "header") |>
  align(align = "center", part = "all") |>
  align(j = "Sensitivity specification", align = "left", part = "all") |>
  set_caption(
    "Supplementary Table. DLNM time-adjustment sensitivity analyses for the association between daily maximum temperature and OHCA ICU admission."
  ) |>
  add_footer_lines(
    "Random-effects pooled cumulative relative risks compare site-specific hot temperature with the site-specific median warm-season reference temperature. Hot and reference temperatures are n-day-weighted averages of site-specific values. I² indicates between-site heterogeneity."
  )

flextable::save_as_docx(ft, path = file.path(table_dir, "supplement_table_dlnm_time_sensitivity.docx"))

message("Wrote DLNM time-sensitivity supplement table to ", table_dir)

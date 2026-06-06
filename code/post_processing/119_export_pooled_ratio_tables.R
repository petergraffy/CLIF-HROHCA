#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
dir.create(pooled_dir, recursive = TRUE, showWarnings = FALSE)

read_if_exists <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  readr::read_csv(path, show_col_types = FALSE)
}

format_ratio <- function(ratio, low, high) {
  ifelse(
    is.finite(ratio) & is.finite(low) & is.finite(high),
    sprintf("%.2f (%.2f, %.2f)", ratio, low, high),
    NA_character_
  )
}

term_family <- function(term) {
  dplyr::case_when(
    term == "(Intercept)" ~ "intercept",
    grepl("tmax|temperature", term, ignore.case = TRUE) ~ "temperature_spline_basis",
    grepl("rmax|humidity", term, ignore.case = TRUE) ~ "humidity_spline_basis",
    TRUE ~ "covariate"
  )
}

dlnm_primary <- read_if_exists(file.path(pooled_dir, "pooled_dlnm_random_effects_results.csv"))
dlnm_lags <- read_if_exists(file.path(pooled_dir, "pooled_dlnm_random_effects_lag_summaries.csv"))
dlnm_lag_specific <- read_if_exists(file.path(pooled_dir, "pooled_dlnm_random_effects_lag_specific_summaries.csv"))

dlnm_ratio_results <- tibble::tibble()
if (nrow(dlnm_primary) > 0) {
  dlnm_ratio_results <- dlnm_primary |>
    dplyr::transmute(
      analysis_family = "DLNM",
      analysis = .data$export_family,
      model = .data$model,
      stratum = .data$stratum,
      contrast = paste0(.data$reference_type, "-referenced cumulative effect"),
      lag_window = "Modeled lag window",
      effect_measure = "RR",
      log_effect = .data$log_rr,
      standard_error = .data$log_rr_se,
      ratio = .data$rr,
      ratio_low = .data$rr_low,
      ratio_high = .data$rr_high,
      ratio_ci = format_ratio(.data$rr, .data$rr_low, .data$rr_high),
      k_sites = .data$k_sites,
      tau2 = .data$tau2,
      i2 = .data$i2,
      p_value = NA_real_,
      interpretation_note = "Random-effects pooled DLNM relative risk."
    )
}

dlnm_lag_ratio_results <- tibble::tibble()
if (nrow(dlnm_lags) > 0) {
  dlnm_lag_ratio_results <- dlnm_lags |>
    dplyr::transmute(
      analysis_family = "DLNM cumulative lag windows",
      analysis = .data$export_family,
      model = .data$model,
      stratum = .data$stratum,
      contrast = paste0(.data$reference_type, "-referenced cumulative effect"),
      lag_window = .data$lag_label,
      lag_start = .data$lag_start,
      lag_end = .data$lag_end,
      effect_measure = "RR",
      log_effect = .data$log_rr,
      standard_error = .data$log_rr_se,
      ratio = .data$rr,
      ratio_low = .data$rr_low,
      ratio_high = .data$rr_high,
      ratio_ci = format_ratio(.data$rr, .data$rr_low, .data$rr_high),
      k_sites = .data$k_sites,
      tau2 = .data$tau2,
      i2 = .data$i2,
      p_value = NA_real_,
      interpretation_note = "Random-effects pooled DLNM cumulative relative risk for the specified lag window."
    )
}

dlnm_lag_specific_ratio_results <- tibble::tibble()
if (nrow(dlnm_lag_specific) > 0) {
  dlnm_lag_specific_ratio_results <- dlnm_lag_specific |>
    dplyr::transmute(
      analysis_family = "DLNM lag-specific",
      analysis = .data$export_family,
      model = .data$model,
      stratum = .data$stratum,
      contrast = paste0(.data$reference_type, "-referenced lag-specific effect"),
      lag_window = .data$lag_label,
      lag = .data$lag,
      effect_measure = "RR",
      log_effect = .data$log_rr,
      standard_error = .data$log_rr_se,
      ratio = .data$rr,
      ratio_low = .data$rr_low,
      ratio_high = .data$rr_high,
      ratio_ci = format_ratio(.data$rr, .data$rr_low, .data$rr_high),
      k_sites = .data$k_sites,
      tau2 = .data$tau2,
      i2 = .data$i2,
      p_value = NA_real_,
      interpretation_note = "Random-effects pooled DLNM lag-specific relative risk."
    )
}

phenotype_primary <- read_if_exists(file.path(pooled_dir, "pooled_ohca_icu_72h_phenotype_assignment_coefficients.csv"))
phenotype_mechanism <- read_if_exists(file.path(pooled_dir, "pooled_ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted.csv"))
phenotype_lag <- read_if_exists(file.path(pooled_dir, "pooled_ohca_icu_72h_phenotype_assignment_lag_sensitivity_coefficients.csv"))

make_phenotype_ratios <- function(df, adjustment) {
  if (nrow(df) == 0) return(tibble::tibble())
  if (!"exposure_window" %in% names(df)) df$exposure_window <- "lag0"
  df |>
    dplyr::filter(.data$coefficient != "(Intercept)") |>
    dplyr::transmute(
      analysis_family = "72h phenotype multinomial",
      analysis = dplyr::case_when(
        .data$exposure_window == "lag0" ~ adjustment,
        grepl("mechanism", .data$model) ~ paste0("temperature_humidity_demographics_mechanism_", .data$exposure_window),
        TRUE ~ paste0("temperature_humidity_demographics_", .data$exposure_window)
      ),
      exposure_window = .data$exposure_window,
      model = .data$model,
      outcome = .data$outcome_level,
      reference = .data$reference_level,
      term = .data$coefficient,
      term_family = term_family(.data$coefficient),
      effect_measure = "OR",
      log_effect = .data$estimate,
      standard_error = .data$standard_error,
      ratio = .data$odds_ratio,
      ratio_low = .data$odds_ratio_low,
      ratio_high = .data$odds_ratio_high,
      ratio_ci = format_ratio(.data$odds_ratio, .data$odds_ratio_low, .data$odds_ratio_high),
      k_sites = .data$k_sites,
      tau2 = .data$tau2,
      i2 = .data$i2,
      p_value = .data$p_value,
      interpretation_note = dplyr::if_else(
        term_family(.data$coefficient) %in% c("temperature_spline_basis", "humidity_spline_basis"),
        "Spline-basis coefficient OR; use temperature curve plots for exposure-response interpretation.",
        "Random-effects pooled multinomial odds ratio."
      )
    )
}

phenotype_ratio_results <- dplyr::bind_rows(
  make_phenotype_ratios(phenotype_primary, "temperature_humidity_demographics"),
  make_phenotype_ratios(phenotype_mechanism, "temperature_humidity_demographics_mechanism"),
  make_phenotype_ratios(phenotype_lag, "temperature_humidity_lag_sensitivity")
)

fine_gray <- read_if_exists(file.path(pooled_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv"))
fine_gray_lag <- read_if_exists(file.path(pooled_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_coefficients.csv"))
fine_gray_ratio_results <- tibble::tibble()
make_fine_gray_ratios <- function(df) {
  if (nrow(df) == 0) return(tibble::tibble())
  if (!"exposure_window" %in% names(df)) {
    df$exposure_window <- dplyr::case_when(
      grepl("lag0_1", df$model) ~ "lag0_1",
      grepl("lag0_3", df$model) ~ "lag0_3",
      grepl("lag0_5", df$model) ~ "lag0_5",
      TRUE ~ "lag0"
    )
  }
  df |>
    dplyr::transmute(
      analysis_family = "72h competing risk Fine-Gray",
      analysis = ifelse(.data$exposure_window == "lag0", .data$model, paste0(.data$model, "_sensitivity")),
      exposure_window = .data$exposure_window,
      model = .data$model,
      outcome = "Awake/extubated by 72h",
      competing_event = "Death before awake/extubated by 72h",
      term = .data$coefficient,
      term_family = term_family(.data$coefficient),
      effect_measure = "SHR",
      log_effect = .data$estimate,
      standard_error = .data$standard_error,
      ratio = .data$subdistribution_hr,
      ratio_low = .data$subdistribution_hr_low,
      ratio_high = .data$subdistribution_hr_high,
      ratio_ci = format_ratio(.data$subdistribution_hr, .data$subdistribution_hr_low, .data$subdistribution_hr_high),
      k_sites = .data$k_sites,
      tau2 = .data$tau2,
      i2 = .data$i2,
      p_value = .data$p_value,
      interpretation_note = dplyr::if_else(
        term_family(.data$coefficient) %in% c("temperature_spline_basis", "humidity_spline_basis"),
        "Spline-basis coefficient SHR; use CIF plots and model terms together for interpretation.",
        "Random-effects pooled Fine-Gray subdistribution hazard ratio."
      )
    )
}
fine_gray_ratio_results <- dplyr::bind_rows(
  make_fine_gray_ratios(fine_gray),
  make_fine_gray_ratios(fine_gray_lag)
)

ratio_tables <- list(
  dlnm_ratio_results = dlnm_ratio_results,
  dlnm_lag_window_ratio_results = dlnm_lag_ratio_results,
  dlnm_lag_specific_ratio_results = dlnm_lag_specific_ratio_results,
  phenotype_ratio_results = phenotype_ratio_results,
  fine_gray_ratio_results = fine_gray_ratio_results
)

for (name in names(ratio_tables)) {
  if (nrow(ratio_tables[[name]]) > 0) {
    readr::write_csv(ratio_tables[[name]], file.path(pooled_dir, paste0("pooled_", name, ".csv")))
  }
}

combined <- dplyr::bind_rows(
  dlnm_ratio_results,
  dlnm_lag_ratio_results,
  dlnm_lag_specific_ratio_results,
  phenotype_ratio_results,
  fine_gray_ratio_results
)
if (nrow(combined) > 0) {
  readr::write_csv(combined, file.path(pooled_dir, "pooled_all_ratio_results.csv"))
}

message("Wrote pooled ratio result tables to ", pooled_dir)

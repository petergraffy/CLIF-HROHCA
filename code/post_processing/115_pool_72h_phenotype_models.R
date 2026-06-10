#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    ofiles <- vapply(sys.frames(), function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile, character(1))
    ofiles <- stats::na.omit(ofiles)
    if (length(ofiles) > 0) return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = TRUE))
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      active_path <- rstudioapi::getActiveDocumentContext()$path
      if (nzchar(active_path)) return(normalizePath(active_path, winslash = "/", mustWork = TRUE))
    }
    stop("Could not determine script path. Run with Rscript or source the script from RStudio.")
  }
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "output", "final", "federated_exports")
output_dir <- file.path(repo_root, "output", "final", "federated_pooled")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pool_der_simonian_laird <- function(est, se) {
  keep <- is.finite(est) & is.finite(se) & se > 0
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
  data.frame(
    k_sites = k,
    estimate = pooled,
    standard_error = pooled_se,
    ci_low = pooled - 1.96 * pooled_se,
    ci_high = pooled + 1.96 * pooled_se,
    p_value = 2 * stats::pnorm(abs(pooled / pooled_se), lower.tail = FALSE),
    tau2 = tau2,
    q = q,
    i2 = ifelse(k > 1 & q > (k - 1) & q > 0, 100 * (q - (k - 1)) / q, 0),
    stringsAsFactors = FALSE
  )
}

read_bind <- function(pattern) {
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(data.frame())
  dat <- lapply(files, read.csv, stringsAsFactors = FALSE)
  all_names <- unique(unlist(lapply(dat, names)))
  aligned <- lapply(dat, function(x) {
    missing_names <- setdiff(all_names, names(x))
    for (name in missing_names) x[[name]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, aligned)
}

pool_rows <- function(df, group_cols, exponentiate = FALSE, effect_label = "ratio") {
  if (nrow(df) == 0) return(data.frame())
  groups <- unique(df[, group_cols, drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    keep <- rep(TRUE, nrow(df))
    for (col in group_cols) keep <- keep & df[[col]] == g[[col]]
    pooled <- pool_der_simonian_laird(df$estimate[keep], df$standard_error[keep])
    if (is.null(pooled)) next
    rows[[i]] <- cbind(g, pooled, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  if (exponentiate) {
    out[[effect_label]] <- exp(out$estimate)
    out[[paste0(effect_label, "_low")]] <- exp(out$ci_low)
    out[[paste0(effect_label, "_high")]] <- exp(out$ci_high)
  }
  out
}

add_multinomial_standard_errors <- function(coef_df, vcov_df) {
  if (nrow(coef_df) == 0 || nrow(vcov_df) == 0) return(data.frame())
  diag_vcov <- vcov_df[
    vcov_df$outcome_level_row == vcov_df$outcome_level_col &
      vcov_df$coefficient_row == vcov_df$coefficient_col,
    ,
    drop = FALSE
  ]
  diag_vcov$outcome_level <- diag_vcov$outcome_level_row
  diag_vcov$coefficient <- diag_vcov$coefficient_row
  diag_vcov$variance <- diag_vcov$covariance
  join_cols <- intersect(
    c("site_name", "model", "exposure_window", "reference_level", "n", "spline_df", "outcome_level", "coefficient"),
    names(coef_df)
  )
  out <- merge(
    coef_df,
    diag_vcov[, c(join_cols, "variance"), drop = FALSE],
    by = join_cols,
    all.x = TRUE,
    sort = FALSE
  )
  out$standard_error <- sqrt(out$variance)
  out
}

ed_death_summary <- read_bind("^[^_]+_ohca_ed_only_death_never_icu_summary[.]csv$")
if (nrow(ed_death_summary) > 0) {
  write.csv(ed_death_summary, file.path(output_dir, "all_site_ohca_ed_only_death_never_icu_summary.csv"), row.names = FALSE)
  pooled_ed_death <- data.frame(
    site_name = "Pooled",
    diagnosis_source = paste(sort(unique(ed_death_summary$diagnosis_source)), collapse = "; "),
    study_start = min(ed_death_summary$study_start, na.rm = TRUE),
    study_end = max(ed_death_summary$study_end, na.rm = TRUE),
    n_adult_hospitalizations = sum(ed_death_summary$n_adult_hospitalizations, na.rm = TRUE),
    n_ohca_dx_hospitalizations = sum(ed_death_summary$n_ohca_dx_hospitalizations, na.rm = TRUE),
    n_ohca_dx_patients = sum(ed_death_summary$n_ohca_dx_patients, na.rm = TRUE),
    n_ohca_dx_first_location_ed = sum(ed_death_summary$n_ohca_dx_first_location_ed, na.rm = TRUE),
    n_ohca_dx_ed_first_only_never_icu = sum(ed_death_summary$n_ohca_dx_ed_first_only_never_icu, na.rm = TRUE),
    n_ohca_dx_ed_first_only_death_never_icu_hospitalizations = sum(ed_death_summary$n_ohca_dx_ed_first_only_death_never_icu_hospitalizations, na.rm = TRUE),
    n_ohca_dx_ed_first_only_death_never_icu_patients = sum(ed_death_summary$n_ohca_dx_ed_first_only_death_never_icu_patients, na.rm = TRUE),
    definition = paste(sort(unique(ed_death_summary$definition)), collapse = " | "),
    stringsAsFactors = FALSE
  )
  write.csv(pooled_ed_death, file.path(output_dir, "pooled_ohca_ed_only_death_never_icu_summary.csv"), row.names = FALSE)
}

ed_death_pathway <- read_bind("^[^_]+_ohca_ed_only_death_never_icu_pathway_audit[.]csv$")
if (nrow(ed_death_pathway) > 0) {
  write.csv(ed_death_pathway, file.path(output_dir, "all_site_ohca_ed_only_death_never_icu_pathway_audit.csv"), row.names = FALSE)
}

phenotype_summary <- read_bind("^[^_]+_ohca_icu_72h_phenotype_summary[.]csv$")
if (nrow(phenotype_summary) > 0) {
  write.csv(phenotype_summary, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_summary.csv"), row.names = FALSE)
  pooled_phenotypes <- aggregate(n ~ phenotype, data = phenotype_summary, FUN = sum)
  pooled_phenotypes$pct <- 100 * pooled_phenotypes$n / sum(pooled_phenotypes$n)
  pooled_phenotypes <- pooled_phenotypes[, c("phenotype", "n", "pct")]
  write.csv(pooled_phenotypes, file.path(output_dir, "pooled_ohca_icu_72h_phenotype_summary.csv"), row.names = FALSE)
}

for (suffix in c(
  "table1",
  "table2_by_phenotype",
  "consort_flow",
  "gcs_landmark_by_phenotype",
  "phenotype_evidence_summary",
  "admission_temperature_distribution_summary",
  "admission_temperature_density",
  "ohca_mechanism_summary"
)) {
  dat <- read_bind(paste0("^[^_]+_ohca_icu_72h_", suffix, "[.]csv$"))
  if (nrow(dat) > 0) {
    write.csv(dat, file.path(output_dir, paste0("all_site_ohca_icu_72h_", suffix, ".csv")), row.names = FALSE)
  }
}

phenotype_model <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_model[.]csv$")
phenotype_model_mechanism <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted[.]csv$")
if (nrow(phenotype_model) > 0) write.csv(phenotype_model, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_assignment_model.csv"), row.names = FALSE)
if (nrow(phenotype_model_mechanism) > 0) write.csv(phenotype_model_mechanism, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted.csv"), row.names = FALSE)

phenotype_coef <- add_multinomial_standard_errors(
  read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_coefficients[.]csv$"),
  read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_vcov[.]csv$")
)
if (nrow(phenotype_coef) > 0) {
  write.csv(phenotype_coef, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_assignment_coefficients_with_se.csv"), row.names = FALSE)
  pooled <- pool_rows(phenotype_coef, c("model", "reference_level", "outcome_level", "coefficient"), exponentiate = TRUE, effect_label = "odds_ratio")
  write.csv(pooled, file.path(output_dir, "pooled_ohca_icu_72h_phenotype_assignment_coefficients.csv"), row.names = FALSE)
}

phenotype_coef_mechanism <- add_multinomial_standard_errors(
  read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted[.]csv$"),
  read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_vcov_mechanism_adjusted[.]csv$")
)
if (nrow(phenotype_coef_mechanism) > 0) {
  write.csv(phenotype_coef_mechanism, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted_with_se.csv"), row.names = FALSE)
  pooled <- pool_rows(phenotype_coef_mechanism, c("model", "reference_level", "outcome_level", "coefficient"), exponentiate = TRUE, effect_label = "odds_ratio")
  write.csv(pooled, file.path(output_dir, "pooled_ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted.csv"), row.names = FALSE)
}

phenotype_lag_model <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_lag_sensitivity_model[.]csv$")
phenotype_lag_curve <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_lag_sensitivity_temperature_curve[.]csv$")
if (nrow(phenotype_lag_model) > 0) {
  write.csv(phenotype_lag_model, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_assignment_lag_sensitivity_model.csv"), row.names = FALSE)
}
if (nrow(phenotype_lag_curve) > 0) {
  write.csv(phenotype_lag_curve, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_assignment_lag_sensitivity_temperature_curve.csv"), row.names = FALSE)
}

phenotype_lag_coef <- add_multinomial_standard_errors(
  read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_lag_sensitivity_coefficients[.]csv$"),
  read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_lag_sensitivity_vcov[.]csv$")
)
if (nrow(phenotype_lag_coef) > 0) {
  write.csv(phenotype_lag_coef, file.path(output_dir, "all_site_ohca_icu_72h_phenotype_assignment_lag_sensitivity_coefficients_with_se.csv"), row.names = FALSE)
  pooled <- pool_rows(
    phenotype_lag_coef,
    c("model", "exposure_window", "reference_level", "outcome_level", "coefficient"),
    exponentiate = TRUE,
    effect_label = "odds_ratio"
  )
  write.csv(pooled, file.path(output_dir, "pooled_ohca_icu_72h_phenotype_assignment_lag_sensitivity_coefficients.csv"), row.names = FALSE)
}

for (suffix in c("summary", "death_source_summary", "ohca_mechanism_summary", "fine_gray_models")) {
  dat <- read_bind(paste0("^[^_]+_ohca_icu_competing_risk_awake_extubated_72h_", suffix, "[.]csv$"))
  if (nrow(dat) > 0) {
    write.csv(dat, file.path(output_dir, paste0("all_site_ohca_icu_competing_risk_awake_extubated_72h_", suffix, ".csv")), row.names = FALSE)
  }
}
death_source <- read_bind("^[^_]+_ohca_icu_competing_risk_death_source_summary[.]csv$")
if (nrow(death_source) > 0) write.csv(death_source, file.path(output_dir, "all_site_ohca_icu_competing_risk_death_source_summary.csv"), row.names = FALSE)

competing_summary <- read_bind("^[^_]+_ohca_icu_competing_risk_awake_extubated_72h_summary[.]csv$")
if (nrow(competing_summary) > 0) {
  pooled_competing <- aggregate(n ~ status_label, data = competing_summary, FUN = sum)
  pooled_competing$pct <- 100 * pooled_competing$n / sum(pooled_competing$n)
  pooled_competing <- pooled_competing[, c("status_label", "n", "pct")]
  write.csv(pooled_competing, file.path(output_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_summary.csv"), row.names = FALSE)
}

cif_curves <- read_bind("^[^_]+_ohca_icu_competing_risk_awake_extubated_72h_cif_curves[.]csv$")
if (nrow(cif_curves) > 0) {
  if (!"stratification" %in% names(cif_curves)) cif_curves$stratification <- "Overall"
  if (!"stratum" %in% names(cif_curves)) cif_curves$stratum <- "Overall"
  if (!"stratum_order" %in% names(cif_curves)) cif_curves$stratum_order <- 0
  if (!"stratum_min_temp_c" %in% names(cif_curves)) cif_curves$stratum_min_temp_c <- NA_real_
  if (!"stratum_max_temp_c" %in% names(cif_curves)) cif_curves$stratum_max_temp_c <- NA_real_
  write.csv(cif_curves, file.path(output_dir, "all_site_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv"), row.names = FALSE)
  pooled_cif <- aggregate(
    cbind(weighted_cif = cif * n, weighted_variance = variance * n^2, n) ~ model + stratification + stratum + stratum_order + event_code + event_type + time_hours,
    data = cif_curves,
    FUN = sum
  )
  pooled_cif$cif <- pooled_cif$weighted_cif / pooled_cif$n
  pooled_cif$variance <- pooled_cif$weighted_variance / pooled_cif$n^2
  pooled_cif$standard_error <- sqrt(pmax(pooled_cif$variance, 0))
  pooled_cif$cif_low <- pmax(0, pooled_cif$cif - 1.96 * pooled_cif$standard_error)
  pooled_cif$cif_high <- pmin(1, pooled_cif$cif + 1.96 * pooled_cif$standard_error)
  cif_group_cols <- c("model", "stratification", "stratum", "stratum_order", "event_code", "event_type", "time_hours")
  k_sites <- aggregate(
    site_name ~ model + stratification + stratum + stratum_order + event_code + event_type + time_hours,
    data = cif_curves,
    FUN = function(x) length(unique(x))
  )
  names(k_sites)[names(k_sites) == "site_name"] <- "k_sites"
  pooled_cif <- merge(pooled_cif, k_sites, by = cif_group_cols, all.x = TRUE, sort = FALSE)

  temp_ranges <- lapply(seq_len(nrow(pooled_cif)), function(i) {
    keep <- rep(TRUE, nrow(cif_curves))
    for (col in cif_group_cols) keep <- keep & cif_curves[[col]] == pooled_cif[[col]][i]
    temp_min <- cif_curves$stratum_min_temp_c[keep]
    temp_max <- cif_curves$stratum_max_temp_c[keep]
    data.frame(
      stratum_min_temp_c = if (all(is.na(temp_min))) NA_real_ else min(temp_min, na.rm = TRUE),
      stratum_max_temp_c = if (all(is.na(temp_max))) NA_real_ else max(temp_max, na.rm = TRUE)
    )
  })
  temp_ranges <- do.call(rbind, temp_ranges)
  pooled_cif$stratum_min_temp_c <- temp_ranges$stratum_min_temp_c
  pooled_cif$stratum_max_temp_c <- temp_ranges$stratum_max_temp_c
  pooled_cif <- pooled_cif[, c(
    "model", "stratification", "stratum", "stratum_order", "stratum_min_temp_c", "stratum_max_temp_c",
    "event_code", "event_type", "time_hours", "cif", "cif_low", "cif_high", "variance", "standard_error", "n", "k_sites"
  )]
  write.csv(pooled_cif, file.path(output_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv"), row.names = FALSE)
}

fg_coef <- read_bind("^[^_]+_ohca_icu_competing_risk_awake_extubated_72h_coefficients[.]csv$")
if (nrow(fg_coef) > 0) {
  write.csv(fg_coef, file.path(output_dir, "all_site_ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv"), row.names = FALSE)
  pooled_fg <- pool_rows(fg_coef, c("model", "coefficient"), exponentiate = TRUE, effect_label = "subdistribution_hr")
  write.csv(pooled_fg, file.path(output_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv"), row.names = FALSE)
}

fg_lag_models <- read_bind("^[^_]+_ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_fine_gray_models[.]csv$")
if (nrow(fg_lag_models) > 0) {
  write.csv(fg_lag_models, file.path(output_dir, "all_site_ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_fine_gray_models.csv"), row.names = FALSE)
}

fg_lag_coef <- read_bind("^[^_]+_ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_coefficients[.]csv$")
if (nrow(fg_lag_coef) > 0) {
  write.csv(fg_lag_coef, file.path(output_dir, "all_site_ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_coefficients.csv"), row.names = FALSE)
  pooled_fg_lag <- pool_rows(fg_lag_coef, c("model", "coefficient"), exponentiate = TRUE, effect_label = "subdistribution_hr")
  write.csv(pooled_fg_lag, file.path(output_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_lag_sensitivity_coefficients.csv"), row.names = FALSE)
}

imv24_summary <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_summary[.]csv$")
if (nrow(imv24_summary) > 0) {
  write.csv(imv24_summary, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_summary.csv"), row.names = FALSE)
  pooled_imv24 <- aggregate(n ~ status_label, data = imv24_summary, FUN = sum)
  pooled_imv24$pct <- 100 * pooled_imv24$n / sum(pooled_imv24$n)
  pooled_imv24 <- pooled_imv24[, c("status_label", "n", "pct")]
  write.csv(pooled_imv24, file.path(output_dir, "pooled_ohca_icu_imv24_time_to_event_summary.csv"), row.names = FALSE)
}

imv24_flow <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_landmark_flow[.]csv$")
if (nrow(imv24_flow) > 0) write.csv(imv24_flow, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_landmark_flow.csv"), row.names = FALSE)

imv24_evidence <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_evidence_summary[.]csv$")
if (nrow(imv24_evidence) > 0) write.csv(imv24_evidence, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_evidence_summary.csv"), row.names = FALSE)

imv24_table <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_table_by_outcome[.]csv$")
if (nrow(imv24_table) > 0) write.csv(imv24_table, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_table_by_outcome.csv"), row.names = FALSE)

imv24_models <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_fine_gray_models[.]csv$")
if (nrow(imv24_models) > 0) write.csv(imv24_models, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_fine_gray_models.csv"), row.names = FALSE)

imv24_coef <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_coefficients[.]csv$")
if (nrow(imv24_coef) > 0) {
  write.csv(imv24_coef, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_coefficients.csv"), row.names = FALSE)
  pooled_imv24_coef <- pool_rows(imv24_coef, c("model", "event_type", "coefficient"), exponentiate = TRUE, effect_label = "subdistribution_hr")
  write.csv(pooled_imv24_coef, file.path(output_dir, "pooled_ohca_icu_imv24_time_to_event_coefficients.csv"), row.names = FALSE)
}

imv24_vcov <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_vcov[.]csv$")
if (nrow(imv24_vcov) > 0) write.csv(imv24_vcov, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_vcov.csv"), row.names = FALSE)

imv24_cif <- read_bind("^[^_]+_ohca_icu_imv24_time_to_event_cif_curves[.]csv$")
if (nrow(imv24_cif) > 0) {
  write.csv(imv24_cif, file.path(output_dir, "all_site_ohca_icu_imv24_time_to_event_cif_curves.csv"), row.names = FALSE)
  pooled_imv24_cif <- aggregate(
    cbind(weighted_cif = cif * n, weighted_variance = variance * n^2, n) ~ model + stratification + stratum + stratum_order + event_code + event_type + time_days + time_hours,
    data = imv24_cif,
    FUN = sum
  )
  pooled_imv24_cif$cif <- pooled_imv24_cif$weighted_cif / pooled_imv24_cif$n
  pooled_imv24_cif$variance <- pooled_imv24_cif$weighted_variance / pooled_imv24_cif$n^2
  pooled_imv24_cif$standard_error <- sqrt(pmax(pooled_imv24_cif$variance, 0))
  pooled_imv24_cif$cif_low <- pmax(0, pooled_imv24_cif$cif - 1.96 * pooled_imv24_cif$standard_error)
  pooled_imv24_cif$cif_high <- pmin(1, pooled_imv24_cif$cif + 1.96 * pooled_imv24_cif$standard_error)
  imv24_cif_group_cols <- c("model", "stratification", "stratum", "stratum_order", "event_code", "event_type", "time_days", "time_hours")
  k_sites <- aggregate(
    site_name ~ model + stratification + stratum + stratum_order + event_code + event_type + time_days + time_hours,
    data = imv24_cif,
    FUN = function(x) length(unique(x))
  )
  names(k_sites)[names(k_sites) == "site_name"] <- "k_sites"
  pooled_imv24_cif <- merge(pooled_imv24_cif, k_sites, by = imv24_cif_group_cols, all.x = TRUE, sort = FALSE)
  pooled_imv24_cif <- pooled_imv24_cif[, c(
    "model", "stratification", "stratum", "stratum_order", "event_code", "event_type",
    "time_days", "time_hours", "cif", "cif_low", "cif_high", "variance", "standard_error", "n", "k_sites"
  )]
  write.csv(pooled_imv24_cif, file.path(output_dir, "pooled_ohca_icu_imv24_time_to_event_cif_curves.csv"), row.names = FALSE)
}

imv12_discharge_summary <- read_bind("^[^_]+_ohca_icu_imv12_discharge_outcome_summary[.]csv$")
if (nrow(imv12_discharge_summary) > 0) {
  write.csv(imv12_discharge_summary, file.path(output_dir, "all_site_ohca_icu_imv12_discharge_outcome_summary.csv"), row.names = FALSE)
  pooled_imv12_summary <- aggregate(n ~ discharge_outcome, data = imv12_discharge_summary, FUN = sum)
  pooled_imv12_summary$pct_among_imv12 <- 100 * pooled_imv12_summary$n / sum(pooled_imv12_summary$n)
  pooled_imv12_summary <- pooled_imv12_summary[, c("discharge_outcome", "n", "pct_among_imv12")]
  write.csv(pooled_imv12_summary, file.path(output_dir, "pooled_ohca_icu_imv12_discharge_outcome_summary.csv"), row.names = FALSE)
}

imv12_discharge_flow <- read_bind("^[^_]+_ohca_icu_imv12_discharge_outcome_flow[.]csv$")
if (nrow(imv12_discharge_flow) > 0) write.csv(imv12_discharge_flow, file.path(output_dir, "all_site_ohca_icu_imv12_discharge_outcome_flow.csv"), row.names = FALSE)

imv12_discharge_table <- read_bind("^[^_]+_ohca_icu_imv12_discharge_outcome_table_by_outcome[.]csv$")
if (nrow(imv12_discharge_table) > 0) write.csv(imv12_discharge_table, file.path(output_dir, "all_site_ohca_icu_imv12_discharge_outcome_table_by_outcome.csv"), row.names = FALSE)

imv12_discharge_model <- read_bind("^[^_]+_ohca_icu_imv12_discharge_outcome_model[.]csv$")
if (nrow(imv12_discharge_model) > 0) write.csv(imv12_discharge_model, file.path(output_dir, "all_site_ohca_icu_imv12_discharge_outcome_model.csv"), row.names = FALSE)

imv12_discharge_coef <- read_bind("^[^_]+_ohca_icu_imv12_discharge_outcome_coefficients[.]csv$")
if (nrow(imv12_discharge_coef) > 0) {
  write.csv(imv12_discharge_coef, file.path(output_dir, "all_site_ohca_icu_imv12_discharge_outcome_coefficients.csv"), row.names = FALSE)
  pooled_imv12_coef <- pool_rows(
    imv12_discharge_coef,
    c("model", "outcome_level", "reference_level", "coefficient"),
    exponentiate = TRUE,
    effect_label = "odds_ratio"
  )
  write.csv(pooled_imv12_coef, file.path(output_dir, "pooled_ohca_icu_imv12_discharge_outcome_coefficients.csv"), row.names = FALSE)
}

imv12_discharge_vcov <- read_bind("^[^_]+_ohca_icu_imv12_discharge_outcome_vcov[.]csv$")
if (nrow(imv12_discharge_vcov) > 0) write.csv(imv12_discharge_vcov, file.path(output_dir, "all_site_ohca_icu_imv12_discharge_outcome_vcov.csv"), row.names = FALSE)

message("Wrote pooled 72-hour phenotype and competing-risk outputs to ", output_dir)

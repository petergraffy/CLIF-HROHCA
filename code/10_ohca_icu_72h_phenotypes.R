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

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(lubridate)
  library(readr)
  library(splines)
  library(stringr)
  library(tidyr)
})

OHCA_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
OHCA_SUMMARY_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024_summary.csv")
OHCA_RESTRICTION_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_pathway_timing_restriction_summary.csv")
TMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
RMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_rmax_2018_2024.parquet")
OUTPUT_DIR <- file.path(repo_root, "output", "final", "ohca_icu_phenotypes")
FIGURE_DIR <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(OUTPUT_DIR, c(
  "ohca_icu_72h_phenotype_heat_thresholds.csv",
  "ohca_icu_72h_heat_by_phenotype_summary.csv",
  "ohca_icu_72h_phenotype_heat_models.csv",
  "ohca_icu_72h_phenotype_multinomial_heat_models.csv",
  "ohca_icu_72h_phenotype_spline_heat_models.csv",
  "ohca_icu_72h_phenotype_spline_heat_curves.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_model.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_curves.csv",
  "ohca_icu_72h_phenotype_heat_models_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_multinomial_heat_models_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_spline_heat_models_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_spline_heat_curves_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_model_mechanism_adjusted.csv",
  "ohca_icu_72h_phenotype_multinomial_spline_heat_curves_mechanism_adjusted.csv",
  "ohca_icu_72h_middle_imv_mental_state_summary.csv",
  "ohca_icu_72h_middle_imv_mental_state_heat_summary.csv",
  "ohca_icu_72h_middle_imv_mental_state_binary_models.csv",
  "ohca_icu_72h_middle_imv_mental_state_continuous_models.csv",
  "ohca_icu_72h_middle_imv_mental_state_spline_models.csv",
  "ohca_icu_72h_middle_imv_mental_state_spline_curves.csv",
  "ohca_icu_72h_middle_imv_mental_state_patient_audit.csv"
)))

WINDOW_HOURS <- 72
LATE_START_HOUR <- 24
EXTUBATION_ASSESSMENT_START_HOUR <- 48

config <- load_project_config(repo_root)
validate_site_name(config, repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)
site_name <- config$site_name %||% Sys.getenv("CLIF_SITE_NAME", unset = "unknown_site")

norm_code <- function(x) {
  x |>
    tidyr::replace_na("") |>
    as.character() |>
    stringr::str_to_upper() |>
    stringr::str_replace_all("[^A-Z0-9]", "") |>
    trimws()
}

infer_diagnosis_format <- function(dx_clean, diagnosis_code_format) {
  explicit <- diagnosis_code_format |>
    tidyr::replace_na("") |>
    as.character() |>
    stringr::str_to_upper()
  inferred <- ifelse(
    stringr::str_detect(dx_clean, "^[A-Z]"),
    "ICD10",
    ifelse(stringr::str_detect(dx_clean, "^[0-9]"), "ICD9", "")
  )
  ifelse(nchar(explicit) > 0, explicit, inferred)
}

safe_last <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  tail(x, 1)
}

yesno <- function(x) {
  ifelse(is.na(x), FALSE, x)
}

clean_binary_text <- function(x) {
  stringr::str_to_lower(stringr::str_squish(tidyr::replace_na(as.character(x), "")))
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(stringr::str_squish(x))] <- "Missing"
  stringr::str_squish(x)
}

fmt_n_pct <- function(n, denom) {
  if (is.na(denom) || denom <= 0) return(NA_character_)
  sprintf("%s (%.1f%%)", format(as.integer(n), big.mark = ","), 100 * n / denom)
}

fmt_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  sprintf(
    "%.1f [%.1f, %.1f]",
    stats::median(x, na.rm = TRUE),
    stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  )
}

fmt_yes <- function(x) {
  x <- ifelse(is.na(x), NA, as.integer(x %in% c(TRUE, 1L, "1", "TRUE", "true", "Yes", "yes")))
  fmt_n_pct(sum(x == 1L, na.rm = TRUE), sum(!is.na(x)))
}

add_table1_continuous <- function(df, section, characteristic, variable) {
  if (!variable %in% names(df)) return(tibble::tibble())
  tibble::tibble(section = section, characteristic = characteristic, level = "", value = fmt_median_iqr(df[[variable]]))
}

add_table1_binary <- function(df, section, characteristic, variable) {
  if (!variable %in% names(df)) return(tibble::tibble())
  tibble::tibble(section = section, characteristic = characteristic, level = "Yes", value = fmt_yes(df[[variable]]))
}

add_table1_categorical <- function(df, section, characteristic, variable) {
  if (!variable %in% names(df)) return(tibble::tibble())
  counts <- sort(table(safe_chr(df[[variable]])), decreasing = TRUE)
  tibble::tibble(
    section = section,
    characteristic = characteristic,
    level = names(counts),
    value = vapply(as.integer(counts), fmt_n_pct, character(1), denom = nrow(df))
  )
}

add_table2_continuous <- function(df, section, characteristic, variable, phenotype_levels) {
  if (!variable %in% names(df)) return(tibble::tibble())
  vals <- vapply(phenotype_levels, function(ph) fmt_median_iqr(df[[variable]][df$phenotype == ph]), character(1))
  tibble::tibble(section = section, characteristic = characteristic, level = "") |>
    bind_cols(as_tibble(as.list(vals)))
}

add_table2_binary <- function(df, section, characteristic, variable, phenotype_levels) {
  if (!variable %in% names(df)) return(tibble::tibble())
  vals <- vapply(phenotype_levels, function(ph) fmt_yes(df[[variable]][df$phenotype == ph]), character(1))
  tibble::tibble(section = section, characteristic = characteristic, level = "Yes") |>
    bind_cols(as_tibble(as.list(vals)))
}

add_table2_categorical <- function(df, section, characteristic, variable, phenotype_levels) {
  if (!variable %in% names(df)) return(tibble::tibble())
  levels_all <- names(sort(table(safe_chr(df[[variable]])), decreasing = TRUE))
  rows <- lapply(levels_all, function(lvl) {
    vals <- vapply(phenotype_levels, function(ph) {
      sub <- safe_chr(df[[variable]][df$phenotype == ph])
      fmt_n_pct(sum(sub == lvl, na.rm = TRUE), length(sub))
    }, character(1))
    tibble::tibble(section = section, characteristic = characteristic, level = lvl) |>
      bind_cols(as_tibble(as.list(vals)))
  })
  bind_rows(rows)
}

fmt_or <- function(beta, se) {
  if (!is.finite(beta) || !is.finite(se)) return(c(NA_real_, NA_real_, NA_real_))
  c(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
}

choose_terms <- function(df, terms) {
  keep <- character()
  omitted <- character()
  for (term in terms) {
    if (!term %in% names(df)) {
      omitted <- c(omitted, paste0(term, "_missing"))
      next
    }
    x <- df[[term]]
    if (is.numeric(x)) {
      if (sum(is.finite(x)) > 0 && length(unique(x[is.finite(x)])) >= 2) keep <- c(keep, term) else omitted <- c(omitted, paste0(term, "_constant"))
    } else {
      if (length(unique(stats::na.omit(x))) >= 2) keep <- c(keep, term) else omitted <- c(omitted, paste0(term, "_one_level"))
    }
  }
  list(keep = keep, omitted = omitted)
}

run_logistic <- function(df, outcome, exposure, label, adjust_terms = character()) {
  base <- df |>
    filter(!is.na(.data[[outcome]]), !is.na(.data[[exposure]]))
  events <- sum(base[[outcome]] == 1, na.rm = TRUE)
  nonevents <- sum(base[[outcome]] == 0, na.rm = TRUE)
  if (nrow(base) < 25 || events < 5 || nonevents < 5 || length(unique(base[[exposure]])) < 2) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      events = events,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Insufficient events, non-events, sample size, or exposure variation"
    ))
  }
  term_info <- choose_terms(base, adjust_terms)
  formula_txt <- paste(outcome, "~", paste(c(exposure, term_info$keep), collapse = " + "))
  fit <- tryCatch(
    stats::glm(stats::as.formula(formula_txt), data = base, family = stats::binomial()),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      events = events,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = conditionMessage(fit)
    ))
  }
  coef_tbl <- summary(fit)$coefficients
  if (!exposure %in% rownames(coef_tbl)) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(stats::model.frame(fit)),
      events = sum(stats::model.frame(fit)[[outcome]] == 1, na.rm = TRUE),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Exposure term unavailable after model fit"
    ))
  }
  if (
    !is.finite(coef_tbl[exposure, "Std. Error"]) ||
      coef_tbl[exposure, "Std. Error"] <= 1e-6 ||
      coef_tbl[exposure, "Std. Error"] > 10 ||
      !is.finite(coef_tbl[exposure, "Estimate"])
  ) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(stats::model.frame(fit)),
      events = sum(stats::model.frame(fit)[[outcome]] == 1, na.rm = TRUE),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Potential separation or numerical instability in logistic exposure estimate"
    ))
  }
  ors <- fmt_or(coef_tbl[exposure, "Estimate"], coef_tbl[exposure, "Std. Error"])
  tibble::tibble(
    outcome = outcome,
    model = label,
    exposure = exposure,
    n = nrow(stats::model.frame(fit)),
    events = sum(stats::model.frame(fit)[[outcome]] == 1, na.rm = TRUE),
    odds_ratio = ors[[1]],
    ci_low = ors[[2]],
    ci_high = ors[[3]],
    p_value = coef_tbl[exposure, "Pr(>|z|)"],
    covariates = paste(term_info$keep, collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
}

run_linear <- function(df, outcome, exposure, label, adjust_terms = character()) {
  base <- df |>
    filter(!is.na(.data[[outcome]]), is.finite(.data[[outcome]]), !is.na(.data[[exposure]]), is.finite(.data[[exposure]]))
  if (nrow(base) < 25 || length(unique(base[[outcome]])) < 2 || length(unique(base[[exposure]])) < 2) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      outcome_median = median(base[[outcome]], na.rm = TRUE),
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Insufficient sample size or outcome/exposure variation"
    ))
  }
  term_info <- choose_terms(base, adjust_terms)
  formula_txt <- paste(outcome, "~", paste(c(exposure, term_info$keep), collapse = " + "))
  fit <- tryCatch(
    stats::lm(stats::as.formula(formula_txt), data = base),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(base),
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      outcome_median = median(base[[outcome]], na.rm = TRUE),
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = conditionMessage(fit)
    ))
  }
  coef_tbl <- summary(fit)$coefficients
  if (!exposure %in% rownames(coef_tbl)) {
    return(tibble::tibble(
      outcome = outcome,
      model = label,
      exposure = exposure,
      n = nrow(stats::model.frame(fit)),
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      outcome_median = median(stats::model.frame(fit)[[outcome]], na.rm = TRUE),
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Exposure term unavailable after model fit"
    ))
  }
  beta <- coef_tbl[exposure, "Estimate"]
  se <- coef_tbl[exposure, "Std. Error"]
  tibble::tibble(
    outcome = outcome,
    model = label,
    exposure = exposure,
    n = nrow(stats::model.frame(fit)),
    estimate = beta,
    ci_low = beta - 1.96 * se,
    ci_high = beta + 1.96 * se,
    p_value = coef_tbl[exposure, "Pr(>|t|)"],
    outcome_median = median(stats::model.frame(fit)[[outcome]], na.rm = TRUE),
    covariates = paste(term_info$keep, collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
}

run_multinomial <- function(df, exposure, label, adjust_terms = character()) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = NA_character_,
      reference_level = "regained_consciousness_extubated",
      n = 0L,
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Package nnet unavailable"
    ))
  }
  keep_levels <- c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury")
  base <- df |>
    filter(.data$phenotype %in% keep_levels, !is.na(.data[[exposure]])) |>
    mutate(phenotype = stats::relevel(factor(.data$phenotype, levels = keep_levels), ref = "regained_consciousness_extubated"))
  if (nrow(base) < 50 || length(unique(base$phenotype)) < 3 || any(table(base$phenotype) < 5) || length(unique(base[[exposure]])) < 2) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = setdiff(keep_levels, "regained_consciousness_extubated"),
      reference_level = "regained_consciousness_extubated",
      n = nrow(base),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(adjust_terms, collapse = " + "),
      omitted_terms = NA_character_,
      estimable = FALSE,
      skip_reason = "Insufficient phenotype counts, sample size, or exposure variation"
    ))
  }
  term_info <- choose_terms(base, adjust_terms)
  formula_txt <- paste("phenotype ~", paste(c(exposure, term_info$keep), collapse = " + "))
  fit <- tryCatch(
    nnet::multinom(stats::as.formula(formula_txt), data = base, trace = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = setdiff(keep_levels, "regained_consciousness_extubated"),
      reference_level = "regained_consciousness_extubated",
      n = nrow(base),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = conditionMessage(fit)
    ))
  }
  s <- summary(fit)
  coef_mat <- s$coefficients
  se_mat <- s$standard.errors
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, nrow = 1, dimnames = list(names(coef_mat)[[1]], names(coef_mat)))
    se_mat <- matrix(se_mat, nrow = 1, dimnames = dimnames(coef_mat))
  }
  if (
    !exposure %in% colnames(coef_mat) ||
      any(!is.finite(se_mat[, exposure])) ||
      any(se_mat[, exposure] <= 1e-6) ||
      any(abs(coef_mat[, exposure] / se_mat[, exposure]) > 8)
  ) {
    return(tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = rownames(coef_mat),
      reference_level = "regained_consciousness_extubated",
      n = nrow(stats::model.frame(fit)),
      odds_ratio = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = FALSE,
      skip_reason = "Potential separation or numerical instability in multinomial exposure estimate"
    ))
  }
  rows <- rownames(coef_mat)
  out <- lapply(rows, function(level) {
    if (!exposure %in% colnames(coef_mat)) return(NULL)
    beta <- coef_mat[level, exposure]
    se <- se_mat[level, exposure]
    ors <- fmt_or(beta, se)
    z <- beta / se
    tibble::tibble(
      model = label,
      exposure = exposure,
      outcome_level = level,
      reference_level = "regained_consciousness_extubated",
      n = nrow(stats::model.frame(fit)),
      odds_ratio = ors[[1]],
      ci_low = ors[[2]],
      ci_high = ors[[3]],
      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = TRUE,
      skip_reason = NA_character_
    )
  })
  dplyr::bind_rows(out)
}

mode_value <- function(x) {
  x <- stats::na.omit(x)
  if (length(x) == 0L) return(NA)
  names(sort(table(x), decreasing = TRUE))[[1]]
}

make_reference_newdata <- function(df, temp_grid, adjust_terms = character()) {
  out <- tibble::tibble(tmax_mean_c = temp_grid, tmax_per5c = temp_grid / 5)
  for (term in adjust_terms) {
    if (!term %in% names(df)) next
    if (is.numeric(df[[term]])) {
      out[[term]] <- mean(df[[term]], na.rm = TRUE)
    } else {
      out[[term]] <- mode_value(df[[term]])
    }
  }
  out
}

run_spline_logistic <- function(df, outcome, label, adjust_terms = character(), spline_df = 3L) {
  base <- df |>
    filter(!is.na(.data[[outcome]]), !is.na(.data$tmax_mean_c), is.finite(.data$tmax_mean_c))
  events <- sum(base[[outcome]] == 1, na.rm = TRUE)
  nonevents <- sum(base[[outcome]] == 0, na.rm = TRUE)
  empty_summary <- tibble::tibble(
    outcome = outcome,
    model = label,
    n = nrow(base),
    events = events,
    spline_df = spline_df,
    reference_temp_c = NA_real_,
    min_temp_c = NA_real_,
    max_temp_c = NA_real_,
    linear_lrt_p_value = NA_real_,
    spline_lrt_p_value = NA_real_,
    nonlinear_lrt_p_value = NA_real_,
    aic = NA_real_,
    covariates = paste(adjust_terms, collapse = " + "),
    omitted_terms = NA_character_,
    estimable = FALSE,
    skip_reason = "Insufficient events, non-events, sample size, or temperature variation"
  )
  if (nrow(base) < 50 || events < 5 || nonevents < 5 || length(unique(base$tmax_mean_c)) < 5) {
    return(list(summary = empty_summary, curve = tibble::tibble()))
  }

  term_info <- choose_terms(base, adjust_terms)
  adjust_txt <- if (length(term_info$keep) > 0) paste("+", paste(term_info$keep, collapse = " + ")) else ""
  null_formula <- stats::as.formula(paste(outcome, "~", if (length(term_info$keep) > 0) paste(term_info$keep, collapse = " + ") else "1"))
  linear_formula <- stats::as.formula(paste(outcome, "~ tmax_per5c", adjust_txt))
  spline_formula <- stats::as.formula(paste0(outcome, " ~ splines::ns(tmax_mean_c, df = ", spline_df, ") ", adjust_txt))

  null_fit <- tryCatch(stats::glm(null_formula, data = base, family = stats::binomial()), error = function(e) e)
  linear_fit <- tryCatch(stats::glm(linear_formula, data = base, family = stats::binomial()), error = function(e) e)
  spline_fit <- tryCatch(stats::glm(spline_formula, data = base, family = stats::binomial()), error = function(e) e)
  if (inherits(null_fit, "error") || inherits(linear_fit, "error") || inherits(spline_fit, "error")) {
    reason <- paste(
      stats::na.omit(c(
        if (inherits(null_fit, "error")) paste0("null: ", conditionMessage(null_fit)),
        if (inherits(linear_fit, "error")) paste0("linear: ", conditionMessage(linear_fit)),
        if (inherits(spline_fit, "error")) paste0("spline: ", conditionMessage(spline_fit))
      )),
      collapse = "; "
    )
    empty_summary$skip_reason <- reason
    empty_summary$omitted_terms <- paste(term_info$omitted, collapse = "; ")
    return(list(summary = empty_summary, curve = tibble::tibble()))
  }

  temp_grid <- seq(
    stats::quantile(base$tmax_mean_c, 0.01, na.rm = TRUE, names = FALSE),
    stats::quantile(base$tmax_mean_c, 0.99, na.rm = TRUE, names = FALSE),
    length.out = 100
  )
  ref_temp <- stats::median(base$tmax_mean_c, na.rm = TRUE)
  newdata <- make_reference_newdata(base, temp_grid, term_info$keep)
  pred <- predict(spline_fit, newdata = newdata, type = "link", se.fit = TRUE)
  curve <- tibble::tibble(
    model = label,
    outcome = outcome,
    tmax_mean_c = temp_grid,
    predicted_probability = stats::plogis(pred$fit),
    predicted_probability_low = stats::plogis(pred$fit - 1.96 * pred$se.fit),
    predicted_probability_high = stats::plogis(pred$fit + 1.96 * pred$se.fit),
    reference_temp_c = ref_temp,
    spline_df = spline_df
  )

  linear_lrt <- tryCatch(stats::anova(null_fit, linear_fit, test = "LRT"), error = function(e) NULL)
  spline_lrt <- tryCatch(stats::anova(null_fit, spline_fit, test = "LRT"), error = function(e) NULL)
  nonlinear_lrt <- tryCatch(stats::anova(linear_fit, spline_fit, test = "LRT"), error = function(e) NULL)

  summary <- tibble::tibble(
    outcome = outcome,
    model = label,
    n = nrow(stats::model.frame(spline_fit)),
    events = sum(stats::model.frame(spline_fit)[[outcome]] == 1, na.rm = TRUE),
    spline_df = spline_df,
    reference_temp_c = ref_temp,
    min_temp_c = min(base$tmax_mean_c, na.rm = TRUE),
    max_temp_c = max(base$tmax_mean_c, na.rm = TRUE),
    linear_lrt_p_value = if (!is.null(linear_lrt)) linear_lrt$`Pr(>Chi)`[[2]] else NA_real_,
    spline_lrt_p_value = if (!is.null(spline_lrt)) spline_lrt$`Pr(>Chi)`[[2]] else NA_real_,
    nonlinear_lrt_p_value = if (!is.null(nonlinear_lrt)) nonlinear_lrt$`Pr(>Chi)`[[2]] else NA_real_,
    aic = stats::AIC(spline_fit),
    covariates = paste(term_info$keep, collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
  list(summary = summary, curve = curve)
}

run_phenotype_assignment_model <- function(df, label, adjust_terms = character(), spline_df = 3L) {
  keep_levels <- c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury")
  base <- df |>
    filter(
      .data$phenotype %in% keep_levels,
      !is.na(.data$tmax_mean_c), is.finite(.data$tmax_mean_c),
      !is.na(.data$rmax_mean_pct), is.finite(.data$rmax_mean_pct)
    ) |>
    mutate(phenotype = stats::relevel(factor(.data$phenotype, levels = keep_levels), ref = "regained_consciousness_extubated"))

  empty_summary <- tibble::tibble(
    model = label,
    n = nrow(base),
    spline_df = spline_df,
    reference_level = "regained_consciousness_extubated",
    reference_temp_c = NA_real_,
    reference_humidity_pct = NA_real_,
    min_temp_c = NA_real_,
    max_temp_c = NA_real_,
    min_humidity_pct = NA_real_,
    max_humidity_pct = NA_real_,
    temperature_lrt_p_value = NA_real_,
    humidity_lrt_p_value = NA_real_,
    temperature_nonlinear_lrt_p_value = NA_real_,
    aic = NA_real_,
    covariates = paste(c("splines::ns(tmax_mean_c, df = 3)", "splines::ns(rmax_mean_pct, df = 3)", adjust_terms), collapse = " + "),
    omitted_terms = NA_character_,
    estimable = FALSE,
    skip_reason = "Insufficient phenotype counts, sample size, temperature variation, or humidity variation"
  )
  if (!requireNamespace("nnet", quietly = TRUE)) {
    empty_summary$skip_reason <- "Package nnet unavailable"
    return(list(summary = empty_summary, curve = tibble::tibble()))
  }
  if (
    nrow(base) < 50 ||
      length(unique(base$phenotype)) < 3 ||
      any(table(base$phenotype) < 5) ||
      length(unique(base$tmax_mean_c)) < 5 ||
      length(unique(base$rmax_mean_pct)) < 5
  ) {
    return(list(summary = empty_summary, curve = tibble::tibble()))
  }

  term_info <- choose_terms(base, adjust_terms)
  adjust_txt <- if (length(term_info$keep) > 0) paste("+", paste(term_info$keep, collapse = " + ")) else ""
  covar_rhs <- if (length(term_info$keep) > 0) paste(term_info$keep, collapse = " + ") else "1"
  humidity_term <- paste0("splines::ns(rmax_mean_pct, df = ", spline_df, ")")
  temp_term <- paste0("splines::ns(tmax_mean_c, df = ", spline_df, ")")
  full_formula <- stats::as.formula(paste("phenotype ~", paste(c(temp_term, humidity_term, term_info$keep), collapse = " + ")))
  no_temp_formula <- stats::as.formula(paste("phenotype ~", paste(c(humidity_term, term_info$keep), collapse = " + ")))
  no_humidity_formula <- stats::as.formula(paste("phenotype ~", paste(c(temp_term, term_info$keep), collapse = " + ")))
  linear_temp_formula <- stats::as.formula(paste("phenotype ~", paste(c("tmax_per5c", humidity_term, term_info$keep), collapse = " + ")))

  full_fit <- tryCatch(nnet::multinom(full_formula, data = base, trace = FALSE), error = function(e) e)
  no_temp_fit <- tryCatch(nnet::multinom(no_temp_formula, data = base, trace = FALSE), error = function(e) e)
  no_humidity_fit <- tryCatch(nnet::multinom(no_humidity_formula, data = base, trace = FALSE), error = function(e) e)
  linear_temp_fit <- tryCatch(nnet::multinom(linear_temp_formula, data = base, trace = FALSE), error = function(e) e)
  if (inherits(full_fit, "error") || inherits(no_temp_fit, "error") || inherits(no_humidity_fit, "error") || inherits(linear_temp_fit, "error")) {
    reason <- paste(
      stats::na.omit(c(
        if (inherits(full_fit, "error")) paste0("full: ", conditionMessage(full_fit)),
        if (inherits(no_temp_fit, "error")) paste0("no_temp: ", conditionMessage(no_temp_fit)),
        if (inherits(no_humidity_fit, "error")) paste0("no_humidity: ", conditionMessage(no_humidity_fit)),
        if (inherits(linear_temp_fit, "error")) paste0("linear_temp: ", conditionMessage(linear_temp_fit))
      )),
      collapse = "; "
    )
    empty_summary$skip_reason <- reason
    empty_summary$omitted_terms <- paste(term_info$omitted, collapse = "; ")
    return(list(summary = empty_summary, curve = tibble::tibble()))
  }

  lrt_p <- function(reduced, full) {
    ll_reduced <- stats::logLik(reduced)
    ll_full <- stats::logLik(full)
    chi <- 2 * (as.numeric(ll_full) - as.numeric(ll_reduced))
    df_diff <- attr(ll_full, "df") - attr(ll_reduced, "df")
    if (!is.finite(chi) || !is.finite(df_diff) || df_diff <= 0) return(NA_real_)
    stats::pchisq(chi, df = df_diff, lower.tail = FALSE)
  }

  temp_grid <- seq(
    stats::quantile(base$tmax_mean_c, 0.01, na.rm = TRUE, names = FALSE),
    stats::quantile(base$tmax_mean_c, 0.99, na.rm = TRUE, names = FALSE),
    length.out = 100
  )
  ref_temp <- stats::median(base$tmax_mean_c, na.rm = TRUE)
  ref_humidity <- stats::median(base$rmax_mean_pct, na.rm = TRUE)
  newdata <- make_reference_newdata(base, temp_grid, term_info$keep) |>
    mutate(rmax_mean_pct = ref_humidity)
  pred <- as.data.frame(predict(full_fit, newdata = newdata, type = "probs"), stringsAsFactors = FALSE)
  if (is.null(dim(pred))) pred <- as.data.frame(t(pred), stringsAsFactors = FALSE)
  pred$tmax_mean_c <- temp_grid
  pred$rmax_mean_pct <- ref_humidity
  curve <- pred |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(keep_levels),
      names_to = "phenotype",
      values_to = "predicted_probability"
    ) |>
    mutate(
      model = label,
      reference_temp_c = ref_temp,
      reference_humidity_pct = ref_humidity,
      spline_df = spline_df,
      .before = 1
    )

  summary <- tibble::tibble(
    model = label,
    n = nrow(stats::model.frame(full_fit)),
    spline_df = spline_df,
    reference_level = "regained_consciousness_extubated",
    reference_temp_c = ref_temp,
    reference_humidity_pct = ref_humidity,
    min_temp_c = min(base$tmax_mean_c, na.rm = TRUE),
    max_temp_c = max(base$tmax_mean_c, na.rm = TRUE),
    min_humidity_pct = min(base$rmax_mean_pct, na.rm = TRUE),
    max_humidity_pct = max(base$rmax_mean_pct, na.rm = TRUE),
    temperature_lrt_p_value = lrt_p(no_temp_fit, full_fit),
    humidity_lrt_p_value = lrt_p(no_humidity_fit, full_fit),
    temperature_nonlinear_lrt_p_value = lrt_p(linear_temp_fit, full_fit),
    aic = stats::AIC(full_fit),
    covariates = paste(c(temp_term, humidity_term, term_info$keep), collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = TRUE,
    skip_reason = NA_character_
  )
  list(summary = summary, curve = curve)
}

ohca <- readr::read_csv(OHCA_PATH, show_col_types = FALSE) |>
  mutate(
    hospitalization_id = as.character(.data$hospitalization_id),
    patient_id = as.character(.data$patient_id),
    county_fips = normalize_county_fips(.data$county_fips),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    first_icu_in = as_utc_datetime(.data$first_icu_in),
    last_icu_out = as_utc_datetime(.data$last_icu_out),
    admission_date = as.Date(.data$admission_date),
    admit_year = lubridate::year(.data$admission_date),
    age_group = ifelse(as.numeric(.data$age_at_admission) >= 65, ">=65", "<65"),
    race_group = ifelse(is_black_race(.data$race_category), "Black", "Non-Black"),
    sex_group = dplyr::case_when(
      is_male(.data$sex_category) ~ "Male",
      is_female(.data$sex_category) ~ "Female",
      TRUE ~ "Other/Unknown"
    ),
    hospital_death = ifelse(.data$hospital_death == 1L | is_expired_discharge(.data$discharge_category), 1L, 0L)
  )

tmax_raw <- arrow::read_parquet(TMAX_PATH)
tmax_county_col <- if ("county_fips" %in% names(tmax_raw)) "county_fips" else if ("geoid" %in% names(tmax_raw)) "geoid" else NA_character_
if (is.na(tmax_county_col)) stop("Tmax file must contain county_fips or geoid.", call. = FALSE)
tmax <- tmax_raw |>
  transmute(
    county_fips = normalize_county_fips(.data[[tmax_county_col]]),
    admission_date = as.Date(.data$date),
    tmax_mean_c = suppressWarnings(as.numeric(.data$tmax_mean_c))
  )

rmax_raw <- arrow::read_parquet(RMAX_PATH)
rmax_county_col <- if ("county_fips" %in% names(rmax_raw)) "county_fips" else if ("geoid" %in% names(rmax_raw)) "geoid" else NA_character_
if (is.na(rmax_county_col)) stop("Rmax file must contain county_fips or geoid.", call. = FALSE)
rmax <- rmax_raw |>
  transmute(
    county_fips = normalize_county_fips(.data[[rmax_county_col]]),
    admission_date = as.Date(.data$date),
    rmax_mean_pct = suppressWarnings(as.numeric(.data$rmax_mean_pct))
  )

cohort <- ohca |>
  left_join(tmax, by = c("county_fips", "admission_date")) |>
  left_join(rmax, by = c("county_fips", "admission_date")) |>
  filter(!is.na(.data$tmax_mean_c), !is.na(.data$rmax_mean_pct), !is.na(.data$first_icu_in)) |>
  mutate(
    tmax_per5c = .data$tmax_mean_c / 5,
    rmax_per10pct = .data$rmax_mean_pct / 10,
    death_within_72h = !is.na(.data$discharge_dttm) &
      .data$hospital_death == 1L &
      as.numeric(difftime(.data$discharge_dttm, .data$first_icu_in, units = "hours")) <= WINDOW_HOURS
  )

cohort_ids <- unique(cohort$hospitalization_id)

diagnosis_source <- NULL
for (candidate in c("hospital_diagnosis", "admission_diagnosis")) {
  if (file.exists(file.path(tables_path, sprintf("clif_%s.%s", candidate, file_type)))) {
    diagnosis_source <- candidate
    break
  }
}
if (is.null(diagnosis_source)) stop("Could not find clif_hospital_diagnosis or clif_admission_diagnosis.")

diagnosis_raw <- read_clif_table(tables_path, file_type, diagnosis_source, required = TRUE)
ohca_mechanism <- derive_ohca_mechanism(diagnosis_raw, cohort_ids, diagnosis_source)
diagnosis <- diagnosis_raw |>
  transmute(
    hospitalization_id = as.character(.data$hospitalization_id),
    diagnosis_code = if ("diagnosis_code" %in% names(diagnosis_raw)) as.character(.data$diagnosis_code) else NA_character_,
    diagnosis_code_format = if ("diagnosis_code_format" %in% names(diagnosis_raw)) as.character(.data$diagnosis_code_format) else NA_character_
  ) |>
  filter(.data$hospitalization_id %in% cohort_ids) |>
  mutate(
    diagnosis_code_clean = norm_code(.data$diagnosis_code),
    diagnosis_code_format = infer_diagnosis_format(.data$diagnosis_code_clean, .data$diagnosis_code_format),
    anoxic_brain_injury_dx = (
      stringr::str_detect(.data$diagnosis_code_format, "10") & stringr::str_starts(.data$diagnosis_code_clean, "G931")
    ) | (
      stringr::str_detect(.data$diagnosis_code_format, "9") & stringr::str_starts(.data$diagnosis_code_clean, "3481")
    ),
    brain_death_dx = stringr::str_detect(.data$diagnosis_code_format, "10") & stringr::str_starts(.data$diagnosis_code_clean, "G9382")
  )

neuro_dx <- diagnosis |>
  group_by(.data$hospitalization_id) |>
  summarise(
    anoxic_brain_injury_dx = any(.data$anoxic_brain_injury_dx, na.rm = TRUE),
    brain_death_dx = any(.data$brain_death_dx, na.rm = TRUE),
    anoxic_brain_injury_codes = paste(sort(unique(.data$diagnosis_code_clean[.data$anoxic_brain_injury_dx | .data$brain_death_dx])), collapse = " | "),
    .groups = "drop"
  )

respiratory <- read_clif_table(tables_path, file_type, "respiratory_support", required = FALSE)
if (is.null(respiratory) || nrow(respiratory) == 0) {
  resp_summary <- tibble::tibble(hospitalization_id = cohort_ids)
} else {
  resp_summary <- respiratory |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      device_category = stringr::str_to_upper(tidyr::replace_na(as.character(.data$device_category), "")),
      artificial_airway = stringr::str_to_upper(tidyr::replace_na(as.character(.data$artificial_airway), "")),
      tracheostomy = stringr::str_to_upper(tidyr::replace_na(as.character(.data$tracheostomy), ""))
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    inner_join(cohort |> select("hospitalization_id", "first_icu_in"), by = "hospitalization_id") |>
    mutate(
      icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")),
      imv = .data$device_category %in% c("IMV", "VENT") |
        .data$artificial_airway %in% c("ETT", "NASAL ETT", "TRACH") |
        .data$tracheostomy %in% c("1", "TRUE", "YES")
    ) |>
    filter(.data$icu_hour >= 0, .data$icu_hour <= WINDOW_HOURS) |>
    arrange(.data$hospitalization_id, .data$recorded_dttm) |>
    group_by(.data$hospitalization_id) |>
    summarise(
      any_imv_0_72h = any(.data$imv, na.rm = TRUE),
      any_imv_48_72h = any(.data$imv & .data$icu_hour >= EXTUBATION_ASSESSMENT_START_HOUR, na.rm = TRUE),
      any_non_imv_after_imv = {
        imv_times <- .data$icu_hour[.data$imv]
        if (length(imv_times) == 0L) FALSE else any(!.data$imv & .data$icu_hour > max(imv_times, na.rm = TRUE), na.rm = TRUE)
      },
      last_resp_hour_0_72h = max(.data$icu_hour, na.rm = TRUE),
      last_resp_device_0_72h = safe_last(.data$device_category),
      last_resp_imv_0_72h = safe_last(.data$imv),
      .groups = "drop"
    )
}

assessments <- read_clif_table(tables_path, file_type, "patient_assessments", required = FALSE)
if (is.null(assessments) || nrow(assessments) == 0) {
  neuro_summary <- tibble::tibble(hospitalization_id = cohort_ids)
} else {
  neuro_summary <- assessments |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      assessment_category = as.character(.data$assessment_category),
      numerical_value = suppressWarnings(as.numeric(.data$numerical_value)),
      categorical_value = clean_binary_text(.data$categorical_value),
      text_value = clean_binary_text(.data$text_value)
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    inner_join(cohort |> select("hospitalization_id", "first_icu_in"), by = "hospitalization_id") |>
    mutate(
      assessment_category_clean = stringr::str_to_lower(.data$assessment_category),
      icu_hour = as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours"))
    ) |>
    filter(.data$icu_hour >= 0, .data$icu_hour <= WINDOW_HOURS) |>
    group_by(.data$hospitalization_id) |>
    summarise(
      n_neuro_assessments_0_72h = sum(.data$assessment_category_clean %in% c("gcs_total", "gcs_motor", "rass", "avpu"), na.rm = TRUE),
      min_gcs_0_72h = suppressWarnings(min(.data$numerical_value[.data$assessment_category_clean == "gcs_total" & .data$numerical_value >= 3 & .data$numerical_value <= 15], na.rm = TRUE)),
      best_gcs_24_72h = suppressWarnings(max(.data$numerical_value[.data$assessment_category_clean == "gcs_total" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 3 & .data$numerical_value <= 15], na.rm = TRUE)),
      last_gcs_24_72h = safe_last(.data$numerical_value[.data$assessment_category_clean == "gcs_total" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 3 & .data$numerical_value <= 15]),
      min_gcs_motor_0_72h = suppressWarnings(min(.data$numerical_value[.data$assessment_category_clean == "gcs_motor" & .data$numerical_value >= 1 & .data$numerical_value <= 6], na.rm = TRUE)),
      best_gcs_motor_24_72h = suppressWarnings(max(.data$numerical_value[.data$assessment_category_clean == "gcs_motor" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 1 & .data$numerical_value <= 6], na.rm = TRUE)),
      last_gcs_motor_24_72h = safe_last(.data$numerical_value[.data$assessment_category_clean == "gcs_motor" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= 1 & .data$numerical_value <= 6]),
      min_rass_0_72h = suppressWarnings(min(.data$numerical_value[.data$assessment_category_clean == "rass" & .data$numerical_value >= -5 & .data$numerical_value <= 4], na.rm = TRUE)),
      best_rass_24_72h = suppressWarnings(max(.data$numerical_value[.data$assessment_category_clean == "rass" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= -5 & .data$numerical_value <= 4], na.rm = TRUE)),
      last_rass_24_72h = safe_last(.data$numerical_value[.data$assessment_category_clean == "rass" & .data$icu_hour >= LATE_START_HOUR & .data$numerical_value >= -5 & .data$numerical_value <= 4]),
      any_avpu_alert_0_72h = any(.data$assessment_category_clean == "avpu" & stringr::str_detect(paste(.data$categorical_value, .data$text_value), "alert"), na.rm = TRUE),
      any_avpu_unresponsive_0_72h = any(.data$assessment_category_clean == "avpu" & stringr::str_detect(paste(.data$categorical_value, .data$text_value), "unresponsive|pain"), na.rm = TRUE),
      sat_pass_0_72h = any(.data$assessment_category_clean == "sat_delivery_pass_fail" & .data$categorical_value == "pass", na.rm = TRUE),
      sbt_pass_0_72h = any(.data$assessment_category_clean == "sbt_delivery_pass_fail" & .data$categorical_value == "pass", na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      across(c(min_gcs_0_72h, best_gcs_24_72h, min_gcs_motor_0_72h, best_gcs_motor_24_72h, min_rass_0_72h, best_rass_24_72h), ~ ifelse(is.infinite(.x), NA_real_, .x))
    )
}

phenotype_cohort <- cohort |>
  left_join(ohca_mechanism, by = "hospitalization_id") |>
  left_join(neuro_dx, by = "hospitalization_id") |>
  left_join(resp_summary, by = "hospitalization_id") |>
  left_join(neuro_summary, by = "hospitalization_id") |>
  mutate(
    ohca_mechanism = tidyr::replace_na(.data$ohca_mechanism, "unclear_other"),
    anoxic_brain_injury_dx = yesno(.data$anoxic_brain_injury_dx),
    brain_death_dx = yesno(.data$brain_death_dx),
    any_imv_0_72h = yesno(.data$any_imv_0_72h),
    any_imv_48_72h = yesno(.data$any_imv_48_72h),
    any_non_imv_after_imv = yesno(.data$any_non_imv_after_imv),
    last_resp_imv_0_72h = yesno(.data$last_resp_imv_0_72h),
    n_neuro_assessments_0_72h = tidyr::replace_na(.data$n_neuro_assessments_0_72h, 0L),
    severe_neuro_signal = (
      (!is.na(.data$min_gcs_0_72h) & .data$min_gcs_0_72h <= 5) |
        (!is.na(.data$min_gcs_motor_0_72h) & .data$min_gcs_motor_0_72h <= 2) |
        (!is.na(.data$min_rass_0_72h) & .data$min_rass_0_72h <= -4)
    ),
    awake_signal = (
      (!is.na(.data$best_gcs_24_72h) & .data$best_gcs_24_72h >= 13) |
        (!is.na(.data$last_gcs_24_72h) & .data$last_gcs_24_72h >= 13) |
        (!is.na(.data$best_gcs_motor_24_72h) & .data$best_gcs_motor_24_72h >= 6) |
        (!is.na(.data$best_rass_24_72h) & .data$best_rass_24_72h >= -1) |
        yesno(.data$any_avpu_alert_0_72h) |
        yesno(.data$sat_pass_0_72h) |
        yesno(.data$sbt_pass_0_72h)
    ),
    impaired_neuro_signal = (
      (!is.na(.data$last_gcs_24_72h) & .data$last_gcs_24_72h < 13) |
        (!is.na(.data$min_gcs_0_72h) & .data$min_gcs_0_72h <= 8) |
        (!is.na(.data$last_gcs_motor_24_72h) & .data$last_gcs_motor_24_72h < 6) |
        (!is.na(.data$last_rass_24_72h) & .data$last_rass_24_72h <= -2) |
        yesno(.data$any_avpu_unresponsive_0_72h)
    ),
    extubated_by_72h = .data$any_imv_0_72h & (
      .data$any_non_imv_after_imv |
        (!.data$any_imv_48_72h & !.data$last_resp_imv_0_72h)
    ),
    phenotype = dplyr::case_when(
      .data$death_within_72h ~ "anoxic_brain_injury",
      .data$extubated_by_72h & .data$awake_signal ~ "regained_consciousness_extubated",
      .data$severe_neuro_signal | .data$impaired_neuro_signal | .data$any_imv_48_72h | (.data$any_imv_0_72h & !.data$extubated_by_72h) ~ "limited_brain_function",
      TRUE ~ "unclassified"
    ),
    phenotype = factor(.data$phenotype, levels = c(
      "regained_consciousness_extubated",
      "limited_brain_function",
      "anoxic_brain_injury",
      "unclassified"
    )),
    anoxic_outcome = as.integer(.data$phenotype == "anoxic_brain_injury"),
    limited_vs_regained = dplyr::case_when(
      .data$phenotype == "limited_brain_function" ~ 1L,
      .data$phenotype == "regained_consciousness_extubated" ~ 0L,
      TRUE ~ NA_integer_
    ),
    extubated_awake_outcome = as.integer(.data$phenotype == "regained_consciousness_extubated")
  )

phenotype_summary <- phenotype_cohort |>
  count(.data$phenotype, name = "n") |>
  mutate(
    site_name = site_name,
    pct = 100 * .data$n / sum(.data$n),
    .before = 1
  )

mechanism_summary <- phenotype_cohort |>
  count(.data$ohca_mechanism, .data$phenotype, name = "n") |>
  group_by(.data$ohca_mechanism) |>
  mutate(pct_within_mechanism = 100 * .data$n / sum(.data$n)) |>
  ungroup() |>
  group_by(.data$phenotype) |>
  mutate(pct_within_phenotype = 100 * .data$n / sum(.data$n)) |>
  ungroup() |>
  mutate(site_name = site_name, .before = 1)

evidence_summary <- phenotype_cohort |>
  group_by(.data$phenotype) |>
  summarise(
    n = n(),
    anoxic_dx_pct = safe_pct(sum(.data$anoxic_brain_injury_dx, na.rm = TRUE), n()),
    brain_death_dx_pct = safe_pct(sum(.data$brain_death_dx, na.rm = TRUE), n()),
    severe_neuro_signal_pct = safe_pct(sum(.data$severe_neuro_signal, na.rm = TRUE), n()),
    awake_signal_pct = safe_pct(sum(.data$awake_signal, na.rm = TRUE), n()),
    impaired_neuro_signal_pct = safe_pct(sum(.data$impaired_neuro_signal, na.rm = TRUE), n()),
    extubated_by_72h_pct = safe_pct(sum(.data$extubated_by_72h, na.rm = TRUE), n()),
    imv_48_72h_pct = safe_pct(sum(.data$any_imv_48_72h, na.rm = TRUE), n()),
    hospital_death_pct = safe_pct(sum(.data$hospital_death == 1, na.rm = TRUE), n()),
    median_tmax_c = median(.data$tmax_mean_c, na.rm = TRUE),
    median_rmax_pct = median(.data$rmax_mean_pct, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(site_name = site_name, .before = 1)

base_adjust <- c("age_at_admission", "sex_group", "race_group", "admit_year")
mechanism_adjust <- c(base_adjust, "ohca_mechanism")

phenotype_assignment_primary <- run_phenotype_assignment_model(
  phenotype_cohort,
  "phenotype_assignment_temperature_humidity_demographics",
  base_adjust
)
phenotype_assignment_model <- phenotype_assignment_primary$summary |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_curve <- phenotype_assignment_primary$curve |>
  mutate(site_name = site_name, .before = 1)

phenotype_assignment_mechanism <- run_phenotype_assignment_model(
  phenotype_cohort,
  "phenotype_assignment_temperature_humidity_demographics_mechanism",
  mechanism_adjust
)
phenotype_assignment_model_mechanism <- phenotype_assignment_mechanism$summary |>
  mutate(site_name = site_name, .before = 1)
phenotype_assignment_curve_mechanism <- phenotype_assignment_mechanism$curve |>
  mutate(site_name = site_name, .before = 1)

phenotype_definitions <- tibble::tibble(
  phenotype = c(
    "anoxic_brain_injury",
    "regained_consciousness_extubated",
    "limited_brain_function",
    "unclassified"
  ),
  hierarchy_order = c(1L, 2L, 3L, 4L),
  definition = c(
    "Death within 72 hours of ICU entry. No diagnosis-code criterion is applied.",
    "Not anoxic phenotype; invasive ventilation occurred in first 72 ICU hours; evidence of extubation/no ongoing IMV by 72h; and awake signal from GCS >=13, GCS motor 6, RASS >=-1, AVPU alert, SAT pass, or SBT pass.",
    "Not anoxic or extubated-awake phenotype; severe or impaired neurologic signal in first 72 ICU hours, ongoing IMV at 48-72h, or IMV in first 72h without evidence of extubation.",
    "Does not meet structured-data criteria for the three primary phenotypes."
  )
) |>
  mutate(site_name = site_name, .before = 1)

primary_phenotype_levels <- c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury")
phenotype_analysis_cohort <- phenotype_cohort |>
  filter(.data$phenotype %in% primary_phenotype_levels) |>
  mutate(phenotype = factor(as.character(.data$phenotype), levels = primary_phenotype_levels))

ohca_summary_raw <- if (file.exists(OHCA_SUMMARY_PATH)) readr::read_csv(OHCA_SUMMARY_PATH, show_col_types = FALSE) else tibble::tibble()
summary_value <- function(col) {
  if (col %in% names(ohca_summary_raw) && nrow(ohca_summary_raw) > 0) suppressWarnings(as.numeric(ohca_summary_raw[[col]][[1]])) else NA_real_
}

consort_counts <- tibble::tibble(
  step_order = seq_len(8),
  step = c(
    "All ICU admissions in CLIF, 2018-2024",
    "OHCA present-on-admission ICU admissions before pathway/timing restriction",
    "OHCA ICU admissions with allowed ED/procedure/direct ICU pathway and ICU entry <24h",
    "OHCA ICU admissions with admission-day county Tmax and humidity",
    "Classified into one of three primary 72h phenotypes",
    "Regained consciousness and extubated",
    "Limited brain function",
    "Anoxic/severe group: death within 72h"
  ),
  n = c(
    summary_value("n_all_icu_admissions"),
    summary_value("n_ohca_admissions_before_pathway_timing_restriction"),
    summary_value("n_ohca_admissions"),
    nrow(phenotype_cohort),
    nrow(phenotype_analysis_cohort),
    sum(phenotype_analysis_cohort$phenotype == "regained_consciousness_extubated", na.rm = TRUE),
    sum(phenotype_analysis_cohort$phenotype == "limited_brain_function", na.rm = TRUE),
    sum(phenotype_analysis_cohort$phenotype == "anoxic_brain_injury", na.rm = TRUE)
  ),
  notes = c(
    "Source: clif_hospitalization/ADT ICU cohort summary",
    "POA cardiac arrest/OHCA-proxy diagnosis plus ICU admission",
    "Cohort definition used for all OHCA ICU heat analyses",
    "Required exposure fields for 72h phenotype and competing-risk models",
    "Excludes unclassified structured-data phenotype",
    "Not death within 72h; IMV in first 72h; extubated/no ongoing IMV by 72h; awake signal",
    "Not death within 72h or regained/extubated; impaired neurologic signal and/or ongoing IMV evidence",
    "Death within 72h of ICU entry"
  )
) |>
  mutate(
    site_name = site_name,
    n_excluded_since_previous = dplyr::lag(.data$n) - .data$n,
    n_excluded_since_previous = ifelse(.data$step_order == 1L | .data$step_order >= 6L, NA_real_, .data$n_excluded_since_previous),
    .before = 1
  )

table1 <- bind_rows(
  tibble::tibble(section = "Cohort", characteristic = "Three primary phenotype cohort", level = "", value = format(nrow(phenotype_analysis_cohort), big.mark = ",")),
  add_table1_continuous(phenotype_analysis_cohort, "Demographics", "Age, years, median [IQR]", "age_at_admission"),
  add_table1_categorical(phenotype_analysis_cohort, "Demographics", "Sex", "sex_category"),
  add_table1_categorical(phenotype_analysis_cohort, "Demographics", "Race", "race_category"),
  add_table1_categorical(phenotype_analysis_cohort, "Demographics", "Ethnicity", "ethnicity_category"),
  add_table1_continuous(phenotype_analysis_cohort, "Admission exposure", "Admission-day Tmax, C, median [IQR]", "tmax_mean_c"),
  add_table1_continuous(phenotype_analysis_cohort, "Admission exposure", "Admission-day relative humidity, %, median [IQR]", "rmax_mean_pct"),
  add_table1_categorical(phenotype_analysis_cohort, "OHCA mechanism", "POA diagnosis mechanism", "ohca_mechanism"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Any IMV in first 72h", "any_imv_0_72h"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "IMV evidence at ICU hours 48-72", "any_imv_48_72h"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Extubated/no ongoing IMV by 72h", "extubated_by_72h"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Awake signal", "awake_signal"),
  add_table1_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Impaired neurologic signal", "impaired_neuro_signal"),
  add_table1_binary(phenotype_analysis_cohort, "Outcome", "Hospital death", "hospital_death"),
  add_table1_binary(phenotype_analysis_cohort, "Outcome", "Death within 72h", "death_within_72h")
) |>
  mutate(site_name = site_name, .before = 1)

table2 <- bind_rows(
  tibble::tibble(section = "Cohort", characteristic = "Phenotype cohort size", level = "") |>
    bind_cols(as_tibble(as.list(setNames(vapply(primary_phenotype_levels, function(ph) format(sum(phenotype_analysis_cohort$phenotype == ph), big.mark = ","), character(1)), primary_phenotype_levels)))),
  add_table2_continuous(phenotype_analysis_cohort, "Demographics", "Age, years, median [IQR]", "age_at_admission", primary_phenotype_levels),
  add_table2_categorical(phenotype_analysis_cohort, "Demographics", "Sex", "sex_category", primary_phenotype_levels),
  add_table2_categorical(phenotype_analysis_cohort, "Demographics", "Race", "race_category", primary_phenotype_levels),
  add_table2_categorical(phenotype_analysis_cohort, "Demographics", "Ethnicity", "ethnicity_category", primary_phenotype_levels),
  add_table2_continuous(phenotype_analysis_cohort, "Admission exposure", "Admission-day Tmax, C, median [IQR]", "tmax_mean_c", primary_phenotype_levels),
  add_table2_continuous(phenotype_analysis_cohort, "Admission exposure", "Admission-day relative humidity, %, median [IQR]", "rmax_mean_pct", primary_phenotype_levels),
  add_table2_categorical(phenotype_analysis_cohort, "OHCA mechanism", "POA diagnosis mechanism", "ohca_mechanism", primary_phenotype_levels),
  add_table2_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Any IMV in first 72h", "any_imv_0_72h", primary_phenotype_levels),
  add_table2_binary(phenotype_analysis_cohort, "72h phenotype evidence", "IMV evidence at ICU hours 48-72", "any_imv_48_72h", primary_phenotype_levels),
  add_table2_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Extubated/no ongoing IMV by 72h", "extubated_by_72h", primary_phenotype_levels),
  add_table2_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Awake signal", "awake_signal", primary_phenotype_levels),
  add_table2_binary(phenotype_analysis_cohort, "72h phenotype evidence", "Impaired neurologic signal", "impaired_neuro_signal", primary_phenotype_levels),
  add_table2_binary(phenotype_analysis_cohort, "Outcome", "Hospital death", "hospital_death", primary_phenotype_levels),
  add_table2_binary(phenotype_analysis_cohort, "Outcome", "Death within 72h", "death_within_72h", primary_phenotype_levels)
) |>
  mutate(site_name = site_name, .before = 1)

readr::write_csv(
  phenotype_cohort |> select(
    "hospitalization_id", "patient_id", "admission_date", "tmax_mean_c", "rmax_mean_pct", "phenotype",
    "ohca_mechanism", "ohca_mechanism_detail",
    "anoxic_brain_injury_dx", "brain_death_dx", "severe_neuro_signal", "awake_signal",
    "impaired_neuro_signal", "extubated_by_72h", "any_imv_0_72h", "any_imv_48_72h",
    "min_gcs_0_72h", "best_gcs_24_72h", "last_gcs_24_72h", "best_rass_24_72h",
    "last_rass_24_72h", "hospital_death", "death_within_72h"
  ),
  file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_patient_audit.csv")
)
readr::write_csv(phenotype_summary, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_summary.csv"))
readr::write_csv(consort_counts, file.path(OUTPUT_DIR, "ohca_icu_72h_consort_flow.csv"))
readr::write_csv(table1, file.path(OUTPUT_DIR, "ohca_icu_72h_table1.csv"))
readr::write_csv(table2, file.path(OUTPUT_DIR, "ohca_icu_72h_table2_by_phenotype.csv"))
readr::write_csv(mechanism_summary, file.path(OUTPUT_DIR, "ohca_icu_72h_ohca_mechanism_summary.csv"))
readr::write_csv(evidence_summary, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_evidence_summary.csv"))
readr::write_csv(phenotype_assignment_model, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_model.csv"))
readr::write_csv(phenotype_assignment_curve, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_temperature_curve.csv"))
readr::write_csv(phenotype_assignment_model_mechanism, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted.csv"))
readr::write_csv(phenotype_assignment_curve_mechanism, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted.csv"))
readr::write_csv(phenotype_definitions, file.path(OUTPUT_DIR, "ohca_icu_72h_phenotype_definitions.csv"))

plot_df <- phenotype_summary |>
  filter(.data$phenotype != "unclassified") |>
  mutate(phenotype = factor(.data$phenotype, levels = rev(c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury"))))

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = .data$n, y = .data$phenotype, fill = .data$phenotype)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%s (%.1f%%)", .data$n, .data$pct)), hjust = -0.05, size = 3.5) +
    scale_fill_manual(values = c(
      "regained_consciousness_extubated" = "#2A9D8F",
      "limited_brain_function" = "#E9C46A",
      "anoxic_brain_injury" = "#C44536"
    ), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "OHCA ICU 72-hour Structured Phenotypes",
      x = "Patients",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIGURE_DIR, "figure_ohca_icu_72h_phenotype_counts.png"), p, width = 8.5, height = 4.8, dpi = 300)
}

if (nrow(phenotype_assignment_curve) > 0 && isTRUE(phenotype_assignment_model$estimable[[1]])) {
  curve_plot <- phenotype_assignment_curve |>
    mutate(
      phenotype = factor(
        .data$phenotype,
        levels = c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury"),
        labels = c("Regained/extubated", "Limited brain function", "Death within 72h")
      )
    )

  p_curve <- ggplot(curve_plot, aes(x = .data$tmax_mean_c, y = .data$predicted_probability, color = .data$phenotype)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(
      "Regained/extubated" = "#2A9D8F",
      "Limited brain function" = "#E9C46A",
      "Death within 72h" = "#C44536"
    )) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(
      title = "Adjusted OHCA ICU Phenotype Probability by Admission Temperature",
      x = "Admission-day assigned-county Tmax, C",
      y = "Predicted phenotype probability",
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  ggsave(file.path(FIGURE_DIR, "figure_ohca_icu_72h_phenotype_temperature_curves.png"), p_curve, width = 8.5, height = 5.2, dpi = 300)
}

print(phenotype_summary)
print(phenotype_assignment_model)
print(phenotype_assignment_model_mechanism)
message("Wrote OHCA ICU 72-hour phenotype outputs to ", OUTPUT_DIR)

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
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

OHCA_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
TMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
OUTPUT_DIR <- file.path(repo_root, "output", "final", "ohca_icu_phenotypes")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

unlink(file.path(OUTPUT_DIR, c(
  "ohca_icu_imv12_discharge_outcome_summary.csv",
  "ohca_icu_imv12_discharge_outcome_flow.csv",
  "ohca_icu_imv12_discharge_outcome_table_by_outcome.csv",
  "ohca_icu_imv12_discharge_outcome_model.csv",
  "ohca_icu_imv12_discharge_outcome_coefficients.csv",
  "ohca_icu_imv12_discharge_outcome_vcov.csv"
)))

config <- load_project_config(repo_root)
validate_site_name(config, repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)
site_name <- config$site_name %||% Sys.getenv("CLIF_SITE_NAME", unset = "unknown_site")

ensure_columns <- function(df, columns) {
  if (is.null(df)) df <- tibble::tibble()
  for (column in columns) {
    if (!column %in% names(df)) df[[column]] <- NA
  }
  df
}

is_hospice_discharge <- function(x) {
  stringr::str_detect(normalize_category(as.character(x)), "hospice|comfort")
}

is_home_discharge <- function(x) {
  normalize_category(as.character(x)) == "home"
}

is_snf_discharge <- function(x) {
  stringr::str_detect(normalize_category(as.character(x)), "skilled nursing|\\bsnf\\b")
}

is_ltach_discharge <- function(x) {
  stringr::str_detect(normalize_category(as.character(x)), "long term care hospital|ltach")
}

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
    ifelse(nchar(dx_clean) >= 3, "ICD9", NA_character_)
  )
  ifelse(nzchar(explicit), explicit, inferred)
}

format_n_pct <- function(x, denom) {
  fmt_n_pct(sum(x, na.rm = TRUE), denom)
}

format_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  fmt_median_iqr(x)
}

format_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  sprintf("%.1f [%.1f]", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
}

add_outcome_continuous <- function(df, section, characteristic, variable, outcome_levels) {
  row <- tibble::tibble(section = section, characteristic = characteristic, level = "Median [IQR]")
  for (outcome in outcome_levels) {
    x <- df[[variable]][df$descriptive_outcome == outcome]
    row[[outcome]] <- format_median_iqr(x)
  }
  row
}

add_outcome_mean_sd <- function(df, section, characteristic, variable, outcome_levels) {
  row <- tibble::tibble(section = section, characteristic = characteristic, level = "Mean [SD]")
  for (outcome in outcome_levels) {
    x <- df[[variable]][df$descriptive_outcome == outcome]
    row[[outcome]] <- format_mean_sd(x)
  }
  row
}

add_outcome_binary <- function(df, section, characteristic, variable, outcome_levels) {
  row <- tibble::tibble(section = section, characteristic = characteristic, level = "Yes")
  for (outcome in outcome_levels) {
    keep <- df$descriptive_outcome == outcome
    row[[outcome]] <- format_n_pct(df[[variable]][keep] %in% TRUE | df[[variable]][keep] == 1L, sum(keep, na.rm = TRUE))
  }
  row
}

add_outcome_categorical <- function(df, section, characteristic, variable, outcome_levels) {
  levels <- sort(unique(stats::na.omit(as.character(df[[variable]]))))
  if (length(levels) == 0L) return(tibble::tibble())
  rows <- lapply(levels, function(level_value) {
    row <- tibble::tibble(section = section, characteristic = characteristic, level = level_value)
    for (outcome in outcome_levels) {
      keep <- df$descriptive_outcome == outcome
      row[[outcome]] <- fmt_n_pct(sum(as.character(df[[variable]][keep]) == level_value, na.rm = TRUE), sum(keep, na.rm = TRUE))
    }
    row
  })
  dplyr::bind_rows(rows)
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
      if (sum(is.finite(x)) > 0 && length(unique(x[is.finite(x)])) >= 2) {
        keep <- c(keep, term)
      } else {
        omitted <- c(omitted, paste0(term, "_constant"))
      }
    } else {
      if (length(unique(stats::na.omit(x))) >= 2) {
        keep <- c(keep, term)
      } else {
        omitted <- c(omitted, paste0(term, "_one_level"))
      }
    }
  }
  list(keep = keep, omitted = omitted)
}

extract_multinom_coefficients <- function(fit, model, n, outcome_counts, covariates, omitted_terms) {
  if (inherits(fit, "error") || is.null(fit)) return(tibble::tibble())
  s <- summary(fit)
  coef_mat <- s$coefficients
  se_mat <- s$standard.errors
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, nrow = 1, dimnames = list(names(coef_mat)[[1]], names(coef_mat)))
    se_mat <- matrix(se_mat, nrow = 1, dimnames = dimnames(coef_mat))
  }
  purrr::map_dfr(rownames(coef_mat), function(outcome_level) {
    purrr::map_dfr(colnames(coef_mat), function(coefficient) {
      beta <- coef_mat[outcome_level, coefficient]
      se <- se_mat[outcome_level, coefficient]
      tibble::tibble(
        model = model,
        outcome_level = outcome_level,
        reference_level = "alive_home",
        n = n,
        outcome_n = unname(outcome_counts[[outcome_level]] %||% NA_integer_),
        coefficient = coefficient,
        estimate = as.numeric(beta),
        standard_error = as.numeric(se),
        odds_ratio = exp(as.numeric(beta)),
        ci_low = exp(as.numeric(beta) - 1.96 * as.numeric(se)),
        ci_high = exp(as.numeric(beta) + 1.96 * as.numeric(se)),
        p_value = 2 * stats::pnorm(abs(as.numeric(beta) / as.numeric(se)), lower.tail = FALSE),
        covariates = paste(covariates, collapse = " + "),
        omitted_terms = paste(omitted_terms, collapse = "; ")
      )
    })
  })
}

extract_multinom_vcov <- function(fit, model, n, covariates, omitted_terms) {
  if (inherits(fit, "error") || is.null(fit)) return(tibble::tibble())
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(vc) || length(vc) == 0L) return(tibble::tibble())

  parse_name <- function(x) {
    pieces <- strsplit(x, ":", fixed = TRUE)[[1]]
    if (length(pieces) < 2L) return(c(outcome = NA_character_, coefficient = x))
    c(outcome = pieces[[1]], coefficient = paste(pieces[-1], collapse = ":"))
  }
  row_info <- purrr::map_dfr(rownames(vc), function(x) {
    parsed <- parse_name(x)
    tibble::tibble(vcov_row = x, outcome_level_row = parsed[["outcome"]], coefficient_row = parsed[["coefficient"]])
  })
  col_info <- purrr::map_dfr(colnames(vc), function(x) {
    parsed <- parse_name(x)
    tibble::tibble(vcov_col = x, outcome_level_col = parsed[["outcome"]], coefficient_col = parsed[["coefficient"]])
  })

  as.data.frame(as.table(vc), stringsAsFactors = FALSE) |>
    transmute(
      vcov_row = as.character(.data$Var1),
      vcov_col = as.character(.data$Var2),
      covariance = as.numeric(.data$Freq)
    ) |>
    left_join(row_info, by = "vcov_row") |>
    left_join(col_info, by = "vcov_col") |>
    transmute(
      model = model,
      n = n,
      outcome_level_row,
      coefficient_row,
      outcome_level_col,
      coefficient_col,
      covariance,
      covariates = paste(covariates, collapse = " + "),
      omitted_terms = paste(omitted_terms, collapse = "; ")
    )
}

run_discharge_multinomial <- function(df) {
  model <- "imv12_discharge_outcome_temperature_demographics"
  outcome_levels <- c("alive_home", "death_hospice", "ltach_trach", "alive_snf")
  reference_level <- "alive_home"
  adjust_terms <- c("age_at_admission", "sex_group", "race_group")
  base <- df |>
    filter(.data$discharge_outcome %in% outcome_levels, !is.na(.data$tmax_per5c), is.finite(.data$tmax_per5c)) |>
    mutate(discharge_outcome = stats::relevel(factor(.data$discharge_outcome, levels = outcome_levels), ref = reference_level))
  counts <- table(base$discharge_outcome)
  term_info <- choose_terms(base, c("tmax_per5c", adjust_terms))
  empty_summary <- tibble::tibble(
    site_name = site_name,
    model = model,
    exposure = "Admission Tmax per 5 C",
    n = nrow(base),
    reference_level = reference_level,
    modeled_outcomes = paste(outcome_levels, collapse = "; "),
    covariates = paste(term_info$keep, collapse = " + "),
    omitted_terms = paste(term_info$omitted, collapse = "; "),
    estimable = FALSE,
    skip_reason = "Insufficient outcome counts, sample size, exposure variation, or nnet unavailable"
  )
  if (!requireNamespace("nnet", quietly = TRUE) ||
      nrow(base) < 50 ||
      length(unique(base$discharge_outcome)) < 4 ||
      any(counts[outcome_levels] < 5) ||
      !"tmax_per5c" %in% term_info$keep) {
    return(list(summary = empty_summary, coefficients = tibble::tibble(), vcov = tibble::tibble()))
  }

  formula_txt <- paste("discharge_outcome ~", paste(term_info$keep, collapse = " + "))
  fit <- tryCatch(
    nnet::multinom(stats::as.formula(formula_txt), data = base, trace = FALSE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    empty_summary$skip_reason <- conditionMessage(fit)
    return(list(summary = empty_summary, coefficients = tibble::tibble(), vcov = tibble::tibble()))
  }
  coefficients <- extract_multinom_coefficients(
    fit,
    model,
    nrow(stats::model.frame(fit)),
    table(stats::model.frame(fit)$discharge_outcome),
    term_info$keep,
    term_info$omitted
  )
  temp_rows <- coefficients |>
    filter(.data$coefficient == "tmax_per5c") |>
    transmute(
      site_name = site_name,
      model = model,
      exposure = "Admission Tmax per 5 C",
      outcome_level,
      reference_level,
      n,
      outcome_n,
      odds_ratio,
      ci_low,
      ci_high,
      p_value,
      covariates = paste(term_info$keep, collapse = " + "),
      omitted_terms = paste(term_info$omitted, collapse = "; "),
      estimable = TRUE,
      skip_reason = NA_character_
    )
  list(
    summary = temp_rows,
    coefficients = coefficients |> mutate(site_name = site_name, .before = 1),
    vcov = extract_multinom_vcov(fit, model, nrow(stats::model.frame(fit)), term_info$keep, term_info$omitted) |>
      mutate(site_name = site_name, .before = 1)
  )
}

ohca <- readr::read_csv(OHCA_PATH, show_col_types = FALSE) |>
  ensure_columns(c(
    "discharge_name", "imv_duration_hours", "imv_any", "imv_record_count",
    "first_imv_dttm", "last_imv_dttm", "ethnicity_category",
    "admission_dttm", "discharge_dttm", "first_icu_in", "last_icu_out"
  )) |>
  mutate(
    hospitalization_id = as.character(.data$hospitalization_id),
    county_fips = normalize_county_fips(.data$county_fips),
    admission_date = as.Date(.data$admission_date),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    first_icu_in = as_utc_datetime(.data$first_icu_in),
    last_icu_out = as_utc_datetime(.data$last_icu_out),
    discharge_category = as.character(.data$discharge_category),
    discharge_name = as.character(.data$discharge_name),
    imv_duration_hours = suppressWarnings(as.numeric(.data$imv_duration_hours)),
    imv_record_count = suppressWarnings(as.numeric(.data$imv_record_count)),
    first_imv_dttm = as_utc_datetime(.data$first_imv_dttm),
    last_imv_dttm = as_utc_datetime(.data$last_imv_dttm),
    race_group = ifelse(is_black_race(.data$race_category), "Black", "Non-Black"),
    sex_group = dplyr::case_when(
      is_male(.data$sex_category) ~ "Male",
      is_female(.data$sex_category) ~ "Female",
      TRUE ~ "Other/Unknown"
    ),
    hospital_death = ifelse(.data$hospital_death == 1L | is_expired_discharge(.data$discharge_category), 1L, 0L),
    hospice_discharge = is_hospice_discharge(.data$discharge_category) | is_hospice_discharge(.data$discharge_name),
    death_hospice_discharge = .data$hospital_death == 1L | .data$hospice_discharge
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

cohort <- ohca |>
  left_join(tmax, by = c("county_fips", "admission_date")) |>
  filter(!is.na(.data$tmax_mean_c)) |>
  mutate(tmax_per5c = .data$tmax_mean_c / 5)

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
    .groups = "drop"
  )

respiratory <- read_clif_table(
  tables_path,
  file_type,
  "respiratory_support",
  columns = c("hospitalization_id", "recorded_dttm", "tracheostomy"),
  required = FALSE
)

if (is.null(respiratory) || nrow(respiratory) == 0 || !all(c("hospitalization_id", "recorded_dttm") %in% names(respiratory))) {
  trach_summary <- tibble::tibble(hospitalization_id = cohort_ids, any_trach = FALSE, first_trach_dttm = as.POSIXct(NA, tz = "UTC"))
} else {
  respiratory <- ensure_columns(respiratory, c("tracheostomy"))
  trach_summary <- respiratory |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$recorded_dttm),
      tracheostomy = stringr::str_to_upper(tidyr::replace_na(as.character(.data$tracheostomy), ""))
    ) |>
    filter(.data$hospitalization_id %in% cohort_ids, !is.na(.data$recorded_dttm)) |>
    mutate(
      trach = .data$tracheostomy %in% c("1", "TRUE", "YES")
    ) |>
    group_by(.data$hospitalization_id) |>
    summarise(
      any_trach = any(.data$trach, na.rm = TRUE),
      first_trach_dttm = {
        trach_times <- .data$recorded_dttm[.data$trach %in% TRUE]
        if (length(trach_times) == 0L) as.POSIXct(NA, tz = "UTC") else min(trach_times, na.rm = TRUE)
      },
      .groups = "drop"
    )
}

analysis_base <- cohort |>
  left_join(neuro_dx, by = "hospitalization_id") |>
  left_join(trach_summary, by = "hospitalization_id") |>
  mutate(
    any_trach = tidyr::replace_na(.data$any_trach, FALSE),
    anoxic_brain_injury_dx = tidyr::replace_na(as.logical(.data$anoxic_brain_injury_dx), FALSE),
    brain_death_dx = tidyr::replace_na(as.logical(.data$brain_death_dx), FALSE),
    imv_12h_cohort = !is.na(.data$imv_duration_hours) & .data$imv_duration_hours >= 12,
    discharge_outcome = dplyr::case_when(
      .data$death_hospice_discharge ~ "death_hospice",
      is_ltach_discharge(.data$discharge_category) & .data$any_trach ~ "ltach_trach",
      !.data$death_hospice_discharge & is_home_discharge(.data$discharge_category) ~ "alive_home",
      !.data$death_hospice_discharge & is_snf_discharge(.data$discharge_category) ~ "alive_snf",
      TRUE ~ NA_character_
    ),
    hospital_los_days = as.numeric(difftime(.data$discharge_dttm, .data$admission_dttm, units = "days")),
    icu_entry_to_discharge_days = as.numeric(difftime(.data$discharge_dttm, .data$first_icu_in, units = "days")),
    icu_los_days = as.numeric(difftime(.data$last_icu_out, .data$first_icu_in, units = "days")),
    imv_start_to_discharge_days = as.numeric(difftime(.data$discharge_dttm, .data$first_imv_dttm, units = "days")),
    imv_stop_to_discharge_days = as.numeric(difftime(.data$discharge_dttm, .data$last_imv_dttm, units = "days")),
    time_to_trach_days = as.numeric(difftime(.data$first_trach_dttm, .data$first_icu_in, units = "days"))
  )

analysis <- analysis_base |>
  filter(.data$imv_12h_cohort) |>
  mutate(
    discharge_outcome = factor(
      .data$discharge_outcome,
      levels = c("alive_home", "death_hospice", "ltach_trach", "alive_snf")
    ),
    descriptive_outcome = factor(
      tidyr::replace_na(as.character(.data$discharge_outcome), "excluded_other_unknown"),
      levels = c("alive_home", "death_hospice", "ltach_trach", "alive_snf", "excluded_other_unknown")
    )
  )

flow <- tibble::tibble(
  site_name = site_name,
  step = c(
    "OHCA ICU admissions with admission Tmax",
    "IMV duration >=12 hours",
    "Classified as death/hospice",
    "Classified as LTACH with tracheostomy",
    "Classified as alive discharge home",
    "Classified as alive discharge SNF",
    "Excluded other/unknown discharge outcome",
    "Final modeled cohort"
  ),
  n = c(
    nrow(cohort),
    nrow(analysis),
    sum(analysis$discharge_outcome == "death_hospice", na.rm = TRUE),
    sum(analysis$discharge_outcome == "ltach_trach", na.rm = TRUE),
    sum(analysis$discharge_outcome == "alive_home", na.rm = TRUE),
    sum(analysis$discharge_outcome == "alive_snf", na.rm = TRUE),
    sum(is.na(analysis$discharge_outcome)),
    sum(!is.na(analysis$discharge_outcome))
  ),
  definition = "IMV-12 discharge-outcome cohort includes OHCA ICU admissions with imv_duration_hours >=12 based on first-to-last IMV respiratory-support timestamps from the OHCA cohort build. Final outcomes are death/hospice, LTACH discharge with any tracheostomy evidence, alive discharge home, or alive discharge SNF. Other discharge dispositions are excluded from the multinomial model."
)

summary_tbl <- analysis |>
  mutate(discharge_outcome = as.character(.data$discharge_outcome), discharge_outcome = tidyr::replace_na(.data$discharge_outcome, "excluded_other_unknown")) |>
  count(.data$discharge_outcome, name = "n") |>
  mutate(
    site_name = site_name,
    pct_among_imv12 = safe_pct(.data$n, nrow(analysis)),
    definition = "Death/hospice includes expired/death discharge or hospice/comfort discharge. LTACH with trach requires LTACH discharge category and any tracheostomy evidence in respiratory support. Home and SNF outcomes exclude death/hospice.",
    .before = 1
  )

descriptive_outcome_levels <- levels(analysis$descriptive_outcome)
descriptive_table <- bind_rows(
  add_outcome_continuous(analysis, "Demographics", "Age, years", "age_at_admission", descriptive_outcome_levels),
  add_outcome_categorical(analysis, "Demographics", "Sex", "sex_group", descriptive_outcome_levels),
  add_outcome_categorical(analysis, "Demographics", "Race", "race_group", descriptive_outcome_levels),
  add_outcome_categorical(analysis, "Demographics", "Ethnicity", "ethnicity_category", descriptive_outcome_levels),
  add_outcome_continuous(analysis, "Admission exposure", "Admission-day Tmax, C", "tmax_mean_c", descriptive_outcome_levels),
  add_outcome_continuous(analysis, "IMV exposure", "IMV duration, hours", "imv_duration_hours", descriptive_outcome_levels),
  add_outcome_continuous(analysis, "IMV exposure", "IMV record count", "imv_record_count", descriptive_outcome_levels),
  add_outcome_mean_sd(analysis, "Analytic metrics", "Hospital admission to discharge/death, days", "hospital_los_days", descriptive_outcome_levels),
  add_outcome_mean_sd(analysis, "Analytic metrics", "ICU entry to discharge/death, days", "icu_entry_to_discharge_days", descriptive_outcome_levels),
  add_outcome_mean_sd(analysis, "Analytic metrics", "ICU length of stay, days", "icu_los_days", descriptive_outcome_levels),
  add_outcome_mean_sd(analysis, "Analytic metrics", "First IMV to discharge/death, days", "imv_start_to_discharge_days", descriptive_outcome_levels),
  add_outcome_mean_sd(analysis, "Analytic metrics", "Final IMV record to discharge/death, days", "imv_stop_to_discharge_days", descriptive_outcome_levels),
  add_outcome_mean_sd(analysis, "Analytic metrics", "ICU entry to first tracheostomy evidence, days", "time_to_trach_days", descriptive_outcome_levels),
  add_outcome_binary(analysis, "Tracheostomy evidence", "Any tracheostomy evidence", "any_trach", descriptive_outcome_levels),
  add_outcome_binary(analysis, "ICD evidence", "Anoxic brain injury ICD", "anoxic_brain_injury_dx", descriptive_outcome_levels),
  add_outcome_binary(analysis, "ICD evidence", "Brain death ICD", "brain_death_dx", descriptive_outcome_levels),
  add_outcome_binary(analysis, "Outcome evidence", "Hospital death", "hospital_death", descriptive_outcome_levels),
  add_outcome_binary(analysis, "Outcome evidence", "Hospice/comfort discharge", "hospice_discharge", descriptive_outcome_levels),
  add_outcome_categorical(analysis, "Discharge disposition", "Original discharge category", "discharge_category", descriptive_outcome_levels)
) |>
  mutate(site_name = site_name, .before = 1)

model_results <- run_discharge_multinomial(analysis)

readr::write_csv(summary_tbl, file.path(OUTPUT_DIR, "ohca_icu_imv12_discharge_outcome_summary.csv"))
readr::write_csv(flow, file.path(OUTPUT_DIR, "ohca_icu_imv12_discharge_outcome_flow.csv"))
readr::write_csv(descriptive_table, file.path(OUTPUT_DIR, "ohca_icu_imv12_discharge_outcome_table_by_outcome.csv"))
readr::write_csv(model_results$summary, file.path(OUTPUT_DIR, "ohca_icu_imv12_discharge_outcome_model.csv"))
readr::write_csv(model_results$coefficients, file.path(OUTPUT_DIR, "ohca_icu_imv12_discharge_outcome_coefficients.csv"))
readr::write_csv(model_results$vcov, file.path(OUTPUT_DIR, "ohca_icu_imv12_discharge_outcome_vcov.csv"))

print(flow)
print(summary_tbl)
print(model_results$summary)

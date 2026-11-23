# First-day measured-component SOFA, using the same clinical cut points as
# the lung cancer CLIF analysis. Missing organ domains stay missing.
aim2_sofa_window <- function(tables_path, file_type, table, cohort, time, fields,
                             category = NULL, categories = NULL) {
  columns <- c("hospitalization_id", time, fields)
  path <- file.path(tables_path, sprintf("clif_%s.%s", table, file_type))
  if (!file.exists(path) || !nrow(cohort))
    return(tibble::as_tibble(setNames(rep(list(character()), length(columns)), columns)))

  ids <- cohort$hospitalization_id
  if (file_type == "parquet") {
    message("Reading first-day SOFA inputs from ", basename(path))
    source <- arrow::open_dataset(path, format = "parquet")
    if (all(c("hospitalization_id", time) %in% names(source))) {
      rows <- source |>
        filter(as.character(.data$hospitalization_id) %in% !!ids) |>
        select(any_of(columns)) |>
        collect()
      for (field in setdiff(columns, names(rows))) rows[[field]] <- rep(NA, nrow(rows))
    } else {
      rows <- tibble::as_tibble(setNames(rep(list(character()), length(columns)), columns))
    }
  } else {
    rows <- read_clif_table(tables_path, file_type, table, columns, required = FALSE)
  }
  if (!is.null(category))
    rows <- rows[normalize_category(rows[[category]]) %in% categories, , drop = FALSE]
  rows |>
    mutate(hospitalization_id = as.character(hospitalization_id),
           across(all_of(time), as_utc_datetime)) |>
    inner_join(select(cohort, hospitalization_id, first_icu_in), by = "hospitalization_id") |>
    filter(.data[[time]] >= first_icu_in, .data[[time]] < first_icu_in + 86400) |>
    select(all_of(columns))
}

aim2_sofa_dose <- function(dose, unit, weight) {
  unit <- gsub(" ", "", normalize_category(unit), fixed = TRUE)
  unit <- gsub("µg|μg|^ug", "mcg", unit)
  unit <- sub("/(hour|hours)$", "/hr", unit)
  case_when(
    unit == "mcg/kg/min" ~ dose,
    unit == "mcg/kg/hr" ~ dose / 60,
    unit == "mg/kg/min" ~ dose * 1000,
    unit == "mg/kg/hr" ~ dose * 1000 / 60,
    unit == "mcg/min" ~ dose / weight,
    unit == "mcg/hr" ~ dose / weight / 60,
    unit == "mg/min" ~ dose * 1000 / weight,
    unit == "mg/hr" ~ dose * 1000 / weight / 60
  )
}

aim2_sofa_respiratory <- function(labs, support) {
  device <- normalize_category(support$device_category)
  fio2 <- suppressWarnings(as.numeric(support$fio2_set))
  fio2 <- ifelse(fio2 > 1, fio2 / 100, fio2)
  fio2[which(is.na(fio2) & device %in% c("room air", "ra"))] <- 0.21
  eligible <- which(is.finite(fio2) & fio2 >= 0.21 & fio2 <= 1 &
                      !device %in% c("", "unknown", "missing", "other"))
  pao2 <- labs$lab_value_numeric
  blood_gas <- which(normalize_category(labs$lab_category) == "po2_arterial" &
                       is.finite(pao2) & pao2 >= 20 & pao2 <= 800)
  if (!length(eligible) || !length(blood_gas)) return(NA_real_)

  # Pair each gas to the closest usable FiO2 within one hour. Both rows are
  # already restricted to the pre-landmark window. Resolve ties toward earlier data.
  eligible <- eligible[order(support$recorded_dttm[eligible])]
  support_time <- as.numeric(support$recorded_dttm[eligible])
  gas_time <- as.numeric(labs$lab_result_dttm[blood_gas])
  before <- findInterval(gas_time, support_time)
  left <- pmax(1L, before)
  right <- pmin(length(eligible), before + 1L)
  nearest <- ifelse(abs(gas_time - support_time[left]) <= abs(gas_time - support_time[right]), left, right)
  paired <- abs(gas_time - support_time[nearest]) <= 3600
  use <- eligible[nearest[paired]]
  pf <- pao2[blood_gas[paired]] / fio2[use]
  ventilation <- device[use] %in% c("imv", "vent", "invasive mechanical ventilation", "nippv", "niv", "bipap", "cpap")
  aim2_sofa_max(case_when(pf < 100 & ventilation ~ 4, pf < 200 & ventilation ~ 3,
                          pf < 300 ~ 2, pf < 400 ~ 1, TRUE ~ 0))
}

aim2_sofa_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_

aim2_sofa_scores <- function(cohort, tables_path, file_type) {
  pressors <- c("dopamine", "dobutamine", "epinephrine", "norepinephrine")
  read_window <- function(table, time, fields, category = NULL, categories = NULL)
    aim2_sofa_window(tables_path, file_type, table, cohort, time, fields, category, categories)
  labs <- read_window("labs", "lab_result_dttm", c("lab_category", "lab_value_numeric"),
                      "lab_category", c("po2_arterial", "platelet_count", "bilirubin_total", "creatinine")) |>
    mutate(lab_value_numeric = suppressWarnings(as.numeric(lab_value_numeric)))
  vitals <- read_window("vitals", "recorded_dttm", c("vital_category", "vital_value"),
                        "vital_category", c("map", "weight_kg"))
  support <- read_window("respiratory_support", "recorded_dttm", c("device_category", "fio2_set"))
  meds <- read_window("medication_admin_continuous", "admin_dttm",
                      c("med_category", "med_dose", "med_dose_unit", "mar_action_group", "mar_action_category"),
                      "med_category", pressors) |>
    filter(normalize_category(mar_action_group) != "not_administered",
           normalize_category(mar_action_category) != "stop")
  gcs <- read_window("patient_assessments", "recorded_dttm", c("assessment_category", "numerical_value"),
                     "assessment_category", c("gcs_total", "gcs_eye", "gcs_motor", "gcs_verbal"))

  tables <- lapply(list(labs = labs, vitals = vitals, support = support, meds = meds, gcs = gcs),
                   function(x) split(x, x$hospitalization_id))
  score <- function(id) {
    l <- tables$labs[[id]]
    v <- tables$vitals[[id]]
    m <- tables$meds[[id]]
    g <- tables$gcs[[id]]
    values <- function(value, category, name, lower, upper) {
      x <- suppressWarnings(as.numeric(value[normalize_category(category) == name]))
      x[is.finite(x) & x >= lower & x <= upper]
    }
    platelet <- values(l$lab_value_numeric, l$lab_category, "platelet_count", 0, 3000)
    bili <- values(l$lab_value_numeric, l$lab_category, "bilirubin_total", 0, 100)
    creat <- values(l$lab_value_numeric, l$lab_category, "creatinine", 0, 25)
    map <- values(v$vital_value, v$vital_category, "map", 20, 250)
    weight <- values(v$vital_value, v$vital_category, "weight_kg", 10, 500)
    raw <- suppressWarnings(as.numeric(m$med_dose))
    dose <- aim2_sofa_dose(raw, m$med_dose_unit, if (length(weight)) median(weight) else NA_real_)
    vaso <- function(med) aim2_sofa_max(dose[normalize_category(m$med_category) == med & dose > 0])
    dopamine <- vaso("dopamine")
    epinephrine <- vaso("epinephrine")
    norepinephrine <- vaso("norepinephrine")
    dobutamine <- vaso("dobutamine")
    pressor_score <- case_when(epinephrine > 0.1 | norepinephrine > 0.1 | dopamine > 15 ~ 4,
                              epinephrine > 0 | norepinephrine > 0 | dopamine > 5 ~ 3,
                              dopamine > 0 | dobutamine > 0 ~ 2)
    # An unconvertible positive dose leaves the cardiovascular score unknown
    # unless another observation establishes a score of 4.
    cv <- if (any(raw > 0 & !is.finite(dose), na.rm = TRUE) && !isTRUE(pressor_score == 4))
      NA_real_ else aim2_sofa_max(c(as.numeric(map < 70), pressor_score))

    total_gcs <- values(g$numerical_value, g$assessment_category, "gcs_total", 3, 15)
    if (!is.null(g)) {
      # Sum components only when all three were charted at the same timestamp.
      component_gcs <- g |>
        mutate(category = normalize_category(assessment_category),
               value = suppressWarnings(as.numeric(numerical_value)),
               upper = c(gcs_eye = 4, gcs_motor = 6, gcs_verbal = 5)[category]) |>
        filter(is.finite(value), value >= 1, value <= upper, value == floor(value)) |>
        group_by(recorded_dttm, category) |>
        summarise(value = min(value), .groups = "drop_last") |>
        summarise(components = n(), score = sum(value), .groups = "drop") |>
        filter(components == 3L)
      total_gcs <- c(total_gcs, component_gcs$score)
    }
    c(resp = aim2_sofa_respiratory(l, tables$support[[id]]),
      cv = cv,
      coag = aim2_sofa_max(case_when(platelet < 20 ~ 4, platelet < 50 ~ 3, platelet < 100 ~ 2,
                                   platelet < 150 ~ 1, TRUE ~ 0)),
      liver = aim2_sofa_max(findInterval(bili, c(1.2, 2, 6, 12))),
      renal = aim2_sofa_max(findInterval(creat, c(1.2, 2, 3.5, 5))),
      cns = aim2_sofa_max(case_when(total_gcs < 6 ~ 4, total_gcs < 10 ~ 3, total_gcs < 13 ~ 2,
                                  total_gcs < 15 ~ 1, TRUE ~ 0)))
  }
  domains <- t(vapply(cohort$hospitalization_id, score, setNames(numeric(6), aim2_sofa_domains)))
  colnames(domains) <- aim2_sofa_domains
  bind_cols(tibble(hospitalization_id = cohort$hospitalization_id, sofa_24h = rowSums(domains),
                   sofa_domains_observed = rowSums(is.finite(domains))), as_tibble(domains))
}

aim2_sofa_domains <- c("resp", "cv", "coag", "liver", "renal", "cns")

aim2_sofa_specification <- c(
  "SOFA landmark: patients with observed follow-up beyond hour 24, with time reset at hour 24 and administrative censoring at ICU day 28.",
  "Pre-landmark physiology uses index-stay records from first ICU entry inclusive to hour 24 exclusive.",
  "The measured SOFA proxy sums the worst in-window scores across six organ domains. It uses creatinine without urine output, observed GCS without sedation correction, and paired arterial PaO2 and FiO2.",
  "Landmark models compare the standard adjustment with and without linear SOFA on identical records with all six domains observed. An all-eligible landmark model also uses standard adjustment.",
  "Landmark VFD uses the 27 days from hour 24 through ICU day 28 and permits respiratory support intervals to carry across the landmark within the usual six-hour limit."
)

aim2_run_sofa_landmark <- function(d, respiratory, tables_path, file_type) {
  landmark <- droplevels(filter(d, time > 1))
  scores <- aim2_sofa_scores(landmark, tables_path, file_type)
  original_vfd_fields <- c("n_support_records", "imv_hours", "last_support_imv", "vfd", "vfd_status")
  landmark <- landmark |>
    left_join(select(scores, hospitalization_id, sofa_24h, sofa_domains_observed), by = "hospitalization_id") |>
    mutate(time = time - 1) |>
    select(-any_of(original_vfd_fields)) |>
    aim2_vfd(respiratory, start_day = 1)
  complete <- droplevels(filter(landmark, sofa_domains_observed == length(aim2_sofa_domains)))

  outcome_count <- function(data, endpoint) {
    if (endpoint %in% c("death_or_hospice", "ventilator_free_days"))
      return(sum(data$event %in% if (endpoint == "death_or_hospice") c("death", "hospice") else "death"))
    sum(data$event == endpoint)
  }
  status_row <- function(data, endpoint, model, method, outcome, initial = outcome) {
    describe <- function(x) if (inherits(x, "error")) conditionMessage(x) else "ok"
    tibble(model, endpoint, method, n = nrow(data), events = outcome_count(data, endpoint),
           unpenalized_status = describe(initial), status = describe(outcome))
  }
  fit_models <- function(data, baseline_model, adjusted_model = NULL) {
    terms <- setNames(list(aim2_adjustment(data)), baseline_model)
    if (!is.null(adjusted_model)) terms[[adjusted_model]] <- c(terms[[1]], "sofa_24h")
    cox <- aim2_exit_models(data, terms)
    estimates <- status <- list()
    for (model in names(terms)) {
      for (endpoint in names(cox$fits[[model]])) {
        outcome <- cox$fits[[model]][[endpoint]]
        status[[length(status) + 1L]] <- status_row(data, endpoint, model, cox$method, outcome,
                                                  cox$unpenalized[[model]][[endpoint]])
        if (!inherits(outcome, "error")) {
          estimates[[length(estimates) + 1L]] <- bind_rows(lapply(aim2_hot_temperatures, function(hot) {
            cbind(aim2_contrast(outcome$fit, vcov(outcome$fit), nrow(data), outcome$n_events,
                                endpoint, model, data, hot), method = cox$method)
          }))
        }
      }
      vfd_data <- droplevels(filter(data, !is.na(vfd)))
      vfd_terms <- c(aim2_adjustment(vfd_data), intersect(terms[[model]], "sofa_24h"))
      vfd <- tryCatch(aim2_vfd_estimates(vfd_data, vfd_terms), error = identity)
      status[[length(status) + 1L]] <- status_row(vfd_data, "ventilator_free_days", model, "quasi_poisson", vfd)
      if (!inherits(vfd, "error")) {
        vfd$model <- model
        vfd$method <- "quasi_poisson"
        estimates[[length(estimates) + 1L]] <- vfd
      }
    }
    list(estimates = bind_rows(estimates), status = bind_rows(status))
  }

  all <- fit_models(landmark, "landmark_all")
  matched <- fit_models(complete, "landmark_matched", "landmark_sofa_24h")
  coverage <- tibble(
    primary_cox_n = nrow(d), excluded_at_or_before_landmark = sum(d$time <= 1),
    death_at_or_before_landmark = sum(d$time <= 1 & d$event == "death"),
    hospice_at_or_before_landmark = sum(d$time <= 1 & d$event == "hospice"),
    live_exit_at_or_before_landmark = sum(d$time <= 1 & d$event == "live_discharge"),
    censored_at_or_before_landmark = sum(d$time <= 1 & d$event == "censored"),
    landmark_n = nrow(landmark), complete_sofa_n = nrow(complete),
    landmark_vfd_n = sum(!is.na(landmark$vfd)), complete_sofa_vfd_n = sum(!is.na(complete$vfd)),
    complete_sofa_median = if (nrow(complete)) median(complete$sofa_24h) else NA_real_
  )
  domains <- tibble(domain = aim2_sofa_domains,
                    observed_n = colSums(!is.na(scores[aim2_sofa_domains])))
  distribution <- tibble(observed_domains = 0:6,
                         n = tabulate(scores$sofa_domains_observed + 1L, nbins = 7L))
  list(sofa_landmark_estimates = bind_rows(all$estimates, matched$estimates),
       sofa_landmark_status = bind_rows(all$status, matched$status),
       sofa_landmark_coverage = coverage, sofa_landmark_domains = domains,
       sofa_landmark_domain_counts = distribution)
}

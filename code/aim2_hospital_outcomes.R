# Sourced by the shared cohort builder. Outputs are aggregate only.
aim2_horizon <- 28
aim2_temp_basis <- list(knots = c(10, 25), Boundary.knots = c(-40, 45))
aim2_reference_c <- 20
aim2_hot_temperatures <- 30:35

aim2_charlson <- function(cohort, hospitalizations, diagnosis, diagnosis_source) {
  index <- transmute(cohort, index_id = hospitalization_id, patient_id, index_start = admission_dttm,
                     lookback = admission_dttm - 365 * 86400)
  prior <- hospitalizations |>
    transmute(dx_id = hospitalization_id, patient_id, prior_end = discharge_dttm) |>
    inner_join(index, by = join_by(patient_id, prior_end <= index_start, prior_end >= lookback),
               relationship = "many-to-many", na_matches = "never") |>
    filter(dx_id != index_id) |>
    distinct(index_id, dx_id)
  maps <- bind_rows(mutate(prior, is_index = FALSE),
                   transmute(index, index_id, dx_id = index_id, is_index = TRUE))
  codes <- transmute(diagnosis, dx_id = as.character(hospitalization_id),
                     code = norm_code(diagnosis_code),
                     version = infer_diagnosis_format(code, diagnosis_code_format),
                     poa = derive_poa_flag(diagnosis, diagnosis_source)) |>
    inner_join(maps, by = "dx_id", relationship = "many-to-many") |>
    filter(nzchar(code), !is_index | poa) |>
    distinct(index_id, code, version, is_index)
  weights <- c(mi = 1, chf = 1, pvd = 1, cevd = 1, dementia = 1, cpd = 1,
               rheumd = 1, pud = 1, mld = 1, diab = 1, diabwc = 2, hp = 2, rend = 2,
               canc = 2, msld = 3, metacanc = 6, aids = 6)
  domains <- list()
  for (version in c("9", "10")) {
    for (index_codes in c(FALSE, TRUE)) {
      subset <- filter(codes, str_detect(.data$version, .env$version), is_index == index_codes) |>
        # An ICD-10 label must not also match the ICD-9 branch.
        filter(if (.env$version == "9") !str_detect(.data$version, "10") else TRUE) |>
        select(id = index_id, code)
      if (!nrow(subset)) next
      flags <- comorbidity::comorbidity(subset, id = "id", code = "code",
                                        map = paste0("charlson_icd", version, "_quan"), assign0 = FALSE)
      if (index_codes) flags[c("mi", "cevd", "hp")] <- 0L
      domains[[length(domains) + 1L]] <- flags
    }
  }
  score <- tibble(hospitalization_id = character(), charlson = numeric())
  if (length(domains)) {
    flags <- bind_rows(domains) |>
      group_by(id) |>
      summarise(across(all_of(names(weights)), max), .groups = "drop")
    flags$diab[flags$diabwc == 1L] <- 0L
    flags$mld[flags$msld == 1L] <- 0L
    flags$canc[flags$metacanc == 1L] <- 0L
    score <- tibble(hospitalization_id = flags$id,
                    charlson = drop(as.matrix(flags[names(weights)]) %*% weights))
  }
  cohort |>
    left_join(score, by = "hospitalization_id") |>
    mutate(charlson = coalesce(charlson, 0),
           prior_admission_found = hospitalization_id %in% prior$index_id)
}

aim2_exposures <- function(cohort, root) {
  counties <- unique(stats::na.omit(cohort$county_fips))
  read_county_weather <- function(file, column) {
    out <- open_dataset(file.path(root, "exposome", file)) |>
      filter(geoid %in% counties) |>
      select(geoid, date, all_of(column)) |>
      collect() |>
      rename(county_fips = geoid, exposure_date = date)
    stopifnot(!anyDuplicated(select(out, county_fips, exposure_date)))
    out
  }
  cohort |>
    mutate(exposure_date = admission_date - 1L) |>
    left_join(read_county_weather("daymet_county_tmax_2018_2024_conus.parquet", "tmax_mean_c"),
              by = c("county_fips", "exposure_date")) |>
    left_join(read_county_weather("daymet_county_rmax_2018_2024.parquet", "rmax_mean_pct"),
              by = c("county_fips", "exposure_date"))
}

aim2_exits <- function(cohort, patient, adt) {
  death <- tibble(patient_id = as.character(patient$patient_id),
                  patient_death = if ("death_dttm" %in% names(patient))
                    as_utc_datetime(patient$death_dttm) else as_utc_datetime(rep(NA_character_, nrow(patient))))
  stopifnot(!anyDuplicated(death$patient_id))
  observed <- adt |>
    semi_join(cohort, by = "hospitalization_id") |>
    group_by(hospitalization_id) |>
    summarise(last_adt = as_utc_datetime(if (all(is.na(c(in_dttm, out_dttm)))) NA_real_ else
      max(as.numeric(c(in_dttm, out_dttm)), na.rm = TRUE)), .groups = "drop")
  cohort |>
    left_join(death, by = "patient_id", na_matches = "never") |>
    left_join(observed, by = "hospitalization_id") |>
    mutate(
      disposition = normalize_category(gsub("_", " ", discharge_category)),
      acute_transfer = disposition %in% c("acute care hospital", "other hospital"),
      ltach = disposition %in% c("long term care hospital (ltach)", "ltach"),
      outcome = case_when(
        str_detect(disposition, "expired|death|died|deceased") ~ "death",
        str_detect(disposition, "hospice") ~ "hospice",
        disposition %in% c("home", "home health", "skilled nursing facility (snf)", "skilled nursing facility",
                           "acute inpatient rehab facility", "rehab", "long term care hospital (ltach)",
                           "ltach", "acute care hospital", "other hospital", "group home", "chemical dependency",
                           "against medical advice (ama)", "assisted living", "psychiatric hospital",
                           "psychiatric unit", "shelter", "jail", "alive", "discharged alive") ~ "live_discharge",
        TRUE ~ "unknown"
      ),
      validated_death = outcome == "death" & !is.na(patient_death) & patient_death >= first_icu_in &
        (is.na(discharge_dttm) | patient_death <= discharge_dttm),
      event_time = if_else(validated_death, patient_death, discharge_dttm),
      end_time = case_when(
        !is.na(event_time) ~ event_time,
        outcome == "unknown" ~ last_adt,
        TRUE ~ as_utc_datetime(NA_real_)
      ),
      follow_days = as.numeric(difftime(end_time, first_icu_in, units = "days")),
      valid_followup = !is.na(follow_days) & follow_days >= 0 & first_icu_in >= admission_dttm,
      time = if_else(valid_followup, pmin(pmax(follow_days, 1 / 1440), aim2_horizon), NA_real_),
      event = if_else(valid_followup & follow_days <= aim2_horizon & outcome != "unknown",
                      outcome, "censored")
    )
}

aim2_vfd <- function(cohort, respiratory, gap_hours = 6) {
  support <- tibble(hospitalization_id = character(), n_support_records = integer(),
                    imv_hours = numeric(), last_support_imv = logical())
  if (!is.null(respiratory) && nrow(respiratory)) {
    support <- respiratory |>
      transmute(hospitalization_id = as.character(hospitalization_id),
                recorded = as_utc_datetime(recorded_dttm), device = normalize_category(device_category)) |>
      filter(!device %in% c("", "missing", "unknown"), !is.na(recorded)) |>
      inner_join(select(cohort, hospitalization_id, first_icu_in, end_time, valid_followup),
                 by = "hospitalization_id") |>
      filter(valid_followup) |>
      mutate(observed_end = pmin(end_time, first_icu_in + aim2_horizon * 86400),
             imv = device %in% c("imv", "vent", "invasive mechanical ventilation")) |>
      filter(recorded >= first_icu_in, recorded <= observed_end) |>
      group_by(hospitalization_id, recorded) |>
      summarise(imv = any(imv), observed_end = first(observed_end), .groups = "drop_last") |>
      arrange(recorded, .by_group = TRUE) |>
      mutate(interval_end = pmin(recorded + gap_hours * 3600, lead(recorded), observed_end, na.rm = TRUE)) |>
      summarise(n_support_records = n(),
                imv_hours = sum(as.numeric(difftime(interval_end, recorded, units = "hours"))[imv]),
                last_support_imv = last(imv), .groups = "drop")
  }
  cohort |>
    left_join(support, by = "hospitalization_id") |>
    mutate(
      vfd_status = case_when(
        !valid_followup ~ "invalid_followup",
        event == "death" ~ "hospital_death_28d",
        follow_days < aim2_horizon & (acute_transfer | ltach) ~ "early_acute_transfer_or_ltach",
        follow_days < aim2_horizon & event == "hospice" ~ "early_hospice",
        follow_days < aim2_horizon & event != "live_discharge" ~ "incomplete_followup",
        is.na(n_support_records) ~ "no_observed_respiratory_support",
        follow_days < aim2_horizon & last_support_imv ~ "last_record_on_imv_before_live_discharge",
        TRUE ~ "vfd_observable_with_discharge_assumption"
      ),
      vfd = case_when(
        vfd_status == "hospital_death_28d" ~ 0,
        vfd_status == "vfd_observable_with_discharge_assumption" ~ pmax(0, aim2_horizon - imv_hours / 24),
        TRUE ~ NA_real_
      )
    )
}

aim2_temperature <- function(values) {
  basis <- do.call(splines::ns, c(list(x = values), aim2_temp_basis))
  as.data.frame(stats::setNames(as.data.frame(basis), paste0("temp_", seq_len(ncol(basis)))))
}

aim2_model_data <- function(cohort) {
  d <- cohort |>
    filter(valid_followup, is.finite(tmax_mean_c), is.finite(rmax_mean_pct)) |>
    mutate(
      model_patient = if_else(is.na(patient_id) | patient_id == "",
                              paste0("index:", hospitalization_id), paste0("patient:", patient_id)),
      sex = factor(case_when(normalize_category(sex_category) %in% c("male", "m") ~ "Male",
                             normalize_category(sex_category) %in% c("female", "f") ~ "Female",
                             TRUE ~ "Other or unknown")),
      race = factor(case_when(is_black_race(race_category) ~ "Black",
                              normalize_category(race_category) == "white" ~ "White",
                              normalize_category(race_category) %in% c("", "unknown") ~ "Unknown",
                              TRUE ~ "Other")),
      ethnicity = factor(case_when(normalize_category(ethnicity_category) == "hispanic" ~ "Hispanic",
                                   normalize_category(ethnicity_category) %in% c("non-hispanic", "non hispanic") ~ "Non-Hispanic",
                                   TRUE ~ "Unknown")),
      year = factor(format(admission_date, "%Y")),
      season_sin = sin(2 * pi * (lubridate::yday(admission_date) - 1) / 365.25),
      season_cos = cos(2 * pi * (lubridate::yday(admission_date) - 1) / 365.25)
    )
  if (!nrow(d)) stop("No eligible Aim 2 model rows")
  bind_cols(d, aim2_temperature(d$tmax_mean_c))
}

aim2_cause_specific <- function(d, endpoint, rhs) {
  ridge <- survival::ridge
  exits <- if (endpoint == "death_or_hospice") c("death", "hospice") else endpoint
  d$target_event <- as.integer(d$event %in% exits)
  if (sum(d$target_event) < length(rhs) + 2L) stop("Too few events for prespecified covariates")
  form <- reformulate(rhs, response = "survival::Surv(time, target_event)")
  # Efron ties and a patient-level sandwich accommodate repeat hospitalizations.
  fit <- survival::coxph(form, data = d, cluster = model_patient,
                        x = TRUE, model = TRUE, ties = "efron", singular.ok = FALSE)
  if (any(!is.finite(coef(fit)))) stop("Non-estimable Cox coefficients")
  list(fit = fit, n_events = sum(d$target_event))
}

aim2_contrast <- function(fit, covariance, n, events, endpoint, model, d, hot) {
  temps <- c(aim2_reference_c, hot)
  predicted_basis <- aim2_temperature(temps)
  b <- coef(fit)
  delta <- unlist(predicted_basis[2, ] - predicted_basis[1, ], use.names = TRUE)
  estimate <- sum(delta * b[names(delta)])
  temp_covariance <- covariance[names(delta), names(delta)]
  se <- sqrt(drop(t(delta) %*% temp_covariance %*% delta))
  temp_coefficients <- b[names(delta)]
  wald <- tryCatch(drop(crossprod(temp_coefficients, solve(temp_covariance, temp_coefficients))),
                   error = function(e) NA_real_)
  data.frame(model, endpoint, n, events, reference_c = temps[1], hot_c = temps[2],
             reference_support = sum(abs(d$tmax_mean_c - temps[1]) <= 1),
             hot_support = sum(abs(d$tmax_mean_c - temps[2]) <= 1),
             in_exposure_range = min(d$tmax_mean_c) <= temps[1] && max(d$tmax_mean_c) >= temps[2],
             log_ratio = estimate, se, ratio = exp(estimate),
             low = exp(estimate - 1.96 * se), high = exp(estimate + 1.96 * se),
             temperature_wald_p = pchisq(wald, df = length(delta), lower.tail = FALSE),
             formula = paste(deparse(formula(fit)), collapse = " "))
}

aim2_set_temperature <- function(d, temperature) {
  basis <- aim2_temperature(rep(temperature, nrow(d)))
  d[names(basis)] <- basis
  d
}

aim2_standardized_exits <- function(d, fits, endpoints) {
  hazards <- lapply(fits[endpoints], function(fit) survival::basehaz(fit, centered = TRUE))
  report_days <- 0:aim2_horizon
  times <- sort(unique(c(report_days, unlist(lapply(hazards, function(h) h$time[h$time <= aim2_horizon])))))
  baseline_steps <- vapply(hazards, function(h) {
    cumulative <- c(0, h$hazard)[1L + findInterval(times, h$time)]
    diff(c(0, cumulative))
  }, numeric(length(times)))
  bind_rows(lapply(c(aim2_reference_c, aim2_hot_temperatures), function(temp) {
    scenario <- aim2_set_temperature(d, temp)
    risk <- vapply(fits[endpoints], function(fit)
      predict(fit, newdata = scenario, type = "risk", reference = "sample"), numeric(nrow(d)))
    survival_prob <- rep(1, nrow(d))
    cif <- matrix(0, nrow(d), length(endpoints))
    probabilities <- matrix(NA_real_, length(report_days), length(endpoints) + 1L)
    for (time_index in seq_along(times)) {
      if (any(baseline_steps[time_index, ] > 0)) {
        increment <- sweep(risk, 2, baseline_steps[time_index, ], "*")
        total <- rowSums(increment)
        fail <- -expm1(-total)
        cif <- cif + increment * ifelse(total > 0, survival_prob * fail / total, 0)
        survival_prob <- survival_prob * exp(-total)
      }
      report_index <- match(times[time_index], report_days)
      if (!is.na(report_index)) probabilities[report_index, ] <- c(colMeans(cif), mean(survival_prob))
    }
    stopifnot(max(abs(rowSums(cif) + survival_prob - 1)) < 1e-8)
    support <- sum(abs(d$tmax_mean_c - temp) <= 1)
    bind_rows(lapply(seq_len(ncol(probabilities)), function(i) {
      data.frame(temperature_c = temp, endpoint = c(endpoints, "still_hospitalized")[[i]],
                 day = report_days, probability = probabilities[, i], n = nrow(d),
                 support_within_1c = support,
                 in_exposure_range = temp >= min(d$tmax_mean_c) && temp <= max(d$tmax_mean_c))
    }))
  }))
}

aim2_vfd_estimates <- function(d, rhs) {
  fit <- glm(reformulate(rhs, response = "vfd"), data = d,
             family = quasipoisson(link = "log"), na.action = na.fail)
  if (!fit$converged || anyNA(coef(fit))) stop("VFD fit failed or contains aliased covariates")
  covariance <- sandwich::vcovCL(fit, cluster = d$model_patient)
  predict_mean <- function(temp) {
    design <- model.matrix(delete.response(terms(fit)), data = aim2_set_temperature(d, temp),
                           contrasts.arg = fit$contrasts, xlev = fit$xlevels)
    fitted_mean <- drop(exp(design %*% coef(fit)))
    list(mean = mean(fitted_mean), gradient = colMeans(design * fitted_mean))
  }
  reference <- predict_mean(aim2_reference_c)
  bind_rows(lapply(aim2_hot_temperatures, function(hot) {
    result <- aim2_contrast(fit, covariance, nrow(d), sum(d$event == "death"),
                           "ventilator_free_days", "quasi_poisson", d, hot)
    prediction <- predict_mean(hot)
    difference <- prediction$mean - reference$mean
    slope <- prediction$gradient - reference$gradient
    diff_se <- sqrt(drop(t(slope) %*% covariance %*% slope))
    cbind(result, mean_20c = reference$mean, mean_hot = prediction$mean,
          difference_days = difference, difference_low = difference - 1.96 * diff_se,
          difference_high = difference + 1.96 * diff_se, difference_se = diff_se)
  }))
}

run_aim2_site <- function(output_dir, cohort, hospitalizations, patient, diagnosis, adt,
                          diagnosis_source, tables_path, file_type, root, site) {
  for (package in c("comorbidity", "survival", "sandwich")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install the Aim 2 dependency ", package)
  }
  stopifnot(!anyDuplicated(cohort$hospitalization_id))
  cohort <- aim2_charlson(cohort, hospitalizations, diagnosis, diagnosis_source) |>
    aim2_exposures(root) |>
    aim2_exits(patient, adt)
  respiratory <- read_clif_table(tables_path, file_type, "respiratory_support",
                                  columns = c("hospitalization_id", "recorded_dttm", "device_category"),
                                  required = FALSE)
  cohort <- aim2_vfd(cohort, respiratory)
  d <- aim2_model_data(cohort)
  baseline <- c("age_at_admission", "sex", "race", "ethnicity", "charlson", "year", "season_sin", "season_cos")
  choose_rhs <- function(x) {
    kept <- baseline[vapply(x[baseline], n_distinct, integer(1)) > 1L]
    c("temp_1", "temp_2", "temp_3", "splines::ns(rmax_mean_pct, df = 3)", kept)
  }
  endpoints <- c("death", "hospice", "live_discharge")
  all_endpoints <- c(endpoints, "death_or_hospice")
  fit_exits <- function(data, rhs) {
    results <- lapply(all_endpoints, function(endpoint)
      tryCatch(withCallingHandlers(aim2_cause_specific(data, endpoint, rhs),
                                   warning = function(w) stop(conditionMessage(w))), error = identity))
    setNames(results, all_endpoints)
  }
  rhs <- choose_rhs(d)
  unpenalized <- fit_exits(d, rhs)
  separation <- vapply(unpenalized, function(x)
    inherits(x, "error") && grepl("coefficient may be infinite", conditionMessage(x), fixed = TRUE), logical(1))
  demographic_terms <- intersect(c("sex", "race", "ethnicity"), rhs)
  cox_method <- "unpenalized"
  cox_data <- d
  all_fits <- unpenalized
  if (any(separation) && length(demographic_terms)) {
    # Use the same fallback for every exit so the hazards can be combined.
    # The temperature and other continuous covariates stay unpenalized.
    demographics <- model.matrix(reformulate(demographic_terms), data = d)[, -1L, drop = FALSE]
    colnames(demographics) <- paste0("demographic_", seq_len(ncol(demographics)))
    cox_data <- bind_cols(d, as.data.frame(demographics))
    penalty <- paste0("ridge(", paste(colnames(demographics), collapse = ", "),
                      ", theta = 1, scale = FALSE)")
    all_fits <- fit_exits(cox_data, c(setdiff(rhs, demographic_terms), penalty))
    cox_method <- "ridge_demographic_nuisance"
  }
  fits <- list()
  estimates <- list()
  status <- list()
  ph <- list()
  for (endpoint in all_endpoints) {
    outcome <- all_fits[[endpoint]]
    status[[endpoint]] <- tibble(
      endpoint, method = cox_method,
      unpenalized_status = if (inherits(unpenalized[[endpoint]], "error"))
        conditionMessage(unpenalized[[endpoint]]) else "ok",
      status = if (inherits(outcome, "error")) conditionMessage(outcome) else "ok"
    )
    if (inherits(outcome, "error")) next
    fits[[endpoint]] <- outcome$fit
    estimates[[endpoint]] <- bind_rows(lapply(aim2_hot_temperatures, function(hot) {
      cbind(aim2_contrast(outcome$fit, vcov(outcome$fit), nrow(d), outcome$n_events,
                          endpoint, "cause_specific_cox", d, hot), method = cox_method)
    }))
    diagnostics <- tryCatch(survival::cox.zph(outcome$fit)$table, error = function(e) NULL)
    if (!is.null(diagnostics)) {
      valid <- is.finite(diagnostics[, "chisq"]) & diagnostics[, "chisq"] >= 0
      ph[[endpoint]] <- data.frame(
        endpoint, term = rownames(diagnostics), statistic = diagnostics[, "chisq"],
        df = diagnostics[, "df"], p_value = ifelse(valid, diagnostics[, "p"], NA_real_),
        status = ifelse(valid, "ok", "non_estimable"), row.names = NULL)
    }
  }

  d_vfd <- droplevels(filter(d, !is.na(vfd)))
  vfd <- tryCatch(aim2_vfd_estimates(d_vfd, choose_rhs(d_vfd)), error = identity)
  status$vfd <- tibble(endpoint = "ventilator_free_days", method = "quasi_poisson",
                       status = if (inherits(vfd, "error")) conditionMessage(vfd) else "ok")
  if (!inherits(vfd, "error")) estimates$vfd <- vfd

  incidence <- data.frame()
  if (all(endpoints %in% names(fits))) {
    incidence <- aim2_standardized_exits(cox_data, fits, endpoints)
  }
  d$aj_event <- factor(d$event, levels = c("censored", endpoints))
  aj <- survival::survfit(survival::Surv(time, aj_event) ~ 1, data = d,
                          id = hospitalization_id, cluster = model_patient, robust = TRUE)
  aj_summary <- summary(aj, times = c(0, 1, 3, 7, 14, 21, 28), extend = TRUE)
  observed <- bind_rows(lapply(seq_along(aj$states), function(i) {
    data.frame(day = aj_summary$time, endpoint = aj$states[[i]], probability = aj_summary$pstate[, i],
               se = aj_summary$std.err[, i], at_risk = aj_summary$n.risk[, 1])
  }))

  audit <- cohort |>
    summarise(eligible_ohca = n(), valid_followup = sum(valid_followup),
              missing_temperature_or_humidity = sum(!is.finite(tmax_mean_c) | !is.finite(rmax_mean_pct)),
              model_cohort = nrow(d), vfd_cohort = nrow(d_vfd),
              hospital_deaths_28d = sum(d$event == "death"), hospice_28d = sum(d$event == "hospice"),
              live_discharge_28d = sum(d$event == "live_discharge"),
              acute_transfer_28d = sum(d$event == "live_discharge" & d$acute_transfer),
              censored_before_28d = sum(d$event == "censored" & d$time < aim2_horizon),
              prior_admission_found = sum(d$prior_admission_found),
              vfd_deaths_zero = sum(d_vfd$event == "death"),
              vfd_mean = mean(d_vfd$vfd), mean_imv_days_measured = mean(d_vfd$imv_hours / 24, na.rm = TRUE))
  tables <- list(estimates = bind_rows(estimates), model_status = bind_rows(status),
                 coverage = audit, vfd_status = count(cohort, vfd_status, name = "n"),
                 proportional_hazards = bind_rows(ph), standardized_cif = incidence,
                 standardized_day28 = if (nrow(incidence)) filter(incidence, day == aim2_horizon) else incidence,
                 aalen_johansen = observed)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    table <- tables[[name]]
    table$site <- rep(site, nrow(table))
    readr::write_csv(table, file.path(output_dir, paste0(name, ".csv")))
  }
  writeLines(c("Time origin: first ICU admission; horizon: 28 days; all events are from the index hospital stay.",
               "Exposure: preceding-day assigned-county Tmax and maximum relative humidity.",
               paste("Diagnosis table:", diagnosis_source),
               paste("Model adjustment:", paste(choose_rhs(d), collapse = " + ")),
               paste("Cox method:", cox_method),
               "VFD: zero for index-stay death by day 28; otherwise 28 minus documented cumulative IMV days.",
               "Six-hour respiratory carry-forward, clipped at the next documentation, observed end, and day 28.",
               "Early eligible live discharges are assumed alive and free from IMV to day 28.",
               "Adjusted day-28 event probabilities are exploratory point estimates.",
               capture.output(sessionInfo())), file.path(output_dir, "specification.txt"))
}

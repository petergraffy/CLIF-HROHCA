# Fit joint temperature models using county-day subgroup returns.
read_aim1_subgroups <- function(path, panel) {
  site <- panel$site_name[[1L]]
  required <- c("site_name", "county_fips", "admission_date", "subgroup_variable",
                "subgroup_level", "n_icu_admissions", "ohca_admissions")
  x <- read.csv(path, colClasses = c(county_fips = "character"))
  if (!all(required %in% names(x))) stop("Subgroup file is missing required columns")
  x <- x |>
    filter(!is.na(county_fips), nzchar(county_fips)) |>
    transmute(site_name, county_fips, date = as.Date(admission_date),
              subgroup_variable, subgroup_level, n_icu_admissions, ohca_admissions)
  keys <- c("county_fips", "date")
  stopifnot(nrow(x) > 0L, !anyNA(x), all(x$site_name == site),
            !anyDuplicated(x[c(keys, "subgroup_variable", "subgroup_level")]),
            is.numeric(x$n_icu_admissions), is.numeric(x$ohca_admissions),
            all(is.finite(x$n_icu_admissions) & is.finite(x$ohca_admissions)),
            all(x$n_icu_admissions == floor(x$n_icu_admissions)),
            all(x$ohca_admissions == floor(x$ohca_admissions)),
            all(x$n_icu_admissions >= x$ohca_admissions & x$ohca_admissions >= 0))
  expected <- panel |>
    filter(n_icu_admissions > 0) |>
    select(all_of(keys), n_icu_admissions, ohca_admissions)
  for (variable in c("age_group", "sex", "race")) {
    totals <- x |>
      filter(subgroup_variable == variable) |>
      group_by(county_fips, date) |>
      summarise(n_icu_admissions = sum(n_icu_admissions), ohca_admissions = sum(ohca_admissions),
                .groups = "drop")
    columns <- c(keys, "n_icu_admissions", "ohca_admissions")
    if (nrow(anti_join(expected, totals, by = columns)) ||
        nrow(anti_join(totals, expected, by = columns))) {
      stop("Subgroup counts do not reconcile with the county-day export for ", variable)
    }
  }
  x
}

aim1_subgroup_specifications <- list(
  age = list(variable = "age_group", target = "65+", reference = "18-64",
             target_levels = "65+", reference_levels = "18-64"),
  sex = list(variable = "sex", target = "Female", reference = "Male",
             target_levels = "Female", reference_levels = "Male"),
  race = list(variable = "race", target = "Black", reference = "Recorded non-Black",
              target_levels = "Black or African American",
              reference_levels = c("White", "Asian", "American Indian or Alaska Native",
                                   "Native Hawaiian or Other Pacific Islander", "Other")),
  race_black_white = list(variable = "race", target = "Black", reference = "White",
                          target_levels = "Black or African American", reference_levels = "White")
)

subgroup_date_covariance <- function(fit, contrasts, dates) {
  terms <- colnames(contrasts)
  direction <- summary(fit)$cov.unscaled[terms, terms, drop = FALSE] %*% t(contrasts)
  scores <- (model.matrix(fit)[, terms, drop = FALSE] %*% direction) *
    (model.response(model.frame(fit)) - fitted(fit))
  calendar <- seq(min(dates), max(dates), by = "day")
  daily_scores <- rowsum(scores, match(dates, calendar))
  daily <- matrix(0, nrow = length(calendar), ncol = nrow(contrasts))
  daily[as.integer(rownames(daily_scores)), ] <- daily_scores
  bandwidth <- min(14L, nrow(daily) - 1L)
  covariance <- crossprod(daily)
  for (lag in seq_len(bandwidth)) {
    serial <- crossprod(head(daily, -lag), tail(daily, -lag))
    covariance <- covariance + (1 - lag / (bandwidth + 1)) * (serial + t(serial))
  }
  covariance
}

fit_aim1_subgroup <- function(panel, returned, spec, model, temp_basis, hot, cold) {
  keys <- c("county_fips", "date")
  lag_basis <- list(fun = "ns", knots = 2, Boundary.knots = c(0, 5))
  cb <- crossbasis(panel$tmax_mean_c, lag = 5, group = factor(panel$county_fips),
                   argvar = temp_basis, arglag = lag_basis)
  cb_names <- paste0("temperature_", seq_len(ncol(cb)))
  eligible <- panel$observed & complete.cases(cb, panel$rmax_mean_pct)
  base <- panel[eligible, ]
  base[cb_names] <- as.data.frame(unclass(cb)[eligible, , drop = FALSE])
  counts <- returned |>
    filter(subgroup_variable == spec$variable,
           subgroup_level %in% c(spec$target_levels, spec$reference_levels)) |>
    mutate(subgroup = if_else(subgroup_level %in% spec$target_levels, "target", "reference")) |>
    group_by(county_fips, date, subgroup) |>
    summarise(y = sum(ohca_admissions), subgroup_icu = sum(n_icu_admissions), .groups = "drop")
  x <- tidyr::expand_grid(panel_row = seq_len(nrow(base)), subgroup = c("reference", "target"))
  x <- bind_cols(select(x, -panel_row), base[x$panel_row, ]) |>
    left_join(counts, by = c(keys, "subgroup")) |>
    mutate(subgroup = factor(subgroup, levels = c("reference", "target")),
           y = coalesce(y, 0L), subgroup_icu = coalesce(subgroup_icu, 0L))
  if (model == "icu_offset") x <- filter(x, subgroup_icu > 0)
  sample <- x |>
    group_by(subgroup) |>
    summarise(n_days = n(), n_ohca = sum(y), n_icu = sum(subgroup_icu),
              near_hot = sum(abs(tmax_mean_c - hot) <= 1),
              near_ref = sum(abs(tmax_mean_c - cold) <= 1),
              in_range = min(tmax_mean_c) <= cold && max(tmax_mean_c) >= hot,
              .groups = "drop")
  if (nrow(sample) != 2L || any(sample$n_ohca == 0)) stop("No outcome events in one or both subgroups")
  nuisance <- c("ns(rmax_mean_pct, df = 3)", "dow")
  years <- nlevels(droplevels(x$year))
  nuisance <- c(nuisance, sprintf("ns(time_index, df = %d)", 4L * years))
  if (years > 1L) nuisance <- c(nuisance, "year")
  if (n_distinct(x$county_fips) > 1L) nuisance <- c(nuisance, "factor(county_fips)")
  # Let each subgroup have its own temperature and nuisance coefficients.
  rhs <- paste("subgroup * (", paste(c(cb_names, nuisance), collapse = " + "), ")")
  if (model == "icu_offset") rhs <- paste(rhs, "+ offset(log(subgroup_icu))")
  formula <- reformulate(rhs, response = "y")
  fit <- glm(formula, data = x, family = quasipoisson(), na.action = na.fail)
  if (!fit$converged) stop("Model failed to converge")
  coefficients <- coef(fit)
  interaction_names <- paste0("subgrouptarget:", cb_names)
  if (!all(is.finite(coefficients[c(cb_names, interaction_names)]))) {
    stop("Temperature interaction is aliased")
  }
  histories <- rbind(rep(hot, 6L), rep(cold, 6L))
  contrast_cb <- crossbasis(histories, lag = 5, argvar = temp_basis, arglag = lag_basis)
  delta <- contrast_cb[1, ] - contrast_cb[2, ]
  contrasts <- matrix(0, nrow = 3L, ncol = length(coefficients),
                      dimnames = list(NULL, names(coefficients)))
  contrasts[1:2, cb_names] <- rep(delta, each = 2L)
  contrasts[2:3, interaction_names] <- rep(delta, each = 2L)
  # Some nuisance coefficients can be aliased when a subgroup has no ICU days
  # in a county. A zero-weight alias can be removed from the contrast.
  retained <- is.finite(coefficients)
  estimates <- drop(contrasts[, retained, drop = FALSE] %*% coefficients[retained])
  se <- sqrt(diag(subgroup_date_covariance(fit, contrasts[, retained, drop = FALSE], x$date)))
  target_sample <- filter(sample, subgroup == "target")
  ref_sample <- filter(sample, subgroup == "reference")
  data.frame(
    site = panel$site_name[[1L]], model, target = spec$target, reference = spec$reference,
    estimate_type = c("reference_rr", "target_rr", "interaction_ratio"),
    reference_c = cold, hot_c = hot, lag_min = 0L, lag_max = 5L,
    n_days_target = target_sample$n_days, n_days_reference = ref_sample$n_days,
    n_ohca_target = target_sample$n_ohca, n_ohca_reference = ref_sample$n_ohca,
    n_icu_target = target_sample$n_icu, n_icu_reference = ref_sample$n_icu,
    near_hot_target = target_sample$near_hot, near_hot_reference = ref_sample$near_hot,
    near_ref_target = target_sample$near_ref, near_ref_reference = ref_sample$near_ref,
    supported = all(sample$in_range & sample$near_hot >= 10 & sample$near_ref >= 10),
    log_rr = estimates, se, rr = exp(estimates),
    low = exp(estimates - 1.96 * se), high = exp(estimates + 1.96 * se),
    formula = paste(deparse(formula), collapse = " ")
  )
}

run_aim1_subgroups <- function(inputs, panels, temp_basis, hot, cold, pool) {
  results <- list()
  status <- list()
  for (i in seq_along(inputs)) {
    panel <- panels[[i]]
    site <- panel$site_name[[1L]]
    path <- paste0(inputs[[i]], ".subgroups.csv")
    if (!file.exists(path)) {
      status[[site]] <- data.frame(site, comparison = NA_character_, model = NA_character_,
                                   status = "No subgroup export supplied")
      next
    }
    returned <- tryCatch(read_aim1_subgroups(path, panel), error = conditionMessage)
    for (comparison in names(aim1_subgroup_specifications)) {
      spec <- aim1_subgroup_specifications[[comparison]]
      for (model in c("ohca_count", "icu_offset")) {
        key <- paste(site, comparison, model)
        fit <- if (is.character(returned)) returned else tryCatch(
          fit_aim1_subgroup(panel, returned, spec, model, temp_basis, hot, cold),
          error = conditionMessage
        )
        status[[key]] <- data.frame(site, comparison, model,
                                    status = if (is.character(fit)) fit else "ok")
        if (!is.character(fit)) results[[key]] <- cbind(comparison, fit)
      }
    }
  }
  estimates <- bind_rows(results)
  pooled <- data.frame()
  if (nrow(estimates)) {
    pooled <- estimates |>
      filter(supported, is.finite(se), se > 0) |>
      group_by(comparison, model, target, reference, estimate_type, reference_c, hot_c, lag_min, lag_max) |>
      group_modify(function(d, key) if (nrow(d) >= 2L) pool(d) else tibble()) |>
      ungroup()
  }
  list(subgroup_site_estimates = estimates, subgroup_pooled_estimates = pooled,
       subgroup_model_status = bind_rows(status))
}

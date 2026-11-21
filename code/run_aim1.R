#!/usr/bin/env Rscript
# Usage from any directory: Rscript code/run_aim1.R OUTPUT_DIR INPUT_CSV ...
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(dlnm)
  library(splines)
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript code/run_aim1.R OUTPUT_DIR INPUT_CSV ...")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(dirname(normalizePath(script)))
out <- args[[1L]]
inputs <- args[-1L]
geography <- read.csv(file.path(root, "reference/clif_hospital_geography.csv"),
                      colClasses = "character")
pm_path <- Sys.getenv("AIM1_PM25_CSV")
historical_dir <- Sys.getenv("AIM1_HISTORICAL_DIR")
proxy_sites <- strsplit(Sys.getenv("AIM1_PROXY_SITES"), ",", fixed = TRUE)[[1L]]
cold <- 20
hot <- suppressWarnings(as.numeric(Sys.getenv("AIM1_HOT_C", unset = "30")))
if (!is.finite(hot) || hot <= cold) stop("AIM1_HOT_C must be a temperature above 20 C.")
temp_basis <- list(fun = "ns", knots = c(0, 15, 25), Boundary.knots = c(-40, 45))

read_panel <- function(path) {
  fields <- names(read.csv(path, nrows = 0))
  d <- read.csv(path, colClasses = if ("county_fips" %in% fields) c(county_fips = "character") else NA)
  stopifnot(length(unique(d$site_name)) == 1L, !anyNA(d$site_name))
  site <- d$site_name[[1L]]
  mapping <- filter(geography, site_name == site)
  counties <- unique(mapping$hospital_county_fips)
  if (!length(counties)) stop("Missing crosswalk for ", site)
  legacy <- !"county_fips" %in% names(d)
  if (legacy && length(counties) > 1L && !site %in% proxy_sites) {
    stop(site, " requires a county-day export. A flagship-county pilot can be enabled via AIM1_PROXY_SITES.")
  }
  d$source <- if (!legacy) "county_export" else if (length(counties) == 1L)
    "archived_single_county" else "archived_flagship_proxy"
  if (legacy) d$county_fips <- counties[[1L]]
  stopifnot(all(d$county_fips %in% counties))
  d$date <- as.Date(d$admission_date, format = "%Y-%m-%d")
  stopifnot(!anyNA(d$date), !anyDuplicated(d[c("county_fips", "date")]))
  d$historical_tmax <- if ("tmax_mean_c" %in% names(d)) d$tmax_mean_c else NA_real_
  d$historical_rmax <- if ("rmax_mean_pct" %in% names(d)) d$rmax_mean_pct else NA_real_
  if (!"ohca_ed_admissions" %in% names(d)) d$ohca_ed_admissions <- NA_real_
  for (col in c("ohca_admissions", "n_icu_admissions", "ohca_ed_admissions")) {
    stopifnot(is.numeric(d[[col]]) || all(is.na(d[[col]])),
              all(is.na(d[[col]]) | (is.finite(d[[col]]) & d[[col]] >= 0 & d[[col]] == floor(d[[col]]))))
  }
  stopifnot(!any(d$ohca_admissions > d$n_icu_admissions, na.rm = TRUE),
            !any(d$ohca_ed_admissions > d$ohca_admissions, na.rm = TRUE))
  system_dates <- unique(d$date[!is.na(d$n_icu_admissions) & d$n_icu_admissions > 0])
  panel <- d |>
    group_by(county_fips) |>
    group_modify(function(x, key) {
      calendar <- data.frame(date = seq(min(x$date), max(x$date), by = "day"))
      left_join(calendar, x, by = "date") |>
        mutate(site_name = site, source = d$source[[1L]],
               inferred_zero = !legacy & is.na(admission_date) & date %in% system_dates,
               across(c(n_icu_admissions, ohca_admissions), ~ ifelse(inferred_zero, 0, .x)),
               ohca_ed_admissions = ifelse(inferred_zero, 0, ohca_ed_admissions))
    }) |>
    ungroup() |>
    select(site_name, county_fips, date, source, inferred_zero, n_icu_admissions, ohca_admissions,
           ohca_ed_admissions, historical_tmax, historical_rmax)
  historical_path <- file.path(historical_dir, paste0(site, "_ohca_icu_admission_rate_daily_timeseries.csv"))
  if (!legacy && nzchar(historical_dir) && file.exists(historical_path)) {
    old <- read.csv(historical_path) |>
      transmute(date = as.Date(admission_date), historical_tmax = tmax_mean_c,
                historical_rmax = rmax_mean_pct)
    stopifnot(!anyNA(old$date), !anyDuplicated(old$date))
    panel <- left_join(select(panel, -historical_tmax, -historical_rmax), old, by = "date")
  }
  panel
}

panels <- lapply(inputs, read_panel)
stopifnot(!anyDuplicated(vapply(panels, function(x) x$site_name[[1L]], character(1))))
county_codes <- unique(unlist(lapply(panels, function(x) x$county_fips)))
read_weather <- function(file, column) {
  weather <- open_dataset(file.path(root, "exposome", file)) |>
    filter(geoid %in% county_codes) |>
    select(geoid, date, all_of(column)) |>
    collect() |>
    rename(county_fips = geoid)
  stopifnot(!anyDuplicated(weather[c("county_fips", "date")]))
  weather
}
weather <- full_join(
  read_weather("daymet_county_tmax_2018_2024_conus.parquet", "tmax_mean_c"),
  read_weather("daymet_county_rmax_2018_2024.parquet", "rmax_mean_pct"),
  by = c("county_fips", "date")
)
if (nzchar(pm_path)) {
  pm <- read.csv(pm_path, colClasses = c(county_fips = "character"))
  pm$date <- as.Date(pm$date, format = "%Y-%m-%d")
  stopifnot(is.numeric(pm$pm25_ug_m3), !anyNA(pm$date), !anyNA(pm$county_fips),
            all(grepl("^[0-9]{5}$", pm$county_fips)),
            !anyDuplicated(pm[c("county_fips", "date")]),
            all(is.na(pm$pm25_ug_m3) | (is.finite(pm$pm25_ug_m3) & pm$pm25_ug_m3 >= 0)))
  weather <- left_join(weather, select(pm, county_fips, date, pm25_ug_m3), by = c("county_fips", "date"))
} else {
  weather$pm25_ug_m3 <- NA_real_
}
panels <- lapply(panels, function(d) {
  d |>
    left_join(weather, by = c("county_fips", "date")) |>
    arrange(county_fips, date) |>
    mutate(time_index = as.numeric(date - min(date)), year = factor(format(date, "%Y")),
           month = format(date, "%m"), dow = factor(format(date, "%u")),
           match_id = interaction(county_fips, year, month, dow),
           observed = !is.na(ohca_admissions) & !is.na(n_icu_admissions))
})

specs <- data.frame(
  model = c("ohca_count", "ohca_count_positive_icu", "icu_offset", "all_icu_count", "lag_1_5", "lag_0_10",
            "time_df3", "time_df6", "historical_weather", "local_percentiles",
            "case_crossover", "ed_pathway", "pm25_matched_baseline", "pm25_adjusted"),
  stringsAsFactors = FALSE
)
fit_site <- function(d, model) {
  lag <- switch(model, lag_1_5 = c(1, 5), lag_0_10 = c(0, 10), c(0, 5))
  lag_basis <- list(fun = "ns", knots = lag[1] + diff(lag) * 0.4, Boundary.knots = lag)
  df_year <- switch(model, time_df3 = 3L, time_df6 = 6L, 4L)
  temperature <- if (model == "historical_weather") d$historical_tmax else d$tmax_mean_c
  humidity <- if (model == "historical_weather") d$historical_rmax else d$rmax_mean_pct
  outcome <- switch(model, all_icu_count = d$n_icu_admissions,
                    ed_pathway = d$ohca_ed_admissions, d$ohca_admissions)
  if (startsWith(model, "pm25_") && !nzchar(pm_path)) stop("Awaiting daily county LGHAP PM2.5 file")
  if (!any(is.finite(temperature)) || !any(is.finite(outcome))) stop("Required input unavailable")
  cb <- crossbasis(temperature, lag = lag, group = factor(d$county_fips),
                   argvar = temp_basis, arglag = lag_basis)
  keep <- d$observed & complete.cases(cb, humidity, outcome)
  if (model %in% c("icu_offset", "ohca_count_positive_icu")) keep <- keep & d$n_icu_admissions > 0
  if (model %in% c("pm25_matched_baseline", "pm25_adjusted"))
    keep <- keep & !is.na(d$pm25_ug_m3)
  if (!any(keep) || sum(outcome[keep]) == 0) stop("No eligible event days")
  x <- d[keep, , drop = FALSE]
  x$y <- outcome[keep]
  x$humidity <- humidity[keep]
  cb_names <- paste0("temperature_", seq_len(ncol(cb)))
  x[cb_names] <- as.data.frame(unclass(cb)[keep, , drop = FALSE])
  rhs <- c(cb_names, "ns(humidity, df = 3)")
  if (model == "case_crossover") {
    if (!requireNamespace("fixest", quietly = TRUE)) stop("fixest required for case-crossover")
    f <- as.formula(paste("y ~", paste(rhs, collapse = " + "), "| match_id"))
    fit <- fixest::fepois(f, data = x, vcov = ~match_id, notes = FALSE)
    used <- fixest::obs(fit)
  } else {
    rhs <- c(rhs, sprintf("ns(time_index, df = %d)", df_year * nlevels(droplevels(x$year))))
    if (nlevels(droplevels(x$year)) > 1L) rhs <- c(rhs, "year")
    rhs <- c(rhs, "dow")
    if (length(unique(x$county_fips)) > 1L) rhs <- c(rhs, "factor(county_fips)")
    if (model == "icu_offset") rhs <- c(rhs, "offset(log(n_icu_admissions))")
    if (model == "pm25_adjusted") rhs <- c(rhs, "ns(pm25_ug_m3, df = 3)")
    f <- reformulate(rhs, response = "y")
    fit <- glm(f, data = x, family = quasipoisson(), na.action = na.fail)
    if (!fit$converged || anyNA(coef(fit))) stop("Model failed to converge or has aliased coefficients")
    used <- seq_len(nrow(x))
  }
  if (!all(cb_names %in% names(coef(fit)))) stop("Temperature basis not estimable")
  reduced <- crossreduce(cb, coef = coef(fit)[cb_names],
                         vcov = vcov(fit)[cb_names, cb_names],
                         model.link = "log", cen = cold)
  x <- x[used, , drop = FALSE]
  exposure_index <- match(paste(x$county_fips, x$date - lag[1]),
                          paste(d$county_fips, d$date))
  temp <- temperature[exposure_index]
  ref <- cold
  high <- hot
  if (model == "local_percentiles") {
    z <- temperature[d$observed & is.finite(temperature)]
    ref <- median(z)
    high <- as.numeric(quantile(z, .95))
  }
  contrast <- do.call(ns, c(list(x = c(high, ref)), temp_basis[names(temp_basis) != "fun"]))
  contrast <- contrast[1, , drop = FALSE] - contrast[2, , drop = FALSE]
  estimate <- drop(contrast %*% coef(reduced))
  se <- sqrt(drop(contrast %*% vcov(reduced) %*% t(contrast)))
  near_ref <- sum(abs(temp - ref) <= 1, na.rm = TRUE)
  near_hot <- sum(abs(temp - high) <= 1, na.rm = TRUE)
  support <- near_ref >= 10L && near_hot >= 10L &&
    min(temp, na.rm = TRUE) <= ref && max(temp, na.rm = TRUE) >= high
  lag1_acf <- NA_real_
  daily_influence <- NULL
  if (model != "case_crossover") {
    r <- data.frame(county_fips = x$county_fips, date = x$date, pearson = residuals(fit, type = "pearson"))
    next_day <- transform(r, date = date + 1)
    pairs <- inner_join(r, next_day, by = c("county_fips", "date"))
    if (nrow(pairs) > 2L) lag1_acf <- cor(pairs$pearson.x, pairs$pearson.y)
    if (model %in% c("ohca_count", "icu_offset", "all_icu_count")) {
      # Sum the contrast's estimating-function contributions within date,
      # preserving same-day dependence across counties and health systems.
      histories <- rbind(rep(high, diff(lag) + 1L), rep(ref, diff(lag) + 1L))
      contrast_cb <- crossbasis(histories, lag = lag, argvar = temp_basis, arglag = lag_basis)
      contrast_full <- setNames(numeric(length(coef(fit))), names(coef(fit)))
      contrast_full[cb_names] <- contrast_cb[1, ] - contrast_cb[2, ]
      stopifnot(isTRUE(all.equal(drop(contrast_full %*% coef(fit)), estimate, tolerance = 1e-8)))
      direction <- drop(summary(fit)$cov.unscaled %*% contrast_full)
      score <- drop(model.matrix(fit) %*% direction) * (x$y - fitted(fit))
      daily_influence <- aggregate(score ~ date, data = data.frame(date = x$date, score), FUN = sum)
    }
  }
  list(estimates = data.frame(
    site = d$site_name[[1L]], source = d$source[[1L]], model, lag_min = lag[1], lag_max = lag[2],
    n_days = nrow(x), n_ohca = sum(x$ohca_admissions), n_icu = sum(x$n_icu_admissions),
    n_outcomes = sum(x$y), reference_c = ref, hot_c = high,
    near_ref, near_hot, supported = support, log_rr = estimate, se,
    rr = exp(estimate), low = exp(estimate - 1.96 * se), high = exp(estimate + 1.96 * se),
    dispersion = if (model == "case_crossover") NA_real_ else summary(fit)$dispersion,
    lag1_residual_correlation = lag1_acf, formula = paste(deparse(f), collapse = " ")
  ), coefficients = coef(reduced), covariance = vcov(reduced),
  temperature_support = if (model %in% c("ohca_count", "all_icu_count", "icu_offset")) {
    bins <- floor(temp + 0.5)
    data.frame(site = d$site_name[[1L]], temp_c = sort(unique(bins)),
               observations = as.integer(table(bins)))
  },
  temp_range = range(temp, na.rm = TRUE), daily_influence = daily_influence)
}

fits <- list()
status <- list()
for (d in panels) {
  site <- d$site_name[[1L]]
  for (model in specs$model) {
    key <- paste(site, model, sep = ".")
    result <- tryCatch(fit_site(d, model), error = function(e) conditionMessage(e))
    if (is.character(result)) status[[key]] <- data.frame(site, model, status = result)
    else {
      fits[[key]] <- result
      status[[key]] <- data.frame(site, model, status = "ok")
    }
  }
}
est <- bind_rows(lapply(fits, "[[", "estimates"))
if (!nrow(est)) stop("No models estimable. ", paste(unique(bind_rows(status)$status), collapse = "; "))
pool <- function(x) {
  model <- rma.uni(yi = log_rr, sei = se, data = x, method = "REML", test = "adhoc")
  data.frame(k = nrow(x), sites = paste(x$site, collapse = "+"), rr = exp(as.numeric(coef(model))),
             low = exp(model$ci.lb), high = exp(model$ci.ub), tau2 = model$tau2)
}
pooled <- list()
leaveout <- list()
for (name in unique(est$model)) {
  eligible <- filter(est, model == name, supported)
  verified <- filter(eligible, source != "archived_flagship_proxy")
  if (nrow(verified) >= 2L)
    pooled[[paste0(name, "_county")]] <- cbind(model = name, set = "county_linked", pool(verified))
  if (nrow(eligible) > nrow(verified) && nrow(eligible) >= 2L)
    pooled[[paste0(name, "_proxy")]] <- cbind(model = name, set = "including_flagship_proxy", pool(eligible))
  if (nrow(eligible) >= 3L && name %in% c("ohca_count", "icu_offset")) {
    for (site in eligible$site)
      leaveout[[paste(name, site)]] <- cbind(model = name, excluded = site, pool(filter(eligible, .data$site != .env$site)))
  }
}
audit <- bind_rows(panels) |>
  group_by(site_name, county_fips, source) |>
  summarise(start = min(date), end = max(date), calendar_days = n(), observed_days = sum(observed),
            inferred_zero_icu_days = sum(inferred_zero),
            zero_ohca_days = sum(observed & ohca_admissions == 0),
            missing_tmax = sum(is.na(tmax_mean_c)), missing_rmax = sum(is.na(rmax_mean_pct)),
            missing_pm25 = sum(is.na(pm25_ug_m3)), n_ohca = sum(ohca_admissions, na.rm = TRUE),
            n_icu = sum(n_icu_admissions, na.rm = TRUE), .groups = "drop")

date_hac <- function(selected, model) {
  series <- lapply(selected, "[[", "daily_influence")
  dates <- do.call(c, lapply(series, "[[", "date"))
  calendar <- seq(min(dates), max(dates), by = "day")
  scores <- vapply(series, function(s) {
    values <- s$score[match(calendar, s$date)]
    replace(values, is.na(values), 0)
  }, numeric(length(calendar)))
  V <- crossprod(scores)
  bandwidth <- 14L
  for (lag in seq_len(bandwidth)) {
    prod <- crossprod(head(scores, -lag), tail(scores, -lag))
    V <- V + (1 - lag / (bandwidth + 1L)) * (prod + t(prod))
  }
  site_est <- bind_rows(lapply(selected, "[[", "estimates"))
  # Keep the inverse-variance weights of the working fixed-effect summary.
  w <- 1 / site_est$se^2
  w <- w / sum(w)
  beta <- sum(w * site_est$log_rr)
  se <- sqrt(drop(crossprod(w, V %*% w)))
  data.frame(model, sites = paste(site_est$site, collapse = "+"), bandwidth_days = bandwidth,
             rr = exp(beta), low = exp(beta - 1.96 * se), high = exp(beta + 1.96 * se))
}
hac <- list()
for (model_name in c("ohca_count", "icu_offset", "all_icu_count")) {
  selected <- Filter(function(x) x$estimates$model == model_name && x$estimates$supported &&
                       x$estimates$source != "archived_flagship_proxy", fits)
  if (length(selected) >= 2L) hac[[model_name]] <- date_hac(selected, model_name)
}

# Pool reduced cumulative coefficients using a parsimonious between-site covariance.
make_curve <- function(selected, set) {
  rows <- bind_rows(lapply(seq_along(selected), function(i) {
    b <- selected[[i]]$coefficients
    data.frame(site = selected[[i]]$estimates$site, component = factor(seq_along(b)), beta = b)
  }))
  V <- as.matrix(Matrix::bdiag(lapply(selected, "[[", "covariance")))
  m <- rma.mv(beta, V, mods = ~ component - 1, random = ~ component | site,
              struct = "ID", data = rows, method = "REML")
  ranges <- do.call(rbind, lapply(selected, "[[", "temp_range"))
  grid <- seq(ceiling(max(ranges[, 1])), floor(min(ranges[, 2])), by = 0.5)
  cbind(set, sites = paste(unique(rows$site), collapse = "+"),
        predict_curve(grid, coef(m), vcov(m)), tau2 = m$tau2)
}
predict_curve <- function(grid, beta, covariance) {
  b <- do.call(ns, c(list(x = grid), temp_basis[names(temp_basis) != "fun"]))
  contrasts <- sweep(b, 2, predict(b, cold)[1, ], "-")
  log_rr <- drop(contrasts %*% beta)
  se <- sqrt(rowSums((contrasts %*% covariance) * contrasts))
  data.frame(temp_c = grid, reference_c = cold, rr = exp(log_rr),
             low = exp(log_rr - 1.96 * se), high = exp(log_rr + 1.96 * se))
}
primary <- Filter(function(f) f$estimates$model == "ohca_count" && f$estimates$supported, fits)
county <- Filter(function(f) f$estimates$source != "archived_flagship_proxy", primary)
curves <- list()
if (length(county) >= 2L) curves$county <- make_curve(county, "county_linked")
if (length(primary) > length(county) && length(primary) >= 2L)
  curves$proxy <- make_curve(primary, "including_flagship_proxy")
count_fits <- Filter(function(f) f$estimates$model == "ohca_count", fits)
make_site_curves <- function(selected) bind_rows(lapply(selected, function(f) {
  grid <- seq(ceiling(f$temp_range[1] * 2) / 2, floor(f$temp_range[2] * 2) / 2, by = 0.5)
  cbind(site = f$estimates$site, source = f$estimates$source,
        predict_curve(grid, f$coefficients, f$covariance))
}))
site_curves <- make_site_curves(count_fits)
all_icu_fits <- Filter(function(f) f$estimates$model == "all_icu_count", fits)
all_icu_eligible <- Filter(function(f) f$estimates$supported &&
                            f$estimates$source != "archived_flagship_proxy", all_icu_fits)
all_icu_pooled <- if (length(all_icu_eligible) >= 2L)
  make_curve(all_icu_eligible, "county_linked") else data.frame()
offset_fits <- Filter(function(f) f$estimates$model == "icu_offset", fits)
offset_eligible <- Filter(function(f) f$estimates$supported &&
                           f$estimates$source != "archived_flagship_proxy", offset_fits)
offset_pooled <- if (length(offset_eligible) >= 2L)
  make_curve(offset_eligible, "county_linked") else data.frame()
coefficients <- bind_rows(lapply(fits, function(x) {
  b <- x$coefficients
  data.frame(site = x$estimates$site, model = x$estimates$model, term = seq_along(b), coefficient = b)
}))
covariances <- bind_rows(lapply(fits, function(x)
  cbind(site = x$estimates$site, model = x$estimates$model, as.data.frame(as.table(x$covariance)))))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
outputs <- list(site_estimates = est, model_status = bind_rows(status), coverage = audit,
                pooled_estimates = bind_rows(pooled), leave_one_site_out = bind_rows(leaveout),
                shared_date_sensitivity = bind_rows(hac),
                pooled_curve = bind_rows(curves), site_curves = site_curves,
                temperature_support = bind_rows(lapply(count_fits, "[[", "temperature_support")),
                all_icu_pooled_curve = all_icu_pooled, all_icu_site_curves = make_site_curves(all_icu_fits),
                all_icu_temperature_support = bind_rows(lapply(all_icu_fits, "[[", "temperature_support")),
                offset_pooled_curve = offset_pooled, offset_site_curves = make_site_curves(offset_fits),
                offset_temperature_support = bind_rows(lapply(offset_fits, "[[", "temperature_support")),
                reduced_coefficients = coefficients,
                reduced_covariances = covariances)
for (name in names(outputs)) write.csv(outputs[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE)
sources <- c(inputs, file.path(root, "reference/clif_hospital_geography.csv"),
             Sys.glob(file.path(root, "exposome", "daymet_county_*2018_2024*.parquet")),
             if (nzchar(pm_path)) pm_path)
if (nzchar(historical_dir)) {
  historical_sources <- file.path(historical_dir, paste0(audit$site_name, "_ohca_icu_admission_rate_daily_timeseries.csv"))
  sources <- unique(c(sources, historical_sources[file.exists(historical_sources)]))
}
write.csv(data.frame(path = normalizePath(sources), md5 = unname(tools::md5sum(sources))),
          file.path(out, "input_manifest.csv"), row.names = FALSE)
writeLines(c(sprintf("Fixed temperature contrast: %g versus %g C", hot, cold),
             "Temperature interior knots: 0, 15, 25 C; boundary knots: -40, 45 C",
             "Lag spline interior knot: min lag + 0.4 * lag range; boundaries: lag endpoints",
             "PM2.5: contemporaneous daily county mean, natural spline with 3 df",
             "Scalar pooling: REML with modified Hartung-Knapp; assumes independent sites",
             "Curve pooling: REML with common component variance and zero cross-component covariance",
             capture.output(sessionInfo())), file.path(out, "specification.txt"))
print(filter(est, model %in% c("ohca_count", "icu_offset", "all_icu_count")) |>
        select(site, source, model, n_days, n_ohca, rr, low, high))

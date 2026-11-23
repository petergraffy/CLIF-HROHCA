#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript code/post_processing/pool_aim2.R OUTPUT_DIR SITE_DIR ...")
output_dir <- args[[1L]]
sites <- args[-1L]
read_sites <- function(file) {
  bind_rows(lapply(sites, function(site) readr::read_csv(file.path(site, file), show_col_types = FALSE)))
}
estimates <- read_sites("estimates.csv")
audit <- read_sites("coverage.csv")
status <- read_sites("model_status.csv")
stopifnot(!anyDuplicated(audit$site), !anyDuplicated(select(estimates, site, model, endpoint, reference_c, hot_c)))

pool_effect <- function(x, y, se, ratio) {
  x <- x[is.finite(x[[y]]) & is.finite(x[[se]]) & x[[se]] > 0, ]
  if (nrow(x) < 2L) return(tibble())
  fit <- metafor::rma.uni(yi = x[[y]], sei = x[[se]], method = "REML", test = "adhoc")
  transform_effect <- if (ratio) exp else identity
  tibble(k = nrow(x), sites = paste(x$site, collapse = "+"),
         n = sum(x$n), estimate = transform_effect(as.numeric(coef(fit))),
         low = transform_effect(fit$ci.lb), high = transform_effect(fit$ci.ub), tau2 = fit$tau2)
}
supported <- function(d) {
  filter(d, in_exposure_range, reference_support >= 10, hot_support >= 10)
}

pool_standard <- function(d, scenario = NULL) {
  groups <- c(scenario, "model", "endpoint", "reference_c", "hot_c")
  ratios <- d |>
    group_by(across(all_of(groups))) |>
    group_modify(function(d, key) pool_effect(d, "log_ratio", "se", ratio = TRUE)) |>
    ungroup() |>
    mutate(scale = if_else(model == "cause_specific_cox", "hazard ratio", "ratio of mean VFDs"))
  differences <- d |>
    filter(endpoint == "ventilator_free_days") |>
    group_by(across(all_of(groups))) |>
    group_modify(function(d, key) pool_effect(d, "difference_days", "difference_se", ratio = FALSE)) |>
    ungroup() |>
    mutate(scale = "difference in days")
  bind_rows(ratios, differences)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(estimates, file.path(output_dir, "site_estimates.csv"))
readr::write_csv(audit, file.path(output_dir, "site_coverage.csv"))
readr::write_csv(status, file.path(output_dir, "site_model_status.csv"))
readr::write_csv(pool_standard(supported(estimates)), file.path(output_dir, "pooled_estimates.csv"))

covid_sites <- sites[file.exists(file.path(sites, "covid_era_coverage.csv"))]
if (length(covid_sites)) {
  read_covid <- function(file) {
    bind_rows(lapply(covid_sites, function(site)
      readr::read_csv(file.path(site, paste0("covid_era_", file, ".csv")), show_col_types = FALSE)))
  }
  for (file in c("coverage", "model_status")) {
    readr::write_csv(read_covid(file), file.path(output_dir, paste0("site_covid_era_", file, ".csv")))
  }
  covid <- read_covid("estimates")
  readr::write_csv(covid, file.path(output_dir, "site_covid_era_estimates.csv"))
  pooled_covid <- tibble()
  if (nrow(covid)) {
    keys <- c("site", "model", "endpoint", "reference_c", "hot_c")
    stopifnot(!anyDuplicated(select(covid, all_of(keys))))
    paired <- bind_rows(full_period = estimates, exclude_covid_era = covid, .id = "scenario") |>
      supported() |>
      filter(is.finite(log_ratio), is.finite(se), se > 0) |>
      group_by(across(all_of(keys))) |>
      filter(n_distinct(scenario) == 2L) |>
      mutate(across(any_of(c("difference_days", "difference_se")),
                    ~ if (all(is.finite(difference_days) & is.finite(difference_se) & difference_se > 0))
                      .x else NA_real_)) |>
      ungroup()
    pooled_covid <- pool_standard(paired, "scenario")
  }
  readr::write_csv(pooled_covid, file.path(output_dir, "pooled_covid_era.csv"))
}

landmark_sites <- sites[file.exists(file.path(sites, "sofa_landmark_status.csv"))]
if (length(landmark_sites)) {
  read_landmark <- function(file) {
    bind_rows(lapply(landmark_sites, function(site)
      readr::read_csv(file.path(site, paste0("sofa_landmark_", file, ".csv")), show_col_types = FALSE)))
  }
  landmark <- read_landmark("estimates")
  for (file in c("coverage", "domains", "domain_counts", "status")) {
    readr::write_csv(read_landmark(file), file.path(output_dir, paste0("site_sofa_landmark_", file, ".csv")))
  }
  readr::write_csv(landmark, file.path(output_dir, "site_sofa_landmark_estimates.csv"))

  pool_landmark <- function(standard_error, effect, scale, ratio, vfd_only = FALSE) {
    eligible <- landmark |>
      supported() |>
      filter(is.finite(.data[[effect]]), is.finite(.data[[standard_error]]), .data[[standard_error]] > 0)
    if (vfd_only) eligible <- filter(eligible, endpoint == "ventilator_free_days")
    # Compare SOFA and non-SOFA estimates across the same sites per endpoint
    # and temperature contrast as well as the same patients within each site.
    eligible |>
      group_by(site, endpoint, reference_c, hot_c) |>
      filter(model == "landmark_all" |
               (all(c("landmark_matched", "landmark_sofa_24h") %in% model))) |>
      ungroup() |>
      group_by(model, endpoint, reference_c, hot_c) |>
      group_modify(function(d, key) pool_effect(d, effect, standard_error, ratio)) |>
      ungroup() |>
      mutate(scale = if_else(endpoint != "ventilator_free_days", "hazard ratio", scale))
  }
  pooled_landmark <- if (nrow(landmark)) {
    stopifnot(!anyDuplicated(select(landmark, site, model, endpoint, reference_c, hot_c)))
    bind_rows(pool_landmark("se", "log_ratio", "ratio of mean VFDs", ratio = TRUE),
              if ("difference_se" %in% names(landmark))
                pool_landmark("difference_se", "difference_days", "difference in days", ratio = FALSE, vfd_only = TRUE))
  } else tibble()
  readr::write_csv(pooled_landmark, file.path(output_dir, "pooled_sofa_landmark.csv"))
}

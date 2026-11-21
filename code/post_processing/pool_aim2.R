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
ratios <- estimates |>
  filter(in_exposure_range, reference_support >= 10, hot_support >= 10) |>
  group_by(model, endpoint, reference_c, hot_c) |>
  group_modify(function(d, key) pool_effect(d, "log_ratio", "se", ratio = TRUE)) |>
  ungroup() |>
  mutate(scale = if_else(model == "cause_specific_cox", "hazard ratio", "ratio of mean VFDs"))
differences <- estimates |>
  filter(endpoint == "ventilator_free_days", in_exposure_range, reference_support >= 10, hot_support >= 10) |>
  group_by(model, endpoint, reference_c, hot_c) |>
  group_modify(function(d, key) pool_effect(d, "difference_days", "difference_se", ratio = FALSE)) |>
  ungroup() |>
  mutate(scale = "difference in days")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(estimates, file.path(output_dir, "site_estimates.csv"))
readr::write_csv(audit, file.path(output_dir, "site_coverage.csv"))
readr::write_csv(status, file.path(output_dir, "site_model_status.csv"))
readr::write_csv(bind_rows(ratios, differences), file.path(output_dir, "pooled_estimates.csv"))

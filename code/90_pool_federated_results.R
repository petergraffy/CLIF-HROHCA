#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) stop("Could not determine script path from commandArgs().")
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
input_dir <- file.path(repo_root, "output", "final", "federated_exports")
output_dir <- file.path(repo_root, "output", "final", "federated_pooled")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pool_der_simonian_laird <- function(est, se) {
  keep <- is.finite(est) & is.finite(se) & se >= 0
  est <- est[keep]
  se <- se[keep]
  k <- length(est)
  if (k == 0) stop("No valid site estimates to pool.")
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
    log_rr = pooled,
    log_rr_se = pooled_se,
    rr = exp(pooled),
    rr_low = exp(pooled - 1.96 * pooled_se),
    rr_high = exp(pooled + 1.96 * pooled_se),
    tau2 = tau2,
    q = q,
    i2 = ifelse(k > 1 & q > (k - 1) & q > 0, 100 * (q - (k - 1)) / q, 0),
    stringsAsFactors = FALSE
  )
}

files <- list.files(input_dir, pattern = "_dlnm_site_estimates\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No federated site estimate files found in ", input_dir)
}

site_results <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))

pooled_rows <- list()
groups <- unique(site_results[, c("stratum", "model", "reference_type")])
for (i in seq_len(nrow(groups))) {
  g <- groups[i, , drop = FALSE]
  dat <- site_results[
    site_results$stratum == g$stratum &
      site_results$model == g$model &
      site_results$reference_type == g$reference_type,
    ,
    drop = FALSE
  ]
  pooled <- pool_der_simonian_laird(dat$log_rr, dat$log_rr_se)
  pooled$stratum <- g$stratum
  pooled$model <- g$model
  pooled$reference_type <- g$reference_type
  pooled_rows[[i]] <- pooled
}

pooled_results <- do.call(rbind, pooled_rows)
pooled_results <- pooled_results[, c("stratum", "model", "reference_type", setdiff(names(pooled_results), c("stratum", "model", "reference_type")))]

write.csv(site_results, file.path(output_dir, "all_site_dlnm_estimates.csv"), row.names = FALSE)
write.csv(pooled_results, file.path(output_dir, "pooled_dlnm_random_effects_results.csv"), row.names = FALSE)

curve_files <- list.files(input_dir, pattern = "_dlnm_curves\\.csv$", full.names = TRUE)
if (length(curve_files) > 0) {
  site_curves <- do.call(rbind, lapply(curve_files, read.csv, stringsAsFactors = FALSE))

  if (!all(c("log_rr", "log_rr_se") %in% names(site_curves))) {
    site_curves$log_rr <- log(site_curves$cumulative_rr)
    site_curves$log_rr_se <- (log(site_curves$cumulative_rr_high) - log(site_curves$cumulative_rr_low)) / (2 * 1.96)
  }

  curve_groups <- unique(site_curves[, c("stratum", "model", "reference_type", "tmax_mean_c")])
  curve_rows <- list()
  for (i in seq_len(nrow(curve_groups))) {
    g <- curve_groups[i, , drop = FALSE]
    dat <- site_curves[
      site_curves$stratum == g$stratum &
        site_curves$model == g$model &
        site_curves$reference_type == g$reference_type &
        site_curves$tmax_mean_c == g$tmax_mean_c,
      ,
      drop = FALSE
    ]
    pooled <- pool_der_simonian_laird(dat$log_rr, dat$log_rr_se)
    pooled$stratum <- g$stratum
    pooled$model <- g$model
    pooled$reference_type <- g$reference_type
    pooled$tmax_mean_c <- g$tmax_mean_c
    curve_rows[[i]] <- pooled
  }

  pooled_curves <- do.call(rbind, curve_rows)
  pooled_curves <- pooled_curves[order(pooled_curves$stratum, pooled_curves$model, pooled_curves$reference_type, pooled_curves$tmax_mean_c), ]
  pooled_curves <- pooled_curves[, c("stratum", "model", "reference_type", "tmax_mean_c", setdiff(names(pooled_curves), c("stratum", "model", "reference_type", "tmax_mean_c")))]

  write.csv(site_curves, file.path(output_dir, "all_site_dlnm_curves.csv"), row.names = FALSE)
  write.csv(pooled_curves, file.path(output_dir, "pooled_dlnm_random_effects_curves.csv"), row.names = FALSE)
}

message("Wrote pooled federated DLNM results to ", output_dir)

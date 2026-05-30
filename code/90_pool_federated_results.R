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
input_dir <- file.path(repo_root, "output", "final", "federated_exports")
output_dir <- file.path(repo_root, "output", "final", "federated_pooled")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pool_der_simonian_laird <- function(est, se) {
  keep <- is.finite(est) & is.finite(se) & se >= 0
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
if ("estimable" %in% names(site_results)) {
  site_results_for_pooling <- site_results[site_results$estimable %in% c(TRUE, "TRUE", NA), , drop = FALSE]
} else {
  site_results_for_pooling <- site_results
}

pooled_rows <- list()
groups <- unique(site_results_for_pooling[, c("stratum", "model", "reference_type")])
for (i in seq_len(nrow(groups))) {
  g <- groups[i, , drop = FALSE]
  dat <- site_results_for_pooling[
    site_results_for_pooling$stratum == g$stratum &
      site_results_for_pooling$model == g$model &
      site_results_for_pooling$reference_type == g$reference_type,
    ,
    drop = FALSE
  ]
  pooled <- pool_der_simonian_laird(dat$log_rr, dat$log_rr_se)
  if (is.null(pooled)) next
  pooled$stratum <- g$stratum
  pooled$model <- g$model
  pooled$reference_type <- g$reference_type
  pooled_rows[[i]] <- pooled
}

pooled_results <- do.call(rbind, pooled_rows)
if (is.null(pooled_results)) {
  pooled_results <- data.frame()
} else {
  pooled_results <- pooled_results[, c("stratum", "model", "reference_type", setdiff(names(pooled_results), c("stratum", "model", "reference_type")))]
}

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
    if (is.null(pooled)) next
    pooled$stratum <- g$stratum
    pooled$model <- g$model
    pooled$reference_type <- g$reference_type
    pooled$tmax_mean_c <- g$tmax_mean_c
    curve_rows[[i]] <- pooled
  }

  pooled_curves <- do.call(rbind, curve_rows)
  if (is.null(pooled_curves)) {
    pooled_curves <- data.frame()
  } else {
    pooled_curves <- pooled_curves[order(pooled_curves$stratum, pooled_curves$model, pooled_curves$reference_type, pooled_curves$tmax_mean_c), ]
    pooled_curves <- pooled_curves[, c("stratum", "model", "reference_type", "tmax_mean_c", setdiff(names(pooled_curves), c("stratum", "model", "reference_type", "tmax_mean_c")))]
  }

  write.csv(site_curves, file.path(output_dir, "all_site_dlnm_curves.csv"), row.names = FALSE)
  write.csv(pooled_curves, file.path(output_dir, "pooled_dlnm_random_effects_curves.csv"), row.names = FALSE)
}

lag_files <- list.files(input_dir, pattern = "_dlnm_lag_summaries\\.csv$", full.names = TRUE)
if (length(lag_files) > 0) {
  site_lags <- do.call(rbind, lapply(lag_files, read.csv, stringsAsFactors = FALSE))

  if (!all(c("log_rr", "log_rr_se") %in% names(site_lags))) {
    site_lags$log_rr <- log(site_lags$cumulative_rr)
    site_lags$log_rr_se <- (log(site_lags$cumulative_rr_high) - log(site_lags$cumulative_rr_low)) / (2 * 1.96)
  }

  lag_groups <- unique(site_lags[, c("stratum", "model", "reference_type", "lag_start", "lag_end", "lag_label")])
  lag_rows <- list()
  for (i in seq_len(nrow(lag_groups))) {
    g <- lag_groups[i, , drop = FALSE]
    dat <- site_lags[
      site_lags$stratum == g$stratum &
        site_lags$model == g$model &
        site_lags$reference_type == g$reference_type &
        site_lags$lag_start == g$lag_start &
        site_lags$lag_end == g$lag_end,
      ,
      drop = FALSE
    ]
    pooled <- pool_der_simonian_laird(dat$log_rr, dat$log_rr_se)
    if (is.null(pooled)) next
    pooled$stratum <- g$stratum
    pooled$model <- g$model
    pooled$reference_type <- g$reference_type
    pooled$lag_start <- g$lag_start
    pooled$lag_end <- g$lag_end
    pooled$lag_label <- g$lag_label
    lag_rows[[i]] <- pooled
  }

  pooled_lags <- do.call(rbind, lag_rows)
  if (is.null(pooled_lags)) {
    pooled_lags <- data.frame()
  } else {
    pooled_lags <- pooled_lags[order(pooled_lags$stratum, pooled_lags$model, pooled_lags$reference_type, pooled_lags$lag_end), ]
    pooled_lags <- pooled_lags[, c("stratum", "model", "reference_type", "lag_start", "lag_end", "lag_label", setdiff(names(pooled_lags), c("stratum", "model", "reference_type", "lag_start", "lag_end", "lag_label")))]
  }

  write.csv(site_lags, file.path(output_dir, "all_site_dlnm_lag_summaries.csv"), row.names = FALSE)
  write.csv(pooled_lags, file.path(output_dir, "pooled_dlnm_random_effects_lag_summaries.csv"), row.names = FALSE)
}

lag_specific_files <- list.files(input_dir, pattern = "_dlnm_lag_specific_summaries\\.csv$", full.names = TRUE)
if (length(lag_specific_files) > 0) {
  site_lag_specific <- do.call(rbind, lapply(lag_specific_files, read.csv, stringsAsFactors = FALSE))

  if (!all(c("log_rr", "log_rr_se") %in% names(site_lag_specific))) {
    site_lag_specific$log_rr <- log(site_lag_specific$rr)
    site_lag_specific$log_rr_se <- (log(site_lag_specific$rr_high) - log(site_lag_specific$rr_low)) / (2 * 1.96)
  }

  lag_specific_groups <- unique(site_lag_specific[, c("stratum", "model", "reference_type", "lag", "lag_label")])
  lag_specific_rows <- list()
  for (i in seq_len(nrow(lag_specific_groups))) {
    g <- lag_specific_groups[i, , drop = FALSE]
    dat <- site_lag_specific[
      site_lag_specific$stratum == g$stratum &
        site_lag_specific$model == g$model &
        site_lag_specific$reference_type == g$reference_type &
        site_lag_specific$lag == g$lag,
      ,
      drop = FALSE
    ]
    pooled <- pool_der_simonian_laird(dat$log_rr, dat$log_rr_se)
    if (is.null(pooled)) next
    pooled$stratum <- g$stratum
    pooled$model <- g$model
    pooled$reference_type <- g$reference_type
    pooled$lag <- g$lag
    pooled$lag_label <- g$lag_label
    lag_specific_rows[[i]] <- pooled
  }

  pooled_lag_specific <- do.call(rbind, lag_specific_rows)
  if (is.null(pooled_lag_specific)) {
    pooled_lag_specific <- data.frame()
  } else {
    pooled_lag_specific <- pooled_lag_specific[order(pooled_lag_specific$stratum, pooled_lag_specific$model, pooled_lag_specific$reference_type, pooled_lag_specific$lag), ]
    pooled_lag_specific <- pooled_lag_specific[, c("stratum", "model", "reference_type", "lag", "lag_label", setdiff(names(pooled_lag_specific), c("stratum", "model", "reference_type", "lag", "lag_label")))]
  }

  write.csv(site_lag_specific, file.path(output_dir, "all_site_dlnm_lag_specific_summaries.csv"), row.names = FALSE)
  write.csv(pooled_lag_specific, file.path(output_dir, "pooled_dlnm_random_effects_lag_specific_summaries.csv"), row.names = FALSE)
}

message("Wrote pooled federated DLNM results to ", output_dir)

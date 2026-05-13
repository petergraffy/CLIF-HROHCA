#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    ofiles <- vapply(sys.frames(), function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile, character(1))
    ofiles <- stats::na.omit(ofiles)
    if (length(ofiles) > 0) return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = TRUE))
    stop("Could not determine script path. Run with Rscript.")
  }
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
box_root <- "/Users/saborpete/Library/CloudStorage/Box-Box/CLIF/Projects/CLIF-Heat-Related-OHCA"
output_dir <- file.path(repo_root, "output", "final", "federated_pooled")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_csv_base <- function(path) {
  out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"site_name" %in% names(out)) {
    out$site_name <- sub("_.*$", "", basename(path))
  }
  out$source_file <- path
  out
}

find_exports <- function(suffix) {
  files <- list.files(box_root, pattern = paste0("^[^_]+_", suffix, "$"), recursive = TRUE, full.names = TRUE)
  files[grepl("/federated_exports/", files)]
}

find_nonstandard_exports <- function(suffix) {
  jhu_outcomes_dir <- file.path(box_root, "JHU", "final", "ohca_outcomes")
  path <- switch(
    suffix,
    "adverse_outcome_models.csv" = file.path(jhu_outcomes_dir, "ohca_heat_adverse_outcome_models.csv"),
    "continuous_outcome_models.csv" = file.path(jhu_outcomes_dir, "ohca_heat_continuous_outcome_models.csv"),
    "pollution_12m_binary_outcome_models.csv" = file.path(jhu_outcomes_dir, "ohca_pollution_12m_binary_outcome_models.csv"),
    "pollution_12m_continuous_outcome_models.csv" = file.path(jhu_outcomes_dir, "ohca_pollution_12m_continuous_outcome_models.csv"),
    NA_character_
  )
  path <- path[!is.na(path) & file.exists(path)]
  path
}

bind_export <- function(suffix) {
  files <- c(find_exports(suffix), find_nonstandard_exports(suffix))
  if (length(files) == 0) return(data.frame())
  dfs <- lapply(files, read_csv_base)
  dfs <- Map(function(x, path) {
    if (grepl("/JHU/final/ohca_outcomes/", path)) x$site_name <- "JHU"
    x
  }, dfs, files)
  all_names <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(x) {
    missing_names <- setdiff(all_names, names(x))
    for (nm in missing_names) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, dfs)
}

se_from_ci <- function(est, low, high) {
  (log(high) - log(low)) / (2 * 1.96)
}

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
    log_estimate = pooled,
    log_estimate_se = pooled_se,
    ratio = exp(pooled),
    ci_low = exp(pooled - 1.96 * pooled_se),
    ci_high = exp(pooled + 1.96 * pooled_se),
    p_value = 2 * stats::pnorm(abs(pooled / pooled_se), lower.tail = FALSE),
    tau2 = tau2,
    q = q,
    i2 = ifelse(k > 1 && q > (k - 1) && q > 0, 100 * (q - (k - 1)) / q, 0),
    stringsAsFactors = FALSE
  )
}

pool_ratio_table <- function(dat, group_cols, ratio_col, output_measure, low_col = "ci_low", high_col = "ci_high") {
  if (nrow(dat) == 0) return(data.frame())
  if ("estimable" %in% names(dat)) {
    dat <- dat[dat$estimable %in% c(TRUE, "TRUE", NA), , drop = FALSE]
  }
  dat[[ratio_col]] <- suppressWarnings(as.numeric(dat[[ratio_col]]))
  dat[[low_col]] <- suppressWarnings(as.numeric(dat[[low_col]]))
  dat[[high_col]] <- suppressWarnings(as.numeric(dat[[high_col]]))
  dat$log_estimate <- log(dat[[ratio_col]])
  dat$log_estimate_se <- se_from_ci(dat[[ratio_col]], dat[[low_col]], dat[[high_col]])
  groups <- unique(dat[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(dat))
    for (col in group_cols) idx <- idx & dat[[col]] == groups[[col]][i]
    pooled <- pool_der_simonian_laird(dat$log_estimate[idx], dat$log_estimate_se[idx])
    if (is.null(pooled)) next
    rows[[i]] <- cbind(groups[i, , drop = FALSE], pooled)
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out$measure <- output_measure
  out[, c(group_cols, "measure", setdiff(names(out), c(group_cols, "measure")))]
}

parse_count_pct <- function(x) {
  x <- gsub(",", "", x)
  n <- suppressWarnings(as.numeric(sub("^([0-9.]+).*", "\\1", x)))
  pct <- suppressWarnings(as.numeric(sub("^.*\\(([0-9.]+)%\\).*$", "\\1", x)))
  denom <- ifelse(is.finite(n) & is.finite(pct) & pct > 0, round(n / (pct / 100)), NA_real_)
  data.frame(n = n, pct = pct, denom = denom)
}

fmt_ratio <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
normalize_pollution_exposure <- function(x) sub(" per site IQR \\([^)]*\\)", " per site IQR", x)

all_exports <- c(
  "heat_related_vs_non_heat_related_table.csv",
  "all_year_heat_related_vs_non_heat_related_table.csv",
  "heat_related_vs_non_heat_related_table_all_definitions.csv",
  "all_year_heat_related_vs_non_heat_related_table_all_definitions.csv",
  "heat_related_vs_non_heat_related_outcomes.csv",
  "all_year_heat_related_vs_non_heat_related_outcomes.csv",
  "dlnm_site_estimates.csv",
  "dlnm_curves.csv",
  "adverse_outcome_models.csv",
  "continuous_outcome_models.csv",
  "pollution_12m_binary_outcome_models.csv",
  "pollution_12m_continuous_outcome_models.csv",
  "heat_related_hourly_lab_trajectories.csv",
  "heat_related_hourly_vital_trajectories.csv",
  "heat_related_hourly_support_trajectories.csv",
  "heat_related_renal_metabolic_marker_summary.csv",
  "heat_related_crrt_window_summary.csv",
  "all_year_heat_related_hourly_lab_trajectories.csv",
  "all_year_heat_related_hourly_vital_trajectories.csv",
  "all_year_heat_related_hourly_support_trajectories.csv",
  "all_year_heat_related_renal_metabolic_marker_summary.csv",
  "all_year_heat_related_crrt_window_summary.csv",
  "outcomes.csv"
)

combined <- lapply(all_exports, bind_export)
names(combined) <- sub("\\.csv$", "", all_exports)
for (nm in names(combined)) {
  if (nrow(combined[[nm]]) > 0) {
    write.csv(combined[[nm]], file.path(output_dir, paste0("all_sites_", nm, ".csv")), row.names = FALSE)
  }
}

dlnm <- combined$dlnm_site_estimates
dlnm_pooled <- pool_ratio_table(
  dlnm,
  c("stratum", "model", "reference_type"),
  "cumulative_rr",
  "cumulative incidence RR",
  low_col = "cumulative_rr_low",
  high_col = "cumulative_rr_high"
)
write.csv(dlnm_pooled, file.path(output_dir, "pooled_dlnm_random_effects_results.csv"), row.names = FALSE)

dlnm_curves <- combined$dlnm_curves
if (nrow(dlnm_curves) > 0) {
  if (!all(c("log_rr", "log_rr_se") %in% names(dlnm_curves))) {
    dlnm_curves$log_rr <- log(suppressWarnings(as.numeric(dlnm_curves$cumulative_rr)))
    dlnm_curves$log_rr_se <- se_from_ci(
      suppressWarnings(as.numeric(dlnm_curves$cumulative_rr)),
      suppressWarnings(as.numeric(dlnm_curves$cumulative_rr_low)),
      suppressWarnings(as.numeric(dlnm_curves$cumulative_rr_high))
    )
  }
  dlnm_curves$tmax_mean_c <- suppressWarnings(as.numeric(dlnm_curves$tmax_mean_c))
  curve_groups <- unique(dlnm_curves[, c("stratum", "model", "reference_type", "tmax_mean_c")])
  curve_rows <- vector("list", nrow(curve_groups))
  for (i in seq_len(nrow(curve_groups))) {
    idx <- rep(TRUE, nrow(dlnm_curves))
    for (col in c("stratum", "model", "reference_type", "tmax_mean_c")) {
      idx <- idx & dlnm_curves[[col]] == curve_groups[[col]][i]
    }
    pooled <- pool_der_simonian_laird(
      suppressWarnings(as.numeric(dlnm_curves$log_rr[idx])),
      suppressWarnings(as.numeric(dlnm_curves$log_rr_se[idx]))
    )
    if (is.null(pooled)) next
    curve_rows[[i]] <- cbind(curve_groups[i, , drop = FALSE], pooled)
  }
  pooled_dlnm_curves <- do.call(rbind, curve_rows)
  pooled_dlnm_curves <- pooled_dlnm_curves[order(
    pooled_dlnm_curves$stratum,
    pooled_dlnm_curves$model,
    pooled_dlnm_curves$reference_type,
    pooled_dlnm_curves$tmax_mean_c
  ), ]
  names(pooled_dlnm_curves)[names(pooled_dlnm_curves) == "ratio"] <- "cumulative_rr"
  names(pooled_dlnm_curves)[names(pooled_dlnm_curves) == "ci_low"] <- "cumulative_rr_low"
  names(pooled_dlnm_curves)[names(pooled_dlnm_curves) == "ci_high"] <- "cumulative_rr_high"
  write.csv(dlnm_curves, file.path(output_dir, "all_sites_dlnm_curves.csv"), row.names = FALSE)
  write.csv(pooled_dlnm_curves, file.path(output_dir, "pooled_dlnm_random_effects_curves.csv"), row.names = FALSE)

  plot_curve <- function(dat, model_name, reference_type, out_file, title) {
    x <- dat[dat$stratum == "Overall" & dat$model == model_name & dat$reference_type == reference_type, , drop = FALSE]
    if (nrow(x) == 0) return(invisible(FALSE))
    x <- x[order(x$tmax_mean_c), ]
    png(out_file, width = 1600, height = 1000, res = 160)
    ylim <- range(c(x$cumulative_rr_low, x$cumulative_rr_high), finite = TRUE)
    plot(
      x$tmax_mean_c,
      x$cumulative_rr,
      type = "l",
      lwd = 3,
      col = "#1f77b4",
      ylim = ylim,
      xlab = "Daily mean Tmax (C)",
      ylab = "Cumulative RR",
      main = title
    )
    polygon(
      c(x$tmax_mean_c, rev(x$tmax_mean_c)),
      c(x$cumulative_rr_low, rev(x$cumulative_rr_high)),
      col = grDevices::adjustcolor("#1f77b4", alpha.f = 0.18),
      border = NA
    )
    lines(x$tmax_mean_c, x$cumulative_rr, lwd = 3, col = "#1f77b4")
    abline(h = 1, lty = 2, col = "gray40")
    dev.off()
    invisible(TRUE)
  }
  plot_curve(
    pooled_dlnm_curves,
    "primary_humidity_adjusted",
    "median",
    file.path(output_dir, "pooled_dlnm_curve_overall_primary_median.png"),
    "Pooled DLNM Curve: Overall Primary Model"
  )
  plot_curve(
    pooled_dlnm_curves,
    "sensitivity_mrt_reference",
    "mrt",
    file.path(output_dir, "pooled_dlnm_curve_overall_mrt_reference.png"),
    "Pooled DLNM Curve: Overall MRT Reference Sensitivity"
  )
}

adverse <- combined$adverse_outcome_models
adverse_pooled <- pool_ratio_table(
  adverse,
  c("outcome", "exposure"),
  "odds_ratio",
  "odds ratio"
)
write.csv(adverse_pooled, file.path(output_dir, "pooled_adverse_outcome_models.csv"), row.names = FALSE)

continuous <- combined$continuous_outcome_models
continuous_pooled <- pool_ratio_table(
  continuous,
  c("outcome", "exposure"),
  "geometric_mean_ratio",
  "geometric mean ratio"
)
write.csv(continuous_pooled, file.path(output_dir, "pooled_continuous_outcome_models.csv"), row.names = FALSE)

pollution_binary <- combined$pollution_12m_binary_outcome_models
if (nrow(pollution_binary) > 0) {
  pollution_binary$exposure_family <- normalize_pollution_exposure(pollution_binary$exposure)
}
pollution_binary_pooled <- pool_ratio_table(
  pollution_binary,
  c("outcome", "exposure_family", "adjustment_set"),
  "odds_ratio",
  "odds ratio"
)
write.csv(pollution_binary_pooled, file.path(output_dir, "pooled_pollution_12m_binary_outcome_models.csv"), row.names = FALSE)

pollution_continuous <- combined$pollution_12m_continuous_outcome_models
if (nrow(pollution_continuous) > 0) {
  pollution_continuous$exposure_family <- normalize_pollution_exposure(pollution_continuous$exposure)
}
pollution_continuous_pooled <- pool_ratio_table(
  pollution_continuous,
  c("outcome", "exposure_family", "adjustment_set"),
  "geometric_mean_ratio",
  "geometric mean ratio"
)
write.csv(pollution_continuous_pooled, file.path(output_dir, "pooled_pollution_12m_continuous_outcome_models.csv"), row.names = FALSE)

outcomes <- combined$heat_related_vs_non_heat_related_outcomes
if (nrow(outcomes) > 0) {
  pct_cols <- grep("_pct$", names(outcomes), value = TRUE)
  outcome_long <- do.call(rbind, lapply(pct_cols, function(col) {
    data.frame(
      site_name = outcomes$site_name,
      heat_definition = outcomes$heat_definition,
      heat_related_ohca = outcomes$heat_related_ohca,
      outcome = sub("_pct$", "", col),
      events = round(outcomes$n * suppressWarnings(as.numeric(outcomes[[col]])) / 100),
      n = suppressWarnings(as.numeric(outcomes$n)),
      stringsAsFactors = FALSE
    )
  }))
  outcome_counts <- aggregate(cbind(events, n) ~ outcome + heat_definition + heat_related_ohca + site_name, data = outcome_long, FUN = sum)
  outcome_counts$rate_pct <- 100 * outcome_counts$events / outcome_counts$n
  pooled_outcomes <- aggregate(cbind(events, n) ~ outcome + heat_definition + heat_related_ohca, data = outcome_long, FUN = sum)
  pooled_outcomes$rate_pct <- 100 * pooled_outcomes$events / pooled_outcomes$n
  write.csv(outcome_counts, file.path(output_dir, "site_heat_related_vs_non_heat_related_outcome_counts.csv"), row.names = FALSE)
  write.csv(pooled_outcomes, file.path(output_dir, "pooled_heat_related_vs_non_heat_related_outcome_counts.csv"), row.names = FALSE)
}

table_heat95 <- combined$heat_related_vs_non_heat_related_table
binary_rows <- data.frame()
if (nrow(table_heat95) > 0) {
  has_count <- grepl("\\([0-9.]+%\\)", table_heat95$heat_related_ohca) &
    grepl("\\([0-9.]+%\\)", table_heat95$non_heat_related_ohca)
  binary_rows <- table_heat95[has_count, c("site_name", "characteristic", "heat_related_ohca", "non_heat_related_ohca")]
  heat_parsed <- parse_count_pct(binary_rows$heat_related_ohca)
  non_heat_parsed <- parse_count_pct(binary_rows$non_heat_related_ohca)
  binary_rows$heat_events <- heat_parsed$n
  binary_rows$heat_n <- heat_parsed$denom
  binary_rows$non_heat_events <- non_heat_parsed$n
  binary_rows$non_heat_n <- non_heat_parsed$denom
  binary_rows$characteristic_norm <- tolower(trimws(binary_rows$characteristic))
  pooled_binary <- aggregate(
    cbind(heat_events, heat_n, non_heat_events, non_heat_n) ~ characteristic_norm,
    data = binary_rows,
    FUN = sum,
    na.rm = TRUE
  )
  names(pooled_binary)[names(pooled_binary) == "characteristic_norm"] <- "characteristic"
  pooled_binary$heat_pct <- 100 * pooled_binary$heat_events / pooled_binary$heat_n
  pooled_binary$non_heat_pct <- 100 * pooled_binary$non_heat_events / pooled_binary$non_heat_n
  pooled_binary$absolute_difference_pct <- pooled_binary$heat_pct - pooled_binary$non_heat_pct
  pooled_binary <- pooled_binary[order(abs(pooled_binary$absolute_difference_pct), decreasing = TRUE), ]
  write.csv(pooled_binary, file.path(output_dir, "pooled_heat95_binary_characteristic_differences.csv"), row.names = FALSE)
}

weighted_group_summary <- function(dat, group_cols, value_col, weight_col = "n_patients") {
  if (nrow(dat) == 0) return(data.frame())
  dat[[value_col]] <- suppressWarnings(as.numeric(dat[[value_col]]))
  dat[[weight_col]] <- suppressWarnings(as.numeric(dat[[weight_col]]))
  dat <- dat[is.finite(dat[[value_col]]) & is.finite(dat[[weight_col]]) & dat[[weight_col]] > 0, , drop = FALSE]
  if (nrow(dat) == 0) return(data.frame())
  groups <- unique(dat[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(dat))
    for (col in group_cols) idx <- idx & dat[[col]] == groups[[col]][i]
    x <- dat[idx, , drop = FALSE]
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        n_sites = length(unique(x$site_name)),
        n_patients = sum(x[[weight_col]], na.rm = TRUE),
        weighted_median_value = stats::weighted.mean(x[[value_col]], x[[weight_col]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  do.call(rbind, rows)
}

safe_prop_p <- function(x1, n1, x0, n0) {
  if (!all(is.finite(c(x1, n1, x0, n0))) || n1 <= 0 || n0 <= 0) return(NA_real_)
  suppressWarnings(stats::prop.test(c(x1, x0), c(n1, n0), correct = FALSE)$p.value)
}

safe_norm_p <- function(z) {
  if (!is.finite(z)) return(NA_real_)
  2 * stats::pnorm(abs(z), lower.tail = FALSE)
}

add_approx_site_median_se <- function(dat) {
  dat$median_value <- suppressWarnings(as.numeric(dat$median_value))
  dat$q25_value <- suppressWarnings(as.numeric(dat$q25_value))
  dat$q75_value <- suppressWarnings(as.numeric(dat$q75_value))
  dat$n_patients <- suppressWarnings(as.numeric(dat$n_patients))
  dat$median_se <- 1.57 * (dat$q75_value - dat$q25_value) / sqrt(dat$n_patients)
  dat
}

pool_site_median_differences <- function(dat, group_cols) {
  if (nrow(dat) == 0 || !all(c("median_value", "q25_value", "q75_value", "n_patients") %in% names(dat))) {
    return(data.frame())
  }
  dat <- add_approx_site_median_se(dat)
  key_cols <- c(group_cols, "site_name")
  heat <- dat[dat$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- dat[dat$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  names(heat)[names(heat) == "median_value"] <- "heat_site_median"
  names(heat)[names(heat) == "median_se"] <- "heat_site_median_se"
  names(heat)[names(heat) == "n_patients"] <- "heat_site_n"
  names(non_heat)[names(non_heat) == "median_value"] <- "non_heat_site_median"
  names(non_heat)[names(non_heat) == "median_se"] <- "non_heat_site_median_se"
  names(non_heat)[names(non_heat) == "n_patients"] <- "non_heat_site_n"
  site_diffs <- merge(
    heat[, c(key_cols, "heat_site_median", "heat_site_median_se", "heat_site_n"), drop = FALSE],
    non_heat[, c(key_cols, "non_heat_site_median", "non_heat_site_median_se", "non_heat_site_n"), drop = FALSE],
    by = key_cols,
    all = FALSE
  )
  if (nrow(site_diffs) == 0) return(data.frame())
  site_diffs$site_difference <- site_diffs$heat_site_median - site_diffs$non_heat_site_median
  site_diffs$site_difference_se <- sqrt(site_diffs$heat_site_median_se^2 + site_diffs$non_heat_site_median_se^2)
  site_diffs <- site_diffs[is.finite(site_diffs$site_difference) & is.finite(site_diffs$site_difference_se) & site_diffs$site_difference_se > 0, , drop = FALSE]
  if (nrow(site_diffs) == 0) return(data.frame())
  groups <- unique(site_diffs[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(site_diffs))
    for (col in group_cols) idx <- idx & site_diffs[[col]] == groups[[col]][i]
    pooled <- pool_der_simonian_laird(site_diffs$site_difference[idx], site_diffs$site_difference_se[idx])
    if (is.null(pooled)) next
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        approx_difference = pooled$log_estimate,
        approx_difference_se = pooled$log_estimate_se,
        approx_ci_low = pooled$log_estimate - 1.96 * pooled$log_estimate_se,
        approx_ci_high = pooled$log_estimate + 1.96 * pooled$log_estimate_se,
        approx_p_value = pooled$p_value,
        approx_i2 = pooled$i2,
        approx_tau2 = pooled$tau2,
        approx_k_sites = pooled$k_sites,
        stringsAsFactors = FALSE
      )
    )
  }
  do.call(rbind, rows)
}

make_two_group_difference <- function(summary_dat, group_cols, value_col, weight_col) {
  if (nrow(summary_dat) == 0) return(data.frame())
  pooled <- weighted_group_summary(
    summary_dat,
    c(group_cols, "heat_related_ohca"),
    value_col,
    weight_col
  )
  if (nrow(pooled) == 0) return(data.frame())
  heat <- pooled[pooled$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- pooled[pooled$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  names(heat)[names(heat) %in% c("n_sites", "n_patients", "weighted_median_value")] <-
    c("heat_n_sites", "heat_n_patients", "heat_weighted_median")
  names(non_heat)[names(non_heat) %in% c("n_sites", "n_patients", "weighted_median_value")] <-
    c("non_heat_n_sites", "non_heat_n_patients", "non_heat_weighted_median")
  out <- merge(
    heat[, c(group_cols, "heat_n_sites", "heat_n_patients", "heat_weighted_median"), drop = FALSE],
    non_heat[, c(group_cols, "non_heat_n_sites", "non_heat_n_patients", "non_heat_weighted_median"), drop = FALSE],
    by = group_cols,
    all = FALSE
  )
  out$difference_heat_minus_non_heat <- out$heat_weighted_median - out$non_heat_weighted_median
  approx <- pool_site_median_differences(summary_dat, group_cols)
  if (nrow(approx) > 0) {
    out <- merge(out, approx, by = group_cols, all.x = TRUE)
    out$p_value_method <- "Approximate random-effects meta-analysis of site median differences; site SEs estimated from IQR and n"
  }
  out
}

pool_crrt_differences <- function(dat) {
  if (nrow(dat) == 0) return(data.frame())
  dat$n_group <- suppressWarnings(as.numeric(dat$n_group))
  dat$n_crrt <- suppressWarnings(as.numeric(dat$n_crrt))
  counts <- aggregate(
    cbind(n_group, n_crrt) ~ heat_definition + window + heat_related_ohca,
    data = dat,
    FUN = sum,
    na.rm = TRUE
  )
  counts$crrt_pct <- 100 * counts$n_crrt / counts$n_group
  heat <- counts[counts$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- counts[counts$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  names(heat)[names(heat) %in% c("n_group", "n_crrt", "crrt_pct")] <- c("heat_n", "heat_crrt", "heat_crrt_pct")
  names(non_heat)[names(non_heat) %in% c("n_group", "n_crrt", "crrt_pct")] <- c("non_heat_n", "non_heat_crrt", "non_heat_crrt_pct")
  out <- merge(
    heat[, c("heat_definition", "window", "heat_n", "heat_crrt", "heat_crrt_pct")],
    non_heat[, c("heat_definition", "window", "non_heat_n", "non_heat_crrt", "non_heat_crrt_pct")],
    by = c("heat_definition", "window"),
    all = FALSE
  )
  out$absolute_difference_pct <- out$heat_crrt_pct - out$non_heat_crrt_pct
  out$p_value <- mapply(safe_prop_p, out$heat_crrt, out$heat_n, out$non_heat_crrt, out$non_heat_n)
  out$p_value_method <- "Pooled two-sample test for equality of proportions"
  out
}

lab_diffs <- make_two_group_difference(
  combined$heat_related_hourly_lab_trajectories,
  c("heat_definition", "icu_hour", "variable"),
  "median_value",
  "n_patients"
)
write.csv(lab_diffs, file.path(output_dir, "pooled_heat_related_hourly_lab_median_differences.csv"), row.names = FALSE)

vital_diffs <- make_two_group_difference(
  combined$heat_related_hourly_vital_trajectories,
  c("heat_definition", "icu_hour", "variable"),
  "median_value",
  "n_patients"
)
write.csv(vital_diffs, file.path(output_dir, "pooled_heat_related_hourly_vital_median_differences.csv"), row.names = FALSE)

support <- combined$heat_related_hourly_support_trajectories
if (nrow(support) > 0) {
  support$n_at_risk <- suppressWarnings(as.numeric(support$n_at_risk))
  support$n_event <- suppressWarnings(as.numeric(support$n_event))
  support_counts <- aggregate(
    cbind(n_at_risk, n_event) ~ heat_definition + icu_hour + variable + heat_related_ohca,
    data = support,
    FUN = sum,
    na.rm = TRUE
  )
  support_counts$prevalence_pct <- 100 * support_counts$n_event / support_counts$n_at_risk
  heat <- support_counts[support_counts$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- support_counts[support_counts$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  names(heat)[names(heat) %in% c("n_at_risk", "n_event", "prevalence_pct")] <- c("heat_n_at_risk", "heat_n_event", "heat_prevalence_pct")
  names(non_heat)[names(non_heat) %in% c("n_at_risk", "n_event", "prevalence_pct")] <- c("non_heat_n_at_risk", "non_heat_n_event", "non_heat_prevalence_pct")
  support_diffs <- merge(
    heat[, c("heat_definition", "icu_hour", "variable", "heat_n_at_risk", "heat_n_event", "heat_prevalence_pct")],
    non_heat[, c("heat_definition", "icu_hour", "variable", "non_heat_n_at_risk", "non_heat_n_event", "non_heat_prevalence_pct")],
    by = c("heat_definition", "icu_hour", "variable"),
    all = FALSE
  )
  support_diffs$absolute_difference_pct <- support_diffs$heat_prevalence_pct - support_diffs$non_heat_prevalence_pct
  support_diffs$p_value <- mapply(
    safe_prop_p,
    support_diffs$heat_n_event,
    support_diffs$heat_n_at_risk,
    support_diffs$non_heat_n_event,
    support_diffs$non_heat_n_at_risk
  )
  support_diffs$p_value_method <- "Pooled two-sample test for equality of proportions"
  write.csv(support_diffs, file.path(output_dir, "pooled_heat_related_hourly_support_prevalence_differences.csv"), row.names = FALSE)
}

renal_diffs <- make_two_group_difference(
  combined$heat_related_renal_metabolic_marker_summary,
  c("heat_definition", "window", "marker", "direction"),
  "median_value",
  "n_patients"
)
write.csv(renal_diffs, file.path(output_dir, "pooled_heat_related_renal_metabolic_marker_median_differences.csv"), row.names = FALSE)

crrt_diffs <- pool_crrt_differences(combined$heat_related_crrt_window_summary)
write.csv(crrt_diffs, file.path(output_dir, "pooled_heat_related_crrt_window_differences.csv"), row.names = FALSE)

sites <- sort(unique(unlist(lapply(combined, function(x) if (nrow(x) > 0) x$site_name else character()))))
site_availability <- data.frame(site_name = sites)
for (nm in names(combined)) {
  present <- unique(if (nrow(combined[[nm]]) > 0) combined[[nm]]$site_name else character())
  site_availability[[nm]] <- site_availability$site_name %in% present
}
write.csv(site_availability, file.path(output_dir, "site_export_availability.csv"), row.names = FALSE)

primary_dlnm <- dlnm_pooled[
  dlnm_pooled$stratum == "Overall" &
    dlnm_pooled$model == "primary_humidity_adjusted" &
    dlnm_pooled$reference_type == "median",
  ,
  drop = FALSE
]
heat95_adverse <- adverse_pooled[adverse_pooled$exposure == "Tmax >= warm-season 95th percentile", , drop = FALSE]
tmax_adverse <- adverse_pooled[adverse_pooled$exposure == "Tmax per 5 C", , drop = FALSE]
heat95_continuous <- continuous_pooled[continuous_pooled$exposure == "Tmax >= warm-season 95th percentile", , drop = FALSE]
tmax_continuous <- continuous_pooled[continuous_pooled$exposure == "Tmax per 5 C", , drop = FALSE]
pollution_binary_single <- pollution_binary_pooled[pollution_binary_pooled$adjustment_set == "single_pollutant", , drop = FALSE]
pollution_continuous_single <- pollution_continuous_pooled[pollution_continuous_pooled$adjustment_set == "single_pollutant", , drop = FALSE]

report <- c(
  "# Multisite CLIF Heat-Related OHCA Pooled Results",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "## Included Site Outputs",
  "",
  paste0("- Sites found: ", paste(sites, collapse = ", ")),
  paste0("- Coefficient-based heat outcome model exports were available for ",
         length(unique(adverse$site_name)), " sites; estimable site counts vary by outcome/exposure."),
  "",
  "## Temperature-OHCA Association",
  "",
  if (nrow(primary_dlnm) == 1) {
    paste0(
      "- Overall primary humidity-adjusted DLNM, hot temperature versus site median: RR ",
      fmt_ratio(primary_dlnm$ratio), " (95% CI ",
      fmt_ratio(primary_dlnm$ci_low), "-", fmt_ratio(primary_dlnm$ci_high),
      "), p=", fmt_p(primary_dlnm$p_value), ", I2=",
      fmt_ratio(primary_dlnm$i2, 1), "% across ", primary_dlnm$k_sites, " sites."
    )
  } else {
    "- Overall primary DLNM result was not available."
  },
  "",
  "## Heat-Related Versus Non-Heat-Related OHCA",
  "",
  "- The pooled descriptive files preserve site-level privacy, so continuous characteristics are summarized by site rather than patient-level pooled medians.",
  "- Binary/categorical characteristic counts were summed where `n (%)` cells were available; see `pooled_heat95_binary_characteristic_differences.csv` for the largest absolute differences.",
  "",
  "### Pooled Outcome Model Results",
  "",
  "Heat indicator models use Tmax >= warm-season 95th percentile; per-temperature models use Tmax per 5 C.",
  "",
  "```",
  capture.output(print(heat95_adverse[, c("outcome", "measure", "ratio", "ci_low", "ci_high", "p_value", "k_sites", "i2")], row.names = FALSE)),
  capture.output(print(tmax_adverse[, c("outcome", "measure", "ratio", "ci_low", "ci_high", "p_value", "k_sites", "i2")], row.names = FALSE)),
  capture.output(print(heat95_continuous[, c("outcome", "measure", "ratio", "ci_low", "ci_high", "p_value", "k_sites", "i2")], row.names = FALSE)),
  capture.output(print(tmax_continuous[, c("outcome", "measure", "ratio", "ci_low", "ci_high", "p_value", "k_sites", "i2")], row.names = FALSE)),
  "```",
  "",
  "### Pooled 12-Month Pollution Outcome Models",
  "",
  "Pollution estimates are pooled by pollutant per site-specific IQR; site-specific IQR values are preserved in the all-site exports.",
  "",
  "```",
  capture.output(print(pollution_binary_single[, c("outcome", "exposure_family", "measure", "ratio", "ci_low", "ci_high", "p_value", "k_sites", "i2")], row.names = FALSE)),
  capture.output(print(pollution_continuous_single[, c("outcome", "exposure_family", "measure", "ratio", "ci_low", "ci_high", "p_value", "k_sites", "i2")], row.names = FALSE)),
  "```",
  "",
  "## Output Files",
  "",
  "- `all_sites_*.csv`: harmonized row-bound site exports.",
  "- `pooled_*models.csv`: random-effects pooled coefficient-based associations.",
  "- `pooled_heat_related_vs_non_heat_related_outcome_counts.csv`: summed event counts and rates by heat definition/group.",
  "- `pooled_heat95_binary_characteristic_differences.csv`: summed binary/categorical differences between HROHCA and non-HROHCA for the primary heat95 definition.",
  "- `site_export_availability.csv`: which exports were available by site."
)

writeLines(report, file.path(output_dir, "multisite_hrohca_summary.md"))
message("Wrote multisite pooled outputs to ", output_dir)

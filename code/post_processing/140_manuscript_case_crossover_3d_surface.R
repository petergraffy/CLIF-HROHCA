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

repo_root <- normalizePath(file.path(dirname(get_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

surface_path <- file.path(pooled_dir, "pooled_dlnm_random_effects_lag_temperature_surface.csv")
if (!file.exists(surface_path)) {
  stop("Missing pooled lag-temperature surface. Run code/post_processing/90_pool_federated_results.R first.", call. = FALSE)
}

surface <- read.csv(surface_path, stringsAsFactors = FALSE)
plot_dat <- surface[
  surface$export_family == "case_crossover_dlnm" &
    surface$stratum == "Overall" &
    surface$model == "primary_humidity_adjusted" &
    surface$reference_type == "median" &
    is.finite(surface$tmax_mean_c) &
    is.finite(surface$lag) &
    is.finite(surface$rr),
  ,
  drop = FALSE
]

if (nrow(plot_dat) == 0L) {
  stop("Could not find the pooled overall case-crossover median-reference humidity-adjusted lag-temperature surface.", call. = FALSE)
}

plot_dat$tmax_f <- plot_dat$tmax_mean_c * 9 / 5 + 32
EDGE_TRIM_F <- 2
edge_trim_c <- EDGE_TRIM_F * 5 / 9
plot_dat <- plot_dat[plot_dat$tmax_mean_c >= -20 + edge_trim_c & plot_dat$tmax_mean_c <= 38 - edge_trim_c, , drop = FALSE]
plot_dat <- plot_dat[order(plot_dat$lag, plot_dat$tmax_mean_c), , drop = FALSE]

temps_f <- sort(unique(plot_dat$tmax_f))
lags <- sort(unique(plot_dat$lag))
z <- matrix(NA_real_, nrow = length(temps_f), ncol = length(lags))
z[cbind(match(plot_dat$tmax_f, temps_f), match(plot_dat$lag, lags))] <- plot_dat$rr

if (all(is.na(z))) stop("The pooled case-crossover surface contains no finite relative-risk values.", call. = FALSE)

source_table <- plot_dat[, c(
  "export_family", "analysis_period", "stratum", "model", "reference_type",
  "tmax_mean_c", "tmax_f", "lag", "lag_label", "rr", "rr_low", "rr_high",
  "log_rr", "log_rr_se", "k_sites", "tau2", "i2"
)]
write.csv(source_table, file.path(figure_dir, "figure_case_crossover_3d_surface_source.csv"), row.names = FALSE)

draw_surface <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(
    mar = c(2.8, 2.9, 1.8, 1.1),
    mgp = c(1.8, 0.55, 0),
    cex.axis = 0.78,
    cex.lab = 0.88,
    cex.main = 0.95,
    font.lab = 2
  )

  z_facet <- (
    z[-nrow(z), -ncol(z), drop = FALSE] +
      z[-1L, -ncol(z), drop = FALSE] +
      z[-nrow(z), -1L, drop = FALSE] +
      z[-1L, -1L, drop = FALSE]
  ) / 4
  z_range <- range(z_facet, finite = TRUE)
  palette <- grDevices::colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#08519c"))
  colors <- palette(120)
  color_index <- pmax(1L, pmin(120L, round(1 + 119 * (z_facet - z_range[1]) / diff(z_range))))

  graphics::persp(
    x = temps_f,
    y = lags,
    z = z,
    theta = 38,
    phi = 24,
    expand = 0.72,
    col = colors[color_index],
    border = "grey62",
    ltheta = 120,
    shade = 0.24,
    ticktype = "detailed",
    xlab = "Temperature (°F)",
    ylab = "Lag day",
    zlab = "Relative risk",
    zlim = range(z, finite = TRUE),
    main = "Time-stratified case-crossover DLNM"
  )
}

png_path <- file.path(figure_dir, "figure_case_crossover_3d_surface.png")
pdf_path <- file.path(figure_dir, "figure_case_crossover_3d_surface.pdf")

grDevices::png(png_path, width = 7.4, height = 5.6, units = "in", res = 600, type = "cairo")
draw_surface()
invisible(grDevices::dev.off())

grDevices::pdf(pdf_path, width = 7.4, height = 5.6, useDingbats = FALSE)
draw_surface()
invisible(grDevices::dev.off())

message("Wrote case-crossover 3D surface figure to ", png_path, " and ", pdf_path)

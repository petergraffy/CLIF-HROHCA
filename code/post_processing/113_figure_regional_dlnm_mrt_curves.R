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
input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(ragg)
  library(readr)
  library(scales)
  library(svglite)
})

pool_der_simonian_laird <- function(est, se) {
  keep <- is.finite(est) & is.finite(se) & se >= 0
  est <- est[keep]
  se <- se[keep]
  k <- length(est)
  if (k == 0) return(NULL)
  vi <- pmax(se^2, 1e-6)
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
    cumulative_rr = exp(pooled),
    cumulative_rr_low = exp(pooled - 1.96 * pooled_se),
    cumulative_rr_high = exp(pooled + 1.96 * pooled_se),
    p_value = 2 * stats::pnorm(abs(pooled / pooled_se), lower.tail = FALSE),
    tau2 = tau2,
    q = q,
    i2 = ifelse(k > 1 && q > (k - 1) && q > 0, 100 * (q - (k - 1)) / q, 0),
    stringsAsFactors = FALSE
  )
}

fmt_rr <- function(rr, low, high) {
  sprintf("%.2f (%.2f-%.2f)", rr, low, high)
}

regional_map <- tibble::tribble(
  ~site_name, ~region,
  "OHSU", "West Coast",
  "UCSF", "West Coast",
  "Michigan", "Midwest",
  "NU", "Midwest",
  "RUMC", "Midwest",
  "UCMC", "Midwest",
  "UMN", "Midwest",
  "JHU", "Northeast",
  "Penn", "Northeast",
  "Emory", "South"
) |>
  mutate(
    region = factor(.data$region, levels = c("West Coast", "Midwest", "Northeast", "South"))
  )

curve_site <- readr::read_csv(file.path(input_dir, "all_sites_dlnm_curves.csv"), show_col_types = FALSE) |>
  filter(
    .data$stratum == "Overall",
    .data$model == "sensitivity_mrt_reference",
    .data$reference_type == "mrt"
  ) |>
  inner_join(regional_map, by = "site_name") |>
  mutate(
    tmax_mean_c = as.numeric(.data$tmax_mean_c),
    cumulative_rr = as.numeric(.data$cumulative_rr),
    cumulative_rr_low = as.numeric(.data$cumulative_rr_low),
    cumulative_rr_high = as.numeric(.data$cumulative_rr_high),
    log_rr = as.numeric(.data$log_rr),
    log_rr_se = as.numeric(.data$log_rr_se)
  )

missing_sites <- setdiff(regional_map$site_name, unique(curve_site$site_name))
if (length(missing_sites) > 0) {
  warning("Missing MRT curve exports for: ", paste(missing_sites, collapse = ", "))
}

region_groups <- curve_site |>
  distinct(.data$region, .data$site_name) |>
  arrange(.data$region, .data$site_name) |>
  group_by(.data$region) |>
  summarise(
    n_sites_total = n(),
    sites = paste(.data$site_name, collapse = ", "),
    .groups = "drop"
  )

curve_groups <- curve_site |>
  distinct(.data$region, .data$tmax_mean_c) |>
  arrange(.data$region, .data$tmax_mean_c)

curve_rows <- vector("list", nrow(curve_groups))
for (i in seq_len(nrow(curve_groups))) {
  dat <- curve_site |>
    filter(
      .data$region == curve_groups$region[i],
      .data$tmax_mean_c == curve_groups$tmax_mean_c[i]
    )
  pooled <- pool_der_simonian_laird(dat$log_rr, dat$log_rr_se)
  if (is.null(pooled)) next
  curve_rows[[i]] <- bind_cols(
    curve_groups[i, , drop = FALSE],
    tibble::tibble(sites_contributing = paste(sort(unique(dat$site_name)), collapse = ", ")),
    pooled
  )
}

regional_curves <- bind_rows(curve_rows) |>
  left_join(region_groups, by = "region") |>
  mutate(
    region_label = sprintf("%s\n%s", as.character(.data$region), .data$sites),
    region_label = factor(
      .data$region_label,
      levels = region_groups |>
        mutate(region_label = sprintf("%s\n%s", as.character(.data$region), .data$sites)) |>
        pull(.data$region_label)
    )
  ) |>
  arrange(.data$region, .data$tmax_mean_c)

reference_points <- curve_site |>
  distinct(.data$region, .data$site_name, .data$reference_temp_c) |>
  mutate(reference_temp_c = as.numeric(.data$reference_temp_c)) |>
  group_by(.data$region) |>
  summarise(
    median_mrt_c = median(.data$reference_temp_c, na.rm = TRUE),
    min_mrt_c = min(.data$reference_temp_c, na.rm = TRUE),
    max_mrt_c = max(.data$reference_temp_c, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    regional_curves |>
      distinct(.data$region, .data$region_label),
    by = "region"
  )

summary_temps <- c(30, 32, 34, 36)
regional_summary <- bind_rows(lapply(summary_temps, function(temp_c) {
  regional_curves |>
    group_by(.data$region) |>
    slice_min(abs(.data$tmax_mean_c - temp_c), n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(target_temp_c = temp_c)
})) |>
  left_join(region_groups, by = "region") |>
  transmute(
    region = as.character(.data$region),
    sites = .data$sites.x,
    target_temp_c = .data$target_temp_c,
    plotted_temp_c = .data$tmax_mean_c,
    k_sites = .data$k_sites,
    cumulative_rr = .data$cumulative_rr,
    cumulative_rr_low = .data$cumulative_rr_low,
    cumulative_rr_high = .data$cumulative_rr_high,
    rr_95_ci = fmt_rr(.data$cumulative_rr, .data$cumulative_rr_low, .data$cumulative_rr_high),
    p_value = .data$p_value,
    i2 = .data$i2
  ) |>
  arrange(factor(.data$region, levels = levels(regional_map$region)), .data$target_temp_c)

region_colors <- c(
  "West Coast" = "#007C89",
  "Midwest" = "#7A0019",
  "Northeast" = "#1769AA",
  "South" = "#8B1E3F"
)

regional_plot <- ggplot(regional_curves, aes(x = .data$tmax_mean_c, y = .data$cumulative_rr)) +
  geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dashed", color = "#4B5563") +
  geom_ribbon(
    aes(ymin = .data$cumulative_rr_low, ymax = .data$cumulative_rr_high, fill = .data$region),
    alpha = 0.20,
    linewidth = 0
  ) +
  geom_line(aes(color = .data$region), linewidth = 1.15, lineend = "round") +
  geom_vline(
    data = reference_points,
    aes(xintercept = .data$median_mrt_c, color = .data$region),
    linewidth = 0.45,
    linetype = "dotted",
    show.legend = FALSE
  ) +
  facet_wrap(~region_label, ncol = 2) +
  scale_color_manual(values = region_colors, guide = "none") +
  scale_fill_manual(values = region_colors, guide = "none") +
  scale_x_continuous(
    breaks = seq(10, 40, by = 5),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.015, 0.035))
  ) +
  scale_y_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 5, 10, 20, 50, 100, 200),
    labels = label_number(accuracy = 0.01)
  ) +
  coord_cartesian(ylim = c(0.5, 250)) +
  labs(
    title = "Regional MRT-referenced temperature-OHCA DLNM curves",
    x = "County-level daily maximum temperature, Tmax (°C)",
    y = "Cumulative relative risk of OHCA (log scale)"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#111827", margin = margin(b = 8)),
    axis.title = element_text(face = "bold", size = 13, color = "#111827"),
    axis.text = element_text(size = 11.5, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    axis.ticks.length = unit(0.13, "cm"),
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12.5, color = "#111827", lineheight = 1.05),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.0, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 14, 12, 12)
  )

readr::write_csv(curve_site, file.path(figure_dir, "figure_regional_dlnm_mrt_site_curves_source.csv"))
readr::write_csv(regional_curves, file.path(figure_dir, "figure_regional_dlnm_mrt_curves_source.csv"))
readr::write_csv(regional_summary, file.path(table_dir, "supplement_table_regional_dlnm_mrt_summary.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Regional RR")
openxlsx::writeData(wb, "Regional RR", regional_summary)
openxlsx::setColWidths(wb, "Regional RR", cols = seq_along(regional_summary), widths = "auto")
openxlsx::saveWorkbook(wb, file.path(table_dir, "supplement_table_regional_dlnm_mrt_summary.xlsx"), overwrite = TRUE)

ggsave(
  file.path(figure_dir, "figure_regional_dlnm_mrt_curves.pdf"),
  regional_plot,
  width = 9.2,
  height = 7.4,
  device = cairo_pdf,
  bg = "white"
)

svglite::svglite(file.path(figure_dir, "figure_regional_dlnm_mrt_curves.svg"), width = 9.2, height = 7.4, bg = "white")
print(regional_plot)
dev.off()

ragg::agg_png(
  file.path(figure_dir, "figure_regional_dlnm_mrt_curves.png"),
  width = 9.2,
  height = 7.4,
  units = "in",
  res = 600,
  background = "white"
)
print(regional_plot)
dev.off()

ragg::agg_tiff(
  file.path(figure_dir, "figure_regional_dlnm_mrt_curves.tiff"),
  width = 9.2,
  height = 7.4,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(regional_plot)
dev.off()

message("Wrote regional MRT DLNM curves to ", figure_dir)

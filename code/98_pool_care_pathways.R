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
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
})

find_site_exports <- function(patterns) {
  files <- unlist(lapply(patterns, function(pattern) {
    list.files(box_root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  }))
  files[file.exists(files)]
}

site_from_path <- function(path) {
  parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
  idx <- match("CLIF-Heat-Related-OHCA", parts)
  if (!is.na(idx) && length(parts) >= idx + 1L) return(parts[[idx + 1L]])
  sub("_.*$", "", basename(path))
}

read_export <- function(path) {
  out <- readr::read_csv(path, show_col_types = FALSE)
  if (!"site_name" %in% names(out)) out$site_name <- site_from_path(path)
  out$source_file <- path
  out
}

pool_count_table <- function(dat, group_col) {
  dat |>
    mutate(n = suppressWarnings(as.numeric(.data$n))) |>
    filter(!is.na(.data[[group_col]]), is.finite(.data$n)) |>
    group_by(.data[[group_col]]) |>
    summarise(n = sum(.data$n, na.rm = TRUE), .groups = "drop") |>
    rename(label = 1) |>
    arrange(desc(.data$n)) |>
    mutate(pct = 100 * .data$n / sum(.data$n, na.rm = TRUE))
}

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

fmt_n_pct <- function(n, pct) sprintf("%s (%.1f%%)", format(n, big.mark = ","), pct)
fmt_hours <- function(x) ifelse(is.finite(x), sprintf("%.1f", x), "")
fmt_pathway <- function(x) {
  x |>
    str_replace_all("\\bicu\\b", "ICU") |>
    str_replace_all("\\bed\\b", "ED") |>
    str_replace_all("\\bward\\b", "Ward") |>
    str_replace_all("\\bstepdown\\b", "Stepdown") |>
    str_replace_all("\\bprocedural\\b", "Procedural") |>
    str_replace_all("\\bradiology\\b", "Radiology") |>
    str_replace_all("\\bother\\b", "Other") |>
    str_replace_all(" -> ", " to ")
}

pathway_files <- find_site_exports(c("_care_pathway_summary[.]csv$", "ohca_pre_icu_care_pathway_summary[.]csv$"))
first_location_files <- find_site_exports(c("_first_location_summary[.]csv$", "ohca_first_location_summary[.]csv$"))
timing_files <- find_site_exports(c("_icu_timing_summary[.]csv$", "ohca_admission_to_icu_timing_summary[.]csv$"))
timing_bin_files <- find_site_exports(c("_icu_timing_bins[.]csv$", "ohca_admission_to_icu_timing_bins[.]csv$"))

care_pathway_site <- bind_rows(lapply(pathway_files, read_export))
first_location_site <- bind_rows(lapply(first_location_files, read_export))
timing_site <- bind_rows(lapply(timing_files, read_export))
timing_bins_site <- bind_rows(lapply(timing_bin_files, read_export))

care_pathway_pooled <- pool_count_table(care_pathway_site, "care_pathway") |>
  rename(care_pathway = "label") |>
  mutate(display = fmt_n_pct(.data$n, .data$pct))

first_location_pooled <- pool_count_table(first_location_site, "first_location_category") |>
  rename(first_location_category = "label") |>
  mutate(display = fmt_n_pct(.data$n, .data$pct))

timing_bins_pooled <- pool_count_table(timing_bins_site, "icu_entry_timing") |>
  rename(icu_entry_timing = "label") |>
  mutate(
    icu_entry_timing = factor(
      .data$icu_entry_timing,
      levels = c("<=0 hours", "0-6 hours", "6-12 hours", "12-24 hours", "24-48 hours", "48-72 hours", ">72 hours")
    )
  ) |>
  arrange(.data$icu_entry_timing) |>
  mutate(display = fmt_n_pct(.data$n, .data$pct))

timing_site <- timing_site |>
  mutate(across(
    c("n", "median_hours", "iqr_low_hours", "iqr_high_hours", "p90_hours", "p95_hours",
      "within_6h_pct", "within_12h_pct", "within_24h_pct", "within_48h_pct"),
    ~ suppressWarnings(as.numeric(.x))
  ))

timing_pooled <- tibble(
  contributing_sites = n_distinct(timing_site$site_name),
  n = sum(timing_site$n, na.rm = TRUE),
  median_hours_approx = weighted_mean(timing_site$median_hours, timing_site$n),
  iqr_low_hours_approx = weighted_mean(timing_site$iqr_low_hours, timing_site$n),
  iqr_high_hours_approx = weighted_mean(timing_site$iqr_high_hours, timing_site$n),
  p90_hours_approx = weighted_mean(timing_site$p90_hours, timing_site$n),
  p95_hours_approx = weighted_mean(timing_site$p95_hours, timing_site$n),
  within_6h_pct = weighted_mean(timing_site$within_6h_pct, timing_site$n),
  within_12h_pct = weighted_mean(timing_site$within_12h_pct, timing_site$n),
  within_24h_pct = weighted_mean(timing_site$within_24h_pct, timing_site$n),
  within_48h_pct = weighted_mean(timing_site$within_48h_pct, timing_site$n),
  note = "Quantiles are n-weighted averages of site-level quantiles; timing-bin percentages are pooled from site counts."
)

site_availability <- tibble(
  site_name = sort(unique(c(
    care_pathway_site$site_name,
    first_location_site$site_name,
    timing_site$site_name,
    timing_bins_site$site_name
  )))
) |>
  mutate(
    care_pathway = .data$site_name %in% care_pathway_site$site_name,
    first_location = .data$site_name %in% first_location_site$site_name,
    timing_summary = .data$site_name %in% timing_site$site_name,
    timing_bins = .data$site_name %in% timing_bins_site$site_name
  )

readr::write_csv(care_pathway_site, file.path(output_dir, "all_sites_care_pathway_summary.csv"))
readr::write_csv(first_location_site, file.path(output_dir, "all_sites_first_location_summary.csv"))
readr::write_csv(timing_site, file.path(output_dir, "all_sites_icu_timing_summary.csv"))
readr::write_csv(timing_bins_site, file.path(output_dir, "all_sites_icu_timing_bins.csv"))
readr::write_csv(care_pathway_pooled, file.path(output_dir, "pooled_care_pathway_summary.csv"))
readr::write_csv(first_location_pooled, file.path(output_dir, "pooled_first_location_summary.csv"))
readr::write_csv(timing_pooled, file.path(output_dir, "pooled_icu_timing_summary.csv"))
readr::write_csv(timing_bins_pooled, file.path(output_dir, "pooled_icu_timing_bins.csv"))
readr::write_csv(site_availability, file.path(output_dir, "care_pathway_export_availability.csv"))

top_pathways <- care_pathway_pooled |>
  arrange(desc(.data$n)) |>
  mutate(care_pathway_plot = if_else(row_number() <= 5L, .data$care_pathway, "Other")) |>
  group_by(.data$care_pathway_plot) |>
  summarise(n = sum(.data$n), .groups = "drop") |>
  mutate(
    pct = 100 * .data$n / sum(.data$n),
    pathway_label = if_else(.data$care_pathway_plot == "Other", "Other", fmt_pathway(.data$care_pathway_plot)),
    pathway_label = factor(.data$pathway_label, levels = c(fmt_pathway(care_pathway_pooled$care_pathway[1:5]), "Other")),
    display = sprintf("%s: %s (%.1f%%)", .data$pathway_label, scales::comma(.data$n), .data$pct)
  ) |>
  arrange(.data$pathway_label) |>
  mutate(
    ymax = cumsum(.data$n),
    ymin = dplyr::lag(.data$ymax, default = 0),
    ymid = (.data$ymin + .data$ymax) / 2
  )

timing_center <- sprintf(
  "Time to ICU\nMedian %.1f h\nIQR %.1f-%.1f h",
  timing_pooled$median_hours_approx[1],
  timing_pooled$iqr_low_hours_approx[1],
  timing_pooled$iqr_high_hours_approx[1]
)

pathway_palette <- c(
  "ED to ICU" = "#2F6F73",
  "ICU" = "#7A5C99",
  "ED to Procedural to ICU" = "#C7793D",
  "Procedural to ICU" = "#4F79A7",
  "ED to Ward to ICU" = "#8FAE5D",
  "Other" = "#D6D9DE"
)

pathway_plot <- ggplot(top_pathways, aes(ymax = .data$ymax, ymin = .data$ymin, xmax = 4, xmin = 2.35, fill = .data$pathway_label)) +
  geom_rect(color = "white", linewidth = 1.1) +
  annotate(
    "label",
    x = 0,
    y = 0,
    label = timing_center,
    size = 5.0,
    lineheight = 0.92,
    linewidth = 0,
    fontface = "bold",
    fill = "white",
    color = "#17212B"
  ) +
  coord_polar(theta = "y", clip = "off") +
  scale_fill_manual(
    values = pathway_palette,
    breaks = levels(top_pathways$pathway_label),
    labels = setNames(as.character(top_pathways$display), as.character(top_pathways$pathway_label)),
    drop = FALSE
  ) +
  xlim(c(0, 4.25)) +
  labs(fill = NULL) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11, color = "#17212B", margin = margin(b = 6)),
    legend.key.size = unit(0.34, "in"),
    legend.spacing.y = unit(0.12, "in"),
    plot.margin = margin(18, 18, 18, 18)
  )

readr::write_csv(top_pathways, file.path(figure_dir, "figure_pooled_care_pathways_source.csv"))

ggsave(file.path(figure_dir, "figure_pooled_care_pathways.png"), pathway_plot, width = 8.4, height = 6.2, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "figure_pooled_care_pathways.pdf"), pathway_plot, width = 8.4, height = 6.2, device = cairo_pdf, bg = "white")
ggsave(file.path(figure_dir, "figure_pooled_care_pathways.svg"), pathway_plot, width = 8.4, height = 6.2, bg = "white")

message("Wrote pooled care-pathway outputs to: ", output_dir)

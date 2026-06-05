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
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

HEAT_DEFINITION <- Sys.getenv("HEAT_DEFINITION", unset = "heat95")
OUTPUT_SUFFIX <- if (identical(HEAT_DEFINITION, "heat95")) "" else paste0("_", HEAT_DEFINITION)

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(ragg)
  library(readr)
  library(scales)
  library(svglite)
  library(tidyr)
})

time_bins <- tibble::tribble(
  ~time_bin, ~start_hour, ~end_hour,
  "0-6 h", 0, 6,
  "7-12 h", 7, 12,
  "13-24 h", 13, 24,
  "25-48 h", 25, 48,
  "49-72 h", 49, 72
)

assign_time_bin <- function(hour) {
  dplyr::case_when(
    hour >= 0 & hour <= 6 ~ "0-6 h",
    hour >= 7 & hour <= 12 ~ "7-12 h",
    hour >= 13 & hour <= 24 ~ "13-24 h",
    hour >= 25 & hour <= 48 ~ "25-48 h",
    hour >= 49 & hour <= 72 ~ "49-72 h",
    TRUE ~ NA_character_
  )
}

fmt_num <- function(x, digits = 1, suffix = "") {
  ifelse(is.na(x), "", paste0(sprintf(paste0("%.", digits, "f"), x), suffix))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_group <- function(value, n, digits = 1, suffix = "") {
  ifelse(
    is.na(value),
    "",
    paste0(sprintf(paste0("%.", digits, "f"), value), suffix, " (n=", format(round(n), big.mark = ","), ")")
  )
}

theme_manuscript <- function(base_size = 12) {
  theme_classic(base_family = "Helvetica", base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", color = "#111827", size = 15),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, color = "#111827", size = 13),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(color = "#111827", size = 10.5),
      plot.title = element_text(face = "bold", color = "#111827", size = 22, margin = margin(b = 10)),
      plot.margin = margin(8, 12, 8, 10)
    )
}

save_figure <- function(plot, basename, width, height) {
  ggsave(file.path(figure_dir, paste0(basename, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  svglite::svglite(file.path(figure_dir, paste0(basename, ".svg")), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(figure_dir, paste0(basename, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(figure_dir, paste0(basename, ".tiff")), width = width, height = height, units = "in", res = 600, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

summarise_measure_bins <- function(dat, label_col, domain, value_digits = 1) {
  dat |>
    filter(
      .data$heat_definition == .env$HEAT_DEFINITION,
      .data$icu_hour >= 0,
      .data$icu_hour <= 72
    ) |>
    mutate(
      time_bin = assign_time_bin(.data$icu_hour),
      time_bin = factor(.data$time_bin, levels = time_bins$time_bin),
      measure = .data[[label_col]]
    ) |>
    filter(!is.na(.data$time_bin)) |>
    group_by(.data$measure, .data$time_bin) |>
    summarise(
      domain = domain,
      hrohca_value = mean(.data$heat_weighted_median, na.rm = TRUE),
      non_hrohca_value = mean(.data$non_heat_weighted_median, na.rm = TRUE),
      difference = mean(.data$approx_difference, na.rm = TRUE),
      min_p_value = suppressWarnings(min(.data$approx_p_value, na.rm = TRUE)),
      significant_hours = sum(.data$approx_p_value < 0.05, na.rm = TRUE),
      n_hours = n(),
      hrohca_median_n = stats::median(.data$heat_n_patients, na.rm = TRUE),
      non_hrohca_median_n = stats::median(.data$non_heat_n_patients, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      min_p_value = if_else(is.infinite(.data$min_p_value), NA_real_, .data$min_p_value),
      `HROHCA` = fmt_group(.data$hrohca_value, .data$hrohca_median_n, digits = value_digits),
      `Non-HROHCA` = fmt_group(.data$non_hrohca_value, .data$non_hrohca_median_n, digits = value_digits),
      `Difference` = fmt_num(.data$difference, digits = value_digits),
      `Minimum hourly p value` = fmt_p(.data$min_p_value),
      `Significant hours` = paste0(.data$significant_hours, "/", .data$n_hours)
    ) |>
    select(
      "domain",
      "measure",
      "time_bin",
      "HROHCA",
      "Non-HROHCA",
      "Difference",
      "Minimum hourly p value",
      "Significant hours",
      "hrohca_value",
      "non_hrohca_value",
      "difference",
      "min_p_value",
      "significant_hours",
      "n_hours"
    )
}

labs_raw <- read.csv(
  file.path(input_dir, "pooled_heat_related_hourly_lab_median_differences.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lab_labels <- c(
  "lactate" = "Lactate, mmol/L",
  "bicarbonate" = "Bicarbonate, mEq/L",
  "ph_arterial" = "Arterial pH",
  "creatinine" = "Creatinine, mg/dL",
  "bun" = "BUN, mg/dL",
  "potassium" = "Potassium, mEq/L",
  "magnesium" = "Magnesium, mg/dL",
  "phosphate" = "Phosphate, mg/dL",
  "wbc" = "WBC, K/uL"
)

labs_binned <- labs_raw |>
  filter(.data$variable %in% names(lab_labels)) |>
  mutate(lab = lab_labels[.data$variable]) |>
  summarise_measure_bins("lab", "Laboratory", value_digits = 1)

vitals_raw <- read.csv(
  file.path(input_dir, "pooled_heat_related_hourly_vital_median_differences.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

vital_labels <- c(
  "heart_rate" = "Heart rate, bpm",
  "map" = "Mean arterial pressure, mm Hg",
  "sbp" = "Systolic blood pressure, mm Hg",
  "respiratory_rate" = "Respiratory rate, breaths/min",
  "spo2" = "SpO2, %"
)

vitals_binned <- vitals_raw |>
  filter(.data$variable %in% names(vital_labels)) |>
  mutate(vital = vital_labels[.data$variable]) |>
  summarise_measure_bins("vital", "Vital signs", value_digits = 1)

support_raw <- read.csv(
  file.path(input_dir, "pooled_heat_related_hourly_support_prevalence_differences.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

support_labels <- c(
  "Invasive mechanical ventilation" = "IMV prevalence, %",
  "Vasopressor infusion" = "Vasopressor prevalence, %",
  "CRRT" = "CRRT prevalence, %"
)

support_binned <- support_raw |>
  filter(
    .data$heat_definition == .env$HEAT_DEFINITION,
    .data$variable %in% names(support_labels),
    .data$icu_hour >= 0,
    .data$icu_hour <= 72
  ) |>
  mutate(
    time_bin = assign_time_bin(.data$icu_hour),
    time_bin = factor(.data$time_bin, levels = time_bins$time_bin),
    measure = support_labels[.data$variable]
  ) |>
  filter(!is.na(.data$time_bin)) |>
  group_by(.data$measure, .data$time_bin) |>
  summarise(
    domain = "Organ support prevalence",
    hrohca_value = mean(.data$heat_prevalence_pct, na.rm = TRUE),
    non_hrohca_value = mean(.data$non_heat_prevalence_pct, na.rm = TRUE),
    difference = mean(.data$absolute_difference_pct, na.rm = TRUE),
    min_p_value = suppressWarnings(min(.data$p_value, na.rm = TRUE)),
    significant_hours = sum(.data$p_value < 0.05, na.rm = TRUE),
    n_hours = n(),
    hrohca_median_n = stats::median(.data$heat_n_at_risk, na.rm = TRUE),
    non_hrohca_median_n = stats::median(.data$non_heat_n_at_risk, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    min_p_value = if_else(is.infinite(.data$min_p_value), NA_real_, .data$min_p_value),
    `HROHCA` = fmt_group(.data$hrohca_value, .data$hrohca_median_n, digits = 1, suffix = "%"),
    `Non-HROHCA` = fmt_group(.data$non_hrohca_value, .data$non_hrohca_median_n, digits = 1, suffix = "%"),
    `Difference` = fmt_num(.data$difference, digits = 1, suffix = " pp"),
    `Minimum hourly p value` = fmt_p(.data$min_p_value),
    `Significant hours` = paste0(.data$significant_hours, "/", .data$n_hours)
  ) |>
  select(
    "domain",
    "measure",
    "time_bin",
    "HROHCA",
    "Non-HROHCA",
    "Difference",
    "Minimum hourly p value",
    "Significant hours",
    "hrohca_value",
    "non_hrohca_value",
    "difference",
    "min_p_value",
    "significant_hours",
    "n_hours"
  )

cumulative_path <- file.path(input_dir, paste0("pooled_heat_related_hourly_cumulative_incidence", OUTPUT_SUFFIX, ".csv"))
if (!file.exists(cumulative_path)) {
  cumulative_path <- file.path(input_dir, "pooled_heat_related_hourly_cumulative_incidence.csv")
}

cumulative_raw <- read.csv(cumulative_path, stringsAsFactors = FALSE, check.names = FALSE)
event_labels <- c(
  "Hospital death" = "Hospital death cumulative incidence, %",
  "Death or hospice" = "Death/hospice cumulative incidence, %",
  "Discharged alive without hospice" = "Discharged alive without hospice cumulative incidence, %"
)

cumulative_binned <- cumulative_raw |>
  filter(
    .data$heat_definition == .env$HEAT_DEFINITION,
    .data$event %in% names(event_labels),
    .data$icu_hour >= 0,
    .data$icu_hour <= 72
  ) |>
  mutate(
    time_bin = assign_time_bin(.data$icu_hour),
    time_bin = factor(.data$time_bin, levels = time_bins$time_bin),
    measure = event_labels[.data$event],
    group = if_else(.data$heat_related_ohca == "Heat-related OHCA", "HROHCA", "Non-HROHCA")
  ) |>
  filter(!is.na(.data$time_bin)) |>
  group_by(.data$measure, .data$time_bin, .data$group) |>
  filter(.data$icu_hour == max(.data$icu_hour, na.rm = TRUE)) |>
  ungroup() |>
  select("measure", "time_bin", "group", "n_group", "n_events", "cumulative_pct") |>
  pivot_wider(
    names_from = "group",
    values_from = c("n_group", "n_events", "cumulative_pct"),
    names_sep = "__"
  ) |>
  mutate(
    domain = "Cumulative events",
    hrohca_value = .data$`cumulative_pct__HROHCA`,
    non_hrohca_value = .data$`cumulative_pct__Non-HROHCA`,
    difference = .data$hrohca_value - .data$non_hrohca_value,
    p_value = mapply(
      function(x_heat, x_non, n_heat, n_non) {
        if (!all(is.finite(c(x_heat, x_non, n_heat, n_non)))) return(NA_real_)
        suppressWarnings(stats::prop.test(x = c(x_heat, x_non), n = c(n_heat, n_non))$p.value)
      },
      .data$`n_events__HROHCA`,
      .data$`n_events__Non-HROHCA`,
      .data$`n_group__HROHCA`,
      .data$`n_group__Non-HROHCA`
    ),
    `HROHCA` = fmt_group(.data$hrohca_value, .data$`n_group__HROHCA`, digits = 1, suffix = "%"),
    `Non-HROHCA` = fmt_group(.data$non_hrohca_value, .data$`n_group__Non-HROHCA`, digits = 1, suffix = "%"),
    `Difference` = fmt_num(.data$difference, digits = 1, suffix = " pp"),
    `Minimum hourly p value` = fmt_p(.data$p_value),
    `Significant hours` = if_else(.data$p_value < 0.05, "End of bin", "")
  ) |>
  select(
    "domain",
    "measure",
    "time_bin",
    "HROHCA",
    "Non-HROHCA",
    "Difference",
    "Minimum hourly p value",
    "Significant hours",
    "hrohca_value",
    "non_hrohca_value",
    "difference",
    "p_value"
  )

all_binned_display <- bind_rows(
  labs_binned,
  vitals_binned,
  support_binned,
  cumulative_binned
) |>
  mutate(
    time_bin = factor(.data$time_bin, levels = time_bins$time_bin),
    domain = factor(.data$domain, levels = c("Laboratory", "Vital signs", "Organ support prevalence", "Cumulative events"))
  ) |>
  arrange(.data$domain, .data$measure, .data$time_bin)

display_table <- all_binned_display |>
  transmute(
    Domain = as.character(.data$domain),
    Measure = .data$measure,
    `Time bin` = as.character(.data$time_bin),
    HROHCA = .data$HROHCA,
    `Non-HROHCA` = .data$`Non-HROHCA`,
    `HROHCA - non-HROHCA` = .data$Difference,
    `P value` = .data$`Minimum hourly p value`,
    `Significant hours` = .data$`Significant hours`
  )

source_table <- all_binned_display |>
  mutate(time_bin = as.character(.data$time_bin), domain = as.character(.data$domain))

base_name <- paste0("supplement_table_binned_time_summaries", OUTPUT_SUFFIX)
readr::write_csv(display_table, file.path(table_dir, paste0(base_name, ".csv")))
readr::write_csv(source_table, file.path(table_dir, paste0(base_name, "_source.csv")))

wb <- openxlsx::createWorkbook()
for (domain_name in unique(display_table$Domain)) {
  sheet_name <- substr(domain_name, 1, 31)
  openxlsx::addWorksheet(wb, sheet_name)
  sheet_dat <- display_table |>
    filter(.data$Domain == domain_name) |>
    select(-"Domain")
  openxlsx::writeData(wb, sheet_name, sheet_dat)
  openxlsx::setColWidths(wb, sheet_name, cols = 1:ncol(sheet_dat), widths = "auto")
}
openxlsx::addWorksheet(wb, "All")
openxlsx::writeData(wb, "All", display_table)
openxlsx::setColWidths(wb, "All", cols = 1:ncol(display_table), widths = "auto")
openxlsx::saveWorkbook(wb, file.path(table_dir, paste0(base_name, ".xlsx")), overwrite = TRUE)

docx_tables <- lapply(unique(display_table$Domain), function(domain_name) {
  display_table |>
    filter(.data$Domain == domain_name) |>
    select(-"Domain") |>
    flextable() |>
    autofit() |>
    align(align = "left", part = "all") |>
    theme_booktabs() |>
    bold(part = "header")
})
names(docx_tables) <- paste0("Supplemental Table. ", unique(display_table$Domain), " binned summaries")
docx_tables[[length(docx_tables)]] <- docx_tables[[length(docx_tables)]] |>
  add_footer_lines(
    values = paste(
      "Primary heat definition:",
      HEAT_DEFINITION,
      ". Laboratory and vital-sign values are mean pooled medians within each bin.",
      "Support values are mean hourly prevalence within each bin.",
      "Cumulative event values are end-of-bin cumulative incidence.",
      "CRRT is included as hourly support prevalence."
    )
  )
do.call(flextable::save_as_docx, c(docx_tables, list(path = file.path(table_dir, paste0(base_name, ".docx")))))

all_panel_plot_dat <- source_table |>
  mutate(
    time_bin = factor(.data$time_bin, levels = time_bins$time_bin),
    measure_short = case_when(
      .data$domain == "Organ support prevalence" ~ paste0(gsub(" prevalence, %$", "", .data$measure), ", pp"),
      .data$domain == "Cumulative events" ~ paste0(gsub(" cumulative incidence, %$", "", .data$measure), ", pp"),
      TRUE ~ .data$measure
    ),
    panel = paste(.data$domain, .data$measure_short, sep = ": "),
    panel = factor(
      .data$panel,
      levels = c(
        paste0("Laboratory: ", lab_labels),
        paste0("Vital signs: ", vital_labels),
        "Organ support prevalence: IMV, pp",
        "Organ support prevalence: Vasopressor, pp",
        "Organ support prevalence: CRRT, pp",
        "Cumulative events: Hospital death, pp",
        "Cumulative events: Death/hospice, pp",
        "Cumulative events: Discharged alive without hospice, pp"
      )
    ),
    p_for_star = dplyr::coalesce(.data$min_p_value, .data$p_value)
  ) |>
  filter(!is.na(.data$panel)) |>
  group_by(.data$panel) |>
  mutate(
    y_range = diff(range(c(.data$difference, 0), na.rm = TRUE)),
    y_offset = if_else(is.finite(.data$y_range) & .data$y_range > 0, 0.08 * .data$y_range, 0.25),
    star_y = if_else(.data$difference >= 0, .data$difference + .data$y_offset, .data$difference - .data$y_offset),
    star_label = case_when(
      !is.na(.data$p_for_star) & .data$p_for_star < 0.001 ~ "***",
      !is.na(.data$p_for_star) & .data$p_for_star < 0.01 ~ "**",
      !is.na(.data$p_for_star) & .data$p_for_star < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  ungroup()

bin_plot <- ggplot(all_panel_plot_dat, aes(.data$time_bin, .data$difference, fill = .data$difference)) +
  geom_col(aes(fill = .data$difference > 0), width = 0.72, color = "white", linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "#111827", linewidth = 0.45) +
  geom_text(
    data = filter(all_panel_plot_dat, .data$star_label != ""),
    aes(x = .data$time_bin, y = .data$star_y, label = .data$star_label),
    inherit.aes = FALSE,
    color = "#111827",
    family = "Helvetica",
    fontface = "bold",
    size = 4.6
  ) +
  facet_wrap(~panel, ncol = 4, scales = "free_y", labeller = label_wrap_gen(width = 38)) +
  scale_fill_manual(values = c("TRUE" = "#B23A2E", "FALSE" = "#2F5D8C"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.22))) +
  labs(
    title = "Binned clinical trajectory differences after OHCA",
    x = "ICU time bin",
    y = "HROHCA - non-HROHCA"
  ) +
  theme_manuscript(base_size = 13.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 12.5, margin = margin(b = 5)),
    panel.spacing.x = grid::unit(16, "pt"),
    panel.spacing.y = grid::unit(18, "pt")
  )

save_figure(bin_plot, paste0("figure_binned_clinical_time_summaries", OUTPUT_SUFFIX), width = 18.5, height = 18.0)

message("Wrote binned time summary tables and support/CRRT figure to ", table_dir)

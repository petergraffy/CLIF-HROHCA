#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(file_name) {
  path <- file.path(pooled_dir, file_name)
  if (!file.exists(path)) stop("Missing required pooled file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

save_plot <- function(plot, stem, width = 7.4, height = 5.1) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 600)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height)
  invisible(c(png_path, pdf_path))
}

theme_candidate <- function(base_size = 11.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "grey15"),
      strip.background = ggplot2::element_rect(fill = "grey93", color = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 0.5)
    )
}

fmt_effect <- function(est, low, high) sprintf("%.2f (%.2f, %.2f)", est, low, high)

parse_count <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(sub("^\\s*([0-9,]+).*", "\\1", gsub(",", "", x))))
  out[!grepl("[0-9]", x)] <- NA_real_
  out
}

quartile_levels <- c("Q1 coolest", "Q2", "Q3", "Q4 hottest")
quartile_colors <- c(
  "Q1 coolest" = "#335c67",
  "Q2" = "#8f7a2f",
  "Q3" = "#c06c3e",
  "Q4 hottest" = "#9b2f37"
)

fg_cif <- read_required("pooled_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv")
fg_cif <- fg_cif |>
  dplyr::mutate(
    event_type = dplyr::recode(
      .data$event_type,
      "Awake/extubated" = "Awake/extubated",
      "Death before awake/extubated" = "Death before awake/extubated"
    )
  )

fg_event_colors <- c(
  "Awake/extubated" = "#0f6b78",
  "Death before awake/extubated" = "#9b2f37"
)

fg_overall <- fg_cif |>
  dplyr::filter(.data$stratification == "Overall") |>
  dplyr::mutate(event_type = factor(.data$event_type, levels = names(fg_event_colors)))
readr::write_csv(fg_overall, file.path(figure_dir, "candidate_72h_competing_risk_cif_overall_source.csv"))

fig_fg_overall <- ggplot2::ggplot(fg_overall, ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$event_type, fill = .data$event_type)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$cif_low, ymax = .data$cif_high), alpha = 0.13, color = NA) +
  ggplot2::geom_step(linewidth = 1.2) +
  ggplot2::scale_color_manual(values = fg_event_colors, name = NULL) +
  ggplot2::scale_fill_manual(values = fg_event_colors, guide = "none") +
  ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.36), expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(
    title = "72-hour competing-risk cumulative incidence",
    x = "Hours since ICU admission",
    y = "Cumulative incidence"
  ) +
  theme_candidate() +
  ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
save_plot(fig_fg_overall, "candidate_72h_competing_risk_cif_overall")

fg_quartile <- fg_cif |>
  dplyr::filter(.data$stratification == "Admission Tmax quartile") |>
  dplyr::mutate(
    stratum = factor(.data$stratum, levels = quartile_levels),
    event_type = factor(.data$event_type, levels = names(fg_event_colors))
  )
readr::write_csv(fg_quartile, file.path(figure_dir, "candidate_72h_competing_risk_cif_by_tmax_quartile_source.csv"))

fig_fg_quartile <- ggplot2::ggplot(fg_quartile, ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$stratum)) +
  ggplot2::geom_step(linewidth = 1.05) +
  ggplot2::facet_wrap(ggplot2::vars(.data$event_type), nrow = 1) +
  ggplot2::scale_color_manual(values = quartile_colors, drop = FALSE, name = "Admission Tmax quartile") +
  ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.38), expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(
    title = "72-hour cumulative incidence by admission temperature quartile",
    x = "Hours since ICU admission",
    y = "Cumulative incidence"
  ) +
  theme_candidate() +
  ggplot2::theme(legend.position = "bottom")
save_plot(fig_fg_quartile, "candidate_72h_competing_risk_cif_by_tmax_quartile", width = 8.8, height = 4.8)

imv24_cif <- read_required("pooled_ohca_icu_imv24_time_to_event_cif_curves.csv")
imv24_event_order <- c("Successful extubation", "Death/hospice discharge", "Tracheostomy")
imv24_colors <- c(
  "Successful extubation" = "#0f6b78",
  "Death/hospice discharge" = "#9b2f37",
  "Tracheostomy" = "#8f7a2f"
)

imv24_overall <- imv24_cif |>
  dplyr::filter(.data$stratification == "Overall") |>
  dplyr::filter(.data$time_days <= 30) |>
  dplyr::mutate(event_type = factor(.data$event_type, levels = imv24_event_order))
readr::write_csv(imv24_overall, file.path(figure_dir, "candidate_imv24_time_to_event_cif_overall_source.csv"))

fig_imv24_overall <- ggplot2::ggplot(imv24_overall, ggplot2::aes(x = .data$time_days, y = .data$cif, color = .data$event_type, fill = .data$event_type)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$cif_low, ymax = .data$cif_high), alpha = 0.13, color = NA) +
  ggplot2::geom_step(linewidth = 1.2) +
  ggplot2::scale_color_manual(values = imv24_colors, name = NULL) +
  ggplot2::scale_fill_manual(values = imv24_colors, guide = "none") +
  ggplot2::scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.65), expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(
    title = "IMV-24 landmark competing-risk cumulative incidence",
    x = "Days since 24 hours of IMV",
    y = "Cumulative incidence"
  ) +
  theme_candidate() +
  ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
save_plot(fig_imv24_overall, "candidate_imv24_time_to_event_cif_overall")

imv24_quartile <- imv24_cif |>
  dplyr::filter(.data$stratification == "Admission Tmax quartile") |>
  dplyr::filter(.data$time_days <= 30) |>
  dplyr::mutate(
    stratum = factor(.data$stratum, levels = quartile_levels),
    event_type = factor(.data$event_type, levels = imv24_event_order)
  )
readr::write_csv(imv24_quartile, file.path(figure_dir, "candidate_imv24_time_to_event_cif_by_tmax_quartile_source.csv"))

fig_imv24_quartile <- ggplot2::ggplot(imv24_quartile, ggplot2::aes(x = .data$time_days, y = .data$cif, color = .data$stratum)) +
  ggplot2::geom_step(linewidth = 1.0) +
  ggplot2::facet_wrap(ggplot2::vars(.data$event_type), nrow = 1) +
  ggplot2::scale_color_manual(values = quartile_colors, drop = FALSE, name = "Admission Tmax quartile") +
  ggplot2::scale_x_continuous(breaks = seq(0, 30, by = 10), limits = c(0, 30), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.70), expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(
    title = "IMV-24 cumulative incidence by admission temperature quartile",
    x = "Days since 24 hours of IMV",
    y = "Cumulative incidence"
  ) +
  theme_candidate(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")
save_plot(fig_imv24_quartile, "candidate_imv24_time_to_event_cif_by_tmax_quartile", width = 10.2, height = 4.8)

imv24_coef <- read_required("pooled_ohca_icu_imv24_time_to_event_coefficients.csv") |>
  dplyr::filter(.data$coefficient == "tmax_per5c") |>
  dplyr::mutate(
    outcome = factor(.data$event_type, levels = rev(imv24_event_order)),
    label = fmt_effect(.data$subdistribution_hr, .data$subdistribution_hr_low, .data$subdistribution_hr_high),
    label_x = .data$subdistribution_hr_high * 1.035
  )
readr::write_csv(imv24_coef, file.path(figure_dir, "candidate_imv24_temperature_effect_forest_source.csv"))

fig_imv24_forest <- ggplot2::ggplot(imv24_coef, ggplot2::aes(x = .data$subdistribution_hr, y = .data$outcome)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$subdistribution_hr_low, xmax = .data$subdistribution_hr_high), orientation = "y", width = 0.18, linewidth = 0.85, color = "grey20") +
  ggplot2::geom_point(size = 3, color = "grey10") +
  ggplot2::geom_text(ggplot2::aes(x = .data$label_x, label = .data$label), hjust = 0, size = 3.25, color = "grey12") +
  ggplot2::scale_x_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.25, 1.5), labels = c("0.8", "0.9", "1", "1.1", "1.25", "1.5")) +
  ggplot2::coord_cartesian(xlim = c(0.82, max(imv24_coef$label_x, na.rm = TRUE) * 1.35), clip = "off") +
  ggplot2::labs(
    title = "IMV-24 temperature associations",
    x = "Subdistribution HR per 5°C warmer admission day",
    y = NULL
  ) +
  theme_candidate() +
  ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), plot.margin = ggplot2::margin(12, 36, 12, 12))
save_plot(fig_imv24_forest, "candidate_imv24_temperature_effect_forest", width = 7.2, height = 3.6)

imv12_coef <- read_required("pooled_ohca_icu_imv12_discharge_outcome_coefficients.csv") |>
  dplyr::filter(.data$coefficient == "tmax_per5c") |>
  dplyr::mutate(
    outcome_label = dplyr::recode(
      .data$outcome_level,
      death_hospice = "Death/hospice",
      ltach_trach = "LTACH with trach",
      alive_snf = "Alive discharge SNF"
    ),
    outcome = factor(.data$outcome_label, levels = rev(c("Death/hospice", "LTACH with trach", "Alive discharge SNF"))),
    label = fmt_effect(.data$odds_ratio, .data$odds_ratio_low, .data$odds_ratio_high),
    label_x = .data$odds_ratio_high * 1.035
  )
readr::write_csv(imv12_coef, file.path(figure_dir, "candidate_imv12_discharge_temperature_effect_forest_source.csv"))

fig_imv12_forest <- ggplot2::ggplot(imv12_coef, ggplot2::aes(x = .data$odds_ratio, y = .data$outcome)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$odds_ratio_low, xmax = .data$odds_ratio_high), orientation = "y", width = 0.18, linewidth = 0.85, color = "grey20") +
  ggplot2::geom_point(size = 3, color = "grey10") +
  ggplot2::geom_text(ggplot2::aes(x = .data$label_x, label = .data$label), hjust = 0, size = 3.25, color = "grey12") +
  ggplot2::scale_x_log10(breaks = c(0.75, 0.9, 1, 1.1, 1.25), labels = c("0.75", "0.9", "1", "1.1", "1.25")) +
  ggplot2::coord_cartesian(xlim = c(0.78, max(imv12_coef$label_x, na.rm = TRUE) * 1.35), clip = "off") +
  ggplot2::labs(
    title = "IMV >=12h discharge-outcome temperature associations",
    x = "Odds ratio per 5°C warmer admission day vs alive home",
    y = NULL
  ) +
  theme_candidate() +
  ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), plot.margin = ggplot2::margin(12, 36, 12, 12))
save_plot(fig_imv12_forest, "candidate_imv12_discharge_temperature_effect_forest", width = 7.2, height = 3.6)

imv12_summary <- read_required("pooled_ohca_icu_imv12_discharge_outcome_summary.csv") |>
  dplyr::mutate(
    outcome = dplyr::recode(
      .data$discharge_outcome,
      alive_home = "Alive home",
      death_hospice = "Death/hospice",
      ltach_trach = "LTACH with trach",
      alive_snf = "Alive SNF",
      excluded_other_unknown = "Other/unknown"
    ),
    outcome = factor(.data$outcome, levels = rev(c("Death/hospice", "Alive home", "Other/unknown", "LTACH with trach", "Alive SNF"))),
    pct = .data$pct_among_imv12 / 100,
    label = sprintf("%d (%.1f%%)", .data$n, .data$pct_among_imv12)
  )
readr::write_csv(imv12_summary, file.path(figure_dir, "candidate_imv12_discharge_outcome_composition_source.csv"))

fig_imv12_composition <- ggplot2::ggplot(imv12_summary, ggplot2::aes(x = .data$outcome, y = .data$pct, fill = .data$outcome)) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  ggplot2::geom_text(ggplot2::aes(label = .data$label), hjust = -0.08, size = 3.2, color = "grey15") +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(
    values = c(
      "Death/hospice" = "#9b2f37",
      "Alive home" = "#0f6b78",
      "Alive SNF" = "#8f7a2f",
      "LTACH with trach" = "#6f4c8b",
      "Other/unknown" = "#7a7a7a"
    ),
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.68), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::labs(
    title = "IMV >=12h discharge outcomes",
    x = NULL,
    y = "Proportion of IMV >=12h cohort"
  ) +
  theme_candidate() +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(12, 38, 12, 12)
  )
save_plot(fig_imv12_composition, "candidate_imv12_discharge_outcome_composition", width = 7.2, height = 4.2)

imv12_all_summary <- read_required("all_site_ohca_icu_imv12_discharge_outcome_summary.csv") |>
  dplyr::select("site_name", discharge_outcome = "discharge_outcome", outcome_n = "n")

imv12_table <- read_required("all_site_ohca_icu_imv12_discharge_outcome_table_by_outcome.csv")
evidence_rows <- imv12_table |>
  dplyr::filter(
    .data$level == "Yes",
    .data$characteristic %in% c("Any tracheostomy evidence", "Anoxic brain injury ICD", "Brain death ICD", "Hospital death", "Hospice/comfort discharge")
  ) |>
  dplyr::select(
    "site_name",
    evidence = "characteristic",
    alive_home = "alive_home",
    death_hospice = "death_hospice",
    ltach_trach = "ltach_trach",
    alive_snf = "alive_snf",
    excluded_other_unknown = "excluded_other_unknown"
  )

outcome_columns <- c("alive_home", "death_hospice", "ltach_trach", "alive_snf", "excluded_other_unknown")
evidence_long <- do.call(
  rbind,
  lapply(outcome_columns, function(outcome) {
    data.frame(
      site_name = evidence_rows$site_name,
      evidence = evidence_rows$evidence,
      discharge_outcome = outcome,
      numerator = parse_count(evidence_rows[[outcome]]),
      stringsAsFactors = FALSE
    )
  })
) |>
  dplyr::left_join(imv12_all_summary, by = c("site_name", "discharge_outcome")) |>
  dplyr::group_by(.data$evidence, .data$discharge_outcome) |>
  dplyr::summarise(
    numerator = sum(.data$numerator, na.rm = TRUE),
    denominator = sum(.data$outcome_n, na.rm = TRUE),
    prevalence = .data$numerator / .data$denominator,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    evidence = factor(
      .data$evidence,
      levels = c("Hospital death", "Hospice/comfort discharge", "Anoxic brain injury ICD", "Brain death ICD", "Any tracheostomy evidence")
    ),
    outcome = dplyr::recode(
      .data$discharge_outcome,
      alive_home = "Alive home",
      death_hospice = "Death/hospice",
      ltach_trach = "LTACH with trach",
      alive_snf = "Alive SNF",
      excluded_other_unknown = "Other/unknown"
    ),
    outcome = factor(.data$outcome, levels = c("Alive home", "Alive SNF", "LTACH with trach", "Death/hospice", "Other/unknown"))
  )
readr::write_csv(evidence_long, file.path(figure_dir, "candidate_imv12_evidence_by_discharge_outcome_source.csv"))

fig_imv12_evidence <- ggplot2::ggplot(evidence_long, ggplot2::aes(x = .data$outcome, y = .data$prevalence)) +
  ggplot2::geom_col(width = 0.68, fill = "grey25") +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::percent(.data$prevalence, accuracy = 1)),
    vjust = -0.3,
    size = 2.8,
    color = "grey15"
  ) +
  ggplot2::facet_wrap(ggplot2::vars(.data$evidence), ncol = 1, scales = "free_y") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.14))) +
  ggplot2::labs(
    title = "IMV >=12h outcome evidence by discharge outcome",
    x = NULL,
    y = "Prevalence"
  ) +
  theme_candidate(base_size = 10.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    panel.spacing.y = grid::unit(0.55, "lines")
  )
save_plot(fig_imv12_evidence, "candidate_imv12_evidence_by_discharge_outcome", width = 7.4, height = 8.8)

message("Wrote post-admission candidate figures to ", figure_dir)

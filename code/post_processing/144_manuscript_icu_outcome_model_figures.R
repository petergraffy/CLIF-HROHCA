#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "scales", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_if_exists <- function(file_name) {
  path <- file.path(pooled_dir, file_name)
  if (!file.exists(path)) return(tibble::tibble())
  readr::read_csv(path, show_col_types = FALSE)
}

save_plot <- function(plot, stem, width = 7.4, height = 5.1) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 600)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height)
  invisible(c(png_path, pdf_path))
}

theme_outcome <- function(base_size = 11.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "grey15"),
      strip.background = ggplot2::element_rect(fill = "grey93", color = "grey35", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 0.5),
      legend.position = "bottom"
    )
}

format_effect <- function(est, low, high) {
  sprintf("%.2f (%.2f, %.2f)", est, low, high)
}

outcome_colors <- c(
  "Awake/extubated" = "#0f6b78",
  "Death before awake/extubated" = "#9b2f37",
  "Successful extubation" = "#0f6b78",
  "Death before extubation" = "#9b2f37",
  "Death/hospice discharge" = "#9b2f37",
  "Tracheostomy" = "#8f7a2f",
  "Alive home" = "#0f6b78",
  "Death/hospice" = "#9b2f37",
  "LTACH with trach" = "#6f4c8b",
  "Alive SNF" = "#8f7a2f",
  "No IMV in first 72h" = "#335c67",
  "Extubated by 72h" = "#0f6b78",
  "On IMV at 72h" = "#8f7a2f",
  "Death within 72h" = "#9b2f37"
)

quartile_levels <- c("Q1 coolest", "Q2", "Q3", "Q4 hottest")
quartile_colors <- c(
  "Q1 coolest" = "#335c67",
  "Q2" = "#8f7a2f",
  "Q3" = "#c06c3e",
  "Q4 hottest" = "#9b2f37"
)

phenotype_label <- function(x) {
  dplyr::recode(
    x,
    alive_no_imv = "No IMV in first 72h",
    regained_consciousness_extubated = "Extubated by 72h",
    limited_brain_function = "On IMV at 72h",
    anoxic_brain_injury = "Death within 72h",
    .default = x
  )
}

discharge_label <- function(x) {
  dplyr::recode(
    x,
    alive_home = "Alive home",
    death_hospice = "Death/hospice",
    ltach_trach = "LTACH with trach",
    alive_snf = "Alive SNF",
    excluded_other_unknown = "Other/unknown",
    .default = x
  )
}

phenotype_primary <- read_if_exists("pooled_ohca_icu_72h_phenotype_assignment_temperature_curve.csv")
phenotype_mechanism <- read_if_exists("pooled_ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted.csv")
phenotype_curves <- dplyr::bind_rows(
  dplyr::mutate(phenotype_primary, model_label = "Primary"),
  dplyr::mutate(phenotype_mechanism, model_label = "Mechanism-adjusted")
)
if (nrow(phenotype_curves) > 0) {
  phenotype_plot_df <- phenotype_curves |>
    dplyr::mutate(
      phenotype_label = phenotype_label(.data$phenotype),
      phenotype_label = factor(.data$phenotype_label, levels = c(
        "No IMV in first 72h", "Extubated by 72h", "On IMV at 72h", "Death within 72h"
      )),
      model_label = factor(.data$model_label, levels = c("Primary", "Mechanism-adjusted"))
    ) |>
    dplyr::filter(!is.na(.data$phenotype_label))
  readr::write_csv(phenotype_plot_df, file.path(figure_dir, "pooled_icu_outcome_72h_phenotype_temperature_curves_source.csv"))
  phenotype_plot <- ggplot2::ggplot(
    phenotype_plot_df,
    ggplot2::aes(x = .data$tmax_mean_c, y = .data$predicted_probability, color = .data$phenotype_label, fill = .data$phenotype_label)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$predicted_probability_low, ymax = .data$predicted_probability_high),
      alpha = 0.12,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::facet_wrap(ggplot2::vars(.data$model_label), nrow = 1) +
    ggplot2::scale_color_manual(values = outcome_colors, name = NULL, drop = FALSE) +
    ggplot2::scale_fill_manual(values = outcome_colors, guide = "none", drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "72-hour phenotype probability by admission temperature",
      x = "Admission-day Tmax (C)",
      y = "Predicted probability"
    ) +
    theme_outcome()
  save_plot(phenotype_plot, "pooled_icu_outcome_72h_phenotype_temperature_curves", width = 10.2, height = 4.8)
}

fg_cif <- read_if_exists("pooled_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv")
if (nrow(fg_cif) > 0) {
  fg_cif <- fg_cif |>
    dplyr::mutate(
      event_type = dplyr::recode(
        .data$event_type,
        "Awake/extubated" = "Successful extubation",
        "Death before awake/extubated" = "Death before extubation",
        .default = .data$event_type
      )
    )
  fg_overall <- fg_cif |>
    dplyr::filter(.data$stratification == "Overall") |>
    dplyr::mutate(event_type = factor(.data$event_type, levels = c("Successful extubation", "Death before extubation")))
  readr::write_csv(fg_overall, file.path(figure_dir, "pooled_icu_outcome_72h_fine_gray_cif_overall_source.csv"))
  fg_overall_plot <- ggplot2::ggplot(fg_overall, ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$event_type, fill = .data$event_type)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$cif_low, ymax = .data$cif_high), alpha = 0.13, color = NA) +
    ggplot2::geom_step(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = outcome_colors, name = NULL, drop = FALSE) +
    ggplot2::scale_fill_manual(values = outcome_colors, guide = "none", drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "72-hour Fine-Gray cumulative incidence",
      x = "Hours since ICU admission",
      y = "Cumulative incidence"
    ) +
    theme_outcome()
  save_plot(fg_overall_plot, "pooled_icu_outcome_72h_fine_gray_cif_overall", width = 7.4, height = 4.8)

  fg_quartile <- fg_cif |>
    dplyr::filter(.data$stratification == "Admission Tmax quartile") |>
    dplyr::mutate(
      stratum = factor(.data$stratum, levels = quartile_levels),
      event_type = factor(.data$event_type, levels = c("Successful extubation", "Death before extubation"))
    )
  readr::write_csv(fg_quartile, file.path(figure_dir, "pooled_icu_outcome_72h_fine_gray_cif_by_tmax_quartile_source.csv"))
  fg_quartile_plot <- ggplot2::ggplot(fg_quartile, ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$stratum)) +
    ggplot2::geom_step(linewidth = 1.0) +
    ggplot2::facet_wrap(ggplot2::vars(.data$event_type), nrow = 1) +
    ggplot2::scale_color_manual(values = quartile_colors, name = "Admission Tmax quartile", drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 24), limits = c(0, 72), expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "72-hour Fine-Gray cumulative incidence by admission Tmax quartile",
      x = "Hours since ICU admission",
      y = "Cumulative incidence"
    ) +
    theme_outcome(base_size = 11)
  save_plot(fg_quartile_plot, "pooled_icu_outcome_72h_fine_gray_cif_by_tmax_quartile", width = 10.2, height = 4.8)
}

fg_coef <- read_if_exists("pooled_ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv")
if (nrow(fg_coef) > 0) {
  fg_linear <- fg_coef |>
    dplyr::filter(
      .data$coefficient == "tmax_per1c",
      .data$model %in% c(
        "awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_linear",
        "awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_mechanism_linear"
      )
    ) |>
    dplyr::mutate(
      model_label = dplyr::recode(
        .data$model,
        awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_linear = "Primary",
        awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_mechanism_linear = "Mechanism-adjusted"
      ),
      log_effect_per5c = .data$estimate * 5,
      se_per5c = .data$standard_error * 5,
      ratio = exp(.data$log_effect_per5c),
      ratio_low = exp(.data$log_effect_per5c - 1.96 * .data$se_per5c),
      ratio_high = exp(.data$log_effect_per5c + 1.96 * .data$se_per5c),
      outcome_label = .data$model_label
    )
  readr::write_csv(fg_linear, file.path(figure_dir, "pooled_icu_outcome_72h_fine_gray_linear_temperature_forest_source.csv"))
  if (nrow(fg_linear) > 0) {
    fg_linear_plot <- ggplot2::ggplot(fg_linear, ggplot2::aes(x = .data$ratio, y = .data$outcome_label, color = .data$model_label)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.45) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ratio_low, xmax = .data$ratio_high), orientation = "y", width = 0.16, linewidth = 0.85) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = c("Primary" = "#0f6b78", "Mechanism-adjusted" = "#9b2f37"), guide = "none") +
      ggplot2::scale_x_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.25), labels = c("0.8", "0.9", "1", "1.1", "1.25")) +
      ggplot2::labs(
        title = "72-hour Fine-Gray linear Tmax association",
        x = "Subdistribution HR per 5C warmer admission day",
        y = NULL
      ) +
      theme_outcome()
    save_plot(fg_linear_plot, "pooled_icu_outcome_72h_fine_gray_linear_temperature_forest", width = 6.6, height = 3.2)
  }
}

imv24_cif <- read_if_exists("pooled_ohca_icu_imv24_time_to_event_cif_curves.csv")
if (nrow(imv24_cif) > 0) {
  imv24_event_order <- c("Successful extubation", "Death/hospice discharge", "Tracheostomy")
  imv24_overall <- imv24_cif |>
    dplyr::filter(.data$stratification == "Overall", .data$time_days <= 30) |>
    dplyr::mutate(event_type = factor(.data$event_type, levels = imv24_event_order))
  readr::write_csv(imv24_overall, file.path(figure_dir, "pooled_icu_outcome_imv24_cif_overall_source.csv"))
  imv24_overall_plot <- ggplot2::ggplot(imv24_overall, ggplot2::aes(x = .data$time_days, y = .data$cif, color = .data$event_type, fill = .data$event_type)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$cif_low, ymax = .data$cif_high), alpha = 0.13, color = NA) +
    ggplot2::geom_step(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = outcome_colors, name = NULL, drop = FALSE) +
    ggplot2::scale_fill_manual(values = outcome_colors, guide = "none", drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "IMV-24 competing-risk cumulative incidence",
      x = "Days since 24 hours of IMV",
      y = "Cumulative incidence"
    ) +
    theme_outcome()
  save_plot(imv24_overall_plot, "pooled_icu_outcome_imv24_cif_overall", width = 7.4, height = 4.8)

  imv24_quartile <- imv24_cif |>
    dplyr::filter(.data$stratification == "Admission Tmax quartile", .data$time_days <= 30) |>
    dplyr::mutate(
      stratum = factor(.data$stratum, levels = quartile_levels),
      event_type = factor(.data$event_type, levels = imv24_event_order)
    )
  readr::write_csv(imv24_quartile, file.path(figure_dir, "pooled_icu_outcome_imv24_cif_by_tmax_quartile_source.csv"))
  imv24_quartile_plot <- ggplot2::ggplot(imv24_quartile, ggplot2::aes(x = .data$time_days, y = .data$cif, color = .data$stratum)) +
    ggplot2::geom_step(linewidth = 1.0) +
    ggplot2::facet_wrap(ggplot2::vars(.data$event_type), nrow = 1) +
    ggplot2::scale_color_manual(values = quartile_colors, name = "Admission Tmax quartile", drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq(0, 30, by = 10), limits = c(0, 30), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "IMV-24 cumulative incidence by admission Tmax quartile",
      x = "Days since 24 hours of IMV",
      y = "Cumulative incidence"
    ) +
    theme_outcome(base_size = 10.5)
  save_plot(imv24_quartile_plot, "pooled_icu_outcome_imv24_cif_by_tmax_quartile", width = 10.2, height = 4.8)
}

imv24_coef <- read_if_exists("pooled_ohca_icu_imv24_time_to_event_coefficients.csv")
if (nrow(imv24_coef) > 0) {
  imv24_forest <- imv24_coef |>
    dplyr::filter(.data$coefficient == "tmax_per5c") |>
    dplyr::mutate(
      outcome_label = factor(.data$event_type, levels = rev(c("Successful extubation", "Death/hospice discharge", "Tracheostomy"))),
      ratio = .data$subdistribution_hr,
      ratio_low = .data$subdistribution_hr_low,
      ratio_high = .data$subdistribution_hr_high,
      label = format_effect(.data$ratio, .data$ratio_low, .data$ratio_high),
      label_x = .data$ratio_high * 1.035
    )
  readr::write_csv(imv24_forest, file.path(figure_dir, "pooled_icu_outcome_imv24_temperature_forest_source.csv"))
  if (nrow(imv24_forest) > 0) {
    imv24_forest_plot <- ggplot2::ggplot(imv24_forest, ggplot2::aes(x = .data$ratio, y = .data$outcome_label)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.45) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ratio_low, xmax = .data$ratio_high), orientation = "y", width = 0.18, linewidth = 0.85, color = "grey20") +
      ggplot2::geom_point(size = 3, color = "grey10") +
      ggplot2::geom_text(ggplot2::aes(x = .data$label_x, label = .data$label), hjust = 0, size = 3.2, color = "grey12") +
      ggplot2::scale_x_log10(breaks = c(0.75, 0.9, 1, 1.1, 1.25), labels = c("0.75", "0.9", "1", "1.1", "1.25")) +
      ggplot2::coord_cartesian(xlim = c(0.74, max(imv24_forest$label_x, na.rm = TRUE) * 1.32), clip = "off") +
      ggplot2::labs(
        title = "IMV-24 temperature associations",
        x = "Subdistribution HR per 5C warmer admission day",
        y = NULL
      ) +
      theme_outcome() +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), plot.margin = ggplot2::margin(12, 38, 12, 12))
    save_plot(imv24_forest_plot, "pooled_icu_outcome_imv24_temperature_forest", width = 7.2, height = 3.6)
  }
}

imv12_curve <- read_if_exists("pooled_ohca_icu_imv12_discharge_outcome_temperature_curve.csv")
if (nrow(imv12_curve) > 0) {
  imv12_curve_plot_df <- imv12_curve |>
    dplyr::mutate(
      outcome_label = discharge_label(.data$discharge_outcome),
      outcome_label = factor(.data$outcome_label, levels = c("Alive home", "Alive SNF", "LTACH with trach", "Death/hospice"))
    ) |>
    dplyr::filter(!is.na(.data$outcome_label))
  readr::write_csv(imv12_curve_plot_df, file.path(figure_dir, "pooled_icu_outcome_imv12_discharge_temperature_curves_source.csv"))
  imv12_curve_plot <- ggplot2::ggplot(imv12_curve_plot_df, ggplot2::aes(x = .data$tmax_mean_c, y = .data$predicted_probability, color = .data$outcome_label)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::scale_color_manual(values = outcome_colors, name = NULL, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "IMV >=12h discharge outcome probability by admission temperature",
      x = "Admission-day Tmax (C)",
      y = "Predicted probability"
    ) +
    theme_outcome()
  save_plot(imv12_curve_plot, "pooled_icu_outcome_imv12_discharge_temperature_curves", width = 7.4, height = 4.8)
}

imv12_coef <- read_if_exists("pooled_ohca_icu_imv12_discharge_outcome_coefficients.csv")
if (nrow(imv12_coef) > 0) {
  imv12_forest <- imv12_coef |>
    dplyr::filter(.data$coefficient == "tmax_per5c") |>
    dplyr::mutate(
      outcome_label = discharge_label(.data$outcome_level),
      outcome_label = factor(.data$outcome_label, levels = rev(c("Death/hospice", "LTACH with trach", "Alive SNF"))),
      ratio = .data$odds_ratio,
      ratio_low = .data$odds_ratio_low,
      ratio_high = .data$odds_ratio_high,
      label = format_effect(.data$ratio, .data$ratio_low, .data$ratio_high),
      label_x = .data$ratio_high * 1.035
    )
  readr::write_csv(imv12_forest, file.path(figure_dir, "pooled_icu_outcome_imv12_discharge_temperature_forest_source.csv"))
  if (nrow(imv12_forest) > 0) {
    imv12_forest_plot <- ggplot2::ggplot(imv12_forest, ggplot2::aes(x = .data$ratio, y = .data$outcome_label)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.45) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ratio_low, xmax = .data$ratio_high), orientation = "y", width = 0.18, linewidth = 0.85, color = "grey20") +
      ggplot2::geom_point(size = 3, color = "grey10") +
      ggplot2::geom_text(ggplot2::aes(x = .data$label_x, label = .data$label), hjust = 0, size = 3.2, color = "grey12") +
      ggplot2::scale_x_log10(breaks = c(0.75, 0.9, 1, 1.1, 1.25), labels = c("0.75", "0.9", "1", "1.1", "1.25")) +
      ggplot2::coord_cartesian(xlim = c(0.74, max(imv12_forest$label_x, na.rm = TRUE) * 1.32), clip = "off") +
      ggplot2::labs(
        title = "IMV >=12h discharge-outcome temperature associations",
        x = "Odds ratio per 5C warmer admission day vs alive home",
        y = NULL
      ) +
      theme_outcome() +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), plot.margin = ggplot2::margin(12, 38, 12, 12))
    save_plot(imv12_forest_plot, "pooled_icu_outcome_imv12_discharge_temperature_forest", width = 7.2, height = 3.6)
  }
}

combined_effects <- dplyr::bind_rows(
  if (exists("fg_linear") && nrow(fg_linear) > 0) {
    fg_linear |>
      dplyr::transmute(
        model_family = "72h Fine-Gray",
        outcome_label = paste0("Awake/extubated - ", .data$model_label),
        effect_measure = "SHR",
        ratio = .data$ratio,
        ratio_low = .data$ratio_low,
        ratio_high = .data$ratio_high,
        k_sites = .data$k_sites
      )
  },
  if (exists("imv24_forest") && nrow(imv24_forest) > 0) {
    imv24_forest |>
      dplyr::transmute(
        model_family = "IMV-24 Fine-Gray",
        outcome_label = as.character(.data$outcome_label),
        effect_measure = "SHR",
        ratio = .data$ratio,
        ratio_low = .data$ratio_low,
        ratio_high = .data$ratio_high,
        k_sites = .data$k_sites
      )
  },
  if (exists("imv12_forest") && nrow(imv12_forest) > 0) {
    imv12_forest |>
      dplyr::transmute(
        model_family = "IMV >=12h discharge multinomial",
        outcome_label = paste0(as.character(.data$outcome_label), " vs alive home"),
        effect_measure = "OR",
        ratio = .data$ratio,
        ratio_low = .data$ratio_low,
        ratio_high = .data$ratio_high,
        k_sites = .data$k_sites
      )
  }
)

if (nrow(combined_effects) > 0) {
  combined_effects <- combined_effects |>
    dplyr::mutate(
      row_label = paste0(.data$model_family, ": ", .data$outcome_label),
      row_label = factor(.data$row_label, levels = rev(unique(.data$row_label))),
      label = format_effect(.data$ratio, .data$ratio_low, .data$ratio_high),
      label_x = .data$ratio_high * 1.035
    )
  readr::write_csv(combined_effects, file.path(figure_dir, "pooled_icu_outcome_temperature_effect_forest_source.csv"))
  combined_plot <- ggplot2::ggplot(combined_effects, ggplot2::aes(x = .data$ratio, y = .data$row_label, color = .data$model_family)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.45) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ratio_low, xmax = .data$ratio_high), orientation = "y", width = 0.18, linewidth = 0.85) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(x = .data$label_x, label = .data$label), hjust = 0, size = 3.0, color = "grey12") +
    ggplot2::scale_color_manual(
      values = c(
        "72h Fine-Gray" = "#0f6b78",
        "IMV-24 Fine-Gray" = "#8f7a2f",
        "IMV >=12h discharge multinomial" = "#6f4c8b"
      ),
      name = NULL
    ) +
    ggplot2::scale_x_log10(breaks = c(0.75, 0.9, 1, 1.1, 1.25), labels = c("0.75", "0.9", "1", "1.1", "1.25")) +
    ggplot2::coord_cartesian(xlim = c(0.74, max(combined_effects$label_x, na.rm = TRUE) * 1.34), clip = "off") +
    ggplot2::labs(
      title = "Temperature associations across OHCA ICU outcome models",
      x = "Ratio per 5C warmer admission day",
      y = NULL
    ) +
    theme_outcome(base_size = 11) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(12, 42, 12, 12)
    )
  save_plot(combined_plot, "pooled_icu_outcome_temperature_effect_forest", width = 9.2, height = 5.8)
}

message("Wrote pooled ICU outcome model figures to ", figure_dir)

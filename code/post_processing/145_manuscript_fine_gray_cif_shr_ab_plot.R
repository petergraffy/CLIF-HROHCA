#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

cif_path <- file.path(figure_dir, "pooled_icu_outcome_72h_fine_gray_cif_by_tmax_quartile_source.csv")
coef_path <- file.path(pooled_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv")
if (!file.exists(cif_path)) stop("Missing CIF source file: ", cif_path, call. = FALSE)
if (!file.exists(coef_path)) stop("Missing pooled Fine-Gray coefficient file: ", coef_path, call. = FALSE)

quartile_levels <- c("Q1 coolest", "Q2", "Q3", "Q4 hottest")
quartile_colors <- c(
  "Q1 coolest" = "#335c67",
  "Q2" = "#8f7a2f",
  "Q3" = "#c06c3e",
  "Q4 hottest" = "#9b2f37"
)

cif <- readr::read_csv(cif_path, show_col_types = FALSE) |>
  dplyr::mutate(
    event_type = dplyr::recode(
      .data$event_type,
      "Awake/extubated" = "Successful extubation",
      "Death before awake/extubated" = "Death before extubation",
      .default = .data$event_type
    ),
    event_type = factor(.data$event_type, levels = c("Successful extubation", "Death before extubation")),
    stratum = factor(.data$stratum, levels = quartile_levels)
  ) |>
  dplyr::filter(!is.na(.data$event_type), !is.na(.data$stratum))

shr_source <- readr::read_csv(coef_path, show_col_types = FALSE) |>
  dplyr::filter(
    grepl("spline", .data$coefficient),
    .data$model %in% c(
      "awake_extubated_72h_vs_death_72h_temperature_humidity_demographics",
      "awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_mechanism"
    )
  ) |>
  dplyr::mutate(
    model_label = dplyr::recode(
      .data$model,
      awake_extubated_72h_vs_death_72h_temperature_humidity_demographics = "Primary",
      awake_extubated_72h_vs_death_72h_temperature_humidity_demographics_mechanism = "Mechanism-adjusted"
    ),
    exposure = dplyr::case_when(
      grepl("^tmax", .data$coefficient) ~ "Tmax",
      grepl("^rmax", .data$coefficient) ~ "Rmax",
      TRUE ~ NA_character_
    ),
    spline_num = as.integer(sub(".*spline[^0-9]*([0-9]+)$", "\\1", .data$coefficient)),
    term_label = paste(.data$exposure, "spline", .data$spline_num),
    ci_label = sprintf(
      "%.2f (%.2f, %.2f)",
      .data$subdistribution_hr,
      .data$subdistribution_hr_low,
      .data$subdistribution_hr_high
    )
  ) |>
  dplyr::filter(!is.na(.data$exposure), is.finite(.data$spline_num))

readr::write_csv(
  shr_source,
  file.path(figure_dir, "pooled_fine_gray_primary_spline_shr_forest_source.csv")
)
readr::write_csv(
  shr_source,
  file.path(figure_dir, "figure_fine_gray_spline_term_forest_source.csv")
)

shr <- shr_source |>
  dplyr::filter(.data$model_label == "Primary") |>
  dplyr::mutate(
    term_plot = dplyr::case_when(
      .data$exposure == "Tmax" ~ paste0("T[max]~spline~", .data$spline_num),
      .data$exposure == "Rmax" ~ paste0("R[max]~spline~", .data$spline_num),
      TRUE ~ .data$term_label
    ),
    term_order = factor(
      .data$term_plot,
      levels = rev(c(
        "T[max]~spline~1",
        "T[max]~spline~2",
        "T[max]~spline~3",
        "R[max]~spline~1",
        "R[max]~spline~2",
        "R[max]~spline~3"
      ))
    ),
    estimate_label = sprintf(
      "%.2f (%.2f, %.2f)",
      .data$subdistribution_hr,
      .data$subdistribution_hr_low,
      .data$subdistribution_hr_high
    )
  ) |>
  dplyr::filter(!is.na(.data$term_order))

readr::write_csv(
  cif,
  file.path(figure_dir, "fine_gray_cif_shr_ab_plot_panel_a_source.csv")
)
readr::write_csv(
  shr,
  file.path(figure_dir, "fine_gray_cif_shr_ab_plot_panel_b_source.csv")
)

y_limit <- max(cif$cif, cif$cif_high, na.rm = TRUE)
y_limit <- min(0.65, max(0.35, ceiling((y_limit + 0.035) * 20) / 20))
x_max <- max(shr$subdistribution_hr_high, na.rm = TRUE) * 1.18

theme_panel <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      axis.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.text = ggplot2::element_text(color = "grey15"),
      strip.background = ggplot2::element_rect(fill = "grey93", color = "grey35", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      panel.border = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 0.45),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.title = ggplot2::element_text(size = base_size - 0.5)
    )
}

panel_a <- ggplot2::ggplot(cif, ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$stratum)) +
  ggplot2::geom_step(linewidth = 0.95) +
  ggplot2::facet_wrap(ggplot2::vars(.data$event_type), nrow = 1) +
  ggplot2::scale_color_manual(values = quartile_colors, name = "Admission Tmax quartile", drop = FALSE) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 72, by = 24),
    limits = c(0, 72),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, y_limit),
    expand = ggplot2::expansion(mult = c(0, 0.04))
  ) +
  ggplot2::labs(
    title = "A. Cumulative incidence by admission temperature quartile",
    x = "Hours since ICU admission",
    y = "Cumulative incidence"
  ) +
  theme_panel() +
  ggplot2::theme(
    panel.spacing.x = grid::unit(10, "pt"),
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = ggplot2::margin(5.5, 5.5, 14, 5.5)
  )

panel_b <- ggplot2::ggplot(shr, ggplot2::aes(x = .data$subdistribution_hr, y = .data$term_order)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, color = "grey45") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$subdistribution_hr_low, xmax = .data$subdistribution_hr_high),
    orientation = "y",
    width = 0,
    linewidth = 0.75,
    color = "grey25"
  ) +
  ggplot2::geom_point(size = 2.45, color = "#1F2937") +
  ggplot2::geom_text(
    ggplot2::aes(
      x = pmin(.data$subdistribution_hr_high * 1.08, x_max),
      label = .data$estimate_label
    ),
    hjust = 0,
    size = 3.9,
    color = "grey15"
  ) +
  ggplot2::scale_y_discrete(labels = function(x) parse(text = x)) +
  ggplot2::scale_x_log10(
    limits = c(0.08, x_max),
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
    labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10")
  ) +
  ggplot2::labs(
    title = "B. Subdistribution hazard of extubation vs death in 72h",
    x = "Pooled subdistribution hazard ratio (95% CI)",
    y = NULL
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_panel(base_size = 13) +
  ggplot2::theme(
    legend.position = "none",
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(color = "grey88", linewidth = 0.45),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5.5, 58, 5.5, 5.5)
  )

combined <- panel_a / panel_b +
  patchwork::plot_layout(heights = c(0.57, 0.43), guides = "keep")

png_path <- file.path(figure_dir, "fine_gray_cif_shr_ab_plot.png")
pdf_path <- file.path(figure_dir, "fine_gray_cif_shr_ab_plot.pdf")
ggplot2::ggsave(png_path, combined, width = 10.4, height = 8.4, dpi = 600, bg = "white")
ggplot2::ggsave(pdf_path, combined, width = 10.4, height = 8.4, bg = "white")

message("Wrote ", png_path)
message("Wrote ", pdf_path)

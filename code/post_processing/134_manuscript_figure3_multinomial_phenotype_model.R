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

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

site_path <- file.path(pooled_dir, "all_site_ohca_icu_72h_phenotype_assignment_temperature_curve_interpolated.csv")
pooled_path <- file.path(pooled_dir, "pooled_ohca_icu_72h_phenotype_assignment_temperature_curve.csv")
summary_path <- file.path(pooled_dir, "pooled_ohca_icu_72h_phenotype_summary.csv")

site_curves <- read_required(site_path)
pooled_curves <- read_required(pooled_path)
phenotype_summary <- read_required(summary_path)

target_model <- "phenotype_assignment_temperature_humidity_demographics"
phenotype_order <- c(
  "alive_no_imv",
  "regained_consciousness_extubated",
  "limited_brain_function",
  "anoxic_brain_injury"
)
phenotype_labels <- c(
  alive_no_imv = "No IMV in first 72h",
  regained_consciousness_extubated = "Extubated by 72h",
  limited_brain_function = "On IMV at 72h",
  anoxic_brain_injury = "Death within 72h"
)
phenotype_colors <- c(
  alive_no_imv = "#457B9D",
  regained_consciousness_extubated = "#0f6b78",
  limited_brain_function = "#8f7a2f",
  anoxic_brain_injury = "#9b2f37"
)

modeled_counts <- phenotype_summary |>
  dplyr::filter(.data$phenotype %in% phenotype_order) |>
  dplyr::mutate(
    phenotype = factor(.data$phenotype, levels = phenotype_order),
    phenotype_label = phenotype_labels[as.character(.data$phenotype)]
  )
facet_labels <- phenotype_labels

site_plot <- site_curves |>
  dplyr::filter(.data$model == target_model, .data$phenotype %in% phenotype_order) |>
  dplyr::mutate(
    phenotype = factor(.data$phenotype, levels = phenotype_order),
    phenotype_label = phenotype_labels[as.character(.data$phenotype)]
  )

pooled_plot <- pooled_curves |>
  dplyr::filter(.data$model == target_model, .data$phenotype %in% phenotype_order) |>
  dplyr::mutate(
    phenotype = factor(.data$phenotype, levels = phenotype_order),
    phenotype_label = phenotype_labels[as.character(.data$phenotype)]
  ) |>
  dplyr::left_join(
    modeled_counts |> dplyr::select("phenotype", phenotype_n = "n", phenotype_pct = "pct"),
    by = "phenotype"
  )

if (nrow(site_plot) == 0 || nrow(pooled_plot) == 0) {
  stop("Could not find primary pooled multinomial phenotype temperature curves.", call. = FALSE)
}

source_table <- pooled_plot |>
  dplyr::select(
    "model",
    "curve_type",
    "phenotype",
    "phenotype_label",
    "tmax_mean_c",
    "predicted_probability",
    dplyr::any_of(c("predicted_probability_low", "predicted_probability_high", "probability_ci_method", "probability_ci_simulations")),
    "k_sites",
    "pooled_n",
    "phenotype_n",
    "phenotype_pct"
  )
readr::write_csv(source_table, file.path(figure_dir, "figure3_multinomial_phenotype_model_source.csv"))

x_limits <- range(pooled_plot$tmax_mean_c, na.rm = TRUE)
fahrenheit_breaks <- seq(45, 90, by = 15)
x_breaks_c <- (fahrenheit_breaks - 32) * 5 / 9
y_limits <- range(
  c(site_plot$predicted_probability, pooled_plot$predicted_probability),
  na.rm = TRUE
)
y_padding <- 0.025
y_min <- max(0, floor((y_limits[[1]] - y_padding) * 20) / 20)
y_max <- min(1, ceiling((y_limits[[2]] + y_padding) * 20) / 20)

if (all(c("predicted_probability_low", "predicted_probability_high") %in% names(pooled_plot))) {
  ribbon_plot <- pooled_plot |>
    dplyr::filter(
      is.finite(.data$predicted_probability_low),
      is.finite(.data$predicted_probability_high)
    )
} else {
  ribbon_plot <- pooled_plot[0, , drop = FALSE]
  ribbon_plot$predicted_probability_low <- numeric()
  ribbon_plot$predicted_probability_high <- numeric()
}

figure3 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = ribbon_plot,
    ggplot2::aes(
      x = .data$tmax_mean_c,
      ymin = .data$predicted_probability_low,
      ymax = .data$predicted_probability_high,
      fill = .data$phenotype
    ),
    alpha = 0.14,
    color = NA
  ) +
  ggplot2::geom_line(
    data = site_plot,
    ggplot2::aes(
      x = .data$tmax_mean_c,
      y = .data$predicted_probability,
      group = interaction(.data$site_name, .data$phenotype)
    ),
    color = "grey68",
    linewidth = 0.65,
    alpha = 0.68
  ) +
  ggplot2::geom_line(
    data = pooled_plot,
    ggplot2::aes(
      x = .data$tmax_mean_c,
      y = .data$predicted_probability,
      color = .data$phenotype
    ),
    linewidth = 1.35
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(.data$phenotype),
    nrow = 1,
    labeller = ggplot2::as_labeller(facet_labels)
  ) +
  ggplot2::scale_color_manual(values = phenotype_colors, guide = "none") +
  ggplot2::scale_fill_manual(values = phenotype_colors, guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = x_breaks_c,
    labels = fahrenheit_breaks,
    limits = x_limits,
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(y_min, y_max, by = 0.05),
    limits = c(y_min, y_max),
    expand = ggplot2::expansion(mult = c(0.02, 0.08))
  ) +
  ggplot2::labs(
    title = "Admission temperature and 72-hour OHCA ICU phenotype assignment",
    x = "Admission-day maximum temperature (°F)",
    y = "Predicted probability"
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    strip.background = ggplot2::element_rect(fill = "grey93", color = "grey35", linewidth = 0.55),
    strip.text = ggplot2::element_text(face = "bold", size = 10.5),
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "grey15"),
    panel.border = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 0.55),
    panel.spacing.x = grid::unit(1.8, "lines"),
    plot.margin = ggplot2::margin(12, 16, 12, 12)
  )

png_path <- file.path(figure_dir, "figure3_multinomial_phenotype_model.png")
pdf_path <- file.path(figure_dir, "figure3_multinomial_phenotype_model.pdf")
ggplot2::ggsave(png_path, figure3, width = 10.8, height = 4.9, dpi = 600)
ggplot2::ggsave(pdf_path, figure3, width = 10.8, height = 4.9)

message("Wrote multinomial phenotype model figure to ", png_path, " and ", pdf_path)

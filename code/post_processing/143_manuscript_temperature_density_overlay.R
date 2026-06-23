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

density_path <- file.path(pooled_dir, "all_site_ohca_icu_72h_admission_temperature_density.csv")
if (!file.exists(density_path)) {
  stop("Missing pooled admission-temperature density input. Run code/post_processing/115_pool_72h_phenotype_models.R first.", call. = FALSE)
}

EDGE_TRIM_F <- 2
density <- readr::read_csv(density_path, show_col_types = FALSE) |>
  dplyr::mutate(
    tmax_f = .data$tmax_mean_c * 9 / 5 + 32,
    median_tmax_f = .data$median_tmax_c * 9 / 5 + 32,
    q1_tmax_f = .data$q1_tmax_c * 9 / 5 + 32,
    q3_tmax_f = .data$q3_tmax_c * 9 / 5 + 32
  )

site_count <- dplyr::n_distinct(density$site_name)
if (site_count < 3L) {
  stop(
    "Admission-temperature density overlay requires all 3 site exports; found ",
    site_count,
    " site(s): ",
    paste(sort(unique(density$site_name)), collapse = ", "),
    call. = FALSE
  )
}

x_limits <- range(density$tmax_f, na.rm = TRUE) + c(EDGE_TRIM_F, -EDGE_TRIM_F)
plot_density <- density |>
  dplyr::filter(.data$tmax_f >= x_limits[[1]], .data$tmax_f <= x_limits[[2]])
median_lines <- plot_density |>
  dplyr::distinct(.data$site_name, .data$median_tmax_f)

readr::write_csv(plot_density, file.path(figure_dir, "figure_ohca_icu_admission_temperature_density_overlay_source.csv"))

figure_density <- ggplot2::ggplot(
  plot_density,
  ggplot2::aes(x = .data$tmax_f, y = .data$density, color = .data$site_name)
) +
  ggplot2::geom_line(linewidth = 0.9, alpha = 0.95) +
  ggplot2::geom_vline(
    data = median_lines,
    ggplot2::aes(xintercept = .data$median_tmax_f, color = .data$site_name),
    linetype = "dashed",
    linewidth = 0.55,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(
    name = "Admission Tmax (°F)",
    breaks = seq(-10, 110, 10),
    limits = x_limits,
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    name = "Density",
    labels = scales::label_number(accuracy = 0.001),
    expand = ggplot2::expansion(mult = c(0, 0.06))
  ) +
  ggplot2::labs(color = "Site") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank()
  )

png_path <- file.path(figure_dir, "figure_ohca_icu_admission_temperature_density_overlay.png")
pdf_path <- file.path(figure_dir, "figure_ohca_icu_admission_temperature_density_overlay.pdf")
ggplot2::ggsave(png_path, figure_density, width = 7.4, height = 4.6, dpi = 600)
ggplot2::ggsave(pdf_path, figure_density, width = 7.4, height = 4.6)

label_data <- plot_density |>
  dplyr::group_by(.data$site_name) |>
  dplyr::slice_max(.data$tmax_f, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

figure_density_labeled <- figure_density +
  ggplot2::geom_text(
    data = label_data,
    ggplot2::aes(label = .data$site_name),
    hjust = -0.05,
    size = 3.3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(
    name = "Admission Tmax (°F)",
    breaks = seq(-10, 110, 10),
    limits = x_limits + c(0, 5),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  )

figure_density_facets <- ggplot2::ggplot(
  plot_density,
  ggplot2::aes(x = .data$tmax_f, y = .data$density, color = .data$site_name)
) +
  ggplot2::geom_line(linewidth = 0.9, alpha = 0.95, show.legend = FALSE) +
  ggplot2::geom_vline(
    data = median_lines,
    ggplot2::aes(xintercept = .data$median_tmax_f, color = .data$site_name),
    linetype = "dashed",
    linewidth = 0.55,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(ggplot2::vars(.data$site_name), ncol = 3, scales = "free_y") +
  ggplot2::scale_x_continuous(
    name = "Admission Tmax (°F)",
    breaks = seq(-10, 110, 20),
    limits = x_limits,
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    name = "Density",
    labels = scales::label_number(accuracy = 0.001),
    expand = ggplot2::expansion(mult = c(0, 0.06))
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank()
  )

labeled_png_path <- file.path(figure_dir, "figure_ohca_icu_admission_temperature_density_overlay_labeled.png")
labeled_pdf_path <- file.path(figure_dir, "figure_ohca_icu_admission_temperature_density_overlay_labeled.pdf")
facets_png_path <- file.path(figure_dir, "figure_ohca_icu_admission_temperature_density_by_site_facets.png")
facets_pdf_path <- file.path(figure_dir, "figure_ohca_icu_admission_temperature_density_by_site_facets.pdf")
ggplot2::ggsave(labeled_png_path, figure_density_labeled, width = 8.4, height = 4.6, dpi = 600)
ggplot2::ggsave(labeled_pdf_path, figure_density_labeled, width = 8.4, height = 4.6)
ggplot2::ggsave(facets_png_path, figure_density_facets, width = 8.4, height = 6.8, dpi = 600)
ggplot2::ggsave(facets_pdf_path, figure_density_facets, width = 8.4, height = 6.8)

message(
  "Wrote admission-temperature density figures to ",
  paste(c(png_path, labeled_png_path, facets_png_path), collapse = ", ")
)

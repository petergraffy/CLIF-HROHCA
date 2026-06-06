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

read_required <- function(file_name) {
  path <- file.path(pooled_dir, file_name)
  if (!file.exists(path)) stop("Missing required pooled file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE)
}

all_site <- read_required("all_site_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv")
pooled <- read_required("pooled_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv")

if (!all(c("cif_low", "cif_high") %in% names(pooled))) {
  pooled <- pooled |>
    dplyr::mutate(
      cif_low = pmax(0, .data$cif - 1.96 * .data$standard_error),
      cif_high = pmin(1, .data$cif + 1.96 * .data$standard_error)
    )
}

event_order <- c("Awake/extubated", "Death before awake/extubated")
event_labels <- c(
  "Awake/extubated" = "Awake/extubated",
  "Death before awake/extubated" = "Death before awake/extubated"
)
event_colors <- c(
  "Awake/extubated" = "#0f6b78",
  "Death before awake/extubated" = "#9b2f37"
)
quartile_levels <- c("Q1 coolest", "Q2", "Q3", "Q4 hottest")
quartile_colors <- c(
  "Q1 coolest" = "#335c67",
  "Q2" = "#8f7a2f",
  "Q3" = "#c06c3e",
  "Q4 hottest" = "#9b2f37"
)

overall_site <- all_site |>
  dplyr::filter(.data$stratification == "Overall", .data$event_type %in% event_order) |>
  dplyr::mutate(event_type = factor(.data$event_type, levels = event_order))

overall_pooled <- pooled |>
  dplyr::filter(.data$stratification == "Overall", .data$event_type %in% event_order) |>
  dplyr::mutate(event_type = factor(.data$event_type, levels = event_order))

quartile_site <- all_site |>
  dplyr::filter(.data$stratification == "Admission Tmax quartile", .data$event_type %in% event_order) |>
  dplyr::mutate(
    event_type = factor(.data$event_type, levels = event_order),
    stratum = factor(.data$stratum, levels = quartile_levels)
  )

quartile_pooled <- pooled |>
  dplyr::filter(.data$stratification == "Admission Tmax quartile", .data$event_type %in% event_order) |>
  dplyr::mutate(
    event_type = factor(.data$event_type, levels = event_order),
    stratum = factor(.data$stratum, levels = quartile_levels)
  )

if (nrow(overall_pooled) == 0) {
  stop("Could not find pooled overall Fine-Gray CIF curves.", call. = FALSE)
}

overall_source <- overall_pooled |>
  dplyr::select(
    "model",
    "stratification",
    "stratum",
    "event_type",
    "time_hours",
    "cif",
    "cif_low",
    "cif_high",
    "standard_error",
    "n",
    "k_sites"
  )
quartile_source <- quartile_pooled |>
  dplyr::select(
    "model",
    "stratification",
    "stratum",
    "event_type",
    "time_hours",
    "cif",
    "cif_low",
    "cif_high",
    "standard_error",
    "n",
    "k_sites",
    "stratum_min_temp_c",
    "stratum_max_temp_c"
  )
readr::write_csv(overall_source, file.path(figure_dir, "figure4_fine_gray_cif_overall_source.csv"))
readr::write_csv(quartile_source, file.path(figure_dir, "figure4_fine_gray_cif_by_tmax_quartile_source.csv"))

y_max <- max(c(overall_pooled$cif_high, quartile_pooled$cif_high), na.rm = TRUE)
y_limit <- min(0.65, max(0.35, ceiling((y_max + 0.035) * 20) / 20))

panel_a <- ggplot2::ggplot() +
  ggplot2::geom_step(
    data = overall_site,
    ggplot2::aes(
      x = .data$time_hours,
      y = .data$cif,
      group = interaction(.data$site_name, .data$event_type)
    ),
    color = "grey68",
    linewidth = 0.7,
    alpha = 0.70
  ) +
  ggplot2::geom_ribbon(
    data = overall_pooled,
    ggplot2::aes(
      x = .data$time_hours,
      ymin = .data$cif_low,
      ymax = .data$cif_high,
      fill = .data$event_type
    ),
    alpha = 0.13,
    color = NA
  ) +
  ggplot2::geom_step(
    data = overall_pooled,
    ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$event_type),
    linewidth = 1.25
  ) +
  ggplot2::scale_color_manual(values = event_colors, labels = event_labels, name = NULL) +
  ggplot2::scale_fill_manual(values = event_colors, guide = "none") +
  ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, y_limit), expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(
    title = "A. Overall 72-hour cumulative incidence",
    x = "Hours since ICU admission",
    y = "Cumulative incidence"
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = c(0.18, 0.88),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.key.width = grid::unit(1.2, "lines")
  )

panel_b <- ggplot2::ggplot() +
  ggplot2::geom_step(
    data = quartile_site,
    ggplot2::aes(
      x = .data$time_hours,
      y = .data$cif,
      group = interaction(.data$site_name, .data$stratum)
    ),
    color = "grey72",
    linewidth = 0.55,
    alpha = 0.55
  ) +
  ggplot2::geom_step(
    data = quartile_pooled,
    ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$stratum),
    linewidth = 1.05
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(.data$event_type),
    nrow = 1,
    labeller = ggplot2::as_labeller(event_labels)
  ) +
  ggplot2::scale_color_manual(values = quartile_colors, drop = FALSE) +
  ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72), expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, y_limit), expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(
    title = "B. 72-hour cumulative incidence by admission temperature quartile",
    x = "Hours since ICU admission",
    y = "Cumulative incidence",
    color = "Admission Tmax quartile"
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    axis.title = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "grey93", color = "grey35", linewidth = 0.5),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.border = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 0.5),
    legend.position = "bottom"
  )

figure4 <- panel_a / panel_b +
  patchwork::plot_layout(heights = c(0.47, 0.53))

png_path <- file.path(figure_dir, "figure4_pooled_fine_gray_cifs.png")
pdf_path <- file.path(figure_dir, "figure4_pooled_fine_gray_cifs.pdf")
ggplot2::ggsave(png_path, figure4, width = 10.5, height = 8.2, dpi = 600)
ggplot2::ggsave(pdf_path, figure4, width = 10.5, height = 8.2)

message("Wrote pooled Fine-Gray CIF manuscript figure to ", png_path, " and ", pdf_path)

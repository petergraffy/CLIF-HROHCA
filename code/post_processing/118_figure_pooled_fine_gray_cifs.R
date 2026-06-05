#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(pooled_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

all_site_path <- file.path(pooled_dir, "all_site_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv")
pooled_path <- file.path(pooled_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.csv")

if (!file.exists(all_site_path) || !file.exists(pooled_path)) {
  stop(
    "Pooled CIF inputs are missing. Run site pipelines with the CIF export, then run code/post_processing/115_pool_72h_phenotype_models.R.",
    call. = FALSE
  )
}

all_site <- readr::read_csv(all_site_path, show_col_types = FALSE)
pooled <- readr::read_csv(pooled_path, show_col_types = FALSE)
if (!"stratification" %in% names(all_site)) all_site$stratification <- "Overall"
if (!"stratification" %in% names(pooled)) pooled$stratification <- "Overall"
if (!"stratum" %in% names(all_site)) all_site$stratum <- "Overall"
if (!"stratum" %in% names(pooled)) pooled$stratum <- "Overall"

event_colors <- c(
  "Awake/extubated" = "#0f6b78",
  "Death before awake/extubated" = "#8c3b3b"
)

p <- ggplot2::ggplot() +
  ggplot2::geom_step(
    data = all_site |> dplyr::filter(.data$stratification == "Overall"),
    ggplot2::aes(
      x = .data$time_hours,
      y = .data$cif,
      group = interaction(.data$site_name, .data$event_type)
    ),
    color = "grey70",
    linewidth = 0.7,
    alpha = 0.8
  ) +
  ggplot2::geom_step(
    data = pooled |> dplyr::filter(.data$stratification == "Overall"),
    ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$event_type),
    linewidth = 1.25
  ) +
  ggplot2::scale_color_manual(values = event_colors) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, NA)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72)) +
  ggplot2::labs(
    title = "Pooled 72-hour competing-risk cumulative incidence after OHCA ICU admission",
    subtitle = "Light grey lines are site-specific CIFs; colored lines are n-weighted pooled CIFs.",
    x = "Hours since ICU admission",
    y = "Cumulative incidence",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(color = "grey35"),
    legend.position = "bottom"
  )

output_file <- file.path(figure_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_cif_curves.png")
ggplot2::ggsave(output_file, p, width = 8.5, height = 5, dpi = 300)

quartile_levels <- c("Q1 coolest", "Q2", "Q3", "Q4 hottest")
quartile_colors <- c("Q1 coolest" = "#335c67", "Q2" = "#8f7a2f", "Q3" = "#c06c3e", "Q4 hottest" = "#9b2f37")
all_site_quartile <- all_site |>
  dplyr::filter(.data$stratification == "Admission Tmax quartile") |>
  dplyr::mutate(stratum = factor(.data$stratum, levels = quartile_levels))
pooled_quartile <- pooled |>
  dplyr::filter(.data$stratification == "Admission Tmax quartile") |>
  dplyr::mutate(stratum = factor(.data$stratum, levels = quartile_levels))

if (nrow(pooled_quartile) > 0) {
  pq <- ggplot2::ggplot() +
    ggplot2::geom_step(
      data = all_site_quartile,
      ggplot2::aes(
        x = .data$time_hours,
        y = .data$cif,
        group = interaction(.data$site_name, .data$event_type, .data$stratum)
      ),
      color = "grey72",
      linewidth = 0.55,
      alpha = 0.65
    ) +
    ggplot2::geom_step(
      data = pooled_quartile,
      ggplot2::aes(x = .data$time_hours, y = .data$cif, color = .data$stratum),
      linewidth = 1.15
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$event_type), nrow = 1) +
    ggplot2::scale_color_manual(values = quartile_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, NA)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 72, by = 12), limits = c(0, 72)) +
    ggplot2::labs(
      title = "Pooled 72-hour competing-risk CIFs by admission temperature quartile",
      subtitle = "Quartiles are site-specific admission Tmax quartiles; grey lines are sites and colored lines are n-weighted pooled CIFs.",
      x = "Hours since ICU admission",
      y = "Cumulative incidence",
      color = "Tmax quartile"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey35"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  quartile_output_file <- file.path(figure_dir, "pooled_ohca_icu_competing_risk_awake_extubated_72h_cif_curves_by_tmax_quartile.png")
  ggplot2::ggsave(quartile_output_file, pq, width = 10, height = 5.4, dpi = 300)
  message("Wrote pooled Fine-Gray CIF quartile figure to ", quartile_output_file)
}

message("Wrote pooled Fine-Gray CIF figure to ", output_file)

#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root_dir <- normalizePath(file.path(getwd()), mustWork = TRUE)
export_dir <- file.path(root_dir, "output", "final", "federated_exports")
pooled_dir <- file.path(root_dir, "output", "final", "federated_pooled")
figure_dir <- file.path(pooled_dir, "figures")
dir.create(pooled_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_bind <- function(pattern) {
  files <- list.files(export_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(data.frame())
  dplyr::bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))
}

interpolate_site_curves <- function(curves, weights, curve_label) {
  if (nrow(curves) == 0) return(list(site = data.frame(), pooled = data.frame()))

  site_ranges <- curves |>
    dplyr::group_by(.data$site_name, .data$model) |>
    dplyr::summarise(
      min_temp_c = min(.data$tmax_mean_c, na.rm = TRUE),
      max_temp_c = max(.data$tmax_mean_c, na.rm = TRUE),
      .groups = "drop"
    )

  interpolated <- lapply(split(curves, curves$model), function(model_df) {
    ranges <- site_ranges[site_ranges$model == unique(model_df$model), ]
    grid_min <- max(ranges$min_temp_c, na.rm = TRUE)
    grid_max <- min(ranges$max_temp_c, na.rm = TRUE)
    if (!is.finite(grid_min) || !is.finite(grid_max) || grid_min >= grid_max) return(data.frame())

    temp_grid <- seq(grid_min, grid_max, length.out = 100)
    lapply(split(model_df, list(model_df$site_name, model_df$phenotype), drop = TRUE), function(site_df) {
      site_df <- site_df[order(site_df$tmax_mean_c), ]
      unique_df <- site_df[!duplicated(site_df$tmax_mean_c), ]
      if (nrow(unique_df) < 2) return(data.frame())
      data.frame(
        site_name = unique(site_df$site_name),
        model = unique(site_df$model),
        curve_type = curve_label,
        phenotype = unique(site_df$phenotype),
        tmax_mean_c = temp_grid,
        predicted_probability = stats::approx(
          x = unique_df$tmax_mean_c,
          y = unique_df$predicted_probability,
          xout = temp_grid,
          rule = 1
        )$y
      )
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()

  if (nrow(interpolated) == 0) return(list(site = data.frame(), pooled = data.frame()))

  site_interpolated <- interpolated |>
    dplyr::left_join(
      weights |> dplyr::select("site_name", "model", weight_n = "n"),
      by = c("site_name", "model")
    ) |>
    dplyr::mutate(weight_n = dplyr::if_else(is.na(.data$weight_n), 1, .data$weight_n))

  pooled <- site_interpolated |>
    dplyr::group_by(.data$model, .data$curve_type, .data$phenotype, .data$tmax_mean_c) |>
    dplyr::summarise(
      predicted_probability = stats::weighted.mean(.data$predicted_probability, w = .data$weight_n, na.rm = TRUE),
      k_sites = dplyr::n_distinct(.data$site_name),
      pooled_n = sum(.data$weight_n[!is.na(.data$predicted_probability)], na.rm = TRUE),
      .groups = "drop"
    )

  list(site = site_interpolated, pooled = pooled)
}

plot_multinomial_curve <- function(site_curves, pooled_curves, title, output_file) {
  phenotype_order <- c("regained_consciousness_extubated", "limited_brain_function", "anoxic_brain_injury")
  phenotype_labels <- c(
    regained_consciousness_extubated = "Awake/extubated",
    limited_brain_function = "Limited brain function",
    anoxic_brain_injury = "Death within 72h"
  )

  site_curves$phenotype <- factor(site_curves$phenotype, levels = phenotype_order)
  pooled_curves$phenotype <- factor(pooled_curves$phenotype, levels = phenotype_order)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = site_curves,
      ggplot2::aes(
        x = .data$tmax_mean_c,
        y = .data$predicted_probability,
        group = interaction(.data$site_name, .data$phenotype)
      ),
      color = "grey70",
      linewidth = 0.7,
      alpha = 0.75
    ) +
    ggplot2::geom_line(
      data = pooled_curves,
      ggplot2::aes(x = .data$tmax_mean_c, y = .data$predicted_probability),
      color = "#0f6b78",
      linewidth = 1.3
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$phenotype),
      nrow = 1,
      labeller = ggplot2::as_labeller(phenotype_labels)
    ) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, NA)) +
    ggplot2::labs(
      title = title,
      subtitle = "Light grey lines are site-specific curves; teal line is the n-weighted pooled curve over common temperature support.",
      x = "Admission daily maximum temperature (C)",
      y = "Predicted probability"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey35")
    )

  ggplot2::ggsave(output_file, p, width = 10.5, height = 4.8, dpi = 300)
}

primary_curves <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_temperature_curve[.]csv$")
primary_weights <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_model[.]csv$")
mechanism_curves <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted[.]csv$")
mechanism_weights <- read_bind("^[^_]+_ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted[.]csv$")

primary_pooled <- interpolate_site_curves(primary_curves, primary_weights, "primary")
mechanism_pooled <- interpolate_site_curves(mechanism_curves, mechanism_weights, "mechanism_adjusted")

if (nrow(primary_pooled$pooled) > 0) {
  readr::write_csv(primary_pooled$site, file.path(pooled_dir, "all_site_ohca_icu_72h_phenotype_assignment_temperature_curve_interpolated.csv"))
  readr::write_csv(primary_pooled$pooled, file.path(pooled_dir, "pooled_ohca_icu_72h_phenotype_assignment_temperature_curve.csv"))
  plot_multinomial_curve(
    primary_pooled$site,
    primary_pooled$pooled,
    "Pooled Nonlinear Temperature Curves for 72-hour OHCA ICU Phenotype Assignment",
    file.path(figure_dir, "pooled_ohca_icu_72h_phenotype_assignment_temperature_curves.png")
  )
}

if (nrow(mechanism_pooled$pooled) > 0) {
  readr::write_csv(mechanism_pooled$site, file.path(pooled_dir, "all_site_ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted_interpolated.csv"))
  readr::write_csv(mechanism_pooled$pooled, file.path(pooled_dir, "pooled_ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted.csv"))
  plot_multinomial_curve(
    mechanism_pooled$site,
    mechanism_pooled$pooled,
    "Pooled Mechanism-adjusted Temperature Curves for 72-hour Phenotype Assignment",
    file.path(figure_dir, "pooled_ohca_icu_72h_phenotype_assignment_temperature_curves_mechanism_adjusted.png")
  )
}

fine_gray_curve_exports <- list.files(
  export_dir,
  pattern = "^[^_]+_ohca_icu_competing_risk_awake_extubated_72h_temperature_curve[.]csv$",
  full.names = TRUE
)
if (length(fine_gray_curve_exports) == 0) {
  message(
    "Fine-Gray temperature curve exports were not found. ",
    "Current site exports include Fine-Gray coefficients/vcov, but not prediction curves or spline-basis attributes, ",
    "so exact pooled Fine-Gray curves cannot be reconstructed from these exports alone."
  )
}

message("Wrote pooled 72-hour nonlinear multinomial temperature curves to ", figure_dir)

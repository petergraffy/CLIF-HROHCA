#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(pooled_dir, "pooled_dlnm_ratio_results.csv")
if (!file.exists(input_path)) {
  stop("Missing pooled ratio table. Run code/post_processing/119_export_pooled_ratio_tables.R first.", call. = FALSE)
}

fmt_rr <- function(rr, low, high) sprintf("%.2f (%.2f, %.2f)", rr, low, high)
fmt_i2 <- function(x) ifelse(is.finite(x), paste0(round(x), "%"), "")

pooled <- readr::read_csv(input_path, show_col_types = FALSE)

overall <- pooled |>
  dplyr::filter(
    .data$analysis == "case_crossover_dlnm",
    .data$stratum == "Overall",
    .data$model == "primary_humidity_adjusted",
    .data$contrast == "median-referenced cumulative effect"
  )

strata <- pooled |>
  dplyr::filter(
    .data$analysis == "case_crossover_dlnm",
    .data$stratum != "Overall",
    .data$model == "stratified_humidity_adjusted",
    .data$contrast == "median-referenced cumulative effect"
  )

forest <- dplyr::bind_rows(overall, strata) |>
  dplyr::mutate(
    group = dplyr::case_when(
      .data$stratum == "Overall" ~ "Overall",
      .data$stratum %in% c("Male", "Female") ~ "Sex",
      .data$stratum %in% c("<65", ">=65") ~ "Age",
      .data$stratum %in% c("Black", "Non-Black") ~ "Race",
      .data$stratum %in% c("Hispanic", "Non-Hispanic") ~ "Ethnicity",
      TRUE ~ "Other"
    ),
    label = dplyr::case_when(
      .data$stratum == "<65" ~ "<65 years",
      .data$stratum == ">=65" ~ "≥65 years",
      TRUE ~ .data$stratum
    ),
    row_order = dplyr::case_when(
      .data$stratum == "Overall" ~ 1,
      .data$stratum == "Male" ~ 3,
      .data$stratum == "Female" ~ 4,
      .data$stratum == "<65" ~ 6,
      .data$stratum == ">=65" ~ 7,
      .data$stratum == "Black" ~ 9,
      .data$stratum == "Non-Black" ~ 10,
      .data$stratum == "Hispanic" ~ 12,
      .data$stratum == "Non-Hispanic" ~ 13,
      TRUE ~ 99
    ),
    rr_label = fmt_rr(.data$ratio, .data$ratio_low, .data$ratio_high),
    rr_label_plot = sprintf("%.2f\n(%.2f, %.2f)", .data$ratio, .data$ratio_low, .data$ratio_high),
    i2_label = fmt_i2(.data$i2),
    group = factor(.data$group, levels = c("Overall", "Sex", "Age", "Race", "Ethnicity"))
  ) |>
  dplyr::arrange(.data$row_order)

if (nrow(forest) == 0) {
  stop("Could not find case-crossover humidity-adjusted stratified pooled estimates.", call. = FALSE)
}

group_headers <- data.frame(
  group = c("Sex", "Age", "Race", "Ethnicity"),
  y = c(5.75, 2.75, -0.25, -3.25),
  label = c("Sex", "Age", "Race", "Ethnicity"),
  stringsAsFactors = FALSE
)

forest <- forest |>
  dplyr::mutate(
    label = factor(.data$label, levels = unique(.data$label)),
    section = factor(as.character(.data$group), levels = c("Overall", "Sex", "Age", "Race", "Ethnicity")),
    label_y = .data$ratio_high * 1.08
  )

y_min <- min(0.55, min(forest$ratio_low, na.rm = TRUE) * 0.85)
y_max <- max(3, max(forest$label_y, na.rm = TRUE) * 1.2)

source_table <- forest |>
  dplyr::select(
    group,
    stratum,
    ratio,
    ratio_low,
    ratio_high,
    rr_label,
    k_sites,
    i2,
    i2_label,
    tau2
  )
readr::write_csv(source_table, file.path(figure_dir, "figure2_case_crossover_stratified_forest_source.csv"))

figure2 <- ggplot2::ggplot(forest, ggplot2::aes(x = .data$label, y = .data$ratio)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey35") +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin = .data$ratio_low, ymax = .data$ratio_high),
    linewidth = 0.85,
    color = "grey20"
  ) +
  ggplot2::geom_point(size = 3.1, color = "grey10") +
  ggplot2::geom_text(
    ggplot2::aes(y = .data$label_y, label = .data$rr_label_plot),
    hjust = 0.5,
    vjust = 0,
    color = "grey12",
    size = 2.85,
    lineheight = 0.9,
    show.legend = FALSE
  ) +
  ggplot2::facet_grid(. ~ section, scales = "free_x", space = "free_x") +
  ggplot2::scale_y_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
    labels = c("0.5", "0.75", "1", "1.5", "2", "3"),
    limits = c(y_min, y_max),
    expand = ggplot2::expansion(mult = c(0.03, 0.12))
  ) +
  ggplot2::labs(
    title = "Stratified heat-OHCA associations",
    x = NULL,
    y = "Cumulative relative risk"
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "grey15"),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    strip.background = ggplot2::element_rect(fill = "grey92", color = "grey65", linewidth = 0.5),
    strip.text = ggplot2::element_text(face = "bold", color = "grey12"),
    panel.spacing.x = grid::unit(0.5, "lines"),
    plot.margin = ggplot2::margin(12, 18, 12, 12)
  ) +
  ggplot2::coord_cartesian(clip = "off")

png_path <- file.path(figure_dir, "figure2_case_crossover_stratified_forest.png")
pdf_path <- file.path(figure_dir, "figure2_case_crossover_stratified_forest.pdf")
ggplot2::ggsave(png_path, figure2, width = 8.2, height = 5.1, dpi = 600, device = ragg::agg_png)
ggplot2::ggsave(pdf_path, figure2, width = 8.2, height = 5.1, device = grDevices::cairo_pdf)

message("Wrote Figure 2 to ", png_path, " and ", pdf_path)

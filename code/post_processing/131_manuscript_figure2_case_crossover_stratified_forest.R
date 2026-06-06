#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "ggplot2")
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
      TRUE ~ "Other"
    ),
    label = dplyr::case_when(
      .data$stratum == "<65" ~ "<65 years",
      .data$stratum == ">=65" ~ ">=65 years",
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
      TRUE ~ 99
    ),
    rr_label = fmt_rr(.data$ratio, .data$ratio_low, .data$ratio_high),
    i2_label = fmt_i2(.data$i2),
    group = factor(.data$group, levels = c("Overall", "Sex", "Age", "Race"))
  ) |>
  dplyr::arrange(.data$row_order)

if (nrow(forest) == 0) {
  stop("Could not find case-crossover humidity-adjusted stratified pooled estimates.", call. = FALSE)
}

group_headers <- data.frame(
  group = c("Sex", "Age", "Race"),
  y = c(5.75, 2.75, -0.25),
  label = c("Sex", "Age", "Race"),
  stringsAsFactors = FALSE
)

forest$y <- c(7, 5, 4, 2, 1, -1, -2)[seq_len(nrow(forest))]
label_lookup <- forest[, c("label", "y")]

x_min <- min(0.55, min(forest$ratio_low, na.rm = TRUE) * 0.85)
x_max <- max(4, max(forest$ratio_high, na.rm = TRUE) * 3.6)
x_text <- max(forest$ratio_high, na.rm = TRUE) * 1.25
x_i2 <- max(forest$ratio_high, na.rm = TRUE) * 2.35

colors <- c("Overall" = "#0f6b78", "Sex" = "#335c67", "Age" = "#8f7a2f", "Race" = "#9b2f37")

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

figure2 <- ggplot2::ggplot(forest, ggplot2::aes(x = .data$ratio, y = .data$y, color = .data$group)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey35") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = .data$ratio_low, xmax = .data$ratio_high),
    orientation = "y",
    width = 0.18,
    linewidth = 0.85,
    alpha = 0.95
  ) +
  ggplot2::geom_point(size = 3.1) +
  ggplot2::geom_text(
    ggplot2::aes(x = x_text, label = .data$rr_label),
    hjust = 0,
    color = "grey12",
    size = 3.25,
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = x_i2, label = .data$i2_label),
    hjust = 0,
    color = "grey25",
    size = 3.25,
    show.legend = FALSE
  ) +
  ggplot2::annotate("text", x = x_text, y = 7.65, label = "RR (95% CI)", hjust = 0, fontface = "bold", size = 3.35, color = "grey12") +
  ggplot2::annotate("text", x = x_i2, y = 7.65, label = "I2", hjust = 0, fontface = "bold", size = 3.35, color = "grey12") +
  ggplot2::geom_text(
    data = group_headers,
    ggplot2::aes(x = x_min, y = .data$y, label = .data$label),
    inherit.aes = FALSE,
    hjust = 0,
    fontface = "bold",
    color = "grey12",
    size = 3.35
  ) +
  ggplot2::scale_color_manual(values = colors, guide = "none") +
  ggplot2::scale_y_continuous(
    breaks = label_lookup$y,
    labels = label_lookup$label,
    limits = c(-2.6, 8),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
    labels = c("0.5", "0.75", "1", "1.5", "2", "3"),
    limits = c(x_min, x_max),
    expand = ggplot2::expansion(mult = c(0.015, 0.02))
  ) +
  ggplot2::labs(
    title = "Stratified heat-OHCA associations",
    subtitle = "Time-stratified case-crossover DLNM, median-reference and humidity-adjusted",
    x = "Cumulative relative risk",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "grey30", size = 10.5),
    axis.title.x = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "grey15"),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(12, 22, 12, 12)
  ) +
  ggplot2::coord_cartesian(clip = "off")

png_path <- file.path(figure_dir, "figure2_case_crossover_stratified_forest.png")
pdf_path <- file.path(figure_dir, "figure2_case_crossover_stratified_forest.pdf")
ggplot2::ggsave(png_path, figure2, width = 8.2, height = 5.1, dpi = 600)
ggplot2::ggsave(pdf_path, figure2, width = 8.2, height = 5.1)

message("Wrote Figure 2 to ", png_path, " and ", pdf_path)

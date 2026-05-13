#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    ofiles <- vapply(sys.frames(), function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile, character(1))
    ofiles <- stats::na.omit(ofiles)
    if (length(ofiles) > 0) return(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = TRUE))
    stop("Could not determine script path. Run with Rscript.")
  }
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
input_path <- file.path(repo_root, "output", "final", "federated_pooled", "pooled_dlnm_random_effects_results.csv")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(scales)
  library(svglite)
})

fmt_est <- function(rr, low, high) sprintf("%.2f (%.2f-%.2f)", rr, low, high)
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

raw <- read.csv(input_path, stringsAsFactors = FALSE)

forest <- raw |>
  filter(.data$model == "stratified_humidity_adjusted", .data$reference_type == "median") |>
  mutate(
    domain = case_when(
      .data$stratum %in% c("Male", "Female") ~ "Sex",
      .data$stratum %in% c("<65", ">=65") ~ "Age",
      .data$stratum %in% c("Black", "Non-Black") ~ "Race",
      TRUE ~ "Other"
    ),
    label = case_when(
      .data$stratum == "Male" ~ "Male",
      .data$stratum == "Female" ~ "Female",
      .data$stratum == "<65" ~ "<65 years",
      .data$stratum == ">=65" ~ ">=65 years",
      .data$stratum == "Black" ~ "Black",
      .data$stratum == "Non-Black" ~ "Non-Black",
      TRUE ~ .data$stratum
    ),
    row_order = case_when(
      .data$stratum == "Male" ~ 1,
      .data$stratum == "Female" ~ 2,
      .data$stratum == "<65" ~ 3,
      .data$stratum == ">=65" ~ 4,
      .data$stratum == "Black" ~ 5,
      .data$stratum == "Non-Black" ~ 6,
      TRUE ~ 0
    ),
    label = factor(.data$label, levels = rev(.data$label[order(.data$row_order)])),
    domain = factor(.data$domain, levels = c("Sex", "Age", "Race")),
    est_label = fmt_est(.data$ratio, .data$ci_low, .data$ci_high),
    p_label = fmt_p(.data$p_value),
    i2_label = sprintf("%.0f%%", .data$i2)
  ) |>
  arrange(.data$row_order) |>
  mutate(y = rev(seq_len(n())))

overall_primary <- raw |>
  filter(.data$stratum == "Overall", .data$model == "primary_humidity_adjusted", .data$reference_type == "median") |>
  mutate(est_label = fmt_est(.data$ratio, .data$ci_low, .data$ci_high))

caption_text <- paste(
  "Random-effects pooled site estimates for the stratified humidity-adjusted DLNM.",
  "\nEstimates compare each site's hot-temperature contrast with its median-temperature reference;",
  "\nhorizontal lines show 95% CI. I2 indicates between-site heterogeneity."
)

domain_bands <- forest |>
  group_by(.data$domain) |>
  summarise(y_mid = mean(.data$y), .groups = "drop")

forest_panel <- ggplot(forest, aes(y = .data$y, x = .data$ratio, color = .data$domain)) +
  geom_hline(yintercept = c(2.5, 4.5), linewidth = 0.35, color = "#E5E7EB") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, color = "#6B7280") +
  geom_errorbar(
    aes(xmin = .data$ci_low, xmax = .data$ci_high),
    orientation = "y",
    width = 0.18,
    linewidth = 0.8,
    alpha = 0.9
  ) +
  geom_point(size = 3.0) +
  scale_color_manual(
    values = c("Sex" = "#0B3C5D", "Age" = "#8A4F0F", "Race" = "#3F6B3C"),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = forest$y,
    labels = forest$label,
    limits = c(0.55, 6.45),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.75, 1, 1.25, 1.5, 2),
    labels = label_number(accuracy = 0.01),
    limits = c(0.68, 3.45),
    expand = expansion(mult = c(0.02, 0.01))
  ) +
  annotate("text", x = 0.69, y = domain_bands$y_mid, label = as.character(domain_bands$domain), hjust = 1.28, fontface = "bold", size = 3.4, color = "#111827", family = "Helvetica") +
  labs(
    x = "Cumulative relative risk of OHCA, log scale",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = "Helvetica", base_size = 10.5) +
  theme(
    axis.title.x = element_text(face = "bold", size = 9.8, color = "#111827", margin = margin(t = 8)),
    axis.text.x = element_text(size = 8.8, color = "#374151"),
    axis.text.y = element_text(size = 9.5, color = "#111827"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E5E7EB", linewidth = 0.35),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 4, 4, 38)
  )

table_panel <- ggplot(forest, aes(y = .data$y)) +
  geom_hline(yintercept = c(2.5, 4.5), linewidth = 0.35, color = "#E5E7EB") +
  geom_text(aes(x = 0.00, label = .data$est_label), hjust = 0, size = 3.0, family = "Helvetica", color = "#111827") +
  geom_text(aes(x = 0.69, label = .data$p_label), hjust = 0, size = 3.0, family = "Helvetica", color = "#374151") +
  geom_text(aes(x = 0.91, label = .data$i2_label), hjust = 0, size = 3.0, family = "Helvetica", color = "#374151") +
  annotate("text", x = 0.00, y = 6.55, label = "RR (95% CI)", hjust = 0, size = 3.0, fontface = "bold", family = "Helvetica", color = "#111827") +
  annotate("text", x = 0.69, y = 6.55, label = "p", hjust = 0, size = 3.0, fontface = "bold", family = "Helvetica", color = "#111827") +
  annotate("text", x = 0.91, y = 6.55, label = "I2", hjust = 0, size = 3.0, fontface = "bold", family = "Helvetica", color = "#111827") +
  scale_y_continuous(limits = c(0.55, 6.72), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(0, 1.08), expand = expansion(mult = c(0, 0))) +
  theme_void(base_family = "Helvetica", base_size = 10.5) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 4, 4, 8)
  )

forest_plot <- forest_panel + table_panel +
  patchwork::plot_layout(widths = c(0.66, 0.34)) +
  patchwork::plot_annotation(
    title = "Stratified heat-OHCA associations",
    subtitle = "Pooled stratified DLNM estimates by sex, age group, and race",
    caption = caption_text,
    theme = theme(
      plot.title = element_text(family = "Helvetica", face = "bold", size = 12.4, color = "#111827", margin = margin(b = 4)),
      plot.subtitle = element_text(family = "Helvetica", size = 9.8, color = "#374151", margin = margin(b = 10)),
      plot.caption = element_text(family = "Helvetica", size = 7.5, color = "#4B5563", hjust = 0, lineheight = 1.15, margin = margin(t = 8)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 14, 10, 10)
    )
  )

write.csv(
  forest,
  file.path(output_dir, "figure2_stratified_dlnm_forest_source.csv"),
  row.names = FALSE
)

ggsave(
  file.path(output_dir, "figure2_stratified_dlnm_forest.pdf"),
  forest_plot,
  width = 7.6,
  height = 4.8,
  device = cairo_pdf,
  bg = "white"
)

svglite::svglite(file.path(output_dir, "figure2_stratified_dlnm_forest.svg"), width = 7.6, height = 4.8, bg = "white")
print(forest_plot)
dev.off()

ragg::agg_png(
  file.path(output_dir, "figure2_stratified_dlnm_forest.png"),
  width = 7.6,
  height = 4.8,
  units = "in",
  res = 600,
  background = "white"
)
print(forest_plot)
dev.off()

ragg::agg_tiff(
  file.path(output_dir, "figure2_stratified_dlnm_forest.tiff"),
  width = 7.6,
  height = 4.8,
  units = "in",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(forest_plot)
dev.off()

message("Wrote stratified DLNM forest plot to ", output_dir)

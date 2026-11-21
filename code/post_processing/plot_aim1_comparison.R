#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript code/post_processing/plot_aim1_comparison.R MODEL_DIR OUTPUT_STEM")
panels <- c("OHCA admission count", "Total adult ICU admission count", "OHCA with all-ICU offset")
prefixes <- setNames(c("", "all_icu_", "offset_"), panels)
read_panels <- function(file) {
  bind_rows(lapply(prefixes, function(prefix) read.csv(file.path(args[[1L]], paste0(prefix, file)))),
            .id = "model") |>
    mutate(model = factor(model, levels = panels))
}
site <- read_panels("site_curves.csv")
pooled <- filter(read_panels("pooled_curve.csv"), set == "county_linked")
support <- read_panels("temperature_support.csv")
stopifnot(all(panels %in% pooled$model))

site_labels <- c(UCMC = "UCMC", NU = "Northwestern", RUMC = "Rush")
colors <- c(Pooled = "#172B3A", UCMC = "#674399", Northwestern = "#007B80", Rush = "#A45A00")
site$series <- unname(site_labels[site$site])
support$series <- unname(site_labels[support$site])
stopifnot(!anyNA(site$series), !anyNA(support$series))
line_levels <- c("Pooled", unname(site_labels))
pooled$series <- "Pooled"
lines <- bind_rows(site, pooled) |>
  mutate(series = factor(series, levels = line_levels))

base <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", color = "#172B3A"),
    axis.text = element_text(color = "#465762"),
    panel.spacing.x = grid::unit(1.15, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#DFE6EB", linewidth = 0.4),
    strip.text = element_text(face = "bold", color = "#172B3A", size = 12),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_blank(),
    plot.margin = margin(5, 12, 5, 12)
  )

curve <- ggplot(lines, aes(temp_c, rr, color = series)) +
  annotate("rect", xmin = 30.5, xmax = 33.5, ymin = min(lines$rr, pooled$low) * 0.99,
           ymax = max(lines$rr, pooled$high) * 1.01, fill = "#F8EEDC") +
  geom_ribbon(data = pooled, aes(x = temp_c, ymin = low, ymax = high),
              inherit.aes = FALSE, fill = "#A6BECF", alpha = 0.46) +
  geom_hline(yintercept = 1, color = "#71828E", linewidth = 0.42, linetype = "dashed") +
  geom_vline(xintercept = unique(pooled$reference_c), color = "#71828E", linewidth = 0.38,
             linetype = "dotted") +
  geom_line(aes(linewidth = series), lineend = "round") +
  scale_color_manual(values = colors, breaks = line_levels) +
  scale_linewidth_manual(values = c(Pooled = 1.12, UCMC = 0.61, Northwestern = 0.61, Rush = 0.61),
                         guide = "none") +
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20, 30, 35)) +
  scale_y_log10(breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8), labels = scales::label_number()) +
  facet_grid(cols = vars(model)) +
  labs(x = "Daily maximum hospital-county temperature (°C)", y = "Relative count or rate") +
  base

hot_bins <- seq(25, max(support$temp_c))
warm <- expand.grid(model = levels(support$model), temp_c = hot_bins, series = unname(site_labels)) |>
  left_join(select(support, model, temp_c, series, observations),
            by = c("model", "temp_c", "series")) |>
  mutate(model = factor(model, levels = panels), series = factor(series, levels = line_levels[-1L]),
         observations = coalesce(observations, 0L))
bars <- ggplot(warm, aes(temp_c, observations, fill = series)) +
  annotate("rect", xmin = 30.5, xmax = 33.5, ymin = 0, ymax = max(warm$observations) * 1.05,
           fill = "#F8EEDC") +
  geom_col(position = position_dodge(width = 0.85), width = 0.78) +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_continuous(breaks = intersect(hot_bins, c(25, 27, 29, 31, 32, 33, 35, 37)),
                     expand = expansion(add = 0.5)) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.06))) +
  facet_grid(cols = vars(model)) +
  labs(x = "Temperature bin center (°C)", y = "County-days in model") +
  base +
  theme(strip.text = element_blank())

figure <- (curve / bars) +
  plot_layout(heights = c(1.7, 1)) +
  plot_annotation(
    title = "Temperature and ICU admissions",
    subtitle = "Chicago-area pilot  |  Cumulative association over lag days 0-5  |  Reference 20 °C",
    caption = paste("Blue ribbons show pointwise 95% confidence intervals for the pooled splines. Tan shading marks 31-33 °C.",
                    "Bars count fitted county-days in 1 °C bins. UCMC and Rush share Cook County weather.",
                    "The offset panel estimates OHCA admissions relative to total adult ICU admissions.",
                    sep = "\n"),
    theme = theme(
      plot.title = element_text(size = 25, face = "bold", color = "#172B3A", margin = margin(b = 7)),
      plot.subtitle = element_text(size = 12, color = "#465762"),
      plot.caption = element_text(hjust = 0, size = 10, color = "#465762", lineheight = 1.2),
      plot.margin = margin(17, 24, 14, 19)
    )
  )

ggsave(paste0(args[[2L]], ".png"), figure, width = 14, height = 8.2, dpi = 220, bg = "white")
ggsave(paste0(args[[2L]], ".pdf"), figure, width = 14, height = 8.2, bg = "white")

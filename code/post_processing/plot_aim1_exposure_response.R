#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L) || (length(args) == 3L && args[[3L]] != "all_icu")) {
  stop("Usage: Rscript code/post_processing/plot_aim1_exposure_response.R MODEL_DIR OUTPUT_STEM [all_icu]")
}
all_icu <- length(args) == 3L
prefix <- if (all_icu) "all_icu_" else ""
read_model <- function(file) read.csv(file.path(args[[1L]], paste0(prefix, file)))
site <- read_model("site_curves.csv")
pooled <- subset(read_model("pooled_curve.csv"), set == "county_linked")
support <- read_model("temperature_support.csv")
stopifnot(nrow(pooled) > 0L, nrow(site) > 0L, nrow(support) > 0L)

labels <- c(NU = "Northwestern", RUMC = "Rush")
site_label <- function(x) ifelse(x %in% names(labels), labels[x], x)
site$series <- site_label(site$site)
support$series <- site_label(support$site)
pooled$series <- "Pooled"
series <- c("Pooled", sort(unique(site$series)))
lines <- rbind(site[, c("temp_c", "rr", "series")], pooled[, c("temp_c", "rr", "series")])
lines$series <- factor(lines$series, levels = series)
support$series <- factor(support$series, levels = series[-1L])
pilot_colors <- c(Pooled = "#172B3A", Northwestern = "#007C80", Rush = "#A65B00", UCMC = "#7046A0")
missing_colors <- setdiff(series, names(pilot_colors))
colors <- c(pilot_colors, setNames(scales::hue_pal()(length(series))[seq_along(missing_colors)], missing_colors))

base <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 22, margin = margin(b = 7)),
    plot.subtitle = element_text(color = "#465762", size = 11, margin = margin(b = 12)),
    axis.title = element_text(face = "bold", color = "#172B3A", size = 11),
    axis.text = element_text(color = "#465762"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E0E7EB", linewidth = 0.35),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_blank(),
    plot.margin = margin(10, 22, 7, 15),
    plot.caption = element_text(hjust = 0, size = 9, color = "#465762", lineheight = 1.2)
  )

curve <- ggplot(lines, aes(temp_c, rr, color = series)) +
  annotate("rect", xmin = 30.5, xmax = 33.5, ymin = 0.99 * min(lines$rr, pooled$low),
           ymax = 1.01 * max(lines$rr, pooled$high), fill = "#F8EEDC") +
  geom_ribbon(data = pooled, aes(x = temp_c, ymin = low, ymax = high),
              inherit.aes = FALSE, fill = "#A6BECF", alpha = 0.42) +
  geom_hline(yintercept = 1, color = "#778691", linewidth = 0.45, linetype = "dashed") +
  geom_vline(xintercept = unique(pooled$reference_c), color = "#778691", linewidth = 0.4,
             linetype = "dotted") +
  geom_line(aes(linewidth = series), lineend = "round") +
  scale_color_manual(values = colors, breaks = series) +
  scale_linewidth_manual(values = setNames(c(1.15, rep(0.67, length(series) - 1L)), series),
                         guide = "none") +
  scale_x_continuous(breaks = scales::breaks_width(5)) +
  scale_y_log10(breaks = if (all_icu) scales::breaks_pretty(n = 6) else c(0.125, 0.25, 0.5, 1, 2, 4, 8),
                labels = if (all_icu) scales::label_number(accuracy = 0.01) else scales::label_number()) +
  labs(title = if (all_icu) "Temperature and total adult ICU admissions" else "Temperature and OHCA ICU admissions",
       subtitle = "Chicago-area pilot  |  Cumulative association over lag days 0-5  |  Reference 20 °C",
       x = "Daily maximum hospital-county temperature (°C)",
       y = if (all_icu) "All-ICU admission count ratio" else "OHCA admission count ratio") +
  base

hot_bins <- seq(25, max(33, support$temp_c))
warm <- merge(expand.grid(temp_c = hot_bins, series = series[-1L]),
              support[, c("temp_c", "series", "observations")], all.x = TRUE)
warm$observations[is.na(warm$observations)] <- 0
warm$series <- factor(warm$series, levels = series[-1L])
bars <- ggplot(warm, aes(temp_c, observations, fill = series)) +
  annotate("rect", xmin = 30.5, xmax = 33.5, ymin = -Inf, ymax = Inf, fill = "#F8EEDC") +
  geom_col(position = position_dodge(width = 0.83), width = 0.76) +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_continuous(breaks = hot_bins, expand = expansion(add = 0.5)) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.09))) +
  labs(title = "Observations at warmer temperatures", x = "Temperature bin center (°C)",
       y = "County-days in model",
       caption = paste("Shading marks 31-33 °C. Bars use 1 °C bins centered at each tick.",
                       "UCMC and Rush share Cook County weather. Northwestern can contribute multiple counties per date.",
                       "Blue ribbon shows pointwise 95% confidence intervals for the pooled spline.",
                       if (all_icu) "The vertical scale differs from the OHCA figure.",
                       sep = "\n")) +
  base +
  theme(plot.title = element_text(size = 15), plot.margin = margin(13, 22, 12, 15))

figure <- (curve / bars) + plot_layout(heights = c(1.7, 1))
ggsave(paste0(args[[2L]], ".png"), figure, width = 10, height = 7.4, dpi = 220, bg = "white")
ggsave(paste0(args[[2L]], ".pdf"), figure, width = 10, height = 7.4, bg = "white")

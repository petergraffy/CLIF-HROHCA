#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript code/post_processing/plot_aim1_reference_sensitivity.R POOLED_REFERENCE_CSV OUTPUT_STEM")
}

estimates <- read.csv(args[[1L]])
estimates <- subset(estimates, model == "ohca_count" & lag_min == 0 & lag_max == 5 &
                      hot_c %in% 30:35 & reference_c %in% 15:25)
stopifnot(nrow(estimates) == 66L, !anyDuplicated(estimates[c("hot_c", "reference_c")]),
          all(is.finite(estimates$rr)), all(estimates$low > 0),
          all(estimates$low <= estimates$rr & estimates$rr <= estimates$high))
estimates$hot_temperature <- factor(paste0(estimates$hot_c, " °C"),
                                    levels = paste0(30:35, " °C"))

plot <- ggplot(estimates, aes(reference_c, rr)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "#9DC9CC", alpha = 0.5) +
  geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dashed", color = "#72828D") +
  geom_vline(xintercept = 20, linewidth = 0.45, linetype = "dotted", color = "#516774") +
  geom_line(linewidth = 0.75, color = "#00777D") +
  geom_point(size = 1.45, color = "#00777D") +
  geom_point(data = subset(estimates, reference_c == 20), size = 2.15, color = "#193443") +
  facet_wrap(vars(hot_temperature), ncol = 3) +
  scale_x_continuous(breaks = 15:25, expand = expansion(add = 0.4)) +
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3), labels = scales::label_number()) +
  labs(
    title = "OHCA admission estimates across reference temperatures",
    subtitle = "Three-site pilot  |  Cumulative association over lag days 0-5",
    x = "Reference temperature (°C)", y = "Pooled OHCA admission count ratio",
    caption = paste("Points show estimates at each whole-degree reference. Ribbons show pointwise 95% confidence intervals.",
                    "The dotted vertical line marks the primary 20 °C reference.",
                    sep = "\n")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", color = "#172B3A", size = 21, margin = margin(b = 7)),
    plot.subtitle = element_text(color = "#465762", margin = margin(b = 14)),
    plot.caption = element_text(color = "#465762", hjust = 0, lineheight = 1.16, margin = margin(t = 12)),
    plot.margin = margin(16, 21, 14, 16),
    axis.title = element_text(face = "bold", color = "#172B3A"),
    axis.text = element_text(color = "#465762", size = 9),
    strip.text = element_text(face = "bold", color = "#172B3A", size = 13),
    panel.spacing = grid::unit(1.1, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#DFE7EC", linewidth = 0.4)
  )

ggsave(paste0(args[[2L]], ".png"), plot, width = 11.4, height = 7.4, dpi = 220, bg = "white")
ggsave(paste0(args[[2L]], ".pdf"), plot, width = 11.4, height = 7.4, bg = "white")

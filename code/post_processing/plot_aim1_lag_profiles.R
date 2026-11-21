#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript code/post_processing/plot_aim1_lag_profiles.R RESULTS_DIR OUTPUT_STEM")
}

temperatures <- 30:35
read_lags <- function(file, label) {
  bind_rows(lapply(temperatures, function(temp) {
    model_dir <- file.path(args[[1L]], sprintf("models_%dc", temp))
    result <- read.csv(file.path(model_dir, file)) |>
      filter(model == "ohca_count") |>
      mutate(temp_c = temp, type = label)
    stopifnot(identical(sort(result$lag_max), 0:5), all(is.finite(result$rr)))
    result
  }))
}

daily_label <- "Individual lag day"
cumulative_label <- "Cumulative through lag day"
profile <- bind_rows(
  read_lags("pooled_individual_lags.csv", daily_label),
  read_lags("pooled_partial_lags.csv", cumulative_label)
)
write.csv(profile, paste0(args[[2L]], ".csv"), row.names = FALSE)
profile <- profile |>
  mutate(temperature = factor(paste0(temp_c, " °C"), levels = paste0(temperatures, " °C")),
         type = factor(type, levels = c(daily_label, cumulative_label)))

plot <- ggplot(profile, aes(lag_max, rr, color = type)) +
  geom_hline(yintercept = 1, linewidth = 0.45, linetype = "dashed", color = "#72828D") +
  geom_line(linewidth = 0.65) +
  geom_linerange(aes(ymin = low, ymax = high), linewidth = 0.55) +
  geom_point(size = 1.9) +
  facet_grid(rows = vars(type), cols = vars(temperature), switch = "y") +
  scale_x_continuous(breaks = 0:5, expand = expansion(add = 0.3)) +
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4), labels = scales::label_number()) +
  scale_color_manual(values = setNames(c("#007A80", "#63459A"), c(daily_label, cumulative_label)),
                     guide = "none") +
  labs(
    title = "Temperature and OHCA admissions by lag day",
    subtitle = "Three-site pilot  |  Each temperature compared with 20 °C  |  Full lag 0-5 count model",
    x = "Lag day", y = "Pooled OHCA admission count ratio",
    caption = paste("Points show pooled estimates. Vertical lines show 95% confidence intervals.",
                    "Cumulative estimates combine lag coefficients within each site before pooling.",
                    sep = "\n")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 22, face = "bold", color = "#172B3A", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#465762", margin = margin(b = 13)),
    plot.caption = element_text(hjust = 0, color = "#465762", lineheight = 1.15, margin = margin(t = 10)),
    plot.margin = margin(17, 18, 12, 14),
    axis.title = element_text(face = "bold", color = "#172B3A"),
    axis.text = element_text(color = "#465762"),
    strip.text = element_text(face = "bold", color = "#172B3A"),
    strip.text.y.left = element_text(angle = 90),
    strip.placement = "outside",
    panel.spacing = grid::unit(1, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E1E7EC", linewidth = 0.4)
  )

ggsave(paste0(args[[2L]], ".png"), plot, width = 13.6, height = 6.8, dpi = 220, bg = "white")
ggsave(paste0(args[[2L]], ".pdf"), plot, width = 13.6, height = 6.8, bg = "white")

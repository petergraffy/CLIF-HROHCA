#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript code/post_processing/plot_aim2_cif.R OUTPUT_STEM SITE_DIR ...")
temperatures <- c(20, 30:35)
endpoints <- c(death = "Hospital death", hospice = "Hospice",
               live_discharge = "Other live exits")

site_curves <- bind_rows(lapply(args[-1L], function(site_dir) {
  readr::read_csv(file.path(site_dir, "standardized_cif.csv"), show_col_types = FALSE)
})) |>
  filter(temperature_c %in% temperatures, endpoint %in% names(endpoints))
stopifnot(!anyDuplicated(select(site_curves, site, temperature_c, endpoint, day)))

site_names <- c(UCMC = "UCMC", NU = "Northwestern", RUMC = "Rush")
coverage <- site_curves |>
  distinct(site, temperature_c, in_exposure_range, support_within_1c) |>
  tidyr::complete(site, temperature_c = temperatures,
                  fill = list(in_exposure_range = FALSE, support_within_1c = 0L)) |>
  arrange(site, temperature_c) |>
  group_by(site) |>
  summarise(n_supported = sum(cumall(in_exposure_range & support_within_1c >= 10)),
            .groups = "drop") |>
  filter(n_supported >= 2L)
cutoffs <- sort(unique(coverage$n_supported))
if (!length(cutoffs)) stop("No site supports 20 C and 30 C.")
curves <- bind_rows(lapply(cutoffs, function(cutoff) {
  included <- coverage$site[coverage$n_supported >= cutoff]
  shown <- temperatures[seq_len(cutoff)]
  d <- filter(site_curves, site %in% included, temperature_c %in% shown)
  stopifnot(nrow(d) == length(included) * length(shown) * length(endpoints) * 29L)
  names_shown <- ifelse(included %in% names(site_names), site_names[included], included)
  # Keep the sites and their weights fixed across temperatures in each column.
  d |>
    group_by(temperature_c, endpoint, day) |>
    summarise(probability = weighted.mean(probability, n), n = sum(n),
              n_sites = n_distinct(site), sites = paste(sort(site), collapse = "+"), .groups = "drop") |>
    mutate(max_hot_c = max(shown), cohort = paste0(paste(names_shown, collapse = " + "),
                                                    "\n", format(n, big.mark = ","), " hospitalizations"))
}))
write.csv(curves, paste0(args[[1L]], ".csv"), row.names = FALSE)
members <- curves |>
  distinct(cohort, max_hot_c, n_sites, sites) |>
  mutate(site = strsplit(sites, "+", fixed = TRUE)) |>
  tidyr::unnest(site)
observed_risk <- bind_rows(lapply(args[-1L], function(site_dir) {
  readr::read_csv(file.path(site_dir, "aalen_johansen.csv"), show_col_types = FALSE) |>
    distinct(site, day, at_risk)
}))
stopifnot(!anyDuplicated(select(observed_risk, site, day)))
risk <- inner_join(members, observed_risk, by = "site", relationship = "many-to-many") |>
  group_by(cohort, max_hot_c, sites, day) |>
  summarise(at_risk = sum(at_risk), .groups = "drop")
write.csv(risk, paste0(args[[1L]], "_risk_table.csv"), row.names = FALSE)

day28 <- filter(curves, day == 28)
reference <- filter(day28, temperature_c == 20) |>
  transmute(cohort, endpoint, reference_probability = probability)
differences <- filter(day28, temperature_c != 20) |>
  inner_join(reference, by = c("cohort", "endpoint")) |>
  transmute(cohort, max_hot_c, sites, n, endpoint, reference_c = 20, hot_c = temperature_c,
            reference_probability, hot_probability = probability,
            difference_pp = 100 * (probability - reference_probability))
write.csv(differences, paste0(args[[1L]], "_day28_risk_differences.csv"), row.names = FALSE)
temperature_labels <- paste0(temperatures, " °C")
curves <- curves |>
  mutate(event = factor(endpoints[endpoint], levels = unname(endpoints)),
         cohort = factor(cohort, levels = unique(cohort)),
         temperature = factor(paste0(temperature_c, " °C"), levels = temperature_labels))

colors <- setNames(c("#253644", "#147F89", "#46876C", "#8A8240", "#AC7139", "#B75544", "#922B53"),
                   temperature_labels)
types <- setNames(c("dashed", rep("solid", length(temperatures) - 1L)), temperature_labels)
widths <- setNames(c(0.9, rep(0.65, length(temperatures) - 2L), 1.0), temperature_labels)

figure <- ggplot(curves, aes(day, probability, color = temperature,
                            linetype = temperature, linewidth = temperature)) +
  geom_line(lineend = "round") +
  facet_grid(rows = vars(event), cols = vars(cohort), scales = "free_y", switch = "y") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = types) +
  scale_linewidth_manual(values = widths, guide = "none") +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28), expand = expansion(add = c(0, 0.5))) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "28-day cumulative incidence of hospital outcomes",
    subtitle = "Adjusted predictions from cause-specific Cox models",
    x = NULL, y = "Adjusted cumulative incidence",
    color = "Preceding-day Tmax", linetype = "Preceding-day Tmax"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", color = "#172B3A", size = 22, margin = margin(b = 7)),
    plot.subtitle = element_text(color = "#465762", margin = margin(b = 12)),
    plot.caption = element_text(color = "#465762", hjust = 0, lineheight = 1.15, margin = margin(t = 12)),
    plot.margin = margin(16, 21, 14, 16),
    axis.title = element_text(face = "bold", color = "#172B3A"),
    axis.text = element_text(color = "#465762", size = 10),
    axis.text.x = element_blank(),
    strip.text = element_text(face = "bold", color = "#172B3A", size = 13),
    strip.text.y.left = element_text(angle = 90),
    strip.placement = "outside",
    panel.spacing = grid::unit(1.5, "lines"),
    legend.position = "top",
    legend.justification = "left",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#DFE7EC", linewidth = 0.4)
  )

axis_days <- c(0, 7, 14, 21, 28)
risk_panel <- risk |>
  filter(day %in% axis_days) |>
  mutate(cohort = factor(cohort, levels = levels(curves$cohort))) |>
  ggplot(aes(day, 1, label = scales::label_comma()(at_risk))) +
  geom_text(color = "#172B3A", size = 3.4) +
  facet_grid(cols = vars(cohort)) +
  scale_x_continuous(breaks = axis_days, expand = expansion(add = c(0, 0.5))) +
  coord_cartesian(ylim = c(0.7, 1.3), clip = "off") +
  labs(x = "Days since first ICU admission", y = "Observed\nat risk") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(), strip.text = element_blank(),
    panel.spacing = grid::unit(1.5, "lines"), axis.text.y = element_blank(),
    axis.text.x = element_text(color = "#465762", size = 10),
    axis.title = element_text(face = "bold", color = "#172B3A"),
    axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
    plot.margin = margin(3, 21, 3, 16)
  )
figure <- (figure / risk_panel) +
  plot_layout(heights = c(1, 0.13)) +
  plot_annotation(
    caption = paste("Curves average adjusted probabilities using cohort-size weights. Each column uses a fixed set of sites.",
                    "Outcome rows use different vertical scales. The risk counts are observed and apply to the entire column.",
                    "The curves and the separately reported day-28 risk differences are point estimates.",
                    sep = "\n"),
    theme = theme(plot.caption = element_text(hjust = 0, color = "#465762", margin = margin(l = 16, b = 10)))
  )

ggsave(paste0(args[[1L]], ".png"), figure, width = 11.4, height = 8.7, dpi = 220, bg = "white")
ggsave(paste0(args[[1L]], ".pdf"), figure, width = 11.4, height = 8.7, bg = "white")

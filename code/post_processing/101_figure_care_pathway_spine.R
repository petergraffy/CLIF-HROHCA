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

repo_root <- normalizePath(file.path(dirname(get_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)

input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(readr)
  library(ragg)
  library(scales)
  library(stringr)
  library(svglite)
  library(tidyr)
})

pathways <- readr::read_csv(file.path(input_dir, "pooled_care_pathway_summary.csv"), show_col_types = FALSE) |>
  mutate(n = as.numeric(.data$n), pct = as.numeric(.data$pct)) |>
  arrange(desc(.data$n))

n_total <- sum(pathways$n, na.rm = TRUE)
site_count <- if (file.exists(file.path(input_dir, "care_pathway_export_availability.csv"))) {
  nrow(readr::read_csv(file.path(input_dir, "care_pathway_export_availability.csv"), show_col_types = FALSE))
} else {
  NA_integer_
}

nice_location <- function(x) {
  dplyr::recode(
    x,
    ed = "ED",
    icu = "ICU",
    procedural = "Procedural",
    ward = "Ward",
    stepdown = "Stepdown",
    other = "Other",
    radiology = "Other",
    rehab = "Other",
    .default = str_to_title(x)
  )
}

route_text <- function(pathway) {
  paste(nice_location(str_split(pathway, " -> ", simplify = FALSE)[[1]]), collapse = " -> ")
}

top_routes <- pathways |>
  slice_head(n = 5) |>
  mutate(
    route = vapply(.data$care_pathway, route_text, character(1)),
    lane = c(0, -1.45, 1.25, 2.15, -2.15),
    line_width = rescale(sqrt(.data$n), to = c(1.55, 8.6)),
    alpha = rescale(.data$n, to = c(0.68, 0.94)),
    label = sprintf("%s  |  %s (%.1f%%)", .data$route, comma(.data$n), .data$pct)
  )

n_top_routes <- sum(top_routes$n, na.rm = TRUE)
pct_top_routes <- 100 * n_top_routes / n_total

node_x <- c(
  Admitted = 0,
  ED = 1.25,
  ICU = 4.2,
  Procedural = 2.45,
  Ward = 2.45,
  Stepdown = 2.45,
  Other = 2.45
)

route_points <- top_routes |>
  mutate(tokens = str_split(.data$care_pathway, " -> ", simplify = FALSE)) |>
  select("care_pathway", "route", "n", "pct", "lane", "line_width", "alpha", "label", "tokens") |>
  unnest_longer("tokens", values_to = "raw_location", indices_to = "token_order") |>
  mutate(location = nice_location(.data$raw_location)) |>
  group_by(.data$care_pathway) |>
  mutate(
    token_order = as.integer(.data$token_order),
    n_tokens = n(),
    x = case_when(
      .data$token_order == 1L & .data$location == "ICU" ~ 3.0,
      .data$token_order == 1L ~ node_x[.data$location],
      .data$location == "ICU" ~ node_x["ICU"],
      TRUE ~ 1.25 + 0.85 * (.data$token_order - 1L)
    ),
    y = .data$lane
  ) |>
  ungroup()

start_points <- top_routes |>
  transmute(
    care_pathway,
    route,
    n,
    pct,
    lane,
    line_width,
    alpha,
    label,
    raw_location = "admitted",
    token_order = 0L,
    location = "Admitted",
    x = node_x["Admitted"],
    y = .data$lane
  )

route_points <- bind_rows(start_points, route_points) |>
  arrange(.data$care_pathway, .data$token_order)

segments <- route_points |>
  group_by(.data$care_pathway) |>
  mutate(
    xend = lead(.data$x),
    yend = lead(.data$y),
    next_location = lead(.data$location)
  ) |>
  ungroup() |>
  filter(!is.na(.data$xend)) |>
  mutate(
    color_location = case_when(
      .data$location == "Admitted" & .data$next_location == "ED" ~ "ED",
      .data$location == "Admitted" & .data$next_location == "ICU" ~ "ICU",
      .data$location == "Admitted" ~ .data$next_location,
      TRUE ~ .data$location
    )
  )

nodes <- route_points |>
  group_by(.data$location, .data$x, .data$y) |>
  summarise(n = max(.data$n), .groups = "drop") |>
  mutate(
    fill_location = if_else(.data$location %in% c("Admitted", "ED", "ICU", "Procedural", "Ward", "Stepdown"), .data$location, "Other"),
    label = .data$location
  )

label_df <- top_routes |>
  mutate(
    x = 4.65,
    y = .data$lane,
    label = sprintf(
      "<span style='color:#607D8B'>%02d.</span> <b>%s</b><br><span style='color:#607D8B'>%s admissions (%.1f%%)</span>",
      row_number(), .data$route, comma(.data$n), .data$pct
    )
  )

palette <- c(
  Admitted = "#17212B",
  ED = "#1B6CA8",
  ICU = "#7A3E9D",
  Procedural = "#C85A2E",
  Ward = "#3D8B5A",
  Stepdown = "#C79A2B",
  Other = "#6B7280"
)

spine_plot <- ggplot() +
  geom_curve(
    data = segments,
    aes(
      x = .data$x,
      y = .data$y,
      xend = .data$xend,
      yend = .data$yend,
      color = .data$color_location,
      linewidth = .data$line_width,
      alpha = .data$alpha
    ),
    curvature = 0.16,
    lineend = "round"
  ) +
  geom_point(
    data = nodes,
    aes(x = .data$x, y = .data$y, fill = .data$fill_location),
    shape = 21,
    size = 6.2,
    color = "white",
    stroke = 1.2
  ) +
  geom_text(
    data = nodes |> filter(.data$location %in% c("ED", "ICU")),
    aes(x = .data$x, y = .data$y, label = .data$label),
    color = "white",
    fontface = "bold",
    size = 2.9
  ) +
  geom_richtext(
    data = label_df,
    aes(x = .data$x, y = .data$y, label = .data$label),
    hjust = 0,
    vjust = 0.5,
    fill = NA,
    label.color = NA,
    size = 3.35,
    lineheight = 1.0,
    color = "#17212B"
  ) +
  annotate("text", x = 0, y = 2.95, label = "Admitted", fontface = "bold", size = 3.4, color = "#17212B") +
  annotate("text", x = 1.25, y = 2.95, label = "First location", fontface = "bold", size = 3.4, color = "#17212B") +
  annotate("text", x = 2.45, y = 2.95, label = "Intermediate", fontface = "bold", size = 3.4, color = "#17212B") +
  annotate("text", x = 4.2, y = 2.95, label = "ICU", fontface = "bold", size = 3.4, color = "#17212B") +
  scale_color_manual(values = palette, guide = "none") +
  scale_fill_manual(values = palette, guide = "none") +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  coord_cartesian(xlim = c(-0.25, 7.25), ylim = c(-2.65, 3.12), clip = "off") +
  labs(
    title = "Dominant Care Pathways to First ICU Entry",
    subtitle = sprintf("Top 5 pathways across %d CLIF sites account for %s of %s classified admissions (%.1f%%)", site_count, comma(n_top_routes), comma(n_total), pct_top_routes),
    caption = "Line width scales with pathway count. The central route is ED -> ICU; adjacent routes branch through ICU, procedural, or ward locations before first ICU entry."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15.5, color = "#17212B", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.6, color = "#455A64", margin = margin(b = 12)),
    plot.caption = element_text(size = 8.8, color = "#607D8B", hjust = 0, margin = margin(t = 8)),
    plot.margin = margin(12, 28, 10, 12)
  )

save_figure <- function(plot, basename, width = 10.6, height = 5.3) {
  pdf(file.path(output_dir, paste0(basename, ".pdf")), width = width, height = height, onefile = FALSE)
  print(plot)
  dev.off()
  svglite::svglite(file.path(output_dir, paste0(basename, ".svg")), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(output_dir, paste0(basename, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(output_dir, paste0(basename, ".tiff")), width = width, height = height, units = "in", res = 600, background = "white", compression = "lzw")
  print(plot)
  dev.off()
}

save_figure(spine_plot, "figure_care_pathway_top5_spine")
readr::write_csv(top_routes, file.path(output_dir, "figure_care_pathway_top5_spine_source.csv"))

message("Wrote care pathway spine plot to: ", output_dir)

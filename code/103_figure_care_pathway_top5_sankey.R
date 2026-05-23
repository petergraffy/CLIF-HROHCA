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
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)

input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggalluvial)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(ragg)
  library(scales)
  library(stringr)
  library(svglite)
})

nice_location <- function(x) {
  dplyr::recode(
    x,
    ed = "ED",
    icu = "ICU",
    procedural = "Procedural",
    ward = "Ward",
    stepdown = "Stepdown",
    other = "Other",
    radiology = "Radiology",
    rehab = "Rehab",
    dialysis = "Dialysis",
    .default = str_to_title(x)
  )
}

pathways <- readr::read_csv(file.path(input_dir, "pooled_care_pathway_summary.csv"), show_col_types = FALSE) |>
  mutate(n = as.numeric(.data$n), pct = as.numeric(.data$pct)) |>
  arrange(desc(.data$n))

n_total <- sum(pathways$n, na.rm = TRUE)
site_count <- if (file.exists(file.path(input_dir, "care_pathway_export_availability.csv"))) {
  nrow(readr::read_csv(file.path(input_dir, "care_pathway_export_availability.csv"), show_col_types = FALSE))
} else {
  NA_integer_
}

route_tokens <- str_split(pathways$care_pathway, " -> ", simplify = FALSE)

top5 <- pathways |>
  slice_head(n = 5) |>
  mutate(
    path_id = row_number(),
    route_tokens = route_tokens[seq_len(5)],
    route_pretty = vapply(.data$route_tokens, function(x) paste(nice_location(x), collapse = " -> "), character(1)),
    route_label = sprintf("%s\n%s (%.1f%%)", .data$route_pretty, comma(.data$n), .data$pct),
    flow_group = factor(
      sprintf("%02d. %s", .data$path_id, .data$route_pretty),
      levels = sprintf("%02d. %s", .data$path_id, .data$route_pretty)
    ),
    first_location = vapply(.data$route_tokens, function(x) nice_location(x[[1]]), character(1)),
    next_step = vapply(.data$route_tokens, function(x) {
      if (length(x) == 1 && x[[1]] == "icu") return("Already in ICU")
      if (length(x) == 2 && x[[2]] == "icu") return("Direct to ICU")
      nice_location(x[[2]])
    }, character(1)),
    icu_entry = "ICU"
  )

n_top5 <- sum(top5$n, na.rm = TRUE)
top5_pct <- 100 * n_top5 / n_total
axis_cols <- c("first_location", "next_step", "icu_entry")
axis_labels <- c("First location", "Next location", "First ICU entry")

flow_palette <- c(
  "01. ED -> ICU" = "#1B6CA8",
  "02. ICU" = "#7A3E9D",
  "03. ED -> Procedural -> ICU" = "#2A83B9",
  "04. Procedural -> ICU" = "#C85A2E",
  "05. ED -> Ward -> ICU" = "#3D8B5A"
)

node_palette <- c(
  ED = "#1B6CA8",
  ICU = "#7A3E9D",
  Procedural = "#C85A2E",
  Ward = "#3D8B5A",
  `Direct to ICU` = "#B6C8D8",
  `Already in ICU` = "#CBB2DB"
)

alluvial_data <- top5 |>
  ggalluvial::to_lodes_form(
    axes = axis_cols,
    id = "path_id",
    key = "axis",
    value = "stratum"
  ) |>
  mutate(
    axis = factor(.data$axis, levels = axis_cols, labels = axis_labels),
    stratum = as.character(.data$stratum),
    stratum_fill = if_else(.data$stratum %in% names(node_palette), .data$stratum, "Other")
  )

route_table <- top5 |>
  mutate(
    pathway_label = str_replace_all(.data$route_pretty, " -> ", " to "),
    display = sprintf("n = %s (%.1f%%)", comma(.data$n), .data$pct),
    pathway_label = forcats::fct_reorder(.data$pathway_label, .data$n)
  )

route_table_plot <- ggplot(route_table, aes(x = .data$pathway_label, y = .data$n, fill = .data$flow_group)) +
  geom_col(width = 0.58, alpha = 0.96) +
  geom_text(
    aes(label = .data$display),
    hjust = -0.04,
    size = 3.55,
    color = "#263238",
    fontface = "plain"
  ) +
  scale_fill_manual(values = flow_palette, guide = "none") +
  scale_y_continuous(labels = label_comma(), limits = c(0, max(route_table$n) * 1.42), expand = expansion(mult = c(0, 0.015))) +
  coord_flip(clip = "off") +
  labs(
    title = "Pathway Counts",
    x = NULL,
    y = "Admissions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#17212B"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10.2, color = "#263238", lineheight = 0.95),
    axis.text.x = element_text(size = 9.5, color = "#607D8B"),
    axis.title.x = element_text(size = 10, color = "#455A64", margin = margin(t = 6)),
    plot.margin = margin(8, 34, 8, 8)
  )

sankey_plot <- ggplot(
  alluvial_data,
  aes(
    x = .data$axis,
    stratum = .data$stratum,
    alluvium = .data$path_id,
    y = .data$n
  )
) +
  ggalluvial::geom_alluvium(
    aes(fill = .data$flow_group),
    width = 0.18,
    alpha = 0.78,
    knot.pos = 0.42,
    curve_type = "quintic",
    color = "white",
    linewidth = 0.18
  ) +
  ggalluvial::geom_stratum(
    aes(fill = .data$stratum_fill),
    width = 0.25,
    color = "white",
    linewidth = 0.75,
    alpha = 0.98
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(ifelse(stratum %in% c("ED", "ICU", "Direct to ICU", "Already in ICU"), stratum, ""))),
    color = "white",
    fontface = "bold",
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = c(flow_palette, node_palette), guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.28, 0.18))) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.015, 0.04))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Dominant Care Pathways to First ICU Entry",
    subtitle = sprintf("Top 5 pathways across %d CLIF sites account for %s of %s classified admissions (%.1f%%)", site_count, comma(n_top5), comma(n_total), top5_pct),
    x = NULL,
    y = "Admissions",
    caption = "Flows are color-coded by pathway. Direct ICU routes are represented with direct-to-ICU or already-in-ICU next-location strata."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, color = "#17212B"),
    plot.subtitle = element_text(size = 11.5, color = "#455A64", margin = margin(b = 10)),
    plot.caption = element_text(size = 9.5, color = "#607D8B", hjust = 0, margin = margin(t = 8)),
    axis.text.x = element_text(face = "bold", color = "#263238", size = 11.5),
    axis.text.y = element_text(color = "#607D8B", size = 9.5),
    axis.title.y = element_text(color = "#455A64", size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 14)
  )

figure <- sankey_plot + route_table_plot + plot_layout(widths = c(1.75, 1))

save_figure <- function(plot, basename, width = 13.6, height = 7.2) {
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

save_figure(figure, "figure_care_pathway_top5_sankey")
readr::write_csv(top5, file.path(output_dir, "figure_care_pathway_top5_sankey_source.csv"))

message("Wrote top-5 care-pathway Sankey to: ", output_dir)

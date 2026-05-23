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
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(ragg)
  library(sf)
  library(stringr)
  library(svglite)
  library(tidyr)
  library(tigris)
})

options(tigris_use_cache = TRUE)

conus_excluded_state_fips <- c("02", "15", "60", "66", "69", "72", "78")

hospital_geo <- readr::read_csv(
  file.path(repo_root, "reference", "clif_hospital_geography.csv"),
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
) |>
  mutate(
    site_name = str_squish(.data$site_name),
    hospital_id_name = str_squish(.data$hospital_id_name),
    hospital_full_name = str_squish(.data$hospital_full_name),
    hospital_county_fips = sprintf("%05s", .data$hospital_county_fips),
    adjacent_county_fips = replace_na(.data$adjacent_county_fips, "")
  )

hospital_counties <- hospital_geo |>
  transmute(
    site_name,
    hospital_id_name,
    hospital_full_name,
    geoid = .data$hospital_county_fips,
    attribution = "Hospital county"
  )

adjacent_counties <- hospital_geo |>
  select("site_name", "hospital_id_name", "hospital_full_name", "adjacent_county_fips") |>
  separate_rows(.data$adjacent_county_fips, sep = "[|]") |>
  mutate(adjacent_county_fips = str_squish(.data$adjacent_county_fips)) |>
  filter(.data$adjacent_county_fips != "") |>
  transmute(
    site_name,
    hospital_id_name,
    hospital_full_name,
    geoid = sprintf("%05s", .data$adjacent_county_fips),
    attribution = "Adjacent county"
  )

county_attribution_long <- bind_rows(hospital_counties, adjacent_counties) |>
  distinct() |>
  mutate(state_fips = substr(.data$geoid, 1, 2)) |>
  filter(!.data$state_fips %in% conus_excluded_state_fips)

county_attribution <- county_attribution_long |>
  group_by(.data$geoid) |>
  summarise(
    n_sites = n_distinct(.data$site_name),
    n_hospitals = n_distinct(.data$hospital_id_name),
    sites = paste(sort(unique(.data$site_name)), collapse = "; "),
    hospitals = paste(sort(unique(.data$hospital_id_name)), collapse = "; "),
    has_hospital_county = any(.data$attribution == "Hospital county"),
    map_group = if_else(.data$n_sites > 1, "Multiple sites", first(sort(unique(.data$site_name)))),
    .groups = "drop"
  )

county_shapes <- tigris::counties(cb = TRUE, year = 2023, class = "sf") |>
  transmute(
    geoid = as.character(.data$GEOID),
    county_label = paste0(.data$NAME, ", ", .data$STUSPS),
    state_fips = as.character(.data$STATEFP),
    geometry = .data$geometry
  ) |>
  filter(!.data$state_fips %in% conus_excluded_state_fips) |>
  st_transform(5070)

map_data <- county_shapes |>
  left_join(county_attribution, by = "geoid") |>
  mutate(
    map_group = replace_na(.data$map_group, "Not attributed"),
    map_group = factor(
      .data$map_group,
      levels = c(
        "Emory", "JHU", "Michigan", "NU", "OHSU", "Penn", "RUMC", "UCMC", "UCSF", "UMN",
        "Multiple sites", "Not attributed"
      )
    )
  )

hospital_county_points <- county_shapes |>
  inner_join(
    hospital_geo |>
      transmute(
        site_name,
        hospital_id_name,
        hospital_full_name,
        geoid = .data$hospital_county_fips
      ) |>
      distinct(),
    by = "geoid"
  ) |>
  st_point_on_surface()

site_palette <- c(
  "Emory" = "#8B1E3F",
  "JHU" = "#1769AA",
  "Michigan" = "#FFCB05",
  "NU" = "#4E2A84",
  "OHSU" = "#007C89",
  "Penn" = "#011F5B",
  "RUMC" = "#00A3E0",
  "UCMC" = "#800000",
  "UCSF" = "#052049",
  "UMN" = "#7A0019",
  "Multiple sites" = "#111827",
  "Not attributed" = "#E5E7EB"
)

readr::write_csv(
  county_attribution_long,
  file.path(output_dir, "figure_site_attributed_counties_long_source.csv")
)
readr::write_csv(
  county_attribution |>
    arrange(desc(.data$n_sites), desc(.data$n_hospitals), .data$geoid),
  file.path(output_dir, "figure_site_attributed_counties_source.csv")
)

summary_table <- county_attribution |>
  summarise(
    attributed_counties = n(),
    counties_with_multiple_sites = sum(.data$n_sites > 1),
    counties_with_multiple_hospitals = sum(.data$n_hospitals > 1),
    hospital_counties = sum(.data$has_hospital_county),
    adjacent_only_counties = sum(!.data$has_hospital_county),
    .groups = "drop"
  )
readr::write_csv(summary_table, file.path(table_dir, "supplement_table_site_attributed_county_map_summary.csv"))

theme_map <- function() {
  theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0, color = "#111827"),
      plot.subtitle = element_text(size = 11, hjust = 0, color = "#374151", margin = margin(b = 8)),
      plot.caption = element_text(size = 8.5, hjust = 0, color = "#4B5563", margin = margin(t = 8)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11, color = "#111827"),
      legend.text = element_text(size = 10, color = "#111827"),
      legend.key.height = grid::unit(0.38, "cm"),
      legend.key.width = grid::unit(0.38, "cm"),
      plot.margin = margin(12, 18, 12, 18)
    )
}

county_map <- ggplot() +
  geom_sf(
    data = filter(map_data, .data$map_group == "Not attributed"),
    aes(fill = .data$map_group),
    color = NA
  ) +
  geom_sf(
    data = filter(map_data, .data$map_group != "Not attributed"),
    aes(fill = .data$map_group),
    color = "white",
    linewidth = 0.12
  ) +
  geom_sf(
    data = hospital_county_points,
    shape = 21,
    size = 1.7,
    fill = "white",
    color = "#111827",
    linewidth = 0.35,
    alpha = 0.92
  ) +
  scale_fill_manual(
    values = site_palette,
    name = "County attribution",
    drop = FALSE,
    breaks = names(site_palette)[names(site_palette) != "Not attributed"]
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "CLIF Heat-Related OHCA County Attribution",
    subtitle = "Hospital counties and adjacent counties attributed to each contributing site; black indicates counties attributed to multiple sites.",
    caption = "County attribution uses reference/clif_hospital_geography.csv. Points mark hospital counties. Projection: NAD83 / Conus Albers (EPSG:5070). Alaska, Hawaii, and territories excluded."
  ) +
  theme_map()

save_plot <- function(plot, basename, width = 12, height = 7.6) {
  ggsave(file.path(output_dir, paste0(basename, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  svglite::svglite(file.path(output_dir, paste0(basename, ".svg")), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(output_dir, paste0(basename, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(output_dir, paste0(basename, ".tiff")), width = width, height = height, units = "in", res = 600, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

save_plot(county_map, "figure_site_attributed_county_map")

print(summary_table)
message("Wrote site-attributed county map to ", output_dir)

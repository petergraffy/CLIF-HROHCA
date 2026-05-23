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
weather_dir <- file.path(repo_root, "output", "final", "county_weather_maps")
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(ragg)
  library(sf)
  library(stringr)
  library(svglite)
  library(tidyr)
  library(tigris)
  library(viridis)
})

options(tigris_use_cache = TRUE)

conus_excluded_state_fips <- c("02", "15", "60", "66", "69", "72", "78")

county_weather_path <- file.path(weather_dir, "county_mean_tmax_humidity_2018_2024.csv")
if (!file.exists(county_weather_path)) {
  stop("Missing county weather summary: ", county_weather_path, ". Run code/92_map_county_mean_weather.R first.")
}

county_weather <- readr::read_csv(county_weather_path, show_col_types = FALSE) |>
  mutate(
    geoid = sprintf("%05s", as.character(.data$geoid)),
    state_fips = sprintf("%02s", as.character(.data$state_fips))
  )

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
  separate_rows("adjacent_county_fips", sep = "[|]") |>
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
  left_join(county_weather, by = c("geoid", "state_fips")) |>
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
  "UCMC" = "#F97316",
  "UCSF" = "#052049",
  "UMN" = "#7A0019",
  "Multiple sites" = "#111827",
  "Not attributed" = "#E5E7EB"
)

legend_fill_levels <- factor(
  names(site_palette)[names(site_palette) != "Not attributed"],
  levels = levels(map_data$map_group)
)

legend_seed <- st_sf(
  map_group = legend_fill_levels,
  geometry = st_sfc(
    rep(list(st_geometry(map_data)[[1]]), length(legend_fill_levels)),
    crs = st_crs(map_data)
  )
)

theme_map <- function(base_size = 16) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0, color = "#111827", margin = margin(b = 6)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 16, color = "#111827"),
      legend.text = element_text(size = 15, color = "#111827"),
      legend.key.height = grid::unit(0.48, "cm"),
      legend.key.width = grid::unit(0.48, "cm"),
      plot.margin = margin(8, 14, 8, 14)
    )
}

attribution_plot <- ggplot() +
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
    size = 2.2,
    fill = "white",
    color = "#111827",
    linewidth = 0.45,
    alpha = 0.95
  ) +
  geom_sf(
    data = legend_seed,
    aes(fill = .data$map_group),
    color = NA,
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = site_palette,
    name = "County attribution",
    limits = names(site_palette),
    drop = FALSE,
    breaks = names(site_palette)[names(site_palette) != "Not attributed"],
    guide = guide_legend(override.aes = list(alpha = 1, color = NA))
  ) +
  coord_sf(datum = NA) +
  labs(title = "A. CLIF county catchment (2018-2024)") +
  theme_map()

make_weather_plot <- function(fill_col, title, legend_title, palette) {
  ggplot(map_data) +
    geom_sf(aes(fill = .data[[fill_col]]), color = NA) +
    scale_fill_viridis_c(
      option = palette,
      name = legend_title,
      na.value = "#E5E7EB",
      guide = guide_colorbar(barheight = grid::unit(3.6, "in"), barwidth = grid::unit(0.24, "in"))
    ) +
    coord_sf(datum = NA) +
    labs(title = title) +
    theme_map()
}

tmax_plot <- make_weather_plot(
  "mean_tmax_c",
  "B. Mean daily maximum temperature by county (Daymet, 2018-2024)",
  "Mean Tmax (°C)",
  "magma"
)

rmax_plot <- make_weather_plot(
  "mean_rmax_pct",
  "C. Mean daily maximum relative humidity by county (gridMET, 2018-2024)",
  "Mean Rmax (%)",
  "viridis"
)

combined_plot <- attribution_plot / tmax_plot / rmax_plot +
  patchwork::plot_layout(heights = c(1, 1, 1))

save_plot <- function(plot, basename, width = 13.5, height = 21) {
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

readr::write_csv(
  map_data |>
    st_drop_geometry() |>
    select("geoid", "county_label", "state_fips", "mean_tmax_c", "mean_rmax_pct", "map_group", "n_sites", "n_hospitals", "sites", "hospitals"),
  file.path(output_dir, "figure_county_maps_abc_source.csv")
)

save_plot(combined_plot, "figure_county_maps_abc")

message("Wrote ABC county map figure to ", output_dir)

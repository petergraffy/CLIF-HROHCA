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
output_dir <- file.path(repo_root, "output", "final", "county_weather_maps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(tigris)
  library(viridis)
})

options(tigris_use_cache = TRUE)

study_start <- as.Date("2018-01-01")
study_end <- as.Date("2024-12-31")
conus_excluded_state_fips <- c("02", "15", "60", "66", "69", "72", "78")

summarise_daymet <- function(path, value_col, out_col) {
  arrow::read_parquet(path, as_data_frame = TRUE) |>
    transmute(
      geoid = sprintf("%05s", as.character(.data$geoid)),
      county_name = as.character(.data$county_name),
      state_fips = sprintf("%02s", as.character(.data$state_fips)),
      date = as.Date(.data$date),
      value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) |>
    filter(
      .data$date >= study_start,
      .data$date <= study_end,
      !.data$state_fips %in% conus_excluded_state_fips
    ) |>
    group_by(.data$geoid, .data$county_name, .data$state_fips) |>
    summarise(
      "{out_col}" := mean(.data$value, na.rm = TRUE),
      n_days = sum(!is.na(.data$value)),
      .groups = "drop"
    )
}

tmax_mean <- summarise_daymet(
  file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet"),
  "tmax_mean_c",
  "mean_tmax_c"
)

rmax_mean <- summarise_daymet(
  file.path(repo_root, "exposome", "daymet_county_rmax_2018_2024.parquet"),
  "rmax_mean_pct",
  "mean_rmax_pct"
)

county_weather <- full_join(
  tmax_mean,
  rmax_mean,
  by = c("geoid", "county_name", "state_fips"),
  suffix = c("_tmax", "_rmax")
) |>
  rename(
    n_tmax_days = "n_days_tmax",
    n_rmax_days = "n_days_rmax"
  ) |>
  arrange(.data$geoid)

readr::write_csv(
  county_weather,
  file.path(output_dir, "county_mean_tmax_humidity_2018_2024.csv")
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
  left_join(county_weather, by = c("geoid", "state_fips"))

theme_map <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0, margin = margin(b = 8)),
      plot.caption = element_text(size = 8, color = "gray35", hjust = 0),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(12, 16, 12, 16)
    )
}

make_county_map <- function(data, fill_col, title, legend_title, palette, output_name) {
  p <- ggplot(data) +
    geom_sf(aes(fill = .data[[fill_col]]), color = NA) +
    scale_fill_viridis_c(
      option = palette,
      name = legend_title,
      na.value = "gray90",
      guide = guide_colorbar(barheight = unit(3.8, "in"), barwidth = unit(0.18, "in"))
    ) +
    coord_sf(datum = NA) +
    labs(
      title = title,
      subtitle = "County-level mean across daily Daymet values, 2018-2024",
      caption = "Projection: NAD83 / Conus Albers (EPSG:5070). Alaska, Hawaii, and territories excluded."
    ) +
    theme_map()

  ggsave(
    file.path(output_dir, output_name),
    p,
    width = 11,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  p
}

tmax_plot <- make_county_map(
  map_data,
  "mean_tmax_c",
  "Mean Daily Maximum Temperature by County",
  "Mean Tmax (C)",
  "magma",
  "county_mean_tmax_2018_2024.png"
)

rmax_plot <- make_county_map(
  map_data,
  "mean_rmax_pct",
  "Mean Daily Maximum Relative Humidity by County",
  "Mean max RH (%)",
  "viridis",
  "county_mean_rmax_2018_2024.png"
)

pdf(file.path(output_dir, "county_mean_tmax_humidity_2018_2024.pdf"), width = 11, height = 7)
print(tmax_plot)
print(rmax_plot)
dev.off()

combined_plot <- tmax_plot + rmax_plot + patchwork::plot_layout(ncol = 1)
ggsave(
  file.path(output_dir, "county_mean_tmax_humidity_2018_2024_combined.png"),
  combined_plot,
  width = 11,
  height = 13.5,
  dpi = 300,
  bg = "white"
)

summary_tbl <- tibble::tibble(
  study_start = as.character(study_start),
  study_end = as.character(study_end),
  n_counties_with_tmax = sum(!is.na(county_weather$mean_tmax_c)),
  n_counties_with_rmax = sum(!is.na(county_weather$mean_rmax_pct)),
  tmax_min_c = min(county_weather$mean_tmax_c, na.rm = TRUE),
  tmax_median_c = median(county_weather$mean_tmax_c, na.rm = TRUE),
  tmax_max_c = max(county_weather$mean_tmax_c, na.rm = TRUE),
  rmax_min_pct = min(county_weather$mean_rmax_pct, na.rm = TRUE),
  rmax_median_pct = median(county_weather$mean_rmax_pct, na.rm = TRUE),
  rmax_max_pct = max(county_weather$mean_rmax_pct, na.rm = TRUE)
)

readr::write_csv(summary_tbl, file.path(output_dir, "county_mean_tmax_humidity_2018_2024_summary.csv"))
print(summary_tbl)
message("Wrote county weather maps and summaries to ", output_dir)

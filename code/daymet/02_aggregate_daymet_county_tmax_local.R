#!/usr/bin/env Rscript

# Aggregate manually downloaded Daymet annual Tmax netCDF files to county-day means.
# This script is customized for Peter's 2018-2024 North America Tmax files.

required_packages <- c(
  "dplyr",
  "readr",
  "lubridate",
  "purrr",
  "stringr",
  "tibble",
  "tidyr",
  "sf",
  "terra",
  "tigris",
  "arrow"
)

install_if_missing <- function(packages) {
  missing_packages <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]

  if (length(missing_packages) > 0) {
    install.packages(missing_packages)
  }
}

install_if_missing(required_packages)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(sf)
  library(terra)
  library(tigris)
  library(arrow)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

daymet_config <- list(
  start_year = 2018L,
  end_year = 2024L,
  variable = "tmax",
  region = "na",
  raw_dir = Sys.getenv(
    "DAYMET_RAW_DIR",
    unset = "C:/Users/Peter Graffy/Downloads/Daymet_Daily_V4R1_2129_4.1-20260423_144640"
  ),
  processed_dir = Sys.getenv(
    "DAYMET_PROCESSED_DIR",
    unset = "output/intermediate/daymet/processed"
  ),
  county_boundaries_year = 2024L,
  county_boundaries_path = NULL,
  county_geoid_path = NULL,
  include_state_fips = character(),
  overwrite_existing_years = FALSE
)

dir.create(daymet_config$processed_dir, recursive = TRUE, showWarnings = FALSE)

options(tigris_use_cache = TRUE)
terra::terraOptions(progress = 1)

format_duration <- function(seconds) {
  seconds <- round(as.numeric(seconds))
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  sprintf("%02d:%02d:%02d", hours, minutes, secs)
}

build_nc_path <- function(raw_dir, variable, region, year) {
  file.path(raw_dir, sprintf("daymet_v4_daily_%s_%s_%s.nc", region, variable, year))
}

build_year_output_path <- function(processed_dir, variable, year, extension) {
  file.path(processed_dir, sprintf("daymet_county_%s_%s.%s", variable, year, extension))
}

build_status_path <- function(processed_dir, variable, start_year, end_year) {
  file.path(
    processed_dir,
    sprintf("daymet_county_%s_%s_%s_progress.csv", variable, start_year, end_year)
  )
}

build_final_output_stub <- function(variable, start_year, end_year) {
  sprintf("daymet_county_%s_%s_%s", variable, start_year, end_year)
}

load_counties <- function(daymet_cfg) {
  if (!is.null(daymet_cfg$county_boundaries_path) &&
      nzchar(daymet_cfg$county_boundaries_path)) {
    counties <- sf::read_sf(daymet_cfg$county_boundaries_path)
  } else {
    counties <- tigris::counties(
      cb = TRUE,
      year = daymet_cfg$county_boundaries_year,
      class = "sf"
    )
  }

  counties <- counties %>%
    st_make_valid() %>%
    select(GEOID, NAME, NAMELSAD, STATEFP, geometry) %>%
    filter(!STATEFP %in% c("15", "72", "60", "66", "69", "78"))

  if (length(daymet_cfg$include_state_fips) > 0) {
    counties <- counties %>%
      filter(STATEFP %in% as.character(daymet_cfg$include_state_fips))
  }

  if (!is.null(daymet_cfg$county_geoid_path) &&
      nzchar(daymet_cfg$county_geoid_path)) {
    included_geoids <- readr::read_csv(
      daymet_cfg$county_geoid_path,
      show_col_types = FALSE
    ) %>%
      rename(GEOID = 1) %>%
      mutate(GEOID = stringr::str_pad(GEOID, width = 5, side = "left", pad = "0")) %>%
      pull(GEOID) %>%
      unique()

    counties <- counties %>%
      filter(GEOID %in% included_geoids)
  }

  counties
}

daymet_dates <- function(year, n_layers) {
  dates <- seq.Date(
    from = as.Date(sprintf("%s-01-01", year)),
    to = as.Date(sprintf("%s-12-31", year)),
    by = "day"
  )

  if (leap_year(year)) {
    dates <- dates[dates != as.Date(sprintf("%s-12-31", year))]
  }

  if (length(dates) != n_layers) {
    stop(
      sprintf(
        "Date vector length (%s) does not match raster layers (%s) for %s.",
        length(dates),
        n_layers,
        year
      )
    )
  }

  dates
}

aggregate_year <- function(nc_path, counties_sf, variable, year) {
  message("Aggregating county-level ", variable, " for ", year)

  counties_daymet <- st_transform(
    counties_sf,
    crs = "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  )

  daymet_raster <- terra::rast(nc_path)
  daymet_raster <- terra::crop(daymet_raster, terra::ext(terra::vect(counties_daymet)))

  extracted <- terra::extract(
    daymet_raster,
    terra::vect(counties_daymet),
    fun = mean,
    na.rm = TRUE,
    exact = TRUE,
    ID = FALSE
  )

  names(extracted) <- as.character(daymet_dates(year, terra::nlyr(daymet_raster)))

  bind_cols(
    counties_sf %>%
      st_drop_geometry() %>%
      transmute(
        geoid = GEOID,
        county_name = NAMELSAD,
        state_fips = STATEFP
      ),
    as_tibble(extracted)
  ) %>%
    tidyr::pivot_longer(
      cols = -c(geoid, county_name, state_fips),
      names_to = "date",
      values_to = paste0(variable, "_mean_c")
    ) %>%
    mutate(
      date = as.Date(date),
      year = lubridate::year(date),
      month = lubridate::month(date),
      day = lubridate::day(date),
      yday = lubridate::yday(date)
    ) %>%
    arrange(geoid, date)
}

years <- seq.int(daymet_config$start_year, daymet_config$end_year)
counties_sf <- load_counties(daymet_config)
status_path <- build_status_path(
  processed_dir = daymet_config$processed_dir,
  variable = daymet_config$variable,
  start_year = daymet_config$start_year,
  end_year = daymet_config$end_year
)

if (nrow(counties_sf) == 0) {
  stop("No counties were selected. Check the county filters.")
}

missing_nc_files <- years[!file.exists(vapply(
  years,
  function(year) build_nc_path(daymet_config$raw_dir, daymet_config$variable, daymet_config$region, year),
  FUN.VALUE = character(1)
))]

if (length(missing_nc_files) > 0) {
  stop(
    paste(
      "Missing annual netCDF files for years:",
      paste(missing_nc_files, collapse = ", "),
      "in",
      daymet_config$raw_dir
    )
  )
}

status_tbl <- tibble(
  year = integer(),
  status = character(),
  started_at = character(),
  finished_at = character(),
  elapsed_seconds = double(),
  elapsed_hms = character(),
  estimated_remaining_hms = character(),
  rows_written = integer(),
  output_parquet = character(),
  note = character()
)

overall_start <- Sys.time()
pb <- utils::txtProgressBar(min = 0, max = length(years), style = 3)
yearly_parquet_files <- character(length(years))
completed_durations <- numeric()

for (i in seq_along(years)) {
  year <- years[[i]]
  year_started_at <- Sys.time()
  year_output_parquet <- build_year_output_path(
    processed_dir = daymet_config$processed_dir,
    variable = daymet_config$variable,
    year = year,
    extension = "parquet"
  )

  year_status <- "completed"
  year_note <- ""
  rows_written <- NA_integer_
  estimated_remaining_hms <- NA_character_

  tryCatch({
    if (file.exists(year_output_parquet) && !isTRUE(daymet_config$overwrite_existing_years)) {
      year_status <- "skipped_existing"
      year_note <- "Using existing yearly parquet."
      rows_written <- nrow(arrow::read_parquet(year_output_parquet))
      message(sprintf("[%s/%s] Skipping %s; yearly parquet already exists.", i, length(years), year))
    } else {
      nc_path <- build_nc_path(
        raw_dir = daymet_config$raw_dir,
        variable = daymet_config$variable,
        region = daymet_config$region,
        year = year
      )

      year_result <- aggregate_year(
        nc_path = nc_path,
        counties_sf = counties_sf,
        variable = daymet_config$variable,
        year = year
      )

      rows_written <- nrow(year_result)
      arrow::write_parquet(year_result, sink = year_output_parquet)
      rm(year_result)
      gc(verbose = FALSE)

      message(sprintf("[%s/%s] Wrote %s", i, length(years), basename(year_output_parquet)))
    }
  }, error = function(e) {
    year_status <<- "failed"
    year_note <<- conditionMessage(e)
  })

  year_finished_at <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(year_finished_at, year_started_at, units = "secs"))
  elapsed_hms <- format_duration(elapsed_seconds)

  if (identical(year_status, "completed")) {
    completed_durations <- c(completed_durations, elapsed_seconds)
  }

  if (length(completed_durations) > 0) {
    remaining_years <- length(years) - i
    estimated_remaining_hms <- format_duration(mean(completed_durations) * remaining_years)
  }

  status_tbl <- bind_rows(
    status_tbl,
    tibble(
      year = year,
      status = year_status,
      started_at = format(year_started_at, "%Y-%m-%d %H:%M:%S"),
      finished_at = format(year_finished_at, "%Y-%m-%d %H:%M:%S"),
      elapsed_seconds = elapsed_seconds,
      elapsed_hms = elapsed_hms,
      estimated_remaining_hms = estimated_remaining_hms,
      rows_written = rows_written,
      output_parquet = year_output_parquet,
      note = year_note
    )
  )

  readr::write_csv(status_tbl, status_path)
  yearly_parquet_files[[i]] <- year_output_parquet
  utils::setTxtProgressBar(pb, i)

  message(
    sprintf(
      "[%s/%s] %s finished with status=%s in %s. Estimated remaining: %s",
      i,
      length(years),
      year,
      year_status,
      elapsed_hms,
      estimated_remaining_hms %||% "unknown"
    )
  )

  if (identical(year_status, "failed")) {
    close(pb)
    stop(sprintf("Year %s failed: %s", year, year_note), call. = FALSE)
  }
}

close(pb)

output_stub <- build_final_output_stub(
  variable = daymet_config$variable,
  start_year = daymet_config$start_year,
  end_year = daymet_config$end_year
)

yearly_parquet_files <- yearly_parquet_files[file.exists(yearly_parquet_files)]

message("Combining yearly parquet files into final outputs.")
county_tmax <- purrr::map_dfr(yearly_parquet_files, arrow::read_parquet)

arrow::write_parquet(
  county_tmax,
  sink = file.path(daymet_config$processed_dir, paste0(output_stub, ".parquet"))
)

readr::write_csv(
  county_tmax,
  file = file.path(daymet_config$processed_dir, paste0(output_stub, ".csv"))
)

readr::write_csv(
  counties_sf %>%
    st_drop_geometry() %>%
    transmute(
      geoid = GEOID,
      county_name = NAMELSAD,
      state_fips = STATEFP
    ) %>%
    distinct(),
  file = file.path(daymet_config$processed_dir, "county_lookup.csv")
)

message(
  "Saved county-level Daymet data to ",
  normalizePath(daymet_config$processed_dir),
  ". Total elapsed time: ",
  format_duration(as.numeric(difftime(Sys.time(), overall_start, units = "secs")))
)

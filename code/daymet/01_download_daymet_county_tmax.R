#!/usr/bin/env Rscript

# Build county-level daily Daymet Tmax exposure data for 2018-2024.
# The script downloads annual Daymet netCDF files, computes county-level
# area-weighted mean Tmax, and saves the results for downstream CLIF linkage.

required_packages <- c(
  "jsonlite",
  "dplyr",
  "readr",
  "lubridate",
  "purrr",
  "stringr",
  "tibble",
  "tidyr",
  "httr",
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
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(httr)
  library(sf)
  library(terra)
  library(tigris)
  library(arrow)
})

load_project_config <- function(path = "config/config.json") {
  if (!file.exists(path)) {
    message("config/config.json not found; using default Daymet settings.")
    return(list())
  }

  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

config <- load_project_config()

default_daymet_config <- list(
  start_year = 2018L,
  end_year = 2024L,
  variable = "tmax",
  region = "na",
  dataset_id = "2129",
  collection_id = "C2532426483-ORNL_CLOUD",
  product_prefix = "Daymet_Daily_V4R1",
  direct_download_base = "https://data.ornldaac.earthdata.nasa.gov/protected/daymet/Daymet_Daily_V4R1/data",
  county_boundaries_year = 2024L,
  county_boundaries_path = NULL,
  county_geoid_path = NULL,
  include_state_fips = character(),
  raw_dir = "output/intermediate/daymet/raw",
  processed_dir = "output/intermediate/daymet/processed",
  overwrite_existing_years = FALSE
)

merge_config <- function(defaults, supplied) {
  if (is.null(supplied)) {
    return(defaults)
  }

  defaults[names(supplied)] <- supplied
  defaults
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

daymet_config <- merge_config(default_daymet_config, config$daymet)

apply_env_overrides <- function(daymet_cfg) {
  raw_dir_override <- Sys.getenv("DAYMET_RAW_DIR", unset = "")
  processed_dir_override <- Sys.getenv("DAYMET_PROCESSED_DIR", unset = "")

  if (nzchar(raw_dir_override)) {
    daymet_cfg$raw_dir <- raw_dir_override
  }

  if (nzchar(processed_dir_override)) {
    daymet_cfg$processed_dir <- processed_dir_override
  }

  daymet_cfg
}

daymet_config <- apply_env_overrides(daymet_config)

dir.create(daymet_config$raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(daymet_config$processed_dir, recursive = TRUE, showWarnings = FALSE)

options(tigris_use_cache = TRUE)
terra::terraOptions(progress = 1)

build_daymet_granule_name <- function(year, variable, region, product_prefix) {
  sprintf(
    "%s.daymet_v4_daily_%s_%s_%s.nc",
    product_prefix,
    region,
    variable,
    year
  )
}

build_daymet_opendap_url <- function(year, variable, region, collection_id, product_prefix) {
  granule_name <- build_daymet_granule_name(year, variable, region, product_prefix)
  paste(
    "https://opendap.earthdata.nasa.gov/collections",
    collection_id,
    "granules",
    granule_name,
    sep = "/"
  )
}

build_daymet_direct_url <- function(year, variable, region, direct_download_base) {
  filename <- sprintf("daymet_v4_daily_%s_%s_%s.nc", region, variable, year)
  paste(direct_download_base, filename, sep = "/")
}

earthdata_download_config <- function() {
  username <- Sys.getenv("EARTHDATA_USERNAME", unset = "")
  password <- Sys.getenv("EARTHDATA_PASSWORD", unset = "")

  if (!nzchar(username) || !nzchar(password)) {
    return(NULL)
  }

  netrc_path <- tempfile("earthdata_", fileext = ".netrc")
  cookie_path <- tempfile("earthdata_", fileext = ".cookies")

  writeLines(
    sprintf(
      "machine urs.earthdata.nasa.gov login %s password %s",
      username,
      password
    ),
    con = netrc_path,
    useBytes = TRUE
  )

  cookie_file <- file(cookie_path, open = "wb")
  close(cookie_file)

  list(
    netrc_path = netrc_path,
    cookie_path = cookie_path,
    httr_config = httr::config(
      followlocation = 1L,
      netrc = 1L,
      netrc_file = netrc_path,
      cookiefile = cookie_path,
      cookiejar = cookie_path
    )
  )
}

download_daymet_year <- function(year, variable, region, raw_dir, direct_download_base) {
  filename <- sprintf("daymet_v4_daily_%s_%s_%s.nc", region, variable, year)
  destination <- file.path(raw_dir, filename)

  if (!file.exists(destination)) {
    url <- build_daymet_direct_url(year, variable, region, direct_download_base)
    message("Downloading ", filename)
    auth <- earthdata_download_config()

    tryCatch({
      if (is.null(auth)) {
        download.file(url = url, destfile = destination, mode = "wb", quiet = FALSE)
      } else {
        response <- httr::GET(
          url = url,
          auth$httr_config,
          httr::write_disk(destination, overwrite = TRUE),
          httr::progress()
        )

        content_type <- httr::headers(response)[["content-type"]] %||% ""
        final_url <- response$url %||% url

        if (httr::status_code(response) >= 400) {
          stop(sprintf("HTTP %s from %s", httr::status_code(response), final_url))
        }

        if (grepl("text/html", content_type, ignore.case = TRUE)) {
          stop(
            paste(
              "Received HTML instead of a Daymet netCDF file.",
              "This usually means Earthdata application authorization is still required",
              "for the current ORNL DAAC Earthdata download route."
            )
          )
        }
      }
    }, error = function(e) {
      if (file.exists(destination)) {
        unlink(destination, force = TRUE)
      }

      stop(
        paste(
            "Daymet download failed.",
            "If you are using Earthdata credentials, confirm they are valid and that the ORNL DAAC application has been authorized for your account.",
            "You can also manually place the annual netCDFs into", raw_dir, "and rerun the script.",
            "Original error:", conditionMessage(e)
        ),
        call. = FALSE
      )
    }, finally = {
      if (!is.null(auth)) {
        unlink(c(auth$netrc_path, auth$cookie_path), force = TRUE)
      }
    })
  } else {
    message("Using cached file ", filename)
  }

  destination
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
    select(GEOID, NAME, NAMELSAD, STATEFP, geometry)

  if (identical(daymet_cfg$region, "na")) {
    counties <- counties %>%
      filter(!STATEFP %in% c("15", "72", "60", "66", "69", "78"))
  } else if (identical(daymet_cfg$region, "hi")) {
    counties <- counties %>%
      filter(STATEFP == "15")
  } else if (identical(daymet_cfg$region, "pr")) {
    counties <- counties %>%
      filter(STATEFP == "72")
  }

  if (!is.null(daymet_cfg$include_state_fips) &&
      length(daymet_cfg$include_state_fips) > 0) {
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

format_duration <- function(seconds) {
  seconds <- round(as.numeric(seconds))
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60

  sprintf("%02d:%02d:%02d", hours, minutes, secs)
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

write_status_table <- function(status_tbl, path) {
  readr::write_csv(status_tbl, path)
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
  stop("No counties were selected. Check the county filters in config/config.json.")
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
      rows_written <- arrow::read_parquet(year_output_parquet, as_data_frame = FALSE) %>%
        nrow()
      message(sprintf("[%s/%s] Skipping %s; yearly parquet already exists.", i, length(years), year))
    } else {
      nc_path <- download_daymet_year(
        year = year,
        variable = daymet_config$variable,
        region = daymet_config$region,
        raw_dir = daymet_config$raw_dir,
        direct_download_base = daymet_config$direct_download_base
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

  write_status_table(status_tbl, status_path)
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

required_packages <- c(
  "arrow",
  "dlnm",
  "dplyr",
  "ggplot2",
  "jsonlite",
  "lubridate",
  "MASS",
  "readr",
  "stringr",
  "tibble",
  "tidyr",
  "tsModel"
)

ensure_user_library <- function(repo_root) {
  version_parts <- strsplit(as.character(getRversion()), "\\.")[[1]]
  version_stub <- paste(version_parts[1], version_parts[2], sep = ".")
  user_lib <- file.path(repo_root, ".r-user-lib", version_stub)
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  invisible(user_lib)
}

install_missing_packages <- function(repo_root, packages = required_packages) {
  ensure_user_library(repo_root)
  missing_packages <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_packages) > 0) {
    install.packages(
      missing_packages,
      lib = .libPaths()[1],
      repos = "https://packagemanager.posit.co/cran/latest",
      dependencies = TRUE
    )
  }
  invisible(packages)
}

load_project_config <- function(repo_root) {
  config_path <- file.path(repo_root, "config", "config.json")
  if (!file.exists(config_path)) {
    stop("Missing config/config.json. Copy config/config_template.json and update site_name and tables_path.")
  }
  jsonlite::fromJSON(config_path, simplifyVector = TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

resolve_tables_path <- function(config) {
  env_path <- Sys.getenv("CLIF_TABLES_PATH", unset = "")
  candidate <- if (nzchar(env_path)) env_path else config$tables_path %||% ""
  if (!nzchar(candidate)) {
    stop("Missing CLIF tables path. Set CLIF_TABLES_PATH or config$config.json tables_path.")
  }
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

resolve_file_type <- function(config) {
  file_type <- tolower(trimws(Sys.getenv("CLIF_FILE_TYPE", unset = config$file_type %||% "parquet")))
  if (!file_type %in% c("csv", "parquet")) stop("Unsupported CLIF file_type. Use csv or parquet.")
  file_type
}

read_clif_table <- function(tables_path, file_type, table_name, columns = NULL, required = TRUE) {
  path <- file.path(tables_path, sprintf("clif_%s.%s", table_name, file_type))
  if (!file.exists(path)) {
    if (required) stop("Missing CLIF table: ", path)
    return(NULL)
  }
  message("Reading ", basename(path))
  if (file_type == "parquet") {
    if (is.null(columns)) return(arrow::read_parquet(path))
    return(arrow::read_parquet(path, col_select = dplyr::all_of(columns)))
  }
  out <- readr::read_csv(path, show_col_types = FALSE)
  if (!is.null(columns)) out <- out[, intersect(columns, names(out)), drop = FALSE]
  out
}

as_utc_datetime <- function(x) {
  lubridate::as_datetime(x, tz = "UTC")
}

normalize_county_fips <- function(x) {
  clean <- x |>
    tidyr::replace_na("") |>
    as.character() |>
    stringr::str_replace("\\.0$", "") |>
    stringr::str_replace_all("[^0-9]", "") |>
    stringr::str_trim() |>
    stringr::str_pad(width = 5, side = "left", pad = "0") |>
    stringr::str_sub(-5)
  dplyr::na_if(clean, "00000")
}

normalize_hospital_id <- function(x) {
  clean <- as.character(tidyr::replace_na(x, ""))
  clean <- stringr::str_trim(clean)
  dplyr::na_if(clean, "")
}

parse_pipe_set <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(as.character(x))) return(character())
  out <- unlist(strsplit(as.character(x), "\\|", fixed = FALSE))
  out <- trimws(out)
  out[nzchar(out)]
}

build_hospital_id_lookup <- function(adt) {
  if (!"hospital_id" %in% names(adt)) {
    return(tibble::tibble(hospitalization_id = character(), hospital_id = character()))
  }
  out <- adt |>
    dplyr::transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      hospital_id = normalize_hospital_id(.data$hospital_id),
      in_dttm = if ("in_dttm" %in% names(adt)) as_utc_datetime(.data$in_dttm) else as.POSIXct(NA)
    ) |>
    dplyr::filter(!is.na(.data$hospital_id)) |>
    dplyr::arrange(.data$hospitalization_id, .data$in_dttm, .data$hospital_id) |>
    dplyr::distinct(.data$hospitalization_id, .keep_all = TRUE) |>
    dplyr::select("hospitalization_id", "hospital_id")
  out
}

get_site_geography <- function(repo_root, config) {
  site_name <- stringr::str_trim(as.character(config$site_name %||% ""))
  geography_path <- file.path(repo_root, "reference", "clif_hospital_geography.csv")

  if (!nzchar(site_name)) stop("config$site_name is required for hospital county assignment.")
  if (!file.exists(geography_path)) stop("Missing reference/clif_hospital_geography.csv")

  geo <- readr::read_csv(geography_path, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character())) |>
    dplyr::filter(stringr::str_trim(.data$site_name) == site_name) |>
    dplyr::mutate(
      hospital_id = normalize_hospital_id(.data$hospital_id),
      hospital_county_fips = normalize_county_fips(.data$hospital_county_fips)
    )

  if (nrow(geo) == 0) {
    stop("No rows for site_name = '", site_name, "' in reference/clif_hospital_geography.csv")
  }

  default <- geo[1, , drop = FALSE]
  list(
    site_name = site_name,
    source = "site_hospital_geography",
    default_hospital_id = default$hospital_id[[1]],
    default_hospital_county_fips = default$hospital_county_fips[[1]],
    geography = geo
  )
}

apply_site_county_assignment <- function(df, repo_root, config, hospital_id_col = "hospital_id", county_col = "county_fips") {
  site_geo <- get_site_geography(repo_root, config)
  geo <- site_geo$geography
  out <- df

  if (!hospital_id_col %in% names(out)) out[[hospital_id_col]] <- NA_character_
  out[[hospital_id_col]] <- normalize_hospital_id(out[[hospital_id_col]])
  out[[county_col]] <- normalize_county_fips(out[[county_col]])
  out$home_county_fips <- out[[county_col]]
  out$county_fips_was_overridden <- 0L
  out$county_fips_geocode_reason <- "same_or_adjacent_home_county"
  out$assigned_hospital_id <- out[[hospital_id_col]]
  out$assigned_hospital_county_fips <- NA_character_
  out$hospital_metadata_source <- site_geo$source

  geo_by_id <- split(geo, geo$hospital_id)

  for (i in seq_len(nrow(out))) {
    hospital_id <- out[[hospital_id_col]][[i]]
    hospital_meta <- if (!is.na(hospital_id) && hospital_id %in% names(geo_by_id)) geo_by_id[[hospital_id]][1, , drop = FALSE] else NULL

    if (is.null(hospital_meta)) {
      hospital_meta <- geo[1, , drop = FALSE]
      out$assigned_hospital_id[[i]] <- site_geo$default_hospital_id
    }

    hospital_county <- hospital_meta$hospital_county_fips[[1]]
    adjacent <- parse_pipe_set(hospital_meta$adjacent_county_fips[[1]])
    local_counties <- unique(c(hospital_county, adjacent))
    out$assigned_hospital_county_fips[[i]] <- hospital_county

    home_county <- out$home_county_fips[[i]]
    if (is.na(home_county)) {
      out[[county_col]][[i]] <- hospital_county
      out$county_fips_was_overridden[[i]] <- 1L
      out$county_fips_geocode_reason[[i]] <- "missing_home_county_overridden_to_hospital"
    } else if (!home_county %in% local_counties) {
      out[[county_col]][[i]] <- hospital_county
      out$county_fips_was_overridden[[i]] <- 1L
      out$county_fips_geocode_reason[[i]] <- "nonlocal_home_county_overridden_to_hospital"
    }
  }

  out
}

read_exposome_daymet <- function(repo_root, filename, value_col) {
  path <- file.path(repo_root, "exposome", filename)
  if (!file.exists(path)) stop("Missing exposome file: ", path)
  arrow::read_parquet(path) |>
    dplyr::transmute(
      county_fips = normalize_county_fips(.data$geoid),
      admission_date = as.Date(.data$date),
      !!value_col := suppressWarnings(as.numeric(.data[[value_col]]))
    )
}

read_exposome_pollution <- function(repo_root, filename, value_col) {
  path <- file.path(repo_root, "exposome", filename)
  if (!file.exists(path)) stop("Missing exposome file: ", path)
  readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_guess(), GEOID = readr::col_character())) |>
    dplyr::transmute(
      county_fips = normalize_county_fips(.data$GEOID),
      year = as.integer(.data$year),
      !!value_col := suppressWarnings(as.numeric(.data[[value_col]]))
    )
}

safe_pct <- function(x, denom) {
  ifelse(denom > 0, 100 * x / denom, NA_real_)
}

fmt_n_pct <- function(n, denom) {
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), safe_pct(n, denom))
}

fmt_median_iqr <- function(x) {
  sprintf(
    "%.1f [%.1f, %.1f]",
    stats::median(x, na.rm = TRUE),
    stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  )
}

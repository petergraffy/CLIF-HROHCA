required_packages <- c(
  "arrow",
  "cmprsk",
  "dlnm",
  "dplyr",
  "ggplot2",
  "jsonlite",
  "lubridate",
  "readr",
  "stringr",
  "survival",
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

validate_site_name <- function(config, repo_root) {
  site_name_raw <- if (is.null(config$site_name) || length(config$site_name) == 0) "" else config$site_name
  site_name <- stringr::str_trim(as.character(site_name_raw))
  geography_path <- file.path(repo_root, "reference", "clif_hospital_geography.csv")

  if (!nzchar(site_name)) stop("config$site_name is required.")
  if (!file.exists(geography_path)) stop("Missing reference/clif_hospital_geography.csv")

  geography_sites <- readr::read_csv(
    geography_path,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    dplyr::pull("site_name") |>
    stringr::str_trim()

  if (!site_name %in% geography_sites) {
    stop(
      "config$site_name = '", site_name, "' was not found in reference/clif_hospital_geography.csv. ",
      "Use one of: ", paste(sort(unique(geography_sites)), collapse = ", ")
    )
  }

  invisible(site_name)
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
    available_columns <- names(arrow::open_dataset(path, format = "parquet"))
    missing_columns <- setdiff(columns, available_columns)
    selected_columns <- intersect(columns, available_columns)
    if (length(missing_columns) > 0) {
      warning(
        "Missing requested column(s) in ", basename(path), ": ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    out <- if (length(selected_columns) == 0) {
      tibble::tibble()
    } else {
      arrow::read_parquet(path, col_select = dplyr::all_of(selected_columns))
    }
    for (missing_col in missing_columns) out[[missing_col]] <- NA
    return(out[, columns, drop = FALSE])
  }
  out <- readr::read_csv(path, show_col_types = FALSE)
  if (!is.null(columns)) {
    missing_columns <- setdiff(columns, names(out))
    if (length(missing_columns) > 0) {
      warning(
        "Missing requested column(s) in ", basename(path), ": ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    for (missing_col in missing_columns) out[[missing_col]] <- NA
    out <- out[, columns, drop = FALSE]
  }
  out
}

as_utc_datetime <- function(x) {
  lubridate::as_datetime(x, tz = "UTC")
}

as_clif_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  date_part <- stringr::str_extract(as.character(x), "^\\d{4}-\\d{2}-\\d{2}")
  as.Date(date_part)
}

normalize_category <- function(x) {
  x |>
    tidyr::replace_na("") |>
    as.character() |>
    stringr::str_squish() |>
    stringr::str_to_lower()
}

is_black_race <- function(x) {
  normalize_category(x) %in% c("black", "black or african american", "african american")
}

is_male <- function(x) {
  normalize_category(x) == "male"
}

is_female <- function(x) {
  normalize_category(x) == "female"
}

is_expired_discharge <- function(x) {
  stringr::str_detect(normalize_category(x), "expired|death|dead")
}

derive_poa_flag <- function(df, diagnosis_source = "") {
  if ("poa_present" %in% names(df)) {
    raw <- df$poa_present
    numeric_flag <- suppressWarnings(as.integer(as.character(raw)))
    char_flag <- stringr::str_to_upper(as.character(raw))
    return(ifelse(!is.na(numeric_flag), numeric_flag == 1L, char_flag %in% c("Y", "YES", "TRUE", "T", "1")))
  }
  if (identical(diagnosis_source, "admission_diagnosis")) return(rep(TRUE, nrow(df)))
  rep(FALSE, nrow(df))
}

derive_ohca_mechanism <- function(dx, hospitalization_ids = NULL, diagnosis_source = "") {
  if (is.null(dx) || nrow(dx) == 0) {
    out <- tibble::tibble(
      hospitalization_id = as.character(hospitalization_ids %||% character()),
      ohca_mechanism = "unclear_other",
      ohca_mechanism_detail = NA_character_
    )
    return(out)
  }

  dx2 <- dx |>
    dplyr::mutate(
      hospitalization_id = as.character(.data$hospitalization_id),
      diagnosis_code = if ("diagnosis_code" %in% names(dx)) as.character(.data$diagnosis_code) else NA_character_,
      diagnosis_code_format = if ("diagnosis_code_format" %in% names(dx)) as.character(.data$diagnosis_code_format) else NA_character_,
      diagnosis_code_clean = .data$diagnosis_code |>
        tidyr::replace_na("") |>
        as.character() |>
        stringr::str_to_upper() |>
        stringr::str_replace_all("[^A-Z0-9]", ""),
      diagnosis_code_format = dplyr::if_else(
        nchar(tidyr::replace_na(.data$diagnosis_code_format, "")) > 0,
        stringr::str_to_upper(.data$diagnosis_code_format),
        dplyr::if_else(stringr::str_detect(.data$diagnosis_code_clean, "^[A-Z]"), "ICD10", "ICD9")
      ),
      poa_flag = derive_poa_flag(dx, diagnosis_source),
      icd10 = stringr::str_detect(.data$diagnosis_code_format, "10"),
      icd9 = stringr::str_detect(.data$diagnosis_code_format, "9")
    ) |>
    dplyr::filter(.data$poa_flag, nzchar(.data$diagnosis_code_clean))

  if (!is.null(hospitalization_ids)) {
    dx2 <- dx2 |> dplyr::filter(.data$hospitalization_id %in% as.character(hospitalization_ids))
  }

  flags <- dx2 |>
    dplyr::mutate(
      mechanism = dplyr::case_when(
        .data$icd10 & stringr::str_detect(.data$diagnosis_code_clean, "^(T36|T37|T38|T39|T40|T41|T42|T43|T44|T45|T46|T47|T48|T49|T50|T51|T52|T53|T54|T55|T56|T57|T58|T59|T60|T61|T62|T63|T64|T65)") ~ "overdose_toxicologic",
        .data$icd9 & stringr::str_detect(.data$diagnosis_code_clean, "^(96[0-9]|97[0-9]|98[0-9])") ~ "overdose_toxicologic",
        .data$icd10 & stringr::str_detect(.data$diagnosis_code_clean, "^(S|T0[7-9]|T1[0-4]|T2|T3|T6[6-9]|T7|T8|V|W|X[0-5]|Y0)") ~ "trauma",
        .data$icd9 & stringr::str_detect(.data$diagnosis_code_clean, "^(8[0-9][0-9]|9[0-5][0-9])") ~ "trauma",
        .data$icd10 & stringr::str_detect(.data$diagnosis_code_clean, "^(J|I26|T17|T71|W65|W66|W67|W68|W69|W70|W71|W72|W73|W74|X71)") ~ "respiratory_asphyxial",
        .data$icd9 & stringr::str_detect(.data$diagnosis_code_clean, "^(4[6-9][0-9]|50[0-8]|4151|9941)") ~ "respiratory_asphyxial",
        .data$icd10 & stringr::str_detect(.data$diagnosis_code_clean, "^(I2[0-5]|I3|I4|I50|I7[01])") ~ "cardiac",
        .data$icd9 & stringr::str_detect(.data$diagnosis_code_clean, "^(41[0-4]|42[0-9]|428|441)") ~ "cardiac",
        .data$icd10 & stringr::str_detect(.data$diagnosis_code_clean, "^(I6|G4[0-1]|G93)") ~ "neurologic",
        .data$icd9 & stringr::str_detect(.data$diagnosis_code_clean, "^(43[0-8]|345|348)") ~ "neurologic",
        .data$icd10 & stringr::str_detect(.data$diagnosis_code_clean, "^(A4[0-1]|R65|J1[0-8]|N39)") ~ "sepsis_infection",
        .data$icd9 & stringr::str_detect(.data$diagnosis_code_clean, "^(038|9959|48[0-7]|5990)") ~ "sepsis_infection",
        .data$icd10 & stringr::str_detect(.data$diagnosis_code_clean, "^(E1[0-4]|E8[3-7])") ~ "metabolic_electrolyte",
        .data$icd9 & stringr::str_detect(.data$diagnosis_code_clean, "^(25[0-1]|27[5-6])") ~ "metabolic_electrolyte",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(.data$mechanism))

  priority <- c(
    "trauma",
    "overdose_toxicologic",
    "respiratory_asphyxial",
    "cardiac",
    "neurologic",
    "sepsis_infection",
    "metabolic_electrolyte"
  )

  mech <- flags |>
    dplyr::mutate(priority = match(.data$mechanism, priority)) |>
    dplyr::arrange(.data$hospitalization_id, .data$priority, .data$diagnosis_code_clean) |>
    dplyr::group_by(.data$hospitalization_id) |>
    dplyr::summarise(
      ohca_mechanism = dplyr::first(.data$mechanism),
      ohca_mechanism_detail = paste(sort(unique(paste0(.data$mechanism, ":", .data$diagnosis_code_clean))), collapse = " | "),
      .groups = "drop"
    )

  all_ids <- tibble::tibble(hospitalization_id = as.character(hospitalization_ids %||% unique(dx2$hospitalization_id)))
  all_ids |>
    dplyr::left_join(mech, by = "hospitalization_id") |>
    dplyr::mutate(
      ohca_mechanism = tidyr::replace_na(.data$ohca_mechanism, "unclear_other")
    )
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
  if (length(denom) == 1L) {
    if (is.na(denom) || denom <= 0) return(rep(NA_real_, length(x)))
    return(100 * x / denom)
  }
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

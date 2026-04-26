#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) stop("Could not determine script path from commandArgs().")
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "code", "00_project_functions.R"))
ensure_user_library(repo_root)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(readr)
  library(stringr)
  library(tidyr)
})

OHCA_PATH <- file.path(repo_root, "output", "intermediate", "cohorts", "ohca", "ohca_poa_icu_2018_2024.csv")
TMAX_PATH <- file.path(repo_root, "exposome", "daymet_county_tmax_2018_2024_conus.parquet")
OUTPUT_DIR <- file.path(repo_root, "output", "final", "descriptive")
FIGURE_DIR <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

WARM_MONTHS <- c(5L, 6L, 7L, 8L, 9L)
TRAJECTORY_HOURS <- 0:72
CUMULATIVE_HOURS <- 0:168
SMOOTH_HALF_WINDOW_HOURS <- 3L
HEAT_DEFINITIONS <- tibble::tibble(
  heat_definition = c("heat95", "heat90"),
  heat_percentile = c(0.95, 0.90)
)

config <- load_project_config(repo_root)
tables_path <- resolve_tables_path(config)
file_type <- resolve_file_type(config)

safe_fisher_p <- function(group, value) {
  ok <- !is.na(group) & !is.na(value)
  if (length(unique(group[ok])) < 2 || length(unique(value[ok])) < 2) return(NA_real_)
  stats::fisher.test(table(group[ok], value[ok]), simulate.p.value = TRUE, B = 10000)$p.value
}

safe_wilcox_p <- function(group, value) {
  ok <- !is.na(group) & !is.na(value)
  if (length(unique(group[ok])) < 2 || any(table(group[ok]) == 0)) return(NA_real_)
  suppressWarnings(stats::wilcox.test(value[ok] ~ group[ok])$p.value)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

make_heat_groups <- function(ohca, heat_definition, heat_percentile) {
  threshold <- as.numeric(stats::quantile(ohca$tmax_mean_c, heat_percentile, na.rm = TRUE, names = FALSE))
  out <- ohca |>
    mutate(
      heat_definition = heat_definition,
      heat_percentile = heat_percentile,
      heat_threshold_tmax_c = threshold,
      heat_related_ohca = ifelse(
        !is.na(.data$tmax_mean_c) & .data$tmax_mean_c >= threshold,
        "Heat-related OHCA",
        "Non-heat-related OHCA"
      ),
      heat_related_ohca = factor(.data$heat_related_ohca, levels = c("Heat-related OHCA", "Non-heat-related OHCA"))
    )
  out
}

make_binary_row <- function(df, label, var) {
  heat <- df[df$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- df[df$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  tibble::tibble(
    heat_definition = unique(df$heat_definition),
    characteristic = label,
    heat_related_ohca = fmt_n_pct(sum(heat[[var]] == 1, na.rm = TRUE), sum(!is.na(heat[[var]]))),
    non_heat_related_ohca = fmt_n_pct(sum(non_heat[[var]] == 1, na.rm = TRUE), sum(!is.na(non_heat[[var]]))),
    p_value = fmt_p(safe_fisher_p(df$heat_related_ohca, df[[var]]))
  )
}

make_continuous_row <- function(df, label, var, subset_var = NULL) {
  if (!is.null(subset_var)) df <- df[df[[subset_var]] == 1, , drop = FALSE]
  heat <- df[df$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- df[df$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  tibble::tibble(
    heat_definition = unique(df$heat_definition),
    characteristic = label,
    heat_related_ohca = fmt_median_iqr(heat[[var]]),
    non_heat_related_ohca = fmt_median_iqr(non_heat[[var]]),
    p_value = fmt_p(safe_wilcox_p(df$heat_related_ohca, df[[var]]))
  )
}

make_category_rows <- function(df, label, var) {
  levels_found <- names(sort(table(df[[var]], useNA = "no"), decreasing = TRUE))
  p <- safe_fisher_p(df$heat_related_ohca, df[[var]])
  dplyr::bind_rows(lapply(seq_along(levels_found), function(i) {
    level_value <- levels_found[[i]]
    heat <- df[df$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
    non_heat <- df[df$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
    tibble::tibble(
      heat_definition = unique(df$heat_definition),
      characteristic = paste0(label, ": ", level_value),
      heat_related_ohca = fmt_n_pct(sum(heat[[var]] == level_value, na.rm = TRUE), sum(!is.na(heat[[var]]))),
      non_heat_related_ohca = fmt_n_pct(sum(non_heat[[var]] == level_value, na.rm = TRUE), sum(!is.na(non_heat[[var]]))),
      p_value = ifelse(i == 1, fmt_p(p), "")
    )
  }))
}

make_table2 <- function(df) {
  heat <- df[df$heat_related_ohca == "Heat-related OHCA", , drop = FALSE]
  non_heat <- df[df$heat_related_ohca == "Non-heat-related OHCA", , drop = FALSE]
  dplyr::bind_rows(
    tibble::tibble(
      heat_definition = unique(df$heat_definition),
      characteristic = "N",
      heat_related_ohca = format(nrow(heat), big.mark = ","),
      non_heat_related_ohca = format(nrow(non_heat), big.mark = ","),
      p_value = ""
    ),
    make_continuous_row(df, "Age, years", "age_at_admission"),
    make_category_rows(df, "Age group", "age_group"),
    make_category_rows(df, "Sex", "sex_category"),
    make_category_rows(df, "Race", "race_group"),
    make_category_rows(df, "Ethnicity", "ethnicity_category"),
    make_binary_row(df, "County reassigned to hospital county", "county_fips_was_overridden"),
    make_continuous_row(df, "Assigned-county Tmax, C", "tmax_mean_c"),
    make_binary_row(df, "Invasive mechanical ventilation", "imv_any"),
    make_continuous_row(df, "IMV duration among ventilated, days", "imv_duration_days", subset_var = "imv_any"),
    make_binary_row(df, "Vasopressor use", "vasopressor_any"),
    make_continuous_row(df, "ICU length of stay, days", "icu_los_days"),
    make_continuous_row(df, "Hospital length of stay, days", "hospital_los_days"),
    make_binary_row(df, "Hospital death", "hospital_death"),
    make_binary_row(df, "Death or hospice", "death_or_hospice"),
    make_continuous_row(df, "Time to death among decedents, days", "hospital_los_days", subset_var = "hospital_death")
  )
}

summarize_group <- function(df) {
  df |>
    group_by(.data$heat_definition, .data$heat_percentile, .data$heat_threshold_tmax_c, .data$heat_related_ohca) |>
    summarise(
      n = n(),
      median_age = median(.data$age_at_admission, na.rm = TRUE),
      hospital_mortality_pct = safe_pct(sum(.data$hospital_death == 1, na.rm = TRUE), n()),
      death_or_hospice_pct = safe_pct(sum(.data$death_or_hospice == 1, na.rm = TRUE), n()),
      imv_pct = safe_pct(sum(.data$imv_any == 1, na.rm = TRUE), sum(!is.na(.data$imv_any))),
      median_imv_duration_days = median(.data$imv_duration_days[.data$imv_any == 1], na.rm = TRUE),
      vasopressor_pct = safe_pct(sum(.data$vasopressor_any == 1, na.rm = TRUE), sum(!is.na(.data$vasopressor_any))),
      median_icu_los_days = median(.data$icu_los_days, na.rm = TRUE),
      median_hospital_los_days = median(.data$hospital_los_days, na.rm = TRUE),
      median_time_to_death_days = median(.data$hospital_los_days[.data$hospital_death == 1], na.rm = TRUE),
      .groups = "drop"
    )
}

hourly_denominators <- function(cohort, hours) {
  cohort |>
    select("hospitalization_id", "heat_definition", "heat_related_ohca", "icu_los_hours") |>
    tidyr::crossing(icu_hour = hours) |>
    filter(!is.na(.data$icu_los_hours), .data$icu_los_hours >= .data$icu_hour) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$icu_hour) |>
    summarise(n_at_risk = n_distinct(.data$hospitalization_id), .groups = "drop")
}

summarize_hourly_measurements <- function(measure_df, cohort, datetime_col, category_col, value_col, categories, output_type) {
  if (is.null(measure_df) || nrow(measure_df) == 0) {
    return(tibble::tibble())
  }

  patient_hour <- measure_df |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data[[datetime_col]]),
      variable = as.character(.data[[category_col]]),
      value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) |>
    filter(.data$variable %in% categories, is.finite(.data$value)) |>
    inner_join(
      cohort |> select("hospitalization_id", "first_icu_in", "last_icu_out", "heat_definition", "heat_related_ohca"),
      by = "hospitalization_id",
      relationship = "many-to-many"
    ) |>
    mutate(
      icu_hour = floor(as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")))
    ) |>
    filter(.data$icu_hour %in% TRAJECTORY_HOURS, .data$recorded_dttm >= .data$first_icu_in, .data$recorded_dttm <= .data$last_icu_out) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$hospitalization_id, .data$icu_hour, .data$variable) |>
    summarise(patient_hour_value = median(.data$value, na.rm = TRUE), .groups = "drop")

  patient_hour |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$icu_hour, .data$variable) |>
    summarise(
      n_patients = n_distinct(.data$hospitalization_id),
      median_value = median(.data$patient_hour_value, na.rm = TRUE),
      q25_value = stats::quantile(.data$patient_hour_value, 0.25, na.rm = TRUE, names = FALSE),
      q75_value = stats::quantile(.data$patient_hour_value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) |>
    mutate(output_type = output_type)
}

summarize_smoothed_hourly_measurements <- function(measure_df, cohort, datetime_col, category_col, value_col, categories, output_type, half_window = SMOOTH_HALF_WINDOW_HOURS) {
  if (is.null(measure_df) || nrow(measure_df) == 0) {
    return(tibble::tibble())
  }

  patient_hour <- measure_df |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data[[datetime_col]]),
      variable = as.character(.data[[category_col]]),
      value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) |>
    filter(.data$variable %in% categories, is.finite(.data$value)) |>
    inner_join(
      cohort |> select("hospitalization_id", "first_icu_in", "last_icu_out", "heat_definition", "heat_related_ohca"),
      by = "hospitalization_id",
      relationship = "many-to-many"
    ) |>
    mutate(
      observed_hour = floor(as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")))
    ) |>
    filter(.data$observed_hour %in% TRAJECTORY_HOURS, .data$recorded_dttm >= .data$first_icu_in, .data$recorded_dttm <= .data$last_icu_out) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$hospitalization_id, .data$observed_hour, .data$variable) |>
    summarise(patient_hour_value = median(.data$value, na.rm = TRUE), .groups = "drop")

  if (nrow(patient_hour) == 0) return(tibble::tibble())

  patient_window <- patient_hour |>
    tidyr::crossing(icu_hour = TRAJECTORY_HOURS) |>
    filter(abs(.data$observed_hour - .data$icu_hour) <= half_window) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$hospitalization_id, .data$icu_hour, .data$variable) |>
    summarise(patient_window_value = median(.data$patient_hour_value, na.rm = TRUE), .groups = "drop")

  patient_window |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$icu_hour, .data$variable) |>
    summarise(
      n_patients = n_distinct(.data$hospitalization_id),
      median_value = median(.data$patient_window_value, na.rm = TRUE),
      q25_value = stats::quantile(.data$patient_window_value, 0.25, na.rm = TRUE, names = FALSE),
      q75_value = stats::quantile(.data$patient_window_value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) |>
    mutate(
      output_type = output_type,
      smoothing = paste0("centered_", (2L * half_window + 1L), "h_window")
    )
}

summarize_hourly_support <- function(events, cohort, denominators, event_name) {
  if (is.null(events) || nrow(events) == 0) {
    event_counts <- tibble::tibble(
      heat_definition = character(),
      heat_related_ohca = character(),
      icu_hour = integer(),
      variable = character(),
      n_event = integer()
    )
  } else {
    event_counts <- events |>
      inner_join(
        cohort |> select("hospitalization_id", "first_icu_in", "last_icu_out", "heat_definition", "heat_related_ohca"),
        by = "hospitalization_id",
        relationship = "many-to-many"
      ) |>
      mutate(icu_hour = floor(as.numeric(difftime(.data$event_dttm, .data$first_icu_in, units = "hours")))) |>
      filter(.data$icu_hour %in% TRAJECTORY_HOURS, .data$event_dttm >= .data$first_icu_in, .data$event_dttm <= .data$last_icu_out) |>
      distinct(.data$heat_definition, .data$heat_related_ohca, .data$hospitalization_id, .data$icu_hour) |>
      group_by(.data$heat_definition, .data$heat_related_ohca, .data$icu_hour) |>
      summarise(n_event = n_distinct(.data$hospitalization_id), .groups = "drop") |>
      mutate(variable = event_name)
  }

  denominators |>
    mutate(variable = event_name) |>
    left_join(event_counts, by = c("heat_definition", "heat_related_ohca", "icu_hour", "variable")) |>
    mutate(
      n_event = tidyr::replace_na(.data$n_event, 0L),
      prevalence_pct = 100 * .data$n_event / .data$n_at_risk
    )
}

smooth_support_trajectory <- function(df, half_window = SMOOTH_HALF_WINDOW_HOURS) {
  if (nrow(df) == 0) return(df)
  df |>
    tidyr::crossing(smooth_hour = TRAJECTORY_HOURS) |>
    filter(abs(.data$icu_hour - .data$smooth_hour) <= half_window) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$variable, icu_hour = .data$smooth_hour) |>
    summarise(
      n_at_risk = sum(.data$n_at_risk, na.rm = TRUE),
      n_event = sum(.data$n_event, na.rm = TRUE),
      prevalence_pct = 100 * .data$n_event / .data$n_at_risk,
      .groups = "drop"
    ) |>
    mutate(smoothing = paste0("centered_", (2L * half_window + 1L), "h_window"))
}

make_cumulative_events <- function(cohort, events, event_name, hours = CUMULATIVE_HOURS) {
  base <- cohort |>
    select("hospitalization_id", "heat_definition", "heat_related_ohca") |>
    group_by(.data$heat_definition, .data$heat_related_ohca) |>
    summarise(n_group = n_distinct(.data$hospitalization_id), .groups = "drop") |>
    tidyr::crossing(icu_hour = hours) |>
    mutate(event = event_name)

  if (is.null(events) || nrow(events) == 0) {
    return(base |> mutate(n_events = 0L, cumulative_pct = 0))
  }

  event_times <- events |>
    inner_join(
      cohort |> select("hospitalization_id", "first_icu_in", "heat_definition", "heat_related_ohca"),
      by = "hospitalization_id",
      relationship = "many-to-many"
    ) |>
    mutate(event_hour = pmax(0L, floor(as.numeric(difftime(.data$event_dttm, .data$first_icu_in, units = "hours"))))) |>
    filter(!is.na(.data$event_hour)) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$hospitalization_id) |>
    summarise(event_hour = min(.data$event_hour), .groups = "drop")

  base |>
    left_join(event_times, by = c("heat_definition", "heat_related_ohca"), relationship = "many-to-many") |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$event, .data$icu_hour, .data$n_group) |>
    summarise(n_events = n_distinct(.data$hospitalization_id[!is.na(.data$event_hour) & .data$event_hour <= .data$icu_hour]), .groups = "drop") |>
    mutate(cumulative_pct = 100 * .data$n_events / .data$n_group)
}

summarize_renal_marker <- function(df, label, variable, direction = c("max", "min")) {
  direction <- match.arg(direction)
  empty <- tibble::tibble(
    heat_definition = character(),
    heat_related_ohca = character(),
    window = character(),
    marker = character(),
    direction = character(),
    n_patients = integer(),
    median_value = numeric(),
    q25_value = numeric(),
    q75_value = numeric()
  )
  if (is.null(df) || nrow(df) == 0 || !"variable" %in% names(df)) return(empty)
  marker <- df |>
    filter(.data$variable == !!variable, !is.na(.data$event_hour), !is.na(.data$value)) |>
    mutate(window = dplyr::case_when(
      .data$event_hour >= 0 & .data$event_hour <= 24 ~ "0-24h",
      .data$event_hour > 24 & .data$event_hour <= 72 ~ "24-72h",
      TRUE ~ NA_character_
    )) |>
    filter(!is.na(.data$window)) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$hospitalization_id, .data$window) |>
    summarise(
      value = if (direction == "max") max(.data$value, na.rm = TRUE) else min(.data$value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$window) |>
    summarise(
      marker = label,
      direction = direction,
      n_patients = n_distinct(.data$hospitalization_id),
      median_value = median(.data$value, na.rm = TRUE),
      q25_value = stats::quantile(.data$value, 0.25, na.rm = TRUE, names = FALSE),
      q75_value = stats::quantile(.data$value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
  marker
}

summarize_crrt_windows <- function(cohort, crrt_first_events) {
  base <- cohort |>
    select("hospitalization_id", "heat_definition", "heat_related_ohca") |>
    group_by(.data$heat_definition, .data$heat_related_ohca) |>
    summarise(n_group = n_distinct(.data$hospitalization_id), .groups = "drop")

  if (is.null(crrt_first_events) || nrow(crrt_first_events) == 0) {
    return(base |>
      tidyr::crossing(window = c("0-24h", "0-72h", "0-168h")) |>
      mutate(n_crrt = 0L, crrt_pct = 0))
  }

  event_times <- crrt_first_events |>
    inner_join(
      cohort |> select("hospitalization_id", "first_icu_in", "heat_definition", "heat_related_ohca"),
      by = "hospitalization_id",
      relationship = "many-to-many"
    ) |>
    mutate(event_hour = pmax(0L, floor(as.numeric(difftime(.data$event_dttm, .data$first_icu_in, units = "hours"))))) |>
    filter(!is.na(.data$event_hour)) |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$hospitalization_id) |>
    summarise(event_hour = min(.data$event_hour), .groups = "drop")

  dplyr::bind_rows(
    tibble::tibble(window = "0-24h", max_hour = 24L),
    tibble::tibble(window = "0-72h", max_hour = 72L),
    tibble::tibble(window = "0-168h", max_hour = 168L)
  ) |>
    tidyr::crossing(base) |>
    left_join(event_times, by = c("heat_definition", "heat_related_ohca"), relationship = "many-to-many") |>
    group_by(.data$heat_definition, .data$heat_related_ohca, .data$window, .data$n_group) |>
    summarise(
      n_crrt = n_distinct(.data$hospitalization_id[!is.na(.data$event_hour) & .data$event_hour <= .data$max_hour]),
      .groups = "drop"
    ) |>
    mutate(crrt_pct = 100 * .data$n_crrt / .data$n_group)
}

plot_renal_marker_summary <- function(df, heat_definition, filename, title) {
  plot_df <- df |> filter(.data$heat_definition == heat_definition)
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = .data$window, y = .data$median_value, color = .data$heat_related_ohca, group = .data$heat_related_ohca)) +
    geom_errorbar(aes(ymin = .data$q25_value, ymax = .data$q75_value), position = position_dodge(width = 0.35), width = 0.15) +
    geom_point(position = position_dodge(width = 0.35), size = 2.3) +
    facet_wrap(~ marker, scales = "free_y") +
    scale_color_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(title = title, x = NULL, y = "Median [IQR]", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 10, height = 6.5, dpi = 300)
}

plot_crrt_window_summary <- function(df, heat_definition, filename, title) {
  plot_df <- df |> filter(.data$heat_definition == heat_definition)
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = .data$window, y = .data$crrt_pct, fill = .data$heat_related_ohca)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(
      aes(label = paste0(.data$n_crrt, "/", .data$n_group)),
      position = position_dodge(width = 0.75),
      vjust = -0.25,
      size = 3
    ) +
    scale_fill_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(title = title, x = NULL, y = "CRRT cumulative incidence (%)", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 8, height = 5.5, dpi = 300)
}

plot_crrt_window_comparison <- function(df, filename) {
  if (nrow(df) == 0) return(invisible(NULL))
  plot_df <- df |>
    mutate(
      heat_definition = factor(.data$heat_definition, levels = c("heat90", "heat95")),
      window = factor(.data$window, levels = c("0-24h", "0-72h", "0-168h"))
    )

  p <- ggplot(plot_df, aes(x = .data$window, y = .data$crrt_pct, fill = .data$heat_related_ohca)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(
      aes(label = paste0(.data$n_crrt, "/", .data$n_group)),
      position = position_dodge(width = 0.75),
      vjust = -0.25,
      size = 3
    ) +
    facet_wrap(~ heat_definition) +
    scale_fill_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(
      title = "CRRT Initiation Windows by Heat Definition",
      subtitle = "Heat95 is nested within heat90; labels show n CRRT / group N",
      x = NULL,
      y = "CRRT cumulative incidence (%)",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 10, height = 5.8, dpi = 300)
}

plot_measure_trajectory <- function(df, heat_definition, output_type, filename, title) {
  plot_df <- df |>
    filter(.data$heat_definition == heat_definition, .data$output_type == output_type, !is.na(.data$median_value))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = .data$icu_hour, y = .data$median_value, color = .data$heat_related_ohca, fill = .data$heat_related_ohca)) +
    geom_ribbon(aes(ymin = .data$q25_value, ymax = .data$q75_value), alpha = 0.16, color = NA) +
    geom_line(linewidth = 0.85) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    scale_fill_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(title = title, x = "Hours since ICU entry", y = "Median value", color = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 10, height = 6.5, dpi = 300)
}

plot_smoothed_measure_trajectory <- function(df, heat_definition, output_type, filename, title) {
  plot_df <- df |>
    filter(.data$heat_definition == heat_definition, .data$output_type == output_type, !is.na(.data$median_value))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = .data$icu_hour, y = .data$median_value, color = .data$heat_related_ohca, fill = .data$heat_related_ohca)) +
    geom_ribbon(aes(ymin = .data$q25_value, ymax = .data$q75_value), alpha = 0.14, color = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    scale_fill_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(
      title = title,
      subtitle = "Centered 7-hour rolling windows",
      x = "Hours since ICU entry",
      y = "Rolling median value",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 10, height = 6.5, dpi = 300)
}

plot_support_trajectory <- function(df, heat_definition, filename, title) {
  plot_df <- df |> filter(.data$heat_definition == heat_definition, !is.na(.data$prevalence_pct))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = .data$icu_hour, y = .data$prevalence_pct, color = .data$heat_related_ohca)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ variable) +
    scale_color_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(title = title, x = "Hours since ICU entry", y = "Prevalence among patients still in ICU (%)", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 10, height = 5.8, dpi = 300)
}

plot_smoothed_support_trajectory <- function(df, heat_definition, filename, title) {
  plot_df <- df |> filter(.data$heat_definition == heat_definition, !is.na(.data$prevalence_pct))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = .data$icu_hour, y = .data$prevalence_pct, color = .data$heat_related_ohca)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ variable) +
    scale_color_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(
      title = title,
      subtitle = "Centered 7-hour rolling windows",
      x = "Hours since ICU entry",
      y = "Rolling prevalence among patients still in ICU (%)",
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 10, height = 5.8, dpi = 300)
}

plot_cumulative_incidence <- function(df, heat_definition, filename, title) {
  plot_df <- df |> filter(.data$heat_definition == heat_definition, !is.na(.data$cumulative_pct))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = .data$icu_hour, y = .data$cumulative_pct, color = .data$heat_related_ohca)) +
    geom_step(linewidth = 0.9) +
    facet_wrap(~ event, scales = "free_y") +
    scale_color_manual(values = c("Heat-related OHCA" = "#BC3908", "Non-heat-related OHCA" = "#335C67")) +
    labs(title = title, x = "Hours since ICU entry", y = "Crude cumulative incidence (%)", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURE_DIR, filename), p, width = 10, height = 6.5, dpi = 300)
}

ohca_raw <- readr::read_csv(OHCA_PATH, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_guess()))
if (!"hospital_death" %in% names(ohca_raw)) {
  ohca_raw$hospital_death <- ifelse(ohca_raw$discharge_category == "Expired", 1L, 0L)
}
if (!"death_or_hospice" %in% names(ohca_raw)) {
  ohca_raw$death_or_hospice <- ifelse(
    ohca_raw$discharge_category == "Expired" |
      stringr::str_detect(stringr::str_to_lower(tidyr::replace_na(ohca_raw$discharge_category, "")), "hospice"),
    1L,
    0L
  )
}
if (!"imv_any" %in% names(ohca_raw)) ohca_raw$imv_any <- NA_integer_
if (!"imv_duration_hours" %in% names(ohca_raw)) ohca_raw$imv_duration_hours <- NA_real_
if (!"vasopressor_any" %in% names(ohca_raw)) ohca_raw$vasopressor_any <- NA_integer_

ohca <- ohca_raw |>
  mutate(
    hospitalization_id = as.character(.data$hospitalization_id),
    admission_date = as.Date(.data$admission_date),
    admission_dttm = as_utc_datetime(.data$admission_dttm),
    discharge_dttm = as_utc_datetime(.data$discharge_dttm),
    first_icu_in = as_utc_datetime(.data$first_icu_in),
    last_icu_out = as_utc_datetime(.data$last_icu_out),
    age_at_admission = suppressWarnings(as.numeric(.data$age_at_admission)),
    hospital_los_days = as.numeric(difftime(.data$discharge_dttm, .data$admission_dttm, units = "days")),
    hospital_death = suppressWarnings(as.integer(.data$hospital_death)),
    death_or_hospice = suppressWarnings(as.integer(.data$death_or_hospice)),
    imv_any = suppressWarnings(as.integer(.data$imv_any)),
    imv_duration_hours = suppressWarnings(as.numeric(.data$imv_duration_hours)),
    vasopressor_any = suppressWarnings(as.integer(.data$vasopressor_any)),
    county_fips = normalize_county_fips(.data$county_fips),
    age_group = ifelse(.data$age_at_admission >= 65, ">=65", "<65"),
    race_group = ifelse(.data$race_category == "Black or African American", "Black", "Non-Black"),
    icu_los_days = suppressWarnings(as.numeric(.data$icu_los_hours)) / 24,
    imv_duration_days = suppressWarnings(as.numeric(.data$imv_duration_hours)) / 24
  )

tmax <- arrow::read_parquet(TMAX_PATH) |>
  transmute(
    county_fips = normalize_county_fips(.data$geoid),
    admission_date = as.Date(.data$date),
    tmax_mean_c = suppressWarnings(as.numeric(.data$tmax_mean_c))
  ) |>
  filter(as.integer(format(.data$admission_date, "%m")) %in% WARM_MONTHS)

ohca_exposed <- ohca |>
  left_join(tmax, by = c("county_fips", "admission_date")) |>
  filter(!is.na(.data$tmax_mean_c))

cohort_all_defs <- dplyr::bind_rows(lapply(seq_len(nrow(HEAT_DEFINITIONS)), function(i) {
  make_heat_groups(ohca_exposed, HEAT_DEFINITIONS$heat_definition[[i]], HEAT_DEFINITIONS$heat_percentile[[i]])
}))

thresholds <- cohort_all_defs |>
  distinct(.data$heat_definition, .data$heat_percentile, .data$heat_threshold_tmax_c)

summary_df <- summarize_group(cohort_all_defs)
table2_all <- dplyr::bind_rows(lapply(split(cohort_all_defs, cohort_all_defs$heat_definition), make_table2))
table2_primary <- table2_all |> filter(.data$heat_definition == "heat95") |> select(-"heat_definition")
table2_sensitivity <- table2_all |> filter(.data$heat_definition == "heat90") |> select(-"heat_definition")

by_discharge <- cohort_all_defs |>
  count(.data$heat_definition, .data$heat_related_ohca, .data$discharge_category, name = "n") |>
  group_by(.data$heat_definition, .data$heat_related_ohca) |>
  mutate(pct = safe_pct(.data$n, sum(.data$n))) |>
  ungroup()

denominators <- hourly_denominators(cohort_all_defs, TRAJECTORY_HOURS)

vitals <- read_clif_table(
  tables_path, file_type, "vitals",
  columns = c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value"),
  required = FALSE
)
vital_trajectories <- summarize_hourly_measurements(
  vitals,
  cohort_all_defs,
  datetime_col = "recorded_dttm",
  category_col = "vital_category",
  value_col = "vital_value",
  categories = c("heart_rate", "map", "sbp", "spo2", "respiratory_rate", "temp_c"),
  output_type = "vitals"
)
vital_trajectories_smoothed <- summarize_smoothed_hourly_measurements(
  vitals,
  cohort_all_defs,
  datetime_col = "recorded_dttm",
  category_col = "vital_category",
  value_col = "vital_value",
  categories = c("heart_rate", "map", "sbp", "spo2", "respiratory_rate", "temp_c"),
  output_type = "vitals"
)

labs <- read_clif_table(
  tables_path, file_type, "labs",
  columns = c("hospitalization_id", "lab_collect_dttm", "lab_result_dttm", "lab_category", "lab_value_numeric"),
  required = FALSE
)
if (!is.null(labs) && nrow(labs) > 0) {
  labs$trajectory_dttm <- ifelse(is.na(labs$lab_collect_dttm), labs$lab_result_dttm, labs$lab_collect_dttm)
  labs$trajectory_dttm <- as_utc_datetime(labs$trajectory_dttm)
}
renal_lab_values <- if (!is.null(labs) && nrow(labs) > 0) {
  labs |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      recorded_dttm = as_utc_datetime(.data$trajectory_dttm),
      variable = as.character(.data$lab_category),
      value = suppressWarnings(as.numeric(.data$lab_value_numeric))
    ) |>
    filter(.data$variable %in% c("creatinine", "bun", "potassium", "bicarbonate", "phosphate", "magnesium", "lactate"), is.finite(.data$value)) |>
    inner_join(
      cohort_all_defs |> select("hospitalization_id", "first_icu_in", "last_icu_out", "heat_definition", "heat_related_ohca"),
      by = "hospitalization_id",
      relationship = "many-to-many"
    ) |>
    mutate(event_hour = floor(as.numeric(difftime(.data$recorded_dttm, .data$first_icu_in, units = "hours")))) |>
    filter(.data$event_hour >= 0, .data$event_hour <= 72, .data$recorded_dttm >= .data$first_icu_in, .data$recorded_dttm <= .data$last_icu_out)
} else {
  tibble::tibble()
}
renal_marker_summary <- dplyr::bind_rows(
  summarize_renal_marker(renal_lab_values, "Peak creatinine", "creatinine", "max"),
  summarize_renal_marker(renal_lab_values, "Peak BUN", "bun", "max"),
  summarize_renal_marker(renal_lab_values, "Peak potassium", "potassium", "max"),
  summarize_renal_marker(renal_lab_values, "Lowest bicarbonate", "bicarbonate", "min"),
  summarize_renal_marker(renal_lab_values, "Peak phosphate", "phosphate", "max"),
  summarize_renal_marker(renal_lab_values, "Peak magnesium", "magnesium", "max"),
  summarize_renal_marker(renal_lab_values, "Peak lactate", "lactate", "max")
)
lab_trajectories <- summarize_hourly_measurements(
  labs,
  cohort_all_defs,
  datetime_col = "trajectory_dttm",
  category_col = "lab_category",
  value_col = "lab_value_numeric",
  categories = c("lactate", "ph_arterial", "creatinine", "bun", "bicarbonate", "potassium", "phosphate", "magnesium", "wbc", "troponin_t"),
  output_type = "labs"
)
lab_trajectories_smoothed <- summarize_smoothed_hourly_measurements(
  labs,
  cohort_all_defs,
  datetime_col = "trajectory_dttm",
  category_col = "lab_category",
  value_col = "lab_value_numeric",
  categories = c("lactate", "ph_arterial", "creatinine", "bun", "bicarbonate", "potassium", "phosphate", "magnesium", "wbc", "troponin_t"),
  output_type = "labs"
)

respiratory <- read_clif_table(
  tables_path, file_type, "respiratory_support",
  columns = c("hospitalization_id", "recorded_dttm", "device_category", "artificial_airway", "tracheostomy"),
  required = FALSE
)
imv_events <- NULL
if (!is.null(respiratory) && nrow(respiratory) > 0) {
  imv_events <- respiratory |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      event_dttm = as_utc_datetime(.data$recorded_dttm),
      imv_row =
        stringr::str_to_upper(tidyr::replace_na(as.character(.data$device_category), "")) == "IMV" |
        stringr::str_to_upper(tidyr::replace_na(as.character(.data$artificial_airway), "")) %in% c("ETT", "TRACH") |
        tidyr::replace_na(suppressWarnings(as.numeric(.data$tracheostomy)), 0) == 1
    ) |>
    filter(.data$imv_row, !is.na(.data$event_dttm)) |>
    select("hospitalization_id", "event_dttm")
}

medication <- read_clif_table(
  tables_path, file_type, "medication_admin_continuous",
  columns = c("hospitalization_id", "admin_dttm", "med_category", "med_group", "mar_action_group"),
  required = FALSE
)
vasopressor_events <- NULL
if (!is.null(medication) && nrow(medication) > 0) {
  vasoactive_categories <- c("angiotensin", "dobutamine", "dopamine", "epinephrine", "milrinone", "norepinephrine", "phenylephrine", "vasopressin")
  vasopressor_events <- medication |>
    transmute(
      hospitalization_id = as.character(.data$hospitalization_id),
      event_dttm = as_utc_datetime(.data$admin_dttm),
      med_category = stringr::str_to_lower(as.character(.data$med_category)),
      med_group = stringr::str_to_lower(as.character(.data$med_group)),
      mar_action_group = stringr::str_to_lower(as.character(.data$mar_action_group))
    ) |>
    filter(.data$med_group == "vasoactives", .data$med_category %in% vasoactive_categories, .data$mar_action_group != "not_administered", !is.na(.data$event_dttm)) |>
    select("hospitalization_id", "event_dttm")
}

first_event <- function(events) {
  if (is.null(events) || nrow(events) == 0) return(NULL)
  events |>
    group_by(.data$hospitalization_id) |>
    summarise(event_dttm = min(.data$event_dttm, na.rm = TRUE), .groups = "drop")
}

crrt <- read_clif_table(
  tables_path, file_type, "crrt_therapy",
  columns = c("hospitalization_id", "recorded_dttm"),
  required = FALSE
)
crrt_events <- if (!is.null(crrt) && nrow(crrt) > 0) {
  crrt |> transmute(hospitalization_id = as.character(.data$hospitalization_id), event_dttm = as_utc_datetime(.data$recorded_dttm)) |> filter(!is.na(.data$event_dttm))
} else {
  NULL
}
crrt_first_events <- first_event(crrt_events)
crrt_window_summary <- summarize_crrt_windows(cohort_all_defs, crrt_first_events)

support_trajectories <- dplyr::bind_rows(
  summarize_hourly_support(imv_events, cohort_all_defs, denominators, "Invasive mechanical ventilation"),
  summarize_hourly_support(vasopressor_events, cohort_all_defs, denominators, "Vasopressor infusion"),
  summarize_hourly_support(crrt_events, cohort_all_defs, denominators, "CRRT")
)
support_trajectories_smoothed <- smooth_support_trajectory(support_trajectories)

death_or_hospice_events <- ohca_exposed |>
  filter(.data$death_or_hospice == 1) |>
  transmute(hospitalization_id = .data$hospitalization_id, event_dttm = .data$discharge_dttm)

hospital_death_events <- ohca_exposed |>
  filter(.data$hospital_death == 1) |>
  transmute(hospitalization_id = .data$hospitalization_id, event_dttm = .data$discharge_dttm)

alive_discharge_events <- ohca_exposed |>
  filter(.data$death_or_hospice == 0) |>
  transmute(hospitalization_id = .data$hospitalization_id, event_dttm = .data$discharge_dttm)

cumulative_incidence <- dplyr::bind_rows(
  make_cumulative_events(cohort_all_defs, death_or_hospice_events, "Death or hospice"),
  make_cumulative_events(cohort_all_defs, hospital_death_events, "Hospital death"),
  make_cumulative_events(cohort_all_defs, alive_discharge_events, "Discharged alive without hospice"),
  make_cumulative_events(cohort_all_defs, first_event(imv_events), "First IMV"),
  make_cumulative_events(cohort_all_defs, first_event(vasopressor_events), "First vasopressor"),
  make_cumulative_events(cohort_all_defs, crrt_first_events, "First CRRT")
)

readr::write_csv(thresholds, file.path(OUTPUT_DIR, "heat_related_ohca_thresholds.csv"))
readr::write_csv(thresholds |> filter(.data$heat_definition == "heat95") |> transmute(heat_threshold_tmax_c = .data$heat_threshold_tmax_c), file.path(OUTPUT_DIR, "heat_related_ohca_threshold.csv"))
readr::write_csv(summary_df, file.path(OUTPUT_DIR, "heat_related_vs_non_heat_related_ohca_outcomes.csv"))
readr::write_csv(by_discharge, file.path(OUTPUT_DIR, "heat_related_vs_non_heat_related_discharge_categories.csv"))
readr::write_csv(table2_primary, file.path(OUTPUT_DIR, "table2_heat_related_vs_non_heat_related_ohca.csv"))
readr::write_csv(table2_all, file.path(OUTPUT_DIR, "table2_heat_related_vs_non_heat_related_ohca_all_definitions.csv"))
readr::write_csv(table2_sensitivity, file.path(OUTPUT_DIR, "table2_heat90_vs_non_heat90_ohca.csv"))
readr::write_csv(vital_trajectories, file.path(OUTPUT_DIR, "heat_related_ohca_hourly_vital_trajectories.csv"))
readr::write_csv(lab_trajectories, file.path(OUTPUT_DIR, "heat_related_ohca_hourly_lab_trajectories.csv"))
readr::write_csv(support_trajectories, file.path(OUTPUT_DIR, "heat_related_ohca_hourly_support_trajectories.csv"))
readr::write_csv(vital_trajectories_smoothed, file.path(OUTPUT_DIR, "heat_related_ohca_hourly_vital_trajectories_smoothed.csv"))
readr::write_csv(lab_trajectories_smoothed, file.path(OUTPUT_DIR, "heat_related_ohca_hourly_lab_trajectories_smoothed.csv"))
readr::write_csv(support_trajectories_smoothed, file.path(OUTPUT_DIR, "heat_related_ohca_hourly_support_trajectories_smoothed.csv"))
readr::write_csv(cumulative_incidence, file.path(OUTPUT_DIR, "heat_related_ohca_hourly_cumulative_incidence.csv"))
readr::write_csv(renal_marker_summary, file.path(OUTPUT_DIR, "heat_related_ohca_renal_metabolic_marker_summary.csv"))
readr::write_csv(crrt_window_summary, file.path(OUTPUT_DIR, "heat_related_ohca_crrt_window_summary.csv"))

for (heat_def in HEAT_DEFINITIONS$heat_definition) {
  plot_measure_trajectory(vital_trajectories, heat_def, "vitals", paste0("figure_", heat_def, "_hourly_vital_trajectories.png"), paste0("Hourly Vital Sign Trajectories: ", heat_def))
  plot_measure_trajectory(lab_trajectories, heat_def, "labs", paste0("figure_", heat_def, "_hourly_lab_trajectories.png"), paste0("Hourly Laboratory Trajectories: ", heat_def))
  plot_support_trajectory(support_trajectories, heat_def, paste0("figure_", heat_def, "_hourly_support_trajectories.png"), paste0("Hourly ICU Support Trajectories: ", heat_def))
  plot_smoothed_measure_trajectory(vital_trajectories_smoothed, heat_def, "vitals", paste0("figure_", heat_def, "_hourly_vital_trajectories_smoothed.png"), paste0("Smoothed Hourly Vital Sign Trajectories: ", heat_def))
  plot_smoothed_measure_trajectory(lab_trajectories_smoothed, heat_def, "labs", paste0("figure_", heat_def, "_hourly_lab_trajectories_smoothed.png"), paste0("Smoothed Hourly Laboratory Trajectories: ", heat_def))
  plot_smoothed_support_trajectory(support_trajectories_smoothed, heat_def, paste0("figure_", heat_def, "_hourly_support_trajectories_smoothed.png"), paste0("Smoothed Hourly ICU Support Trajectories: ", heat_def))
  plot_cumulative_incidence(cumulative_incidence, heat_def, paste0("figure_", heat_def, "_cumulative_incidence.png"), paste0("Cumulative Incidence After ICU Entry: ", heat_def))
  plot_renal_marker_summary(renal_marker_summary, heat_def, paste0("figure_", heat_def, "_renal_metabolic_marker_summary.png"), paste0("Renal and Metabolic Marker Summary: ", heat_def))
  plot_crrt_window_summary(crrt_window_summary, heat_def, paste0("figure_", heat_def, "_crrt_window_summary.png"), paste0("CRRT Initiation Windows: ", heat_def))
}
plot_crrt_window_comparison(crrt_window_summary, "figure_crrt_window_heat_definition_comparison.png")

message("Wrote heat-related OHCA clinical phenotype outputs to ", OUTPUT_DIR)

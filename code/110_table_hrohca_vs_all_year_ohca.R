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
input_dir <- file.path(repo_root, "output", "final", "federated_pooled")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(openxlsx)
  library(readr)
  library(stringr)
  library(tidyr)
})

parse_count_pct <- function(x) {
  x <- str_replace_all(as.character(x), ",", "")
  nums <- str_extract_all(x, "-?[0-9]+[.]?[0-9]*")
  count <- suppressWarnings(as.numeric(vapply(nums, function(z) if (length(z) >= 1) z[[1]] else NA_character_, character(1))))
  pct <- suppressWarnings(as.numeric(vapply(seq_along(nums), function(i) {
    if (str_detect(x[[i]], "%") && length(nums[[i]]) >= 2) nums[[i]][[length(nums[[i]])]] else NA_character_
  }, character(1))))
  tibble(count = count, pct = pct)
}

parse_median_iqr <- function(x) {
  x <- str_replace_all(as.character(x), ",", "")
  nums <- str_extract_all(x, "-?[0-9]+[.]?[0-9]*")
  mat <- do.call(rbind, lapply(nums, function(z) {
    if (length(z) >= 3) z[1:3] else rep(NA_character_, 3)
  }))
  tibble(
    median = suppressWarnings(as.numeric(mat[, 1])),
    q1 = suppressWarnings(as.numeric(mat[, 2])),
    q3 = suppressWarnings(as.numeric(mat[, 3]))
  )
}

weighted_quantile <- function(values, weights, probs) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(keep)) return(rep(NA_real_, length(probs)))
  values <- values[keep]
  weights <- weights[keep]
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord]
  cum_w <- cumsum(weights) / sum(weights)
  vapply(probs, function(p) values[which(cum_w >= p)[1]], numeric(1))
}

fmt_count_pct <- function(count, denom) {
  if (!is.finite(count) || !is.finite(denom) || denom <= 0) return("")
  sprintf("%s (%.1f%%)", format(round(count), big.mark = ","), 100 * count / denom)
}

fmt_median_iqr <- function(median, q1, q3, digits = 1) {
  if (!all(is.finite(c(median, q1, q3)))) return("")
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), median, q1, q3)
}

fmt_diff <- function(x, suffix = "") {
  if (!is.finite(x)) return("")
  sprintf("%.1f%s", x, suffix)
}

primary <- read.csv(
  file.path(table_dir, "table1_pooled_hrohca_vs_non_hrohca_source.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all_year <- read.csv(
  file.path(input_dir, "all_sites_all_year_heat_related_vs_non_heat_related_table_all_definitions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  filter(.data$heat_definition == "heat95")

group_n <- all_year |>
  filter(.data$characteristic == "N") |>
  transmute(
    site_name,
    heat_n = parse_count_pct(.data$heat_related_ohca)$count,
    non_heat_n = parse_count_pct(.data$non_heat_related_ohca)$count,
    all_year_n = .data$heat_n + .data$non_heat_n
  )

total_all_year_n <- sum(group_n$all_year_n, na.rm = TRUE)

all_year_count <- function(characteristic) {
  all_year |>
    filter(.data$characteristic == !!characteristic) |>
    left_join(group_n, by = "site_name") |>
    transmute(
      heat_count = parse_count_pct(.data$heat_related_ohca)$count,
      non_heat_count = parse_count_pct(.data$non_heat_related_ohca)$count,
      count = .data$heat_count + .data$non_heat_count,
      denom = .data$all_year_n
    ) |>
    summarise(count = sum(.data$count, na.rm = TRUE), denom = sum(.data$denom, na.rm = TRUE), .groups = "drop")
}

all_year_continuous <- function(characteristic, denom_characteristic = NULL) {
  dat <- all_year |>
    filter(.data$characteristic == !!characteristic) |>
    left_join(group_n, by = "site_name")

  if (!is.null(denom_characteristic)) {
    denom_dat <- all_year |>
      filter(.data$characteristic == !!denom_characteristic) |>
      transmute(
        site_name,
        heat_n_override = parse_count_pct(.data$heat_related_ohca)$count,
        non_heat_n_override = parse_count_pct(.data$non_heat_related_ohca)$count
      )
    dat <- dat |>
      left_join(denom_dat, by = "site_name") |>
      mutate(
        heat_n = if_else(is.finite(.data$heat_n_override), .data$heat_n_override, .data$heat_n),
        non_heat_n = if_else(is.finite(.data$non_heat_n_override), .data$non_heat_n_override, .data$non_heat_n)
      )
  }

  parsed <- bind_cols(
    dat,
    parse_median_iqr(dat$heat_related_ohca) |> rename(heat_median = median, heat_q1 = q1, heat_q3 = q3),
    parse_median_iqr(dat$non_heat_related_ohca) |> rename(non_heat_median = median, non_heat_q1 = q1, non_heat_q3 = q3)
  )

  weighted_quantile(
    values = c(parsed$heat_q1, parsed$heat_median, parsed$heat_q3, parsed$non_heat_q1, parsed$non_heat_median, parsed$non_heat_q3),
    weights = c(0.25 * parsed$heat_n, 0.50 * parsed$heat_n, 0.25 * parsed$heat_n, 0.25 * parsed$non_heat_n, 0.50 * parsed$non_heat_n, 0.25 * parsed$non_heat_n),
    probs = c(0.25, 0.50, 0.75)
  )
}

primary_value <- function(section, characteristic, occurrence = 1L) {
  primary |>
    filter(.data$section == !!section, .data$characteristic == !!characteristic) |>
    slice(occurrence)
}

continuous_row <- function(section, display_characteristic, all_year_characteristic, digits = 1, denom_characteristic = NULL, primary_occurrence = 1L) {
  h <- primary_value(section, display_characteristic, primary_occurrence)
  hq <- parse_median_iqr(h$heat_related_ohca)
  aq <- all_year_continuous(all_year_characteristic, denom_characteristic = denom_characteristic)
  tibble(
    Section = section,
    Characteristic = display_characteristic,
    HROHCA = h$heat_related_ohca,
    `All-year OHCA` = fmt_median_iqr(aq[2], aq[1], aq[3], digits = digits),
    `Difference` = fmt_diff(hq$median[[1]] - aq[2]),
    `Difference unit` = "median difference"
  )
}

categorical_row <- function(section, display_characteristic, all_year_characteristic, primary_occurrence = 1L) {
  h <- primary_value(section, display_characteristic, primary_occurrence)
  hp <- parse_count_pct(h$heat_related_ohca)
  ac <- all_year_count(all_year_characteristic)
  all_pct <- 100 * ac$count / ac$denom
  tibble(
    Section = section,
    Characteristic = display_characteristic,
    HROHCA = h$heat_related_ohca,
    `All-year OHCA` = fmt_count_pct(ac$count, ac$denom),
    `Difference` = fmt_diff(hp$pct[[1]] - all_pct, " pp"),
    `Difference unit` = "percentage-point difference"
  )
}

count_row <- tibble(
  Section = "Cohort",
  Characteristic = "Patients, n",
  HROHCA = primary_value("Cohort", "Patients, n")$heat_related_ohca,
  `All-year OHCA` = format(total_all_year_n, big.mark = ","),
  `Difference` = "",
  `Difference unit` = ""
)

table_out <- bind_rows(
  count_row,
  continuous_row("Demographics", "Age, years", "Age, years"),
  categorical_row("Demographics", "<65 years", "Age group: <65"),
  categorical_row("Demographics", ">=65 years", "Age group: >=65"),
  categorical_row("Demographics", "Male", "Sex: Male"),
  categorical_row("Demographics", "Female", "Sex: Female"),
  categorical_row("Demographics", "Unknown", "Sex: Unknown", primary_occurrence = 1L),
  categorical_row("Demographics", "Black", "Race: Black"),
  categorical_row("Demographics", "Non-Black", "Race: Non-Black"),
  categorical_row("Demographics", "Hispanic", "Ethnicity: Hispanic"),
  categorical_row("Demographics", "Non-Hispanic", "Ethnicity: Non-Hispanic"),
  categorical_row("Demographics", "Unknown", "Ethnicity: Unknown", primary_occurrence = 2L),
  categorical_row("Geography and exposure", "County reassigned to hospital county", "County reassigned to hospital county"),
  continuous_row("Geography and exposure", "Assigned-county Tmax, C", "Assigned-county Tmax, C"),
  categorical_row("ICU therapies and outcomes", "Invasive mechanical ventilation", "Invasive mechanical ventilation"),
  continuous_row(
    "ICU therapies and outcomes",
    "IMV duration among ventilated, days",
    "IMV duration among ventilated, days",
    denom_characteristic = "Invasive mechanical ventilation"
  ),
  categorical_row("ICU therapies and outcomes", "Vasopressor use", "Vasopressor use"),
  continuous_row("ICU therapies and outcomes", "ICU length of stay, days", "ICU length of stay, days"),
  continuous_row("ICU therapies and outcomes", "Hospital length of stay, days", "Hospital length of stay, days"),
  categorical_row("ICU therapies and outcomes", "Hospital death", "Hospital death"),
  categorical_row("ICU therapies and outcomes", "Death or hospice", "Death or hospice"),
  continuous_row(
    "ICU therapies and outcomes",
    "Time to death among decedents, days",
    "Time to death among decedents, days",
    denom_characteristic = "Hospital death"
  )
)

readr::write_csv(table_out, file.path(table_dir, "supplement_table_hrohca_vs_all_year_ohca.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "HROHCA vs all-year OHCA")
openxlsx::writeData(wb, "HROHCA vs all-year OHCA", table_out)
openxlsx::setColWidths(wb, "HROHCA vs all-year OHCA", cols = 1:ncol(table_out), widths = "auto")
openxlsx::saveWorkbook(
  wb,
  file.path(table_dir, "supplement_table_hrohca_vs_all_year_ohca.xlsx"),
  overwrite = TRUE
)

ft <- flextable(table_out) |>
  autofit() |>
  align(align = "left", part = "all") |>
  theme_booktabs() |>
  bold(part = "header") |>
  add_footer_lines(
    values = paste(
      "HROHCA uses the primary warm-season heat95 definition.",
      "The all-year OHCA comparator pools all heat95 heat-related and non-heat-related OHCA across 2018-2024, so it includes HROHCA patients; differences are descriptive.",
      "Continuous all-year summaries are approximate and reconstructed from site-level medians and IQRs."
    )
  )

flextable::save_as_docx(
  `Supplemental Table. Primary HROHCA compared with all-year OHCA` = ft,
  path = file.path(table_dir, "supplement_table_hrohca_vs_all_year_ohca.docx")
)

message("Wrote HROHCA vs all-year OHCA comparison table to ", table_dir)

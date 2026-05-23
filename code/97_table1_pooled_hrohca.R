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
output_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(SASS_CACHE = file.path(tempdir(), "clif_hrohca_sass_cache"))
dir.create(Sys.getenv("SASS_CACHE"), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(gt)
  library(officer)
  library(openxlsx)
  library(readr)
  library(stringr)
})

HEAT_DEFINITION <- Sys.getenv("HEAT_DEFINITION", unset = "heat95")
OUTPUT_SUFFIX <- if (HEAT_DEFINITION == "heat95") "" else paste0("_", HEAT_DEFINITION)

source_file <- file.path(input_dir, "all_sites_heat_related_vs_non_heat_related_table_all_definitions.csv")
if (!file.exists(source_file)) {
  stop("Missing pooled site table file: ", source_file)
}

raw_table <- readr::read_csv(source_file, show_col_types = FALSE) |>
  filter(.data$heat_definition == HEAT_DEFINITION) |>
  mutate(
    characteristic = str_squish(.data$characteristic),
    characteristic = case_when(
      characteristic == "Sex: male" ~ "Sex: Male",
      characteristic == "Sex: female" ~ "Sex: Female",
      characteristic == "Ethnicity: non-hispanic" ~ "Ethnicity: Non-Hispanic",
      characteristic == "Ethnicity: hispanic" ~ "Ethnicity: Hispanic",
      characteristic == "Ethnicity: unknown" ~ "Ethnicity: Unknown",
      TRUE ~ characteristic
    )
  )

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

fmt_p <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

group_n <- raw_table |>
  filter(.data$characteristic == "N") |>
  transmute(
    site_name,
    heat_n = parse_count_pct(.data$heat_related_ohca)$count,
    non_heat_n = parse_count_pct(.data$non_heat_related_ohca)$count
  )

total_heat_n <- sum(group_n$heat_n, na.rm = TRUE)
total_non_heat_n <- sum(group_n$non_heat_n, na.rm = TRUE)
site_count <- n_distinct(raw_table$site_name)

count_row <- function(characteristic) {
  raw_table |>
    filter(.data$characteristic == !!characteristic) |>
    transmute(
      site_name,
      heat_events = parse_count_pct(.data$heat_related_ohca)$count,
      non_heat_events = parse_count_pct(.data$non_heat_related_ohca)$count
    )
}

continuous_row <- function(characteristic, heat_n_override = NULL, non_heat_n_override = NULL) {
  parsed_base <- raw_table |>
    filter(.data$characteristic == !!characteristic) |>
    left_join(group_n, by = "site_name")
  parsed <- bind_cols(
    parsed_base,
    parse_median_iqr(parsed_base$heat_related_ohca) |> rename(heat_median = median, heat_q1 = q1, heat_q3 = q3),
    parse_median_iqr(parsed_base$non_heat_related_ohca) |> rename(non_heat_median = median, non_heat_q1 = q1, non_heat_q3 = q3)
  )

  if (!is.null(heat_n_override)) parsed$heat_n <- heat_n_override(parsed)
  if (!is.null(non_heat_n_override)) parsed$non_heat_n <- non_heat_n_override(parsed)

  approx_quantiles <- function(q1, median, q3, n) {
    weighted_quantile(
      values = c(q1, median, q3),
      weights = c(0.25 * n, 0.50 * n, 0.25 * n),
      probs = c(0.25, 0.50, 0.75)
    )
  }
  hq <- approx_quantiles(parsed$heat_q1, parsed$heat_median, parsed$heat_q3, parsed$heat_n)
  nhq <- approx_quantiles(parsed$non_heat_q1, parsed$non_heat_median, parsed$non_heat_q3, parsed$non_heat_n)

  estimated_mean_sd <- function(q1, median, q3, n) {
    tibble(
      mean = (q1 + median + q3) / 3,
      sd = pmax((q3 - q1) / 1.349, .Machine$double.eps),
      n = n
    ) |>
      filter(is.finite(.data$mean), is.finite(.data$sd), is.finite(.data$n), .data$n > 1)
  }
  combine_mean_sd <- function(dat) {
    n_total <- sum(dat$n)
    mean_total <- sum(dat$n * dat$mean) / n_total
    var_total <- sum((dat$n - 1) * dat$sd^2 + dat$n * (dat$mean - mean_total)^2) / (n_total - 1)
    list(n = n_total, mean = mean_total, sd = sqrt(var_total))
  }

  h <- combine_mean_sd(estimated_mean_sd(parsed$heat_q1, parsed$heat_median, parsed$heat_q3, parsed$heat_n))
  nh <- combine_mean_sd(estimated_mean_sd(parsed$non_heat_q1, parsed$non_heat_median, parsed$non_heat_q3, parsed$non_heat_n))
  se <- sqrt(h$sd^2 / h$n + nh$sd^2 / nh$n)
  df <- (h$sd^2 / h$n + nh$sd^2 / nh$n)^2 /
    ((h$sd^2 / h$n)^2 / (h$n - 1) + (nh$sd^2 / nh$n)^2 / (nh$n - 1))
  p <- 2 * stats::pt(abs((h$mean - nh$mean) / se), df = df, lower.tail = FALSE)

  tibble(
    characteristic = characteristic,
    heat_related_ohca = fmt_median_iqr(hq[2], hq[1], hq[3]),
    non_heat_related_ohca = fmt_median_iqr(nhq[2], nhq[1], nhq[3]),
    p_value = fmt_p(p),
    p_value_numeric = p,
    method = "Continuous summaries and p values are approximate, reconstructed from site medians and IQRs."
  )
}

categorical_rows <- function(section, prefix, levels, labels = levels, p_from_composite = TRUE) {
  rows <- lapply(seq_along(levels), function(i) {
    dat <- count_row(paste0(prefix, ": ", levels[[i]]))
    tibble(
      section = section,
      characteristic = labels[[i]],
      heat_events = sum(dat$heat_events, na.rm = TRUE),
      non_heat_events = sum(dat$non_heat_events, na.rm = TRUE)
    )
  }) |>
    bind_rows()

  if (p_from_composite) {
    mat <- as.matrix(rows[, c("heat_events", "non_heat_events")])
    rownames(mat) <- rows$characteristic
    p <- suppressWarnings(stats::chisq.test(t(mat), simulate.p.value = any(mat < 5), B = 10000)$p.value)
  } else {
    p <- NA_real_
  }

  rows |>
    mutate(
      heat_related_ohca = mapply(fmt_count_pct, .data$heat_events, MoreArgs = list(denom = sum(.data$heat_events, na.rm = TRUE))),
      non_heat_related_ohca = mapply(fmt_count_pct, .data$non_heat_events, MoreArgs = list(denom = sum(.data$non_heat_events, na.rm = TRUE))),
      p_value = if_else(row_number() == 1L, fmt_p(p), ""),
      p_value_numeric = if_else(row_number() == 1L, p, NA_real_),
      method = "Categorical p values are chi-square tests using pooled site counts.",
      .before = .data$characteristic
    ) |>
    select("section", "characteristic", "heat_related_ohca", "non_heat_related_ohca", "p_value", "p_value_numeric", "method")
}

binary_row <- function(section, characteristic, source_characteristic) {
  dat <- count_row(source_characteristic)
  events <- tibble(
    heat_events = sum(dat$heat_events, na.rm = TRUE),
    non_heat_events = sum(dat$non_heat_events, na.rm = TRUE)
  )
  mat <- matrix(
    c(
      events$heat_events,
      total_heat_n - events$heat_events,
      events$non_heat_events,
      total_non_heat_n - events$non_heat_events
    ),
    nrow = 2,
    byrow = TRUE
  )
  p <- if (any(mat < 5)) stats::fisher.test(mat)$p.value else suppressWarnings(stats::chisq.test(mat, correct = FALSE)$p.value)
  tibble(
    section = section,
    characteristic = characteristic,
    heat_related_ohca = fmt_count_pct(events$heat_events, total_heat_n),
    non_heat_related_ohca = fmt_count_pct(events$non_heat_events, total_non_heat_n),
    p_value = fmt_p(p),
    p_value_numeric = p,
    method = "Binary p values are two-sided Fisher exact tests or chi-square tests using pooled site counts."
  )
}

events_for_subset <- function(source_characteristic, group) {
  counts <- count_row(source_characteristic)
  values <- if (group == "heat") counts$heat_events else counts$non_heat_events
  values[match(raw_table |> filter(.data$characteristic == "N") |> pull(.data$site_name), counts$site_name)]
}

table1 <- bind_rows(
  tibble(
    section = "Cohort",
    characteristic = "Patients, n",
    heat_related_ohca = format(total_heat_n, big.mark = ","),
    non_heat_related_ohca = format(total_non_heat_n, big.mark = ","),
    p_value = "",
    p_value_numeric = NA_real_,
    method = "Counts summed across contributing sites."
  ),
  continuous_row("Age, years") |> mutate(section = "Demographics", .before = "characteristic"),
  categorical_rows("Demographics", "Age group", c("<65", ">=65"), c("<65 years", ">=65 years")),
  categorical_rows("Demographics", "Sex", c("Male", "Female", "Unknown")),
  categorical_rows("Demographics", "Race", c("Black", "Non-Black")),
  categorical_rows("Demographics", "Ethnicity", c("Hispanic", "Non-Hispanic", "Unknown")),
  binary_row("Geography and exposure", "County reassigned to hospital county", "County reassigned to hospital county"),
  continuous_row("Assigned-county Tmax, C") |> mutate(section = "Geography and exposure", .before = "characteristic"),
  binary_row("ICU therapies and outcomes", "Invasive mechanical ventilation", "Invasive mechanical ventilation"),
  continuous_row(
    "IMV duration among ventilated, days",
    heat_n_override = function(x) events_for_subset("Invasive mechanical ventilation", "heat"),
    non_heat_n_override = function(x) events_for_subset("Invasive mechanical ventilation", "non_heat")
  ) |> mutate(section = "ICU therapies and outcomes", .before = "characteristic"),
  binary_row("ICU therapies and outcomes", "Vasopressor use", "Vasopressor use"),
  continuous_row("ICU length of stay, days") |> mutate(section = "ICU therapies and outcomes", .before = "characteristic"),
  continuous_row("Hospital length of stay, days") |> mutate(section = "ICU therapies and outcomes", .before = "characteristic"),
  binary_row("ICU therapies and outcomes", "Hospital death", "Hospital death"),
  binary_row("ICU therapies and outcomes", "Death or hospice", "Death or hospice"),
  continuous_row(
    "Time to death among decedents, days",
    heat_n_override = function(x) events_for_subset("Hospital death", "heat"),
    non_heat_n_override = function(x) events_for_subset("Hospital death", "non_heat")
  ) |> mutate(section = "ICU therapies and outcomes", .before = "characteristic")
) |>
  mutate(
    heat_related_ohca = if_else(is.na(.data$heat_related_ohca), "", .data$heat_related_ohca),
    non_heat_related_ohca = if_else(is.na(.data$non_heat_related_ohca), "", .data$non_heat_related_ohca)
  )

display_table <- table1 |>
  transmute(
    Section = .data$section,
    Characteristic = .data$characteristic,
    `Heat-related OHCA` = .data$heat_related_ohca,
    `Non-heat-related OHCA` = .data$non_heat_related_ohca,
    `P value` = .data$p_value
  )

source_table <- table1 |>
  mutate(
    heat_definition = HEAT_DEFINITION,
    heat_related_n = total_heat_n,
    non_heat_related_n = total_non_heat_n,
    contributing_sites = site_count,
    .before = "section"
  )

readr::write_csv(display_table, file.path(output_dir, paste0("table1_pooled_hrohca_vs_non_hrohca", OUTPUT_SUFFIX, ".csv")))
readr::write_csv(source_table, file.path(output_dir, paste0("table1_pooled_hrohca_vs_non_hrohca", OUTPUT_SUFFIX, "_source.csv")))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Table 1")
openxlsx::writeData(wb, "Table 1", display_table)
openxlsx::setColWidths(wb, "Table 1", cols = 1:ncol(display_table), widths = c(24, 34, 26, 30, 12))
openxlsx::addStyle(
  wb,
  "Table 1",
  style = openxlsx::createStyle(textDecoration = "bold", fgFill = "#E9EDF2", border = "Bottom"),
  rows = 1,
  cols = 1:ncol(display_table),
  gridExpand = TRUE
)
openxlsx::saveWorkbook(wb, file.path(output_dir, paste0("table1_pooled_hrohca_vs_non_hrohca", OUTPUT_SUFFIX, ".xlsx")), overwrite = TRUE)

ft <- flextable(display_table) |>
  merge_v(j = "Section") |>
  valign(j = "Section", valign = "top") |>
  bold(part = "header") |>
  bg(part = "header", bg = "#E9EDF2") |>
  align(j = c("Heat-related OHCA", "Non-heat-related OHCA", "P value"), align = "center", part = "all") |>
  width(j = "Section", width = 1.35) |>
  width(j = "Characteristic", width = 2.35) |>
  width(j = c("Heat-related OHCA", "Non-heat-related OHCA"), width = 1.55) |>
  width(j = "P value", width = 0.85) |>
  fontsize(size = 9, part = "all") |>
  padding(padding = 3, part = "all") |>
  theme_booktabs() |>
  set_caption(sprintf(
    "Table 1. Patient characteristics and outcomes among heat-related and non-heat-related OHCA admissions across %d CLIF sites (%s definition).",
    site_count,
    HEAT_DEFINITION
  )) |>
  add_footer_lines(values = c(
    sprintf("Heat-related OHCA is defined as OHCA on days at or above the site-specific %s percentile of warm-season assigned-county daily maximum temperature. Values are n (%%) or median [IQR]. Heat-related OHCA n=%s; non-heat-related OHCA n=%s.",
            str_replace(HEAT_DEFINITION, "heat", ""), format(total_heat_n, big.mark = ","), format(total_non_heat_n, big.mark = ",")),
    "Categorical p values use pooled site counts. Continuous p values are approximate pooled Welch tests from site medians and IQRs because patient-level data are not exchanged in the federated analysis."
  )) |>
  fontsize(size = 8, part = "footer")

flextable::save_as_docx(ft, path = file.path(output_dir, paste0("table1_pooled_hrohca_vs_non_hrohca", OUTPUT_SUFFIX, ".docx")))

gt_table <- display_table |>
  gt(groupname_col = "Section") |>
  tab_header(
    title = "Table 1. Pooled HROHCA vs non-HROHCA cohort characteristics",
    subtitle = sprintf("%d contributing CLIF sites; %s definition", site_count, HEAT_DEFINITION)
  ) |>
  cols_align(align = "center", columns = c(`Heat-related OHCA`, `Non-heat-related OHCA`, `P value`)) |>
  tab_source_note(md(sprintf(
    "Values are n (%%) or median [IQR]. Heat-related OHCA n=%s; non-heat-related OHCA n=%s.",
    format(total_heat_n, big.mark = ","), format(total_non_heat_n, big.mark = ",")
  ))) |>
  tab_source_note(md("Categorical p values use pooled counts; continuous p values are approximate from site medians and IQRs."))

gt::gtsave(gt_table, file.path(output_dir, paste0("table1_pooled_hrohca_vs_non_hrohca", OUTPUT_SUFFIX, ".html")))

message("Wrote pooled Table 1 files to: ", output_dir)

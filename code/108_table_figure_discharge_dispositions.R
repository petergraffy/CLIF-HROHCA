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
box_root <- "/Users/saborpete/Library/CloudStorage/Box-Box/CLIF/Projects/CLIF-Heat-Related-OHCA"
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
table_dir <- file.path(repo_root, "output", "final", "manuscript_tables")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(pooled_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(ragg)
  library(readr)
  library(scales)
  library(stringr)
  library(svglite)
  library(tidyr)
})

fmt_n_pct <- function(n, denom) {
  ifelse(
    is.na(n) | is.na(denom) | denom == 0,
    "",
    sprintf("%s (%.1f%%)", format(round(n), big.mark = ","), 100 * n / denom)
  )
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

theme_manuscript <- function(base_size = 13) {
  theme_classic(base_family = "Helvetica", base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(face = "bold", color = "#111827", size = 12),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(color = "#111827", size = 10.5),
      plot.title = element_text(face = "bold", color = "#111827", size = 16, margin = margin(b = 8)),
      plot.margin = margin(8, 14, 8, 12)
    )
}

save_figure <- function(plot, basename, width, height) {
  ggsave(file.path(figure_dir, paste0(basename, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  svglite::svglite(file.path(figure_dir, paste0(basename, ".svg")), width = width, height = height, bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(figure_dir, paste0(basename, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(figure_dir, paste0(basename, ".tiff")), width = width, height = height, units = "in", res = 600, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

read_export <- function(path) {
  out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"site_name" %in% names(out)) {
    out$site_name <- sub("_.*$", "", basename(path))
  }
  out$n <- as.character(out$n)
  out$source_file <- path
  out
}

files <- list.files(
  box_root,
  pattern = "^[^_]+_heat_related_vs_non_heat_related_discharge_categories[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
files <- files[grepl("/federated_exports/", files)]
if (length(files) == 0) stop("No site discharge category exports found.")

raw_disposition <- bind_rows(lapply(files, read_export)) |>
  mutate(
    n_raw = as.character(.data$n),
    n_numeric = suppressWarnings(as.numeric(.data$n_raw)),
    pct = suppressWarnings(as.numeric(.data$pct)),
    suppressed = is.na(.data$n_numeric) & str_detect(.data$n_raw, "<"),
    discharge_category = str_squish(as.character(.data$discharge_category))
  )

infer_denominator <- function(n_numeric, pct, suppressed) {
  numeric_candidates <- n_numeric[is.finite(n_numeric) & is.finite(pct) & pct > 0] / (pct[is.finite(n_numeric) & is.finite(pct) & pct > 0] / 100)
  if (length(numeric_candidates) > 0) return(round(stats::median(numeric_candidates)))

  pct <- pct[is.finite(pct)]
  if (length(pct) == 0) return(NA_real_)
  candidates <- seq_len(5000)
  rounded_counts <- lapply(candidates, function(denom) round(pct / 100 * denom))
  valid <- vapply(seq_along(candidates), function(i) {
    counts <- rounded_counts[[i]]
    sum(counts) == candidates[[i]] && all(counts > 0) && all(counts < 10)
  }, logical(1))
  if (!any(valid)) return(NA_real_)
  valid_candidates <- candidates[valid]
  errors <- vapply(valid_candidates, function(denom) {
    max(abs(round(pct / 100 * denom) / denom * 100 - pct), na.rm = TRUE)
  }, numeric(1))
  valid_candidates[which.min(errors)]
}

disposition_with_denominator <- raw_disposition |>
  group_by(.data$site_name, .data$heat_definition, .data$heat_related_ohca) |>
  mutate(
    inferred_group_n = infer_denominator(.data$n_numeric, .data$pct, .data$suppressed),
    n_reconstructed = if_else(
      is.na(.data$n_numeric) & is.finite(.data$inferred_group_n) & is.finite(.data$pct),
      round(.data$pct / 100 * .data$inferred_group_n),
      .data$n_numeric
    )
  ) |>
  ungroup()

standardize_disposition <- function(x) {
  x_lower <- str_to_lower(str_squish(x))
  case_when(
    str_detect(x_lower, "expired") ~ "Expired",
    str_detect(x_lower, "hospice") ~ "Hospice",
    str_detect(x_lower, "^home$") ~ "Home",
    str_detect(x_lower, "acute inpatient rehab") ~ "Acute inpatient rehab",
    str_detect(x_lower, "long term care|ltach") ~ "Long-term acute care hospital",
    str_detect(x_lower, "skilled nursing|snf") ~ "Skilled nursing facility",
    str_detect(x_lower, "acute care hospital") ~ "Acute care hospital/transfer",
    str_detect(x_lower, "against|ama") ~ "Against medical advice",
    str_detect(x_lower, "jail|assist living|group home|shelter|psychiatric|chemical dependency|other|missing") ~ "Other/unknown",
    TRUE ~ "Other/unknown"
  )
}

category_levels <- c(
  "Expired",
  "Hospice",
  "Home",
  "Acute inpatient rehab",
  "Skilled nursing facility",
  "Long-term acute care hospital",
  "Acute care hospital/transfer",
  "Against medical advice",
  "Other/unknown"
)

pooled_disposition <- disposition_with_denominator |>
  mutate(disposition = factor(standardize_disposition(.data$discharge_category), levels = category_levels)) |>
  group_by(.data$heat_definition, .data$heat_related_ohca, .data$disposition) |>
  summarise(
    n = sum(.data$n_reconstructed, na.rm = TRUE),
    contributing_sites = n_distinct(.data$site_name),
    suppressed_cells_reconstructed = sum(.data$suppressed, na.rm = TRUE),
    .groups = "drop"
  ) |>
  complete(
    heat_definition,
    heat_related_ohca,
    disposition = factor(category_levels, levels = category_levels),
    fill = list(n = 0, contributing_sites = 0, suppressed_cells_reconstructed = 0)
  ) |>
  group_by(.data$heat_definition, .data$heat_related_ohca) |>
  mutate(total_n = sum(.data$n), pct = 100 * .data$n / .data$total_n) |>
  ungroup()

pooled_table <- pooled_disposition |>
  select("heat_definition", "heat_related_ohca", "disposition", "n", "pct", "total_n", "suppressed_cells_reconstructed")

write.csv(disposition_with_denominator, file.path(pooled_dir, "all_sites_heat_related_vs_non_heat_related_discharge_categories.csv"), row.names = FALSE)
write.csv(pooled_table, file.path(pooled_dir, "pooled_heat_related_vs_non_heat_related_discharge_categories.csv"), row.names = FALSE)

build_table_for_definition <- function(definition = "heat95") {
  dat <- pooled_disposition |>
    filter(.data$heat_definition == definition)

  totals <- dat |>
    distinct(.data$heat_related_ohca, .data$total_n)

  heat_total <- totals$total_n[totals$heat_related_ohca == "Heat-related OHCA"]
  non_heat_total <- totals$total_n[totals$heat_related_ohca == "Non-heat-related OHCA"]

  wide <- dat |>
    mutate(group = if_else(.data$heat_related_ohca == "Heat-related OHCA", "HROHCA", "Non-HROHCA")) |>
    select("disposition", "group", "n", "pct", "total_n") |>
    pivot_wider(
      names_from = "group",
      values_from = c("n", "pct", "total_n"),
      names_sep = "__"
    ) |>
    mutate(
      `HROHCA` = fmt_n_pct(.data$`n__HROHCA`, heat_total),
      `Non-HROHCA` = fmt_n_pct(.data$`n__Non-HROHCA`, non_heat_total),
      `Difference, pp` = sprintf("%.1f", .data$`pct__HROHCA` - .data$`pct__Non-HROHCA`),
      `p value` = fmt_p(mapply(
        function(a, b) {
          if (!all(is.finite(c(a, b, heat_total, non_heat_total)))) return(NA_real_)
          tab <- matrix(c(a, heat_total - a, b, non_heat_total - b), nrow = 2, byrow = TRUE)
          suppressWarnings(stats::fisher.test(tab)$p.value)
        },
        .data$`n__HROHCA`,
        .data$`n__Non-HROHCA`
      ))
    ) |>
    transmute(
      Disposition = as.character(.data$disposition),
      HROHCA = .data$HROHCA,
      `Non-HROHCA` = .data$`Non-HROHCA`,
      `Difference, pp` = .data$`Difference, pp`,
      `p value` = .data$`p value`
    )

  contingency <- dat |>
    select("heat_related_ohca", "disposition", "n") |>
    pivot_wider(names_from = "disposition", values_from = "n", values_fill = 0) |>
    select(-"heat_related_ohca") |>
    as.matrix()

  omnibus_p <- suppressWarnings(stats::chisq.test(contingency)$p.value)
  attr(wide, "omnibus_p") <- omnibus_p
  attr(wide, "heat_total") <- heat_total
  attr(wide, "non_heat_total") <- non_heat_total
  wide
}

table_heat95 <- build_table_for_definition("heat95")
table_heat90 <- build_table_for_definition("heat90")

table_heat95_out <- table_heat95 |>
  mutate(`Heat definition` = "Heat95", .before = 1)
table_heat90_out <- table_heat90 |>
  mutate(`Heat definition` = "Heat90", .before = 1)
combined_table <- bind_rows(table_heat95_out, table_heat90_out)

readr::write_csv(table_heat95_out, file.path(table_dir, "supplement_table_pooled_discharge_dispositions_heat95.csv"))
readr::write_csv(combined_table, file.path(table_dir, "supplement_table_pooled_discharge_dispositions_heat95_heat90.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Heat95")
openxlsx::writeData(wb, "Heat95", table_heat95_out)
openxlsx::setColWidths(wb, "Heat95", cols = 1:ncol(table_heat95_out), widths = "auto")
openxlsx::addWorksheet(wb, "Heat95 and heat90")
openxlsx::writeData(wb, "Heat95 and heat90", combined_table)
openxlsx::setColWidths(wb, "Heat95 and heat90", cols = 1:ncol(combined_table), widths = "auto")
openxlsx::saveWorkbook(wb, file.path(table_dir, "supplement_table_pooled_discharge_dispositions.xlsx"), overwrite = TRUE)

ft <- flextable(table_heat95_out) |>
  autofit() |>
  align(align = "left", part = "all") |>
  theme_booktabs() |>
  bold(part = "header") |>
  add_footer_lines(
    values = sprintf(
      "Primary heat95 definition shown. Omnibus discharge-disposition chi-square p value = %s. Michigan cells reported as <10 were reconstructed from the site-supplied percentages and denominators.",
      fmt_p(attr(table_heat95, "omnibus_p"))
    )
  )

flextable::save_as_docx(
  `Supplemental Table. Pooled discharge dispositions by heat-related OHCA status` = ft,
  path = file.path(table_dir, "supplement_table_pooled_discharge_dispositions.docx")
)

figure_dat <- pooled_disposition |>
  filter(.data$heat_definition == "heat95") |>
  mutate(
    group_label = if_else(
      .data$heat_related_ohca == "Heat-related OHCA",
      paste0("HROHCA\nn=", format(.data$total_n, big.mark = ",")),
      paste0("Non-HROHCA\nn=", format(.data$total_n, big.mark = ","))
    ),
    group = factor(.data$group_label, levels = unique(.data$group_label[order(.data$heat_related_ohca != "Heat-related OHCA")])),
    disposition = factor(.data$disposition, levels = rev(category_levels))
  )

disposition_palette <- c(
  "Expired" = "#7F1D1D",
  "Hospice" = "#C2410C",
  "Home" = "#1D4ED8",
  "Acute inpatient rehab" = "#047857",
  "Skilled nursing facility" = "#65A30D",
  "Long-term acute care hospital" = "#7C3AED",
  "Acute care hospital/transfer" = "#0891B2",
  "Against medical advice" = "#A16207",
  "Other/unknown" = "#6B7280"
)

figure_source <- figure_dat |>
  arrange(.data$group, .data$disposition) |>
  select("heat_definition", "group", "disposition", "n", "pct", "total_n")
write.csv(figure_source, file.path(figure_dir, "figure_pooled_discharge_dispositions_source.csv"), row.names = FALSE)

disposition_plot <- ggplot(figure_dat, aes(x = .data$group, y = .data$pct, fill = .data$disposition)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = if_else(.data$pct >= 5, sprintf("%.0f%%", .data$pct), "")),
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold",
    size = 3.5,
    family = "Helvetica"
  ) +
  scale_fill_manual(values = disposition_palette, breaks = category_levels) +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Pooled hospital discharge disposition after OHCA",
    x = NULL,
    y = "Patients"
  ) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  theme_manuscript(base_size = 13) +
  theme(
    legend.key.height = grid::unit(0.38, "cm"),
    legend.key.width = grid::unit(0.42, "cm"),
    plot.title = element_text(hjust = 0),
    axis.text.x = element_text(face = "bold", size = 12)
  )

save_figure(disposition_plot, "figure_pooled_discharge_dispositions", width = 8.2, height = 5.2)

message("Wrote pooled discharge disposition table and figure.")

#!/usr/bin/env Rscript

get_script_path <- function() {
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) stop("Could not determine script path from commandArgs().")
  normalizePath(sub(file_arg, "", match[[1]]), winslash = "/", mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."), winslash = "/", mustWork = TRUE)

ensure_user_library <- function() {
  version_parts <- strsplit(as.character(getRversion()), "\\.")[[1]]
  version_stub <- paste(version_parts[1], version_parts[2], sep = ".")
  user_lib <- file.path(repo_root, ".r-user-lib", version_stub)
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}

ensure_user_library()

suppressPackageStartupMessages({
  library(ggplot2)
})

manuscript_dir <- file.path(repo_root, "output", "final", "ohca_tmax", "manuscript")
descriptive_dir <- file.path(repo_root, "output", "final", "descriptive")
outcomes_dir <- file.path(repo_root, "output", "final", "ohca_outcomes")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

results <- read.csv(file.path(manuscript_dir, "manuscript_dlnm_results.csv"), stringsAsFactors = FALSE)
curves <- read.csv(file.path(manuscript_dir, "manuscript_dlnm_curves.csv"), stringsAsFactors = FALSE)
table1 <- read.csv(file.path(descriptive_dir, "table1_ohca_cohort_characteristics.csv"), stringsAsFactors = FALSE)
heat_outcomes <- read.csv(file.path(descriptive_dir, "heat_related_vs_non_heat_related_ohca_outcomes.csv"), stringsAsFactors = FALSE)
heat_table2_path <- file.path(descriptive_dir, "table2_heat_related_vs_non_heat_related_ohca.csv")
heat_table2_all_path <- file.path(descriptive_dir, "table2_heat_related_vs_non_heat_related_ohca_all_definitions.csv")
heat90_table2_path <- file.path(descriptive_dir, "table2_heat90_vs_non_heat90_ohca.csv")
adverse_models_path <- file.path(outcomes_dir, "ohca_heat_adverse_outcome_models.csv")
continuous_models_path <- file.path(outcomes_dir, "ohca_heat_continuous_outcome_models.csv")
pollution_binary_models_path <- file.path(outcomes_dir, "ohca_pollution_12m_binary_outcome_models.csv")
pollution_continuous_models_path <- file.path(outcomes_dir, "ohca_pollution_12m_continuous_outcome_models.csv")

fmt_rr <- function(rr, low, high) {
  sprintf("%.2f (%.2f, %.2f)", rr, low, high)
}

final_dlnm_table <- results
final_dlnm_table$rr_95_ci <- fmt_rr(
  final_dlnm_table$cumulative_rr,
  final_dlnm_table$cumulative_rr_low,
  final_dlnm_table$cumulative_rr_high
)
final_dlnm_table <- final_dlnm_table[, c(
  "stratum", "model", "n_ohca", "reference_type", "reference_temp_c",
  "hot_temp_c", "rr_95_ci", "log_rr", "log_rr_se"
)]

write.csv(final_dlnm_table, file.path(manuscript_dir, "table2_manuscript_dlnm_results_formatted.csv"), row.names = FALSE)
write.csv(table1, file.path(manuscript_dir, "table1_ohca_cohort_characteristics_formatted.csv"), row.names = FALSE)
write.csv(heat_outcomes, file.path(manuscript_dir, "table3_heat_related_ohca_outcomes.csv"), row.names = FALSE)
if (file.exists(heat_table2_path)) {
  heat_table2 <- read.csv(heat_table2_path, stringsAsFactors = FALSE)
  write.csv(heat_table2, file.path(manuscript_dir, "table2_heat_related_vs_non_heat_related_ohca.csv"), row.names = FALSE)
}
if (file.exists(heat_table2_all_path)) {
  heat_table2_all <- read.csv(heat_table2_all_path, stringsAsFactors = FALSE)
  write.csv(heat_table2_all, file.path(manuscript_dir, "table2_heat_related_vs_non_heat_related_ohca_all_definitions.csv"), row.names = FALSE)
}
if (file.exists(heat90_table2_path)) {
  heat90_table2 <- read.csv(heat90_table2_path, stringsAsFactors = FALSE)
  write.csv(heat90_table2, file.path(manuscript_dir, "table2_heat90_vs_non_heat90_ohca.csv"), row.names = FALSE)
}
if (file.exists(adverse_models_path)) {
  adverse_models <- read.csv(adverse_models_path, stringsAsFactors = FALSE)
  adverse_models$or_95_ci <- fmt_rr(adverse_models$odds_ratio, adverse_models$ci_low, adverse_models$ci_high)
  adverse_models <- adverse_models[, c("outcome", "exposure", "n", "events", "or_95_ci", "p_value", "estimable", "converged", "note")]
  write.csv(adverse_models, file.path(manuscript_dir, "table4_ohca_heat_adverse_outcome_models.csv"), row.names = FALSE)
}
if (file.exists(continuous_models_path)) {
  continuous_models <- read.csv(continuous_models_path, stringsAsFactors = FALSE)
  continuous_models$ratio_95_ci <- fmt_rr(continuous_models$geometric_mean_ratio, continuous_models$ci_low, continuous_models$ci_high)
  continuous_models <- continuous_models[, c("outcome", "exposure", "n", "ratio_95_ci", "p_value", "estimable", "note")]
  write.csv(continuous_models, file.path(manuscript_dir, "table5_ohca_heat_continuous_outcome_models.csv"), row.names = FALSE)
}
if (file.exists(pollution_binary_models_path)) {
  pollution_binary_models <- read.csv(pollution_binary_models_path, stringsAsFactors = FALSE)
  pollution_binary_models$or_95_ci <- fmt_rr(pollution_binary_models$odds_ratio, pollution_binary_models$ci_low, pollution_binary_models$ci_high)
  pollution_binary_models <- pollution_binary_models[, c("outcome", "exposure", "adjustment_set", "n", "events", "or_95_ci", "p_value", "estimable", "note")]
  write.csv(pollution_binary_models, file.path(manuscript_dir, "table6_ohca_pollution_12m_binary_outcome_models.csv"), row.names = FALSE)
}
if (file.exists(pollution_continuous_models_path)) {
  pollution_continuous_models <- read.csv(pollution_continuous_models_path, stringsAsFactors = FALSE)
  pollution_continuous_models$ratio_95_ci <- fmt_rr(pollution_continuous_models$geometric_mean_ratio, pollution_continuous_models$ci_low, pollution_continuous_models$ci_high)
  pollution_continuous_models <- pollution_continuous_models[, c("outcome", "exposure", "adjustment_set", "n", "ratio_95_ci", "p_value", "estimable", "note")]
  write.csv(pollution_continuous_models, file.path(manuscript_dir, "table7_ohca_pollution_12m_continuous_outcome_models.csv"), row.names = FALSE)
}

primary_curve <- curves[curves$stratum == "Overall" & curves$model == "primary_humidity_adjusted", , drop = FALSE]

p_curve <- ggplot(primary_curve, aes(x = tmax_mean_c, y = cumulative_rr)) +
  geom_ribbon(aes(ymin = cumulative_rr_low, ymax = cumulative_rr_high), fill = "#A7C957", alpha = 0.28) +
  geom_line(color = "#1B4332", linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = unique(primary_curve$reference_temp_c), linetype = "dotted", color = "#1B4332") +
  labs(
    title = "Median-Referenced DLNM Association Between Tmax and OHCA ICU Admissions",
    subtitle = "Primary warm-season quasi-Poisson DLNM adjusted for relative humidity",
    x = "Daily assigned-county Tmax (C)",
    y = "Cumulative relative risk"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(figure_dir, "figure_primary_dlnm_curve.png"), p_curve, width = 8, height = 5.5, dpi = 300)

forest_df <- results[results$model %in% c("primary_humidity_adjusted", "stratified_humidity_adjusted"), , drop = FALSE]
forest_df$stratum <- factor(
  forest_df$stratum,
  levels = rev(c("Overall", "Male", "Female", "<65", ">=65", "Black", "Non-Black"))
)

p_forest <- ggplot(forest_df, aes(y = stratum, x = cumulative_rr)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbar(
    aes(xmin = cumulative_rr_low, xmax = cumulative_rr_high),
    width = 0.18,
    orientation = "y",
    color = "#3D405B"
  ) +
  geom_point(size = 2.4, color = "#3D405B") +
  scale_x_log10() +
  labs(
    title = "Stratified DLNM Heat Associations",
    subtitle = "RR for 95th percentile Tmax versus median Tmax",
    x = "Cumulative relative risk, log scale",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(figure_dir, "figure_stratified_dlnm_forest.png"), p_forest, width = 7.5, height = 5.5, dpi = 300)

message("Wrote manuscript tables and figures to ", manuscript_dir, " and ", figure_dir)

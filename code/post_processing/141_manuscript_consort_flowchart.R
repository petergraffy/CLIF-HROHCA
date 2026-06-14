#!/usr/bin/env Rscript

required_packages <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

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

repo_root <- normalizePath(file.path(dirname(get_script_path()), "..", ".."), winslash = "/", mustWork = TRUE)
pooled_dir <- file.path(repo_root, "output", "final", "federated_pooled")
figure_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

box_project_dir <- path.expand("~/Library/CloudStorage/Box-Box/CLIF/Projects/CLIF-Heat-Related-OHCA")
flow_candidates <- c(
  file.path(pooled_dir, "all_site_ohca_icu_72h_consort_flow.csv"),
  file.path(box_project_dir, "federated_pooled", "all_site_ohca_icu_72h_consort_flow.csv")
)
ed_death_path <- file.path(pooled_dir, "all_site_ohca_ed_only_death_never_icu_summary.csv")
imv12_flow_candidates <- c(
  file.path(pooled_dir, "all_site_ohca_icu_imv12_discharge_outcome_flow.csv"),
  file.path(box_project_dir, "federated_pooled", "all_site_ohca_icu_imv12_discharge_outcome_flow.csv")
)

read_multisite_flow <- function(candidates, source_name) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Missing ", source_name, ". Expected one of: ", paste(candidates, collapse = ", "), call. = FALSE)
  }

  flows <- lapply(existing, function(path) {
    dat <- read.csv(path, stringsAsFactors = FALSE)
    dat$.source_path <- path
    dat
  })
  site_counts <- vapply(flows, function(dat) length(unique(dat$site_name)), integer(1))
  flows[[which.max(site_counts)]]
}

flow_site <- read_multisite_flow(flow_candidates, "pooled CONSORT flow source")
flow <- aggregate(n ~ step_order + step, flow_site, sum)

get_n <- function(step_order) {
  value <- flow$n[flow$step_order == step_order]
  if (length(value) != 1L) stop("Could not find one row for step_order ", step_order, call. = FALSE)
  as.integer(value)
}

fmt_n <- function(x) format(as.integer(round(x)), big.mark = ",")
html_label <- function(...) {
  lines <- unlist(strsplit(paste(..., sep = "\n"), "\n", fixed = TRUE))
  lines <- gsub("&", "&amp;", lines, fixed = TRUE)
  lines <- gsub("<", "&lt;", lines, fixed = TRUE)
  lines <- gsub(">", "&gt;", lines, fixed = TRUE)
  paste0("<", paste(lines, collapse = "<BR/>"), ">")
}

html_label_bold_first <- function(...) {
  lines <- unlist(strsplit(paste(..., sep = "\n"), "\n", fixed = TRUE))
  lines <- gsub("&", "&amp;", lines, fixed = TRUE)
  lines <- gsub("<", "&lt;", lines, fixed = TRUE)
  lines <- gsub(">", "&gt;", lines, fixed = TRUE)
  if (length(lines) > 0L) lines[[1]] <- paste0("<B>", lines[[1]], "</B>")
  paste0("<", paste(lines, collapse = "<BR/>"), ">")
}

n_sites <- length(unique(flow_site$site_name))
n_all_icu <- get_n(1)
n_ohca_poa <- get_n(2)
n_pathway <- get_n(3)
n_exposure <- get_n(4)
n_classified <- get_n(5)
n_alive_no_imv <- get_n(6)
n_awake_extubated <- get_n(7)
n_limited <- get_n(8)
n_anoxic <- get_n(9)

n_excluded_no_ohca <- n_all_icu - n_ohca_poa
n_excluded_pathway <- n_ohca_poa - n_pathway
n_excluded_exposure <- n_pathway - n_exposure
n_unclassified <- n_exposure - n_classified

imv12_flow <- read_multisite_flow(imv12_flow_candidates, "IMV12 pooled flow source")
n_imv12 <- sum(imv12_flow$n[imv12_flow$step == "IMV duration >=12 hours"], na.rm = TRUE)
if (!is.finite(n_imv12) || n_imv12 <= 0) {
  stop("Could not derive IMV12 cohort size from pooled IMV12 flow.", call. = FALSE)
}

ed_hosp <- NA_integer_
ed_patients <- NA_integer_
if (file.exists(ed_death_path)) {
  ed_death <- read.csv(ed_death_path, stringsAsFactors = FALSE)
  ed_hosp <- sum(ed_death$n_ohca_dx_ed_first_only_death_never_icu_hospitalizations, na.rm = TRUE)
  ed_patients <- sum(ed_death$n_ohca_dx_ed_first_only_death_never_icu_patients, na.rm = TRUE)
}

source_table <- data.frame(
  display_order = seq_len(14),
  section = c(rep("Main inclusion flow", 9), rep("72-hour phenotype assignment", 4), "IMV12 cohort"),
  label = c(
    "Contributing sites",
    "All ICU admissions in CLIF, 2018-2024",
    "Excluded: ICU admissions without OHCA present-on-admission diagnosis",
    "OHCA present-on-admission ICU admissions before pathway/timing restriction",
    "Excluded: non-eligible ICU pathway or ICU entry >=24h",
    "Eligible OHCA ICU cohort: ED/procedure/direct ICU pathway and ICU entry <24h",
    "Excluded: missing admission-day county Tmax or humidity",
    "Exposure-complete OHCA ICU cohort",
    "Excluded: unclassified 72-hour structured-data phenotype",
    "No IMV in first 72h",
    "Extubated by 72h",
    "On IMV at 72h",
    "Death within 72h",
    "IMV duration >=12 hours"
  ),
  n = c(
    n_sites,
    n_all_icu,
    n_excluded_no_ohca,
    n_ohca_poa,
    n_excluded_pathway,
    n_pathway,
    n_excluded_exposure,
    n_exposure,
    n_unclassified,
    n_alive_no_imv,
    n_awake_extubated,
    n_limited,
    n_anoxic,
    n_imv12
  ),
  stringsAsFactors = FALSE
)
if (is.finite(ed_hosp)) {
  source_table <- rbind(
    source_table,
    data.frame(
      display_order = 15L,
      section = "Separate ED-only audit",
      label = "OHCA diagnosis, ED only, died before ICU",
      n = ed_hosp,
      stringsAsFactors = FALSE
    )
  )
}

write.csv(source_table, file.path(figure_dir, "figure_consort_flow_source.csv"), row.names = FALSE)
write.csv(flow_site[names(flow_site) != ".source_path"], file.path(figure_dir, "figure_consort_flow_site_source.csv"), row.names = FALSE)

node_line <- function(id, label, fill = "#F7FBFF", color = "#3D3D3D", style = "filled", width = 3.8, height = 0.7, fontsize = 18) {
  sprintf(
    '%s [label=%s, fillcolor="%s", color="%s", style="%s", width=%.2f, height=%.2f, fontsize=%d];',
    id,
    label,
    fill,
    color,
    style,
    width,
    height,
    fontsize
  )
}

edge_line <- function(from, to, style = "solid", weight = 1) {
  sprintf('%s -> %s [style="%s", weight=%d];', from, to, style, weight)
}

edge_line_ports <- function(from, from_port, to, to_port, style = "solid", weight = 1) {
  sprintf('%s:%s -> %s:%s [style="%s", weight=%d];', from, from_port, to, to_port, style, weight)
}

dot <- paste(
  "digraph consort_flow {",
  "graph [layout = dot, rankdir = TB, bgcolor = white, splines = polyline, nodesep = 0.45, ranksep = 0.60, margin = 0.02];",
  "node [shape = box, fontname = Helvetica, margin = 0.12, penwidth = 1.2, color = '#3D3D3D'];",
  "edge [color = '#3D3D3D', arrowsize = 0.75, penwidth = 1.1];",
  node_line("all", html_label("All ICU admissions in CLIF, 2018-2024", paste0("N = ", fmt_n(n_all_icu))), width = 4.4),
  node_line("ex_no_ohca", html_label_bold_first("Excluded", "No OHCA POA ICU diagnosis", paste0("n = ", fmt_n(n_excluded_no_ohca))), fill = "#FFF7F3", width = 3.4, height = 0.65, fontsize = 15),
  node_line("poa", html_label("OHCA present-on-admission ICU admissions", paste0("N = ", fmt_n(n_ohca_poa))), width = 4.4),
  node_line("ex_pathway", html_label_bold_first("Excluded", "Non-eligible pathway or ICU entry >=24h", paste0("n = ", fmt_n(n_excluded_pathway))), fill = "#FFF7F3", width = 3.4, height = 0.65, fontsize = 15),
  node_line("eligible", html_label("Eligible OHCA ICU cohort", "ED/procedure/direct ICU pathway", "ICU entry <24h", paste0("N = ", fmt_n(n_pathway))), width = 4.7, height = 0.9),
  node_line("ex_exposure", html_label_bold_first("Excluded", "Missing admission-day Tmax or humidity", paste0("n = ", fmt_n(n_excluded_exposure))), fill = "#FFF7F3", width = 3.4, height = 0.65, fontsize = 15),
  node_line("exposure", html_label("Admission-day Tmax and humidity available", paste0("N = ", fmt_n(n_exposure))), width = 4.4),
  node_line("ex_unclassified", html_label_bold_first("Excluded", "Unclassified phenotype", paste0("n = ", fmt_n(n_unclassified))), fill = "#FFF7F3", width = 3.4, height = 0.65, fontsize = 15),
  node_line("classified", html_label("Primary analytic cohort of OHCA cases", paste0("N = ", fmt_n(n_classified))), width = 4.4),
  node_line("alive", html_label("No IMV in first 72h", paste0("n = ", fmt_n(n_alive_no_imv))), fill = "#EEF7F2", width = 2.2, height = 0.75, fontsize = 14),
  node_line("awake", html_label("Extubated by 72h", paste0("n = ", fmt_n(n_awake_extubated))), fill = "#EEF7F2", width = 2.2, height = 0.75, fontsize = 14),
  node_line("limited", html_label("On IMV at 72h", paste0("n = ", fmt_n(n_limited))), fill = "#EEF7F2", width = 2.2, height = 0.75, fontsize = 14),
  node_line("anoxic", html_label("Death within 72h", paste0("n = ", fmt_n(n_anoxic))), fill = "#EEF7F2", width = 2.2, height = 0.75, fontsize = 14),
  node_line("imv12", html_label("IMV duration >=12h", paste0("n = ", fmt_n(n_imv12))), fill = "#FFF4DE", width = 3.0, height = 0.75, fontsize = 15),
  "{rank = same; all; ex_no_ohca;}",
  "{rank = same; poa; ex_pathway;}",
  "{rank = same; eligible; ex_exposure;}",
  "{rank = same; exposure; ex_unclassified;}",
  "{rank = same; alive; awake; limited; anoxic;}",
  edge_line("all", "poa", weight = 10),
  edge_line("all", "ex_no_ohca", weight = 1),
  edge_line("poa", "eligible", weight = 10),
  edge_line("poa", "ex_pathway", weight = 1),
  edge_line("eligible", "exposure", weight = 10),
  edge_line("eligible", "ex_exposure", weight = 1),
  edge_line("exposure", "classified", weight = 10),
  edge_line("exposure", "ex_unclassified", weight = 1),
  edge_line_ports("classified", "s", "alive", "n", weight = 2),
  edge_line_ports("classified", "s", "awake", "n", weight = 2),
  edge_line_ports("classified", "s", "limited", "n", weight = 2),
  edge_line_ports("classified", "s", "anoxic", "n", weight = 2),
  edge_line_ports("awake", "s", "imv12", "n", weight = 2),
  edge_line_ports("limited", "s", "imv12", "n", weight = 2),
  edge_line_ports("anoxic", "s", "imv12", "n", weight = 2),
  "}",
  sep = "\n"
)

dot_path <- file.path(figure_dir, "figure_consort_cohort_flow.dot")
writeLines(dot, dot_path)

graph <- DiagrammeR::grViz(dot)
svg <- DiagrammeRsvg::export_svg(graph)

svg_path <- file.path(figure_dir, "figure_consort_cohort_flow.svg")
png_path <- file.path(figure_dir, "figure_consort_cohort_flow.png")
pdf_path <- file.path(figure_dir, "figure_consort_cohort_flow.pdf")
tiff_path <- file.path(figure_dir, "figure_consort_cohort_flow.tiff")

writeLines(svg, svg_path)
rsvg::rsvg_pdf(charToRaw(svg), pdf_path, width = 11, height = 9.0)
png_tmp <- tempfile(fileext = ".png")
rsvg::rsvg_png(charToRaw(svg), png_tmp, width = 6600, height = 5400)
if (requireNamespace("magick", quietly = TRUE)) {
  magick::image_read(png_tmp) |>
    magick::image_background("white", flatten = TRUE) |>
    magick::image_write(path = png_path, format = "png")
  magick::image_read(png_path) |>
    magick::image_write(path = tiff_path, format = "tiff")
} else {
  file.copy(png_tmp, png_path, overwrite = TRUE)
  warning("Package 'magick' is not installed; skipping TIFF export.", call. = FALSE)
}

message("Wrote DiagrammeR CONSORT-style cohort flowchart to ", png_path, ", ", pdf_path, ", ", svg_path, ", and ", tiff_path)

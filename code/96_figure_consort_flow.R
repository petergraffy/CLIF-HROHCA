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
output_dir <- file.path(repo_root, "output", "final", "manuscript_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
})

fmt_n <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

read_flow <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"site_name" %in% names(x)) {
    x$site_name <- if (grepl("/JHU/", path)) "JHU" else sub("_.*$", "", basename(path))
  }
  x$source_file <- path
  x
}

flow_files <- list.files(box_root, pattern = "_cohort_flow[.]csv$", recursive = TRUE, full.names = TRUE)
flow_files <- flow_files[grepl("/federated_exports/", flow_files)]
jhu_flow <- file.path(box_root, "JHU", "final", "quality_checks", "cohort_flow.csv")
if (file.exists(jhu_flow)) flow_files <- c(flow_files, jhu_flow)
flows <- do.call(rbind, lapply(flow_files, read_flow))
site_count <- length(unique(flows$site_name))

pooled_steps <- aggregate(n ~ step_order + step, flows, function(x) sum(x, na.rm = TRUE))
get_step_n <- function(step_order) pooled_steps$n[pooled_steps$step_order == step_order][[1]]

outcome_counts <- read.csv(
  file.path(repo_root, "output", "final", "federated_pooled", "pooled_heat_related_vs_non_heat_related_outcome_counts.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
heat95_counts <- subset(outcome_counts, heat_definition == "heat95" & outcome == "hospital_mortality")
hrohca_n <- heat95_counts$n[heat95_counts$heat_related_ohca == "Heat-related OHCA"]
non_hrohca_n <- heat95_counts$n[heat95_counts$heat_related_ohca == "Non-heat-related OHCA"]
analysis_n <- sum(heat95_counts$n)

n_all <- get_step_n(1)
n_adult <- get_step_n(2)
n_2018 <- get_step_n(3)
n_icu <- get_step_n(4)
n_ohca <- get_step_n(5)
n_geo <- get_step_n(6)
n_warm <- get_step_n(7)
n_complete <- get_step_n(8)

flow_summary <- data.frame(
  measure = c(
    "Contributing sites",
    "Adult ICU hospitalizations, 2018-2024",
    "Excluded: ICU without OHCA POA",
    "ICU admissions with OHCA present on admission",
    "Warm-season OHCA ICU admissions",
    "Heat95 clinical phenotype cohort",
    "Heat-related OHCA",
    "Non-heat-related OHCA"
  ),
  n = c(
    site_count,
    n_icu,
    n_icu - n_ohca,
    n_ohca,
    n_warm,
    analysis_n,
    hrohca_n,
    non_hrohca_n
  )
)
write.csv(flow_summary, file.path(output_dir, "figure_consort_flow_source.csv"), row.names = FALSE)
write.csv(flows, file.path(output_dir, "figure_consort_flow_site_source.csv"), row.names = FALSE)

dot_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\"", "\\\\\"", x)
  x
}

node_text <- function(lines, n = NULL) {
  label <- c(lines, if (!is.null(n)) paste0("n = ", fmt_n(n)) else NULL)
  dot_escape(paste(label, collapse = "\\n"))
}

node_html <- function(lines, n = NULL, bold_first = FALSE) {
  label <- c(lines, if (!is.null(n)) paste0("n = ", fmt_n(n)) else NULL)
  label <- html_escape(label)
  if (bold_first && length(label) > 0) label[1] <- paste0("<B>", label[1], "</B>")
  paste(label, collapse = "<BR/>")
}

flow_nodes <- data.frame(
  id = c("icu", "ohca", "warm", "hrohca", "nonhrohca", "excl_no_ohca", "excl_nonwarm"),
  label = c(
    node_text(c("Adult ICU hospitalizations", "(2018-2024)"), n_icu),
    node_text(c("ICU admissions with OHCA", "present on admission"), n_ohca),
    node_text(c("Warm-season OHCA admissions", "(May-Sept)"), n_warm),
    node_text(c("Heat-Related", "Out-of-Hospital Cardiac Arrest", "(HROHCA)"), hrohca_n),
    node_text(c("Non-Heat-Related", "Out-of-Hospital Cardiac Arrest", "(non-HROHCA)"), non_hrohca_n),
    node_html(c("Excluded", "No OHCA present on admission"), n_icu - n_ohca, bold_first = TRUE),
    node_html(c("Excluded", "Non-warm-season OHCA"), n_ohca - n_warm, bold_first = TRUE)
  ),
  type = c("include", "include", "include", "hrohca", "nonhrohca", "exclude", "exclude"),
  label_type = c("text", "text", "text", "text", "text", "html", "html"),
  stringsAsFactors = FALSE
)

flow_edges <- data.frame(
  from = c("icu", "ohca", "warm", "warm", "icu", "ohca"),
  to = c("ohca", "warm", "hrohca", "nonhrohca", "excl_no_ohca", "excl_nonwarm"),
  edge_type = c("primary", "primary", "primary", "primary", "exclude", "exclude"),
  stringsAsFactors = FALSE
)

node_lines <- vapply(seq_len(nrow(flow_nodes)), function(i) {
  fill <- switch(
    flow_nodes$type[i],
    include = "#FFFFFF",
    hrohca = "#FFF1ED",
    nonhrohca = "#EEF5FF",
    exclude = "#F3F4F6",
    "#FFFFFF"
  )
  color <- if (flow_nodes$type[i] == "exclude") "#6B7280" else "#111827"
  penwidth <- if (flow_nodes$type[i] == "exclude") "1.2" else "1.8"
  label <- if (flow_nodes$label_type[i] == "html") {
    paste0("<", flow_nodes$label[i], ">")
  } else {
    paste0("\"", flow_nodes$label[i], "\"")
  }
  sprintf(
    '%s [label=%s, fillcolor="%s", color="%s", penwidth=%s];',
    flow_nodes$id[i], label, fill, color, penwidth
  )
}, character(1))

edge_lines <- vapply(seq_len(nrow(flow_edges)), function(i) {
  if (flow_edges$edge_type[i] == "exclude") {
    sprintf(
      '%s -> %s [constraint=false, color="#9CA3AF", penwidth=1.2];',
      flow_edges$from[i], flow_edges$to[i]
    )
  } else {
    sprintf('%s -> %s;', flow_edges$from[i], flow_edges$to[i])
  }
}, character(1))

rank_same_blocks <- paste(
  "{rank=same; icu; excl_no_ohca}",
  "{rank=same; ohca; excl_nonwarm}",
  "{rank=same; hrohca; nonhrohca}",
  sep = "\n  "
)

dot <- sprintf('digraph cohort_flow {
  graph [
    rankdir=TB,
    bgcolor="white",
    margin=0.28,
    pad=0.35,
    nodesep=0.45,
    ranksep=0.60,
    splines=ortho,
    outputorder=edgesfirst
  ];
  node [
    shape=box,
    style="rounded,filled",
    fontname="Helvetica",
    fontsize=18,
    margin="0.16,0.11",
    color="#111827",
    fillcolor="white"
  ];
  edge [
    color="#374151",
    penwidth=1.8,
    arrowsize=0.75,
    fontname="Helvetica",
    fontsize=14
  ];

  %s

  %s

  %s
}', paste(node_lines, collapse = "\n  "), paste(edge_lines, collapse = "\n  "), rank_same_blocks)

graph <- DiagrammeR::grViz(dot)
svg <- DiagrammeRsvg::export_svg(graph)

svg_path <- file.path(output_dir, "figure_consort_cohort_flow.svg")
png_path <- file.path(output_dir, "figure_consort_cohort_flow.png")
pdf_path <- file.path(output_dir, "figure_consort_cohort_flow.pdf")
tiff_path <- file.path(output_dir, "figure_consort_cohort_flow.tiff")

writeLines(svg, svg_path, useBytes = TRUE)
rsvg::rsvg_png(charToRaw(svg), png_path, width = 3000)
rsvg::rsvg_pdf(charToRaw(svg), pdf_path)
if (requireNamespace("magick", quietly = TRUE)) {
  img <- magick::image_read(png_path)
  magick::image_write(img, path = tiff_path, format = "tiff", compression = "lzw")
} else {
  warning("Package 'magick' is not available; skipping TIFF export.")
}

message("Wrote DiagrammeR CONSORT-style cohort flow figure to ", output_dir)

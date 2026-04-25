# CLIF-HROHCA

Heat exposure and out-of-hospital cardiac arrest (OHCA) ICU admissions in CLIF.

This repository is organized for a federated CLIF analysis. Each site runs the same R workflow locally against its own CLIF 2.1 tables. Sites should share only the aggregate files created in `output/final/federated_exports/`.

## Study Aim

Estimate whether county-level heat exposure is associated with daily OHCA ICU admissions, using hospital-aware county exposure linkage and distributed lag nonlinear models (DLNMs).

## Required CLIF Tables

- `clif_patient`: demographics for Table 1 and stratified analyses.
- `clif_hospitalization`: admission date, discharge disposition, age, and patient county.
- `clif_adt`: ICU entry, hospital identifier, ICU length of stay, and care pathway quality checks.
- `clif_hospital_diagnosis` or `clif_admission_diagnosis`: present-on-admission cardiac arrest diagnosis codes.
- `clif_respiratory_support`: optional but recommended for IMV rate and IMV duration.
- `clif_medication_admin_continuous`: optional but recommended for vasopressor use.

## Exposome Inputs

The root `exposome/` folder contains county-level environmental files used by all sites:

- `daymet_county_tmax_2018_2024_conus.parquet`
- `daymet_county_rmax_2018_2024.parquet`
- `no2_county_year.csv`
- `pm25_county_year.csv`

These files are county-level only and contain no CLIF patient information.

## Site Setup

1. Copy `config/config_template.json` to `config/config.json`.
2. Set `site_name` to the site label used in `reference/clif_hospital_geography.csv`.
3. Set `tables_path` to your local CLIF 2.1 table directory.
4. Run package setup:

```r
Rscript code/00_install_or_restore_packages.R
```

## Run The Site Analysis

Run the full site workflow:

```r
Rscript code/run_site_analysis.R
```

Or run scripts one at a time from `code/`:

```r
Rscript code/01_build_ohca_cohort.R
Rscript code/02_build_icu_exposure_series.R
Rscript code/03_descriptive_tables.R
Rscript code/04_dlnm_primary_and_sensitivity.R
Rscript code/05_heat_related_vs_non_heat_related_table.R
Rscript code/09_supplementary_ohca_outcome_models.R
Rscript code/08_quality_checks.R
Rscript code/06_manuscript_tables_figures.R
Rscript code/07_export_federated_results.R
```

## Federated Sharing

Share only files in:

```text
output/final/federated_exports/
```

Those files are aggregate-only and designed to be pooled across sites. Do not share `output/intermediate/`, which contains CLIF-derived row-level working files.

The coordinating center can pool returned site DLNM estimates with:

```r
Rscript code/90_pool_federated_results.R
```

This also pools exported aggregate DLNM curve points into `output/final/federated_pooled/pooled_dlnm_random_effects_curves.csv` for a CLIF-wide supplemental DLNM curve.

Sites also export reduced DLNM coefficient and variance-covariance files from `dlnm::crossreduce()`, which can be used for a future multivariate DLNM meta-analysis of curve shape.

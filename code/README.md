# Code Workflow

The site-facing workflow is R-only and uses numbered scripts in this folder.

Most sites should run:

Windows PowerShell:

```powershell
Rscript code\run_site_analysis.R
```

macOS/Linux:

```bash
Rscript code/run_site_analysis.R
```

The numbered scripts can also be run one at a time for debugging.

For example, to install packages only:

Windows PowerShell:

```powershell
Rscript code\00_install_or_restore_packages.R
```

macOS/Linux:

```bash
Rscript code/00_install_or_restore_packages.R
```

## Site Scripts

1. `00_install_or_restore_packages.R`: installs/restores required R packages into the project-local library and snapshots `renv.lock` when `renv` is available.
2. `01_build_ohca_cohort.R`: builds adult ICU OHCA hospitalizations using present-on-admission cardiac arrest diagnosis codes and applies hospital-aware county assignment.
3. `02_build_icu_exposure_series.R`: builds the all-ICU daily patient-address exposure series used as the comparison denominator for daily OHCA models.
4. `04_dlnm_primary_and_sensitivity.R`: fits primary, pollution-adjusted, MRT-reference, time-adjustment sensitivity, and stratified DLNM models, including cumulative lag-window and day-specific hot-temperature contrasts from lag 0 through lag 5.
5. `04c_lag30_diagnostic_plots.R`: creates lag-30 diagnostic summaries and figures for lag-window justification.
6. `04d_ohca_icu_admission_rate_dlnm.R`: fits OHCA-per-ICU-admission-rate DLNMs with the same count-model denominator across strata.
7. `04b_case_crossover_dlnm.R`: fits time-stratified case-crossover DLNMs.
8. `10_ohca_icu_72h_phenotypes.R`: creates structured 72-hour ICU OHCA phenotypes and the multinomial phenotype-assignment model using admission-day temperature and humidity, with a mechanism-adjusted sensitivity.
9. `11_ohca_icu_competing_risks.R`: fits 72-hour Fine-Gray competing-risk models for time to awake-and-extubated recovery versus death before awake/extubated using admission-day temperature and humidity, with a mechanism-adjusted sensitivity.
10. `07_export_federated_results.R`: writes DLNM and 72-hour phenotype aggregate site files plus selected site-level figure PNGs for federated pooling and visual QC.

## Optional Tables

The workflow can run with the required cohort/modeling tables only, but the 72-hour phenotype models are strongest when these optional CLIF tables are available:

- `clif_respiratory_support`: IMV trajectories, duration, 72-hour extubation evidence, and competing-risk awake/extubated timing.
- `clif_patient_assessments`: GCS, RASS, AVPU, SAT, and SBT evidence.
- `clif_vitals`: last-vital fallback for death timing in the competing-risk model.

## Privacy Boundary

Only `output/final/federated_exports/` should be shared. This folder contains aggregate DLNM and 72-hour phenotype CSVs plus selected site-level PNG figures. The scripts intentionally keep row-level CLIF-derived working files under `output/intermediate/`, which is git-ignored and should remain local.

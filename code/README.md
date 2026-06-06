# Code Workflow

The site-facing workflow is R-only and uses numbered scripts in this folder.

Top-level scripts in `code/` are the site-running workflow. Coordinating-center pooling, manuscript figures, and other post-processing helpers live in `code/post_processing/`. Historical scripts that are no longer part of the current analytic flow live in `code/archive/`.

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
3. `01b_ohca_ed_death_never_icu.R`: counts OHCA-diagnosis patients who arrived to the ED, stayed only in ED locations, died there, and never reached ICU.
4. `02_build_icu_exposure_series.R`: builds the all-ICU daily patient-address exposure series used as the comparison denominator for daily OHCA models.
5. `04_dlnm_primary_and_sensitivity.R`: fits primary, pollution-adjusted, MRT-reference, time-adjustment sensitivity, and stratified DLNM models, including cumulative lag-window and day-specific hot-temperature contrasts from lag 0 through lag 5.
6. `04c_lag30_diagnostic_plots.R`: creates count-model lag-30 diagnostic summaries and figures for lag-window justification.
7. `04d_ohca_icu_admission_rate_dlnm.R`: fits OHCA-per-ICU-admission-rate DLNMs with the same count-model denominator across strata.
8. `04e_rate_lag30_diagnostic_plots.R`: creates OHCA-per-ICU-admission-rate lag-30 diagnostic summaries and figures using the total ICU admissions denominator.
9. `04b_case_crossover_dlnm.R`: fits time-stratified case-crossover DLNMs.
10. `10_ohca_icu_72h_phenotypes.R`: creates structured 72-hour ICU OHCA phenotypes and the multinomial phenotype-assignment model using admission-day temperature and humidity, with a mechanism-adjusted sensitivity.
11. `11_ohca_icu_competing_risks.R`: fits 72-hour Fine-Gray competing-risk models for time to awake-and-extubated recovery versus death before awake/extubated using admission-day temperature and humidity, with a mechanism-adjusted sensitivity.
12. `07_export_federated_results.R`: writes DLNM, ED-only death, and 72-hour phenotype aggregate site files plus selected site-level figure PNGs for federated pooling and visual QC.

## Optional Tables

The workflow can run with the required cohort/modeling tables only, but the 72-hour phenotype models are strongest when these optional CLIF tables are available:

- `clif_respiratory_support`: IMV trajectories, duration, 72-hour extubation evidence, and competing-risk awake/extubated timing.
- `clif_patient_assessments`: GCS, RASS, AVPU, SAT, and SBT evidence.
- `clif_medication_admin_continuous`: vasopressor prevalence and cardiovascular SOFA support.
- `clif_labs`: laboratory components of SOFA at 24, 48, and 72 hours.
- `clif_vitals`: last-vital fallback for death timing in the competing-risk model plus MAP, SpO2, and weight inputs for SOFA.

## Privacy Boundary

Only `output/final/federated_exports/` should be shared. This folder contains aggregate DLNM, ED-only death, and 72-hour phenotype CSVs plus selected site-level PNG figures. The scripts intentionally keep row-level CLIF-derived working files under `output/intermediate/`, which is git-ignored and should remain local.

## Post-Processing

After aggregate site exports are available, pooled analyses and manuscript figures can be run from `code/post_processing/`. For example:

```bash
Rscript code/post_processing/90_pool_federated_results.R
```

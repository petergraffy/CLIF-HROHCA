# Code Workflow

The site-facing workflow is R-only and uses numbered scripts in this folder.

Top-level scripts in `code/` are the site-running workflow. Current coordinating-center pooling, manuscript figures, and other post-processing helpers live in `code/post_processing/`.

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

### Aim 1 hospital-county pilot

Run the cohort exporter at each site using its CLIF tables and local configuration.
The optional `CLIF_SITE_NAME` and `CLIF_TABLES_PATH` environment variables allow
an explicit site and tables path for the export. Return the six aggregate CSVs
together. Keep the generated file names so the model runner can find the
coverage calendar next to the county-day file.

```bash
Rscript code/01_build_ohca_cohort.R --aim1-export output/final/aim1/SITE_county_day.csv
Rscript code/run_aim1.R output/final/aim1/models output/final/aim1/*_county_day.csv
```

| Site return | Contents |
|---|---|
| `SITE_county_day.csv` | Total adult ICU, probable OHCA and ED-pathway OHCA counts by hospital county and hospital admission date. |
| `SITE_county_day.csv.audit.csv` | Cohort totals and counts of unresolved or ambiguous county assignments and single-county fallbacks. |
| `SITE_county_day.csv.flow.csv` | Hospitalization and patient counts from the source extract through the study window, adult ICU cohort, cardiac-arrest screen, pathway and timing rules, and county assignment. Includes repeat hospitalizations among known patient IDs, missing patient IDs, ADT ICU segment counts, and the diagnostic table used. Timing anomalies are listed separately. |
| `SITE_county_day.csv.coverage.csv` | A complete study calendar for each county in the crosswalk, with observed ICU counts and a `site_status` column for sites to complete from reporting records. |
| `SITE_county_day.csv.transfers.csv` | Counts by county-day and admission source (`outside_hospital` for CLIF `osh`, `other_known`, or `unknown`). Includes acute-care transfer out, missing discharge information, and multiple or missing hospital IDs in ADT. Multiple ADT hospital IDs are a proxy for a within-system move. |
| `SITE_county_day.csv.subgroups.csv` | Total adult ICU, probable OHCA and ED-pathway OHCA counts by county-day for each age, sex, race and ethnicity variable. Each variable partitions the cohort separately. Unknown categories are retained. |

The export assigns the county of the first ICU hospital using `adt.hospital_id`
and groups admissions by local source hospital admission date. Hospitals in
the same county contribute to one county-day within each system. The transfer
and subgroup files also retain unresolved counties as missing FIPS values so
their site totals can be reconciled against the assignment audit.

For the coverage calendar, enter `reporting` for a confirmed complete day,
`incomplete` for a partial feed, `not_reporting` for an interruption, or
`inactive` if the hospital county was outside its operating or study coverage.
Use a blank or `unknown` when coverage cannot be established. With a calendar,
an empty county-day becomes zero only if explicitly marked `reporting`.
Observed days marked incomplete, not reporting, or inactive are excluded from
model fitting. Older count files without a calendar use the original pilot
rule for inferred zeros. Missing weather can interrupt a lag history.
The count model on positive-ICU days uses the same denominator eligibility as
the offset model. Set `AIM1_HOT_C` to evaluate a different hot temperature
against the 20°C reference. Use a distinct output directory for each contrast.
The support flag requires both temperatures within the fitted exposure range
and at least ten model observations within 1°C of each temperature.
The primary contrast is 30°C versus 20°C over lag days 0–5. Additional fixed
contrasts at 31–35°C use the same 20°C reference. The local-percentile model
uses each site's 95th percentile and median instead. Check the contributing
sites for each pooled result, particularly near the upper temperature tail.

`site_partial_lags.csv` and `pooled_partial_lags.csv` report cumulative
associations through each lag day from the full lag 0–5 fits. The separate
`lag_0_1` and `lag_0_2` models refit the OHCA count model using two and three lag
days, respectively, with a separate coefficient for each lag. Their shorter
history requirement can change the set of eligible model days.

The standard Aim 1 runner includes the daily PM2.5 sensitivity using the
2018–2024 LGHAP county-day Parquet files bundled in `exposome/lghap_pm25_county/`.
It reads the years used by the study and joins PM2.5 by county and date.
`AIM1_LGHAP_DIR` can override the bundled directory. Alternatively, set
`AIM1_PM25_CSV` to a CSV with `county_fips` (five digits), `date` (YYYY-MM-DD),
and `pm25_ug_m3`.

`pm25_sensitivity.csv` reports the OHCA count model with a three-degree-of-freedom
natural spline for same-day PM2.5 alongside its baseline fit on the same complete county-days.
The pooled comparison uses the same sites in both models. Coverage and input
checksums are recorded in `coverage.csv` and `input_manifest.csv`.
Set `AIM1_HISTORICAL_DIR` to the directory with
previous admission-rate exports to run the historical weather comparison.
`fixest` enables the optional time-stratified conditional Poisson sensitivity.
The runner writes site estimates, pooling and shared-date sensitivity results,
coverage, reduced coefficients and covariance matrices, model status, and
input checksums.

1. `00_install_or_restore_packages.R`: installs missing required R packages into the project-local library.
2. `01_build_ohca_cohort.R`: builds adult ICU OHCA hospitalizations using present-on-admission cardiac arrest diagnosis codes and applies hospital-aware county assignment.
3. `01b_ohca_ed_death_never_icu.R`: counts OHCA-diagnosis patients who arrived to the ED, stayed only in ED locations, died there, and never reached ICU.
4. `02_build_icu_exposure_series.R`: builds the all-ICU daily patient-address exposure series used as the comparison denominator for daily OHCA models.
5. `04_dlnm_primary_and_sensitivity.R`: fits primary, pollution-adjusted, MRT-reference, time-adjustment sensitivity, and stratified DLNM models, including cumulative lag-window and day-specific hot-temperature contrasts from lag 0 through lag 5.
6. `04c_lag30_diagnostic_plots.R`: creates count-model lag-30 diagnostic summaries and figures for lag-window justification.
7. `04d_ohca_icu_admission_rate_dlnm.R`: fits OHCA-per-ICU-admission-rate DLNMs with the same count-model denominator across strata.
8. `04e_rate_lag30_diagnostic_plots.R`: creates OHCA-per-ICU-admission-rate lag-30 diagnostic summaries and figures using the total ICU admissions denominator.
9. `04b_case_crossover_dlnm.R`: fits time-stratified case-crossover DLNMs.
10. `04f_case_crossover_lag30_diagnostic_plots.R`: creates time-stratified case-crossover lag-30 diagnostic summaries for lag-specific and cumulative lag-window manuscript figures.
11. `10_ohca_icu_72h_phenotypes.R`: creates structured 72-hour ICU OHCA phenotypes using awake/neuro evidence beginning at ICU hour 12, descriptive tables, admission-temperature density exports, and multinomial phenotype-assignment models using admission-day temperature and humidity plus 0-1, 0-3, and 0-5 day exposure-window sensitivities.
12. `11_ohca_icu_competing_risks.R`: fits 72-hour Fine-Gray competing-risk models for time to awake-and-extubated recovery versus death before awake/extubated, with awake/neuro evidence beginning at ICU hour 12, using admission-day temperature and humidity plus 0-1, 0-3, and 0-5 day exposure-window sensitivities.
13. `12_ohca_icu_imv24_time_to_event.R`: fits an IMV-24 landmark competing-risk analysis for successful extubation, death/hospice discharge, and tracheostomy after 24 hours of IMV, with admission temperature adjusted for age, sex, and race.
14. `13_ohca_icu_imv12_discharge_multinomial.R`: fits an IMV-duration >=12-hour multinomial model for final discharge outcome: death/hospice, LTACH with tracheostomy, alive discharge home, or alive discharge SNF.
15. `07_export_federated_results.R`: writes DLNM, ED-only death, and 72-hour phenotype aggregate site files plus selected site-level figure PNGs for federated pooling and visual QC.

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

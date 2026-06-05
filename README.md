# CLIF-HROHCA

Federated CLIF analysis of heat exposure and out-of-hospital cardiac arrest (OHCA) ICU admissions.

Each site runs the same R workflow locally against its own CLIF 2.1 tables. Sites should share only the aggregate outputs in `output/final/federated_exports/`.

## Study Aim

Estimate whether county-level heat exposure is associated with daily OHCA ICU admissions and, among OHCA patients admitted to the ICU, whether admission-day temperature and humidity are associated with 72-hour neurologic/recovery phenotypes.

The workflow produces:

- Daily OHCA ICU admission counts linked to county-level Tmax, relative humidity, NO2, and PM2.5.
- Primary and sensitivity distributed lag nonlinear models (DLNMs) for heat and OHCA ICU admissions.
- Site-level DLNM curves, reduced coefficient/vector covariance exports, and figure PNGs.
- Structured 72-hour ICU OHCA phenotypes.
- Multinomial phenotype-assignment models and 72-hour competing-risk models using admission-day temperature and humidity.
- Sensitivity versions of the 72-hour models additionally adjusted for POA diagnosis mechanism.

## Required Inputs

### CLIF 2.1 Tables

Required:

- `clif_patient`
- `clif_hospitalization`
- `clif_adt`
- `clif_hospital_diagnosis` or `clif_admission_diagnosis`

Recommended for 72-hour phenotype outputs:

- `clif_respiratory_support`: IMV and extubation evidence.
- `clif_patient_assessments`: GCS, RASS, AVPU, SAT, and SBT evidence.
- `clif_medication_admin_continuous`: vasopressor prevalence and cardiovascular SOFA support.
- `clif_labs`: laboratory components of SOFA at 24, 48, and 72 hours.
- `clif_vitals`: last-vital fallback for death timing in the competing-risk model plus MAP, SpO2, and weight inputs for SOFA.

### Exposome Files

The root `exposome/` folder contains county-level environmental inputs:

- `daymet_county_tmax_2018_2024_conus.parquet`
- `daymet_county_rmax_2018_2024.parquet`
- `no2_county_year.csv`
- `pm25_county_year.csv`

These are county-level files and contain no CLIF patient information.

### Hospital Geography

Hospital county assignment uses:

```text
reference/clif_hospital_geography.csv
```

For each hospitalization, the workflow keeps the patient county if it is the hospital county or an adjacent county. Missing or nonlocal patient counties are assigned to the admitting hospital county.

## Buddy-Test Quickstart

1. Clone the repository.

2. Copy the config template:

Windows PowerShell:

```powershell
Copy-Item config\config_template.json config\config.json
```

macOS/Linux:

```bash
cp config/config_template.json config/config.json
```

3. Edit `config/config.json`:

```json
{
  "site_name": "YOUR_SITE_NAME",
  "tables_path": "C:/path/to/local/CLIF/2.1/tables",
  "file_type": "parquet"
}
```

`site_name` must match a site in `reference/clif_hospital_geography.csv`. `file_type` can be `parquet` or `csv`; `fst` is not currently supported.

4. Install or restore R packages:

Windows PowerShell:

```powershell
Rscript code\00_install_or_restore_packages.R
```

macOS/Linux:

```bash
Rscript code/00_install_or_restore_packages.R
```

5. Run the full site workflow:

Windows PowerShell:

```powershell
Rscript code\run_site_analysis.R
```

macOS/Linux:

```bash
Rscript code/run_site_analysis.R
```

6. Review outputs locally:

```text
output/final/
```

7. Share only:

```text
output/final/federated_exports/
```

Do not share `output/intermediate/`; it contains CLIF-derived row-level working files and is git-ignored.

## Key Site Outputs

The federated export folder includes aggregate CSVs plus site-level figure PNGs.

Core analytic outputs:

- `SITE_dlnm_site_estimates.csv`
- `SITE_dlnm_curves.csv`
- `SITE_dlnm_lag_summaries.csv`
- `SITE_dlnm_lag_specific_summaries.csv`
- `SITE_dlnm_reduced_coefficients.csv`
- `SITE_dlnm_reduced_vcov.csv`
- `SITE_dlnm_time_sensitivity.csv`
- `SITE_dlnm_time_sensitivity_lag_summaries.csv`
- `SITE_dlnm_time_sensitivity_lag_specific_summaries.csv`
- `SITE_lag30_diagnostic_summary.csv`
- `SITE_lag30_temperature_lag_rr_surface.csv`
- `SITE_lag30_hot_temperature_lag_specific_rr.csv`
- `SITE_lag30_hot_temperature_cumulative_rr_by_lag.csv`
- `SITE_lag30_temperature_distribution.csv`
- `SITE_ohca_icu_admission_rate_dlnm_site_estimates.csv`
- `SITE_ohca_icu_admission_rate_dlnm_curves.csv`
- `SITE_ohca_icu_admission_rate_dlnm_lag_summaries.csv`
- `SITE_ohca_icu_admission_rate_dlnm_lag_specific_summaries.csv`
- `SITE_ohca_icu_admission_rate_dlnm_reduced_coefficients.csv`
- `SITE_ohca_icu_admission_rate_dlnm_reduced_vcov.csv`
- `SITE_ohca_icu_admission_rate_denominator_summary.csv`
- `SITE_ohca_icu_admission_rate_daily_timeseries.csv`
- `SITE_case_crossover_dlnm_site_estimates.csv`
- `SITE_case_crossover_dlnm_curves.csv`
- `SITE_case_crossover_dlnm_lag_summaries.csv`
- `SITE_case_crossover_dlnm_lag_specific_summaries.csv`
- `SITE_case_crossover_dlnm_reduced_coefficients.csv`
- `SITE_case_crossover_dlnm_reduced_vcov.csv`
- `SITE_case_crossover_referent_set_summary.csv`
- `SITE_ohca_icu_72h_consort_flow.csv`
- `SITE_ohca_icu_72h_table1.csv`
- `SITE_ohca_icu_72h_table2_by_phenotype.csv`
- `SITE_ohca_icu_72h_gcs_hourly_by_phenotype.csv`
- `SITE_ohca_icu_72h_phenotype_summary.csv`
- `SITE_ohca_icu_72h_ohca_mechanism_summary.csv`
- `SITE_ohca_icu_72h_phenotype_evidence_summary.csv`
- `SITE_ohca_icu_72h_phenotype_definitions.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_model.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_temperature_curve.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_coefficients.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_vcov.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_model_mechanism_adjusted.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_temperature_curve_mechanism_adjusted.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_coefficients_mechanism_adjusted.csv`
- `SITE_ohca_icu_72h_phenotype_assignment_vcov_mechanism_adjusted.csv`
- `SITE_ohca_icu_competing_risk_awake_extubated_72h_summary.csv`
- `SITE_ohca_icu_competing_risk_death_source_summary.csv`
- `SITE_ohca_icu_competing_risk_awake_extubated_72h_ohca_mechanism_summary.csv`
- `SITE_ohca_icu_competing_risk_awake_extubated_72h_fine_gray_models.csv`
- `SITE_ohca_icu_competing_risk_awake_extubated_72h_coefficients.csv`
- `SITE_ohca_icu_competing_risk_awake_extubated_72h_vcov.csv`
- `SITE_ohca_icu_competing_risk_awake_extubated_72h_fine_gray_models.csv`

Visual QC figures:

```text
output/final/federated_exports/figures/SITE_figure_*.png
```

These include site-level DLNM plots, lag-30 diagnostic plots, the OHCA rate time-series plot, and 72-hour phenotype plots.

## Troubleshooting

- If the workflow cannot find CLIF tables, confirm `config/config.json` and `file_type`.
- If `site_name` fails, confirm it exactly matches `reference/clif_hospital_geography.csv`.
- If optional respiratory-support or patient-assessment tables are missing, the 72-hour phenotype evidence may be partially populated.

## Privacy Boundary

Only share `output/final/federated_exports/`. The exported files are aggregate site-level outputs designed for federated pooling and visual QC. Do not share `output/intermediate/` or local CLIF tables.

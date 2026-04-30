# CLIF-HROHCA

Federated CLIF analysis of heat exposure and out-of-hospital cardiac arrest (OHCA) ICU admissions.

Each site runs the same R workflow locally against its own CLIF 2.1 tables. Sites should share only the aggregate outputs in `output/final/federated_exports/`.

## Study Aim

Estimate whether county-level heat exposure is associated with daily OHCA ICU admissions and characterize whether heat-related OHCA differs clinically from non-heat-related OHCA among patients who reach the ICU.

The workflow produces:

- Daily OHCA ICU admission counts linked to county-level Tmax, relative humidity, NO2, and PM2.5.
- Primary and sensitivity distributed lag nonlinear models (DLNMs) for heat and OHCA ICU admissions.
- Site-level DLNM curves, reduced coefficient/vector covariance exports, and figure PNGs.
- OHCA cohort descriptive tables and outcomes.
- Heat-related versus non-heat-related OHCA phenotype tables using 95th percentile heat as primary and 90th percentile heat as sensitivity.
- All-year sensitivity versions of heat-related versus non-heat-related OHCA phenotype tables and trajectories.
- ICU-hour trajectories for vitals, labs, organ support, and cumulative incidence.
- Renal/metabolic phenotype summaries, including CRRT initiation windows and early creatinine/BUN/electrolyte/lactate summaries.
- Supplementary OHCA outcome models for same-day heat and 12-month pollution exposures.
- CONSORT-style aggregate cohort flow counts and a site-level cohort flow diagram.

## Required Inputs

### CLIF 2.1 Tables

Required:

- `clif_patient`
- `clif_hospitalization`
- `clif_adt`
- `clif_hospital_diagnosis` or `clif_admission_diagnosis`

Recommended for phenotype/outcome outputs:

- `clif_respiratory_support`: IMV trajectories and duration.
- `clif_medication_admin_continuous`: vasopressor trajectories.
- `clif_labs`: renal/metabolic and lab trajectories.
- `clif_vitals`: vital sign trajectories.
- `clif_crrt_therapy`: CRRT initiation and trajectory outputs.

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
- `SITE_dlnm_reduced_coefficients.csv`
- `SITE_dlnm_reduced_vcov.csv`
- `SITE_dlnm_time_sensitivity.csv`
- `SITE_table1.csv`
- `SITE_outcomes.csv`
- `SITE_cohort_flow.csv`

Heat-related OHCA phenotype outputs:

- `SITE_heat_related_vs_non_heat_related_table.csv`
- `SITE_heat_related_vs_non_heat_related_table_all_definitions.csv`
- `SITE_heat90_vs_non_heat90_table.csv`
- `SITE_heat_related_hourly_vital_trajectories.csv`
- `SITE_heat_related_hourly_lab_trajectories.csv`
- `SITE_heat_related_hourly_support_trajectories.csv`
- `SITE_heat_related_hourly_vital_trajectories_smoothed.csv`
- `SITE_heat_related_hourly_lab_trajectories_smoothed.csv`
- `SITE_heat_related_hourly_support_trajectories_smoothed.csv`
- `SITE_heat_related_hourly_cumulative_incidence.csv`
- `SITE_heat_related_renal_metabolic_marker_summary.csv`
- `SITE_heat_related_crrt_window_summary.csv`

All-year heat-related OHCA sensitivity outputs use the same naming pattern with an `all_year_` prefix, for example:

- `SITE_all_year_heat_related_vs_non_heat_related_table.csv`
- `SITE_all_year_heat_related_hourly_vital_trajectories_smoothed.csv`
- `SITE_all_year_heat_related_hourly_cumulative_incidence.csv`
- `SITE_all_year_heat_related_crrt_window_summary.csv`

Outcome and sensitivity outputs:

- `SITE_adverse_outcome_models.csv`
- `SITE_continuous_outcome_models.csv`
- `SITE_pollution_12m_binary_outcome_models.csv`
- `SITE_pollution_12m_continuous_outcome_models.csv`
- `SITE_adverse_outcome_rates.csv`

Visual QC figures:

```text
output/final/federated_exports/figures/SITE_figure_*.png
```

These include site-level DLNM plots, trajectory plots, cumulative incidence plots, CRRT-window plots, renal/metabolic marker plots, and the cohort flow diagram.

## Troubleshooting

- If the workflow cannot find CLIF tables, confirm `config/config.json` and `file_type`.
- If `site_name` fails, confirm it exactly matches `reference/clif_hospital_geography.csv`.
- If optional tables are missing, the workflow should skip or partially populate the corresponding phenotype outputs.
- If figure PNGs look duplicated across heat90 and heat95, rerun the latest code; plot filtering was corrected to keep heat definitions separate.

## Privacy Boundary

Only share `output/final/federated_exports/`. The exported files are aggregate site-level outputs designed for federated pooling and visual QC. Do not share `output/intermediate/` or local CLIF tables.

# Code Workflow

The site-facing workflow is R-only and uses numbered scripts in this folder.

## Site Scripts

1. `00_install_or_restore_packages.R`: installs/restores required R packages into the project-local library and snapshots `renv.lock` when `renv` is available.
2. `01_build_ohca_cohort.R`: builds adult ICU OHCA hospitalizations using present-on-admission cardiac arrest diagnosis codes and applies hospital-aware county assignment.
3. `02_build_icu_exposure_series.R`: builds the all-ICU daily patient-address exposure series used as the comparison denominator for daily OHCA models.
4. `03_descriptive_tables.R`: creates aggregate OHCA cohort characteristics and outcomes.
5. `04_dlnm_primary_and_sensitivity.R`: fits primary, pollution-adjusted, MRT-reference, time-adjustment sensitivity, and stratified DLNM models.
6. `05_heat_related_vs_non_heat_related_table.R`: creates heat-related vs non-heat-related OHCA comparison tables, 90th-percentile sensitivity tables, ICU-hour clinical trajectories, cumulative incidence plots, and renal/metabolic phenotype summaries.
7. `09_supplementary_ohca_outcome_models.R`: creates supplementary outcome models for death/hospice, vasopressors, ICU LOS, and IMV duration.
8. `08_quality_checks.R`: creates aggregate denominator, admission-to-ICU timing, and care pathway quality checks.
9. `06_manuscript_tables_figures.R`: formats manuscript tables and figures.
10. `07_export_federated_results.R`: writes aggregate site files and site-level figure PNGs for federated pooling and visual QC.

Run everything with:

```r
Rscript code/run_site_analysis.R
```

## Coordinator Script

`90_pool_federated_results.R` pools site-level DLNM estimates after the coordinating center places all site `*_dlnm_site_estimates.csv` files in `output/final/federated_exports/`.

## Privacy Boundary

Only `output/final/federated_exports/` should be shared. The scripts intentionally keep row-level CLIF-derived working files under `output/intermediate/`, which is git-ignored and should remain local.

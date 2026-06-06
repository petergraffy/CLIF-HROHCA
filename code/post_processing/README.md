# Post-Processing Scripts

These scripts are coordinating-center or manuscript-preparation helpers. They are not required for a site to run the federated CLIF workflow.

Use this folder after site-level aggregate exports have been created under `output/final/federated_exports/`.

Primary pooled-output helpers:

- `90_pool_federated_results.R`: pools site-level DLNM estimates and curves.
- `115_pool_72h_phenotype_models.R`: pools 72-hour phenotype multinomial and Fine-Gray model exports.
- `116_figure_pooled_dlnm_curves.R`: creates pooled DLNM curve figures.
- `117_figure_pooled_72h_temperature_curves.R`: creates pooled 72-hour phenotype temperature-response figures.
- `118_figure_pooled_fine_gray_cifs.R`: creates pooled Fine-Gray cumulative-incidence figures.
- `119_export_pooled_ratio_tables.R`: exports pooled ratio tables for DLNM and 72-hour analyses.
- `120_figure_pooled_ratio_tables.R`: creates pooled ratio-table figures.

Manuscript figures:

- `130_manuscript_figure1_case_crossover_dlnm.R`
- `131_manuscript_figure2_case_crossover_stratified_forest.R`
- `132_manuscript_lag_justification_and_window_forest.R`
- `133_manuscript_rate_dlnm_sensitivity_figures.R`
- `134_manuscript_figure3_multinomial_phenotype_model.R`
- `135_manuscript_figure4_fine_gray_cifs.R`

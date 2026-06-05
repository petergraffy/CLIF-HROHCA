# Post-Processing Scripts

These scripts are coordinating-center or manuscript-preparation helpers. They are not required for a site to run the federated CLIF workflow.

Use this folder after site-level aggregate exports have been created under `output/final/federated_exports/`.

Primary pooled-output helpers:

- `90_pool_federated_results.R`: pools site-level DLNM estimates and curves.
- `91_compile_multisite_results.R`: compiles multisite aggregate exports.
- `92_map_county_mean_weather.R`: creates county-level weather summaries for mapping.

The remaining scripts build manuscript figures, descriptive pooled tables, sensitivity tables, and visual QC outputs from aggregate exports.

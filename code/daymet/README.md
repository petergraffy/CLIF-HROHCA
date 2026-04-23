# Daymet Workflow

This directory contains preprocessing code to build county-level heat exposure files for linkage to CLIF hospitalizations.

## Current script

- `01_download_daymet_county_tmax.R`
- `02_aggregate_daymet_county_tmax_local.R`

## What it does

1. Downloads annual Daymet V4 North America `tmax` files for 2018-2024.
2. Loads U.S. county polygons from `tigris`, unless a custom county boundary file is provided in `config/config.json`.
3. Computes county-level area-weighted daily mean Tmax in degrees C.
4. Writes per-year parquet outputs as each year finishes, then combines them into final outputs in `output/intermediate/daymet/processed`.

## Aggregation-only workflow

- `02_aggregate_daymet_county_tmax_local.R` skips all download/authentication logic.
- It is preconfigured for:
  - years `2018:2024`
  - variable `tmax`
  - region `na`
  - local raw folder `C:/Users/Peter Graffy/Downloads/Daymet_Daily_V4R1_2129_4.1-20260423_144640`
- You can still override the raw folder with `DAYMET_RAW_DIR` and the output folder with `DAYMET_PROCESSED_DIR`.

## Required config fields

The script uses the optional `daymet` section in `config/config.json`:

- `start_year`
- `end_year`
- `variable`
- `region`
- `dataset_id`
- `collection_id`
- `product_prefix`
- `direct_download_base`
- `county_boundaries_year`
- `county_boundaries_path`
- `county_geoid_path`
- `include_state_fips`
- `raw_dir`
- `processed_dir`

## Notes

- The script currently uses county-level mean `tmax` as the heat exposure metric.
- Daymet uses a 365-day calendar. In leap years, February 29 is retained and December 31 is omitted.
- To use Earthdata-authenticated download, set `EARTHDATA_USERNAME` and `EARTHDATA_PASSWORD` in the shell before running the script. The script creates temporary auth files at runtime and removes them after each download attempt.
- As of April 23, 2026, the script uses the Earthdata direct-download host published in the CMR virtual directory for the Daymet V4 R1 collection. The companion OPeNDAP metadata links remain useful for inspection, but the direct file host is the simpler path for full annual netCDF downloads.
- To use local netCDFs without downloading, set `DAYMET_RAW_DIR` to the folder containing the annual `.nc` files before running the script.
- During runtime, watch for:
  - `daymet_county_tmax_<year>.parquet` files appearing one by one
  - `daymet_county_tmax_2018_2024_progress.csv` updating after each year with status, elapsed time, and estimated time remaining

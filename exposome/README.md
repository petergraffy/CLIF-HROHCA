# Exposome Inputs

This folder contains the analysis-ready county-level environmental inputs used by the site workflow.

Required files:

- `daymet_county_tmax_2018_2024_conus.parquet`: county-day Daymet maximum temperature in Celsius.
- `daymet_county_rmax_2018_2024.parquet`: county-day Daymet maximum relative humidity percentage.
- `no2_county_year.csv`: county-year NO2 exposure.
- `pm25_county_year.csv`: county-year PM2.5 exposure.

These files do not contain CLIF patient-level data. Site scripts link local CLIF admissions to these county-level exposures using assigned county FIPS and hospital admission date.

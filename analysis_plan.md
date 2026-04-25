# CLIF-HROHCA Analysis Plan

## Primary Question

Is higher same-day ambient heat associated with increased daily ICU admissions for out-of-hospital cardiac arrest (OHCA)?

## Cohort

- Adult ICU admissions from CLIF 2.1
- Exploratory OHCA phenotype based on present-on-admission cardiac arrest diagnosis codes
- County assignment:
  - keep home county if it matches the admitting hospital county or an adjacent county
  - otherwise assign the admitting hospital county

## Primary Exposure

- Same-day county-linked `tmax`
- Warm season only: May through September

## Primary Model

- Warm-season DLNM
- Reference temperature: warm-season median `tmax`
- Primary adjustment set:
  - calendar time spline
  - day of week
  - year
  - same-day relative humidity max (`rmax`)

## Sensitivity Analyses

- DLNM additionally adjusted for county-year `no2` and `pm25`
- DLNM centered at minimum-risk temperature (MRT)
- Secondary time-stratified case-crossover analysis
- Stratified DLNMs by:
  - sex
  - age group `<65` vs `>=65`
  - race group `Black` vs `non-Black`

## Descriptive Outputs

- Cohort flow and summary
- Table 1 cohort characteristics
- ICU and hospital outcomes
- Hospital and site geography maps
- County exposure maps for `tmax`, `rmax`, `no2`, and `pm25`

## Heat-Related OHCA Comparison

- Exploratory comparison of OHCA admissions occurring on high-heat days versus other OHCA admissions
- High-heat day definition:
  - same-day county-linked `tmax` at or above the warm-season 95th percentile in the OHCA cohort

## Notes

- MRT will be reported as a sensitivity analysis, not the primary reference point.
- Pollution variables are included in time-series models but not the current time-stratified case-crossover models because county-year values are constant within matched strata.

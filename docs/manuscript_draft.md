# Heat Exposure, Out-of-Hospital Cardiac Arrest, and Post-Arrest Critical Illness Across CLIF Hospitals

## Title Page

**Title:** Heat Exposure, Out-of-Hospital Cardiac Arrest, and Post-Arrest Critical Illness Across CLIF Hospitals

**Short title:** Heat and out-of-hospital cardiac arrest

**Authors:** Peter S. [last name], [CLIF Heat-Related OHCA Collaborators]

**Affiliations:** [To be completed]

**Corresponding author:** [To be completed]

**Keywords:** out-of-hospital cardiac arrest; heat; extreme temperature; climate; critical care; distributed lag nonlinear model; environmental epidemiology

## Abstract

**Background:** Extreme heat is increasingly recognized as a threat to cardiovascular health, but its relationship with out-of-hospital cardiac arrest (OHCA) and the clinical phenotype of heat-related OHCA among patients who survive to hospital admission remains incompletely characterized.

**Objective:** To evaluate the association between county-level heat exposure and OHCA across Clinical Information Framework (CLIF) hospitals, and to compare clinical outcomes and early critical illness trajectories between heat-related OHCA (HROHCA) and non-heat-related OHCA.

**Methods:** We conducted a multisite federated analysis of adult OHCA hospitalizations across 8 CLIF hospitals from 2018 through 2024. Daily county-level maximum temperature (Tmax) and maximum relative humidity were linked to each site’s assigned patient county. The primary exposure-response analysis used site-level distributed lag nonlinear models (DLNMs) with cumulative associations expressed relative to the minimum-risk temperature (MRT). Site estimates were pooled using random-effects meta-analysis. HROHCA was defined for descriptive clinical comparisons as OHCA occurring on days with assigned-county Tmax at or above the site-specific warm-season 95th percentile; heat90 definitions were examined secondarily. Federated site summaries were used to compare mortality, discharge disposition, organ support, renal/metabolic markers, and hourly ICU trajectories between HROHCA and non-HROHCA.

**Results:** Across 8 sites, the pooled MRT-referenced DLNM demonstrated a positive heat-OHCA association. At 32.5 C, the cumulative relative risk (RR) was 1.34 (95% CI, 1.06-1.70), and at 35.0 C the RR was 1.98 (95% CI, 1.43-2.73). The pooled site-level hot-temperature contrast relative to MRT was RR 1.54 (95% CI, 1.26-1.89; p<0.001; I2=0%). In contrast, the median-reference primary model was attenuated (RR 1.15; 95% CI, 0.88-1.51). Among 4,913 OHCA hospitalizations contributing to heat95 descriptive comparisons, 251 were classified as HROHCA and 4,662 as non-HROHCA. Hospital mortality was similar between groups (53.8% vs 54.5%), as was death or hospice discharge (58.6% vs 58.3%). Adjusted adverse outcome models showed no clear association between heat95 exposure and death or hospice (OR 1.02; 95% CI, 0.70-1.49) or hospital death (OR 1.01; 95% CI, 0.75-1.38). Early ICU trajectories showed lower invasive mechanical ventilation prevalence among HROHCA patients at hours 6 and 12, but higher vasopressor prevalence by hour 48. CRRT use was numerically higher in HROHCA but not statistically distinct.

**Conclusions:** In this preliminary multisite federated CLIF analysis, high ambient temperature was associated with increased OHCA risk when modeled relative to each site’s MRT. However, among hospitalized OHCA patients, HROHCA and non-HROHCA had broadly similar mortality and organ-support outcomes, suggesting that heat may act more strongly as a trigger for OHCA occurrence than as a determinant of post-arrest hospital prognosis. These findings support further evaluation of heat-related OHCA risk using multicenter critical care data linked to environmental exposure.

## Introduction

Extreme heat is an increasingly important public health hazard as climate change increases the frequency, intensity, and duration of high-temperature events. Heat exposure has been associated with excess all-cause mortality, cardiovascular mortality, emergency medical services utilization, and hospitalizations for acute cardiopulmonary illness. The cardiovascular effects of heat may include dehydration, hemoconcentration, sympathetic activation, increased myocardial oxygen demand, electrolyte disturbance, renal injury, and destabilization of pre-existing cardiovascular disease. These mechanisms make out-of-hospital cardiac arrest (OHCA) a plausible heat-sensitive event.

Despite this biological plausibility, several gaps remain. First, many studies evaluate broad cardiovascular outcomes rather than OHCA specifically. Second, temperature-health associations may depend on local acclimatization, baseline climate, and nonlinear exposure-response relationships. Analyses using a single fixed heat threshold may therefore obscure risk patterns that differ across regions. Third, less is known about the hospitalized clinical phenotype of heat-related OHCA (HROHCA). If heat triggers OHCA through dehydration, metabolic stress, renal injury, or shock physiology, patients with HROHCA might have distinct early ICU trajectories, organ support needs, or outcomes compared with patients whose OHCA occurs outside high-heat conditions.

The CLIF consortium provides a federated critical care data infrastructure that allows multisite analyses while preserving local data governance. In this study, we linked site-level OHCA cohorts to county-level Daymet temperature and humidity data to evaluate two complementary questions. First, we estimated the association between daily heat exposure and OHCA using distributed lag nonlinear models (DLNMs), with the primary analysis referenced to each site’s minimum-risk temperature (MRT). Second, among patients hospitalized after OHCA, we compared clinical outcomes and early ICU trajectories between HROHCA and non-HROHCA.

## Methods

### Study Design and Setting

We conducted a retrospective, multisite federated cohort study across 8 CLIF hospitals contributing OHCA and environmental exposure summaries for calendar years 2018 through 2024. Each site executed a shared analytic pipeline locally and exported aggregate, deidentified results for central pooling. No patient-level data were transferred between institutions.

### Cohort Definition

The cohort included adult hospitalizations with OHCA present on admission who were admitted to an ICU during the index hospitalization. Adult age was defined as age 18 years or older at admission. OHCA identification used the shared CLIF project phenotype, including diagnosis and timing logic implemented at each site. Hospitalizations were linked to an assigned county FIPS code based on available patient county information; when patient county was unavailable or required site-specific reassignment, the analytic pipeline assigned hospital county according to project rules.

**Draft note:** Add exact OHCA phenotype definition, diagnosis code set, inclusion/exclusion criteria, and county reassignment rules from the final protocol.

### Environmental Exposures

Daily county-level mean maximum temperature (Tmax, degrees C) and maximum relative humidity were derived from Daymet for 2018-2024. Each OHCA hospitalization was assigned the daily environmental exposure for the patient’s assigned county and date. Additional county-year air pollution measures, including NO2 and PM2.5, were used in sensitivity and secondary models.

### Primary Heat-OHCA Association

The primary heat-OHCA analysis used site-level DLNMs to estimate the cumulative association between daily Tmax and OHCA counts. The main model adjusted for humidity and calendar/time structure according to the shared site pipeline. The primary exposure-response function was referenced to each site’s MRT, rather than to the median temperature, to better represent risk relative to the modeled nadir of the local temperature-OHCA curve. Site-specific log relative risks and standard errors were pooled using DerSimonian-Laird random-effects meta-analysis.

Pooled DLNM curves were generated by pooling site-specific log relative risks at each temperature grid point. The main manuscript should present the MRT-referenced curve as the primary figure. Median-reference and humidity-plus-pollution adjusted models were retained as sensitivity analyses.

**Draft note:** Verify and standardize wording around “MRT” before submission. The exported site metadata include site-specific reference/hot temperatures that vary across sites; the pooled MRT-referenced curve is currently the cleanest primary estimand.

### HROHCA Definition for Clinical Phenotyping

For clinical comparisons, HROHCA was defined as OHCA occurring on a day when assigned-county Tmax was at or above the site-specific warm-season 95th percentile. Non-HROHCA included OHCA hospitalizations below this threshold. A heat90 definition was evaluated secondarily.

### Clinical Outcomes and Trajectories

Clinical outcomes included hospital mortality, death or hospice discharge, invasive mechanical ventilation (IMV), vasopressor use, ICU length of stay, hospital length of stay, and IMV duration among ventilated patients. Early ICU trajectories included hourly vital signs, laboratory values, respiratory and hemodynamic support, CRRT use, and renal/metabolic marker summaries in 0-24 hour and 24-72 hour windows.

Because patient-level data were not centrally pooled, binary outcomes and support trajectories were summarized using pooled counts. Continuous trajectory summaries were pooled as sample-size-weighted site medians. Approximate p values for median-based trajectory contrasts were estimated from site-level medians, IQRs, and sample sizes and should be interpreted as descriptive rather than definitive patient-level tests.

### Statistical Analysis

DLNM and coefficient-based outcome models were pooled using random-effects meta-analysis. For HROHCA versus non-HROHCA descriptive comparisons, pooled event counts and rates were calculated across sites. Pooled support and CRRT trajectory differences used two-sample tests for equality of proportions. Median-based laboratory, vital sign, and renal/metabolic contrasts used approximate site-level median difference meta-analysis. All analyses were conducted using the shared CLIF project R pipeline.

## Results

### Site Contributions

Eight sites contributed DLNM, descriptive, clinical trajectory, and model outputs: Emory, JHU, NU, OHSU, Penn, RUMC, UCMC, and UMN. DLNM estimates were available from all 8 sites. Coefficient-based heat outcome model exports were also available from all 8 sites, although the number of estimable sites varied by outcome because some site-specific confidence intervals were nonfinite for selected models.

### Primary MRT-Referenced Heat-OHCA Association

The primary MRT-referenced DLNM showed a positive association between high Tmax and OHCA. At 30.0 C, the pooled cumulative RR was 1.01 (95% CI, 0.96-1.05). At 32.5 C, the RR increased to 1.34 (95% CI, 1.06-1.70; p=0.015), and at 35.0 C the RR was 1.98 (95% CI, 1.43-2.73; p<0.001). Heterogeneity was low to modest across the MRT-referenced curve, with I2=15.6% at 32.5 C and I2=0% at 35.0 C.

The pooled site-level MRT-referenced hot-temperature contrast was RR 1.54 (95% CI, 1.26-1.89; p<0.001; I2=0%). Female patients had a numerically larger heat association than male patients in stratified median-reference models (RR 1.32; 95% CI, 0.95-1.83 vs RR 1.06; 95% CI, 0.78-1.44), but confidence intervals overlapped. Associations were also numerically larger among patients aged 65 years or older than among younger patients, though imprecision limited inference.

In sensitivity analyses using median temperature as the reference, the overall heat association was attenuated and not statistically clear (RR 1.15; 95% CI, 0.88-1.51; p=0.304; I2=54.4%). Additional humidity-plus-pollution adjustment produced a similar median-reference estimate (RR 1.15; 95% CI, 0.88-1.51).

### HROHCA Versus Non-HROHCA Clinical Outcomes

Using the heat95 clinical phenotype definition, 251 of 4,913 OHCA hospitalizations were classified as HROHCA and 4,662 as non-HROHCA. Hospital mortality was similar in HROHCA and non-HROHCA (53.8% vs 54.5%), as was death or hospice discharge (58.6% vs 58.3%). IMV use was also similar (88.8% vs 88.6%), and vasopressor use was slightly higher among HROHCA patients (81.3% vs 80.7%).

Adjusted coefficient-based models did not show clear differences in adverse outcomes for heat95 exposure. The pooled odds ratio for death or hospice was 1.02 (95% CI, 0.70-1.49; p=0.915), and the pooled odds ratio for hospital death was 1.01 (95% CI, 0.75-1.38; p=0.929). Heat95 exposure was not clearly associated with IMV (OR 1.43; 95% CI, 0.77-2.68) or vasopressor use (OR 1.10; 95% CI, 0.74-1.65). Continuous outcome models also showed no clear difference in ICU length of stay (geometric mean ratio [GMR] 1.04; 95% CI, 0.84-1.28) or IMV duration (GMR 0.94; 95% CI, 0.80-1.09).

### Early ICU Support Trajectories

Hourly support trajectories suggested that HROHCA patients were less likely to be receiving IMV during the early ICU period. At ICU hour 6, IMV prevalence was 33.2% in HROHCA versus 41.4% in non-HROHCA (difference, -8.2 percentage points; p=0.013). At hour 12, IMV prevalence was 33.0% versus 42.8% (difference, -9.7 percentage points; p=0.004). By hour 24 and later, IMV differences were smaller and less precise.

In contrast, vasopressor infusion prevalence was higher among HROHCA patients by hour 48 (35.1% vs 27.5%; difference, +7.6 percentage points; p=0.041). CRRT prevalence was numerically higher among HROHCA patients at later hours but did not meet conventional statistical thresholds. For example, at hour 48, CRRT prevalence was 11.7% versus 7.6% (difference, +4.0 percentage points; p=0.069).

### Renal and Metabolic Trajectories

Windowed renal/metabolic marker summaries showed broadly similar profiles between HROHCA and non-HROHCA. In the first 24 hours, peak lactate was higher in HROHCA (5.70 vs 4.64), but the approximate site-level p value was not significant. In the 24-72 hour window, lowest bicarbonate was slightly lower among HROHCA patients (21.1 vs 21.7; approximate p=0.058). Peak creatinine, BUN, potassium, magnesium, and phosphate were otherwise similar across groups.

CRRT use over prespecified windows was numerically higher among HROHCA patients but not statistically distinct: 8.8% versus 7.6% at 0-24 hours (p=0.474), 12.9% versus 10.1% at 0-72 hours (p=0.165), and 13.3% versus 11.1% at 0-168 hours (p=0.292).

### Vital Sign and Laboratory Trajectories

Median vital sign trajectories did not reveal large or consistent differences between HROHCA and non-HROHCA. HROHCA patients had numerically higher heart rate at hour 72 (90.3 vs 84.7 bpm), but this difference was imprecise. Blood pressure, oxygen saturation, respiratory rate, and temperature trajectories were otherwise similar.

Hourly laboratory trajectories showed several nominal differences, including higher lactate at hour 24 among HROHCA patients (4.38 vs 2.29) and higher phosphate at hour 72 (4.48 vs 3.20). Because these analyses were based on pooled site-level medians with variable measurement density, these findings should be interpreted as descriptive signals rather than definitive patient-level biomarker associations.

### Secondary Pollution Models

In secondary pollution models, 12-month NO2 exposure per site-specific IQR was associated with death or hospice discharge in single-pollutant models (OR 1.04; 95% CI, 1.01-1.08; p=0.018). PM2.5 was not clearly associated with death or hospice discharge (OR 1.12; 95% CI, 0.91-1.38). Two-pollutant models attenuated these associations. Pollution-related continuous outcome models did not show clear associations with ICU length of stay or IMV duration.

## Discussion

In this multisite federated analysis of OHCA across CLIF hospitals, high ambient temperature was associated with increased OHCA risk when modeled relative to local minimum-risk temperature. The association was apparent on the pooled MRT-referenced exposure-response curve and was strongest at higher temperatures. In contrast, median-reference models were attenuated and less precise, underscoring the importance of local temperature-risk structure when evaluating heat-sensitive acute cardiovascular outcomes across geographically diverse hospitals.

The clinical phenotype analysis provides a complementary message. Although heat was associated with OHCA occurrence, patients classified as HROHCA did not have meaningfully higher hospital mortality, death or hospice discharge, ICU length of stay, or IMV duration than non-HROHCA patients. This suggests that heat may function primarily as an upstream trigger for cardiac arrest rather than as a strong determinant of post-arrest hospital prognosis among patients who survive to ICU admission. This distinction is important for both clinical interpretation and public health framing: heat mitigation may reduce OHCA incidence even if the in-hospital clinical course of HROHCA resembles other OHCA syndromes.

Several trajectory findings may help generate hypotheses about HROHCA pathophysiology. HROHCA patients had lower early IMV prevalence but higher later vasopressor prevalence, and numerically higher CRRT use. Laboratory summaries suggested possible metabolic stress, including higher lactate and modest bicarbonate differences, but these analyses were limited by federated summary data and measurement density. Together, these patterns could reflect heterogeneity in arrest etiology, postarrest shock, thermoregulatory physiology, or site-specific care pathways. They should be considered descriptive and hypothesis-generating.

The strongest methodologic implication is that the choice of temperature reference matters. Median-reference models may be intuitive, but they compare risk to a temperature that is not necessarily the nadir of the local exposure-response curve. In this analysis, MRT-referenced estimates were more coherent and showed lower heterogeneity. This supports using MRT-based curves as the primary estimand in multisite heat-health studies, especially when sites span different climate regions.

This study has several strengths. First, it uses a multisite federated critical care data infrastructure, enabling a geographically diverse analysis while preserving local data governance. Second, it links hospitalization-level OHCA phenotypes to county-level daily environmental exposures over multiple years. Third, it evaluates both population-level heat-OHCA associations and patient-level clinical phenotype summaries, bridging environmental epidemiology and post-arrest critical care.

Several limitations deserve emphasis. First, the current draft is based on aggregate federated outputs rather than centrally pooled patient-level data. This protects privacy but limits adjustment, subgroup analysis, and formal testing of continuous trajectory differences. Second, county-level exposure assignment may misclassify individual heat exposure, especially for patients with mobility, indoor cooling, occupational heat exposure, or county reassignment. Third, HROHCA was operationalized using same-day county Tmax thresholds; true heat-related pathophysiology may depend on multi-day exposure, humidity, nighttime temperature, or individual susceptibility. Fourth, patients who died before hospital arrival or were not admitted to an ICU are not represented in clinical trajectory analyses. Finally, site-specific model convergence and nonfinite confidence intervals limited the number of estimable sites for selected outcome models.

## Conclusions

Across 8 CLIF hospitals, high ambient temperature was associated with increased OHCA risk when modeled relative to local minimum-risk temperature. Among hospitalized OHCA patients, however, HROHCA and non-HROHCA had similar mortality and major clinical outcomes, with only modest differences in early ICU support and metabolic trajectories. These findings suggest that heat may be more important as a trigger of OHCA occurrence than as a determinant of post-arrest hospital prognosis. Future work should refine local heat-risk thresholds, incorporate humidity and multi-day exposure metrics, and evaluate whether individual susceptibility or neighborhood vulnerability modifies heat-related OHCA risk.

## Tables and Figures to Build

**Table 1.** Site contributions and cohort counts.

**Table 2.** Pooled MRT-referenced DLNM estimates and sensitivity models.

**Table 3.** HROHCA versus non-HROHCA outcomes using heat95 and heat90 definitions.

**Table 4.** Clinical trajectory contrasts: support, CRRT, renal/metabolic markers, selected labs/vitals.

**Figure 1.** Map of study counties or site geography with mean Tmax/humidity.

**Figure 2.** Pooled MRT-referenced DLNM curve for heat and OHCA.

**Figure 3.** Stratified DLNM forest plot by sex, age group, and race.

**Figure 4.** HROHCA versus non-HROHCA early ICU support trajectories.

**Supplement.** Median-reference DLNM curve, humidity-plus-pollution sensitivity, heat90 phenotype tables, site-specific estimates, and full trajectory summaries.

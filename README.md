# ENSO and temperature-related mortality: illustrative example (Brazil)

This repository contains an illustrative example of the second-stage analysis (ENSO regressions of attributable fractions and calculation of attributable numbers) for Brazil only. It is not the full pipeline of the paper.

Data file:
https://github.com/nataliashartova/ENSO/releases/tag/v1.0-example-data
contains prepared Brazil input data: `01_prepared_data_allsex_allage_allcause_iSEN.0_BRA_example.RData` 

The release file  was produced by a data-preparation step that is not included in this repository as restricted mortality data cannot be shared.
It was built from:
- first-stage (temperature-mortality) estimates: seasonal attributable fractions (AF) for heat and cold by region, with 1,000 simulations;
- seasonal mean temperature by region (1980-2023);
- the seasonal Nino 3.4 anomaly index (NOAA ERSST v6), with previous and following winter values;
- the region-level adjustment factor used to convert AF to attributable numbers;
- region geometries for Brazil.

After `load()`, the file provides these objects:
`df_enso2`, `polygon`, `nsims`, `sim_cols`, `cf_values`, `group_start`, `group_end`, `ngroups`.

`df_enso2` has one row per region x type (Total Heat / Total Cold) x season x
year, with columns: `REGION`, `TYPE`, `SEASON_REG`, `YEAR`, `att_val`
(point AF), `sim_1`...`sim_1000` (AF simulations), `temp`, `Nino34_0`,
`Nino34_prev`, `Nino34_next`, `continent_new`, `HDI_group`,
`adjusted_factor`.

This example reproduces the second-stage workflow for one country (regressions of AF on the Nino 3.4 index, predicted attributable numbers with confidence intervals for La Nina (-2.0), neutral (0) and El Nino (+2.5)).

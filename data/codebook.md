# Codebook: derived analytical dataset

Farooq Z, Semenza JC. *Spatial diffusion and climatic suitability shape West Nile virus invasion risk across Europe.* Repository release `v2-revision`.

------------------------------------------------------------------------

## Files

| File | Format | Content |
|---|---|---|
| `WNV_SurvivalDataset_withBalkans.csv` | CSV, UTF-8 | Region-year panel used to fit all Cox models |
| `analysis_cases_withBalkans.rds` | R tibble | The same panel (identical columns and values) |
| `spatial_objects_withBalkans.rds` | R list | Spatial objects for all 1088 cohort regions (see below) |

```r
d <- read.csv("data/WNV_SurvivalDataset_withBalkans.csv")
# or
d <- readRDS("data/analysis_cases_withBalkans.rds")
```

**Dimensions:** 16 998 rows × 27 columns\
**Unit of observation:** NUTS3 region (NUTS 2021) × calendar year at risk\
**Years:** 2007–2024 (regions followed from 2006; the 2006 row is not an at-risk year, because distance to the nearest invaded region is undefined before any region is invaded)\
**Regions / countries:** 1084 regions in 21 countries: AL, AT, BG, CY, CZ, DE, EL, ES, FR, HR, HU, IT, ME, MK, NL, PT, RO, RS, SI, SK, TR\
**Events:** 323 first invasions

The cohort contains 1088 regions. Four are absent from the panel: three regions invaded in 2006 (one Hungarian, two Romanian), which have no time at risk but act as sources for other regions throughout, and one Portuguese region without climatic covariates. Five French overseas regions were excluded before the cohort was defined. The 2025 and 2026 seasons, used only for the out-of-sample test, are not in this file.

------------------------------------------------------------------------

## Variables

### Identifiers

| Variable | Type | Description |
|---|---|---|
| `CNTR_CODE` | character | Eurostat country code (`EL` = Greece) |
| `NUTS_ID` | character | NUTS3 code, 2021 release |
| `NUTS_NAME` | character | Region name in the national language and script |

### Time structure

| Variable | Type | Description |
|---|---|---|
| `year` | integer | Calendar year, 2007–2024 |
| `start` | integer | Start of the counting-process interval, `year − 2006` (1–18) |
| `stop` | integer | End of the interval, `start + 1` (2–19) |

A region contributes one row per year until its first invasion or until 2024. A region never invaded by 2024 is right-censored at `stop = 19`.

### Outcome and case counts

| Variable | Type | Description |
|---|---|---|
| `event` | 0/1 | 1 in the year of first recorded locally acquired human WNV infection in the region, 0 otherwise |
| `first_year` | integer / NA | Year of first invasion; NA if not invaded by 2024 |
| `wnf_sum` | integer | Locally acquired human WNV infections reported in the region in that year (ECDC, with the additions documented in `provenance/`) |
| `outbreak` | 0/1 | 1 if `wnf_sum > 0` |
| `cases` | numeric | Case count in the source region `nearest_source` in the previous year; the input to `log_cases_nearest` |

### Diffusion metrics

All distances are between region centroids in EPSG:3035 and use only regions invaded **before** the focal year.

| Variable | Type | Units | Description |
|---|---|---|---|
| `nearest_invaded_km` | numeric | km | Distance to the nearest previously invaded region |
| `nearest_source` | character | — | `NUTS_ID` of that region |
| `distance100` | numeric | 100 km | `nearest_invaded_km / 100`; the form used in the models |
| `invasion_pressure` | numeric | — | Neighbourhood invasion pressure Φ = Σⱼ exp(−dᵢⱼ / λ) over previously invaded regions j, λ = 50 km |
| `pressure_z` | numeric | SD | `invasion_pressure` standardised to mean 0, SD 1 over this panel |
| `log_cases_nearest` | numeric | log scale | Source-region transmission intensity, log(1 + `cases`) |

### Climate

ISIMIP3a observational climate forcing (obsclim, 20crv3-era5, bias-adjusted), aggregated to NUTS3 and re-keyed to NUTS 2021 by `code/10_climate_nuts_crosswalk.R`.

| Variable | Type | Units | Description |
|---|---|---|---|
| `tas` | numeric | °C | Annual mean near-surface air temperature |
| `tasmax` | numeric | °C | Annual mean daily maximum temperature (not in the primary model) |
| `tasmin` | numeric | °C | Annual mean daily minimum temperature (not in the primary model) |
| `hurs` | numeric | % | Annual mean near-surface relative humidity |
| `pr` | numeric | as in the source file (flux basis, kg m⁻² s⁻¹) | Precipitation; enters the models only as `pr_z` |
| `pr_z` | numeric | SD | `pr` standardised to mean 0, SD 1 over this panel. Stored rather than computed with `scale()` in the formula, so that out-of-sample scoring uses the fitting data's mean and SD |

### Landscape

| Variable | Type | Description |
|---|---|---|
| `URBN_TYPE` | integer | Eurostat urban-rural typology: 1 = predominantly urban, 2 = intermediate, 3 = predominantly rural. Entered as `factor(URBN_TYPE)`, reference 1 |
| `COAST_TYPE` | integer | Eurostat coastal typology (sensitivity analyses only, Supplementary Table S13) |
| `MOUNT_TYPE` | integer | Eurostat mountain typology (sensitivity analyses only, Supplementary Table S13) |

### Period

| Variable | Type | Description |
|---|---|---|
| `invasion_phase` | factor | Early (2007–2010), Middle (2011–2015), Pre-COVID (2016–2019), Post-COVID (2020–2024) |

------------------------------------------------------------------------

## `spatial_objects_withBalkans.rds`

A list over all 1088 cohort regions, in the row order of `region_lookup`:

| Element | Class | Content |
|---|---|---|
| `region_lookup` | data.frame, 1088 × 1 | `NUTS_ID` |
| `region_index` | integer, 1088 | Row index of each region |
| `dist_mat` | matrix, 1088 × 1088 | Centroid-to-centroid distance, km, EPSG:3035 |
| `W` | matrix, 1088 × 1088 | Decay weights exp(−d / 50), diagonal set to 0 |
| `first_outbreak` | tibble, 1088 × 2 | `NUTS_ID`, `first_year` (2006–2024; NA if not invaded by 2024) |

------------------------------------------------------------------------

## Missing data

The panel has no missing values in `distance100`, `invasion_pressure`, `tas`, `pr` or `hurs`; rows lacking any of these were removed during preparation (`code/00_data_preparation_withBalkans.R`, section 14). `first_year` is NA by design for regions not invaded by 2024.

------------------------------------------------------------------------

*Codebook version 2.0, September 2026*

# Data Codebook — analysis_cases.rds

## WNV Spatial Diffusion Europe

## Zia Farooq, Umeå University

------------------------------------------------------------------------

## Dataset overview

**File:** `data/analysis_cases.rds`\
**Format:** R data frame (compressed RDS)\
**Dimensions:** 16,318 rows × 25 columns\
**Unit of observation:** NUTS3 administrative region × calendar year\
**Study period:** 2006–2024\
**Countries:** 18 European countries with at least one confirmed locally acquired human WNV infection during the study period\
**Exclusions:** NUTS3 regions outside continental Europe (overseas territories)

Load in R with: `analysis_cases <- readRDS("data/analysis_cases.rds")`

------------------------------------------------------------------------

## Variable definitions

### Region identifiers

| Variable | Type | Description |
|----|----|----|
| `CNTR_CODE` | character | ISO 3166-1 alpha-2 country code (e.g. "DE", "IT", "EL") |
| `NUTS_ID` | character | Eurostat NUTS3 region identifier (e.g. "AT111", "DE212") |
| `NUTS_NAME` | character | NUTS3 region name in local language |

------------------------------------------------------------------------

### Temporal structure

| Variable | Type | Units | Description |
|----|----|----|----|
| `year` | integer | Calendar year | Observation year (2007–2024; baseline t₀ = 2006) |
| `start` | integer | Years since 2006 | Start of counting-process interval (startᵢₜ = t − 2006) |
| `stop` | numeric | Years since 2006 | End of counting-process interval (stopᵢₜ = startᵢₜ + 1) |

------------------------------------------------------------------------

### Outcome variable

| Variable | Type | Values | Description |
|----|----|----|----|
| `event` | numeric (0/1) | 0 = no first invasion, 1 = first invasion | First recorded occurrence of at least one confirmed locally acquired human WNV infection in a previously disease-free region. Set to 1 only in the calendar year of first invasion. Regions remaining disease-free throughout are right-censored (event = 0 for all years). |
| `wnf_sum` | integer | ≥ 0 | Raw count of confirmed locally acquired human WNV cases in the NUTS3 region in the focal year |
| `outbreak` | numeric (0/1) | 0 = no outbreak, 1 = outbreak | Binary indicator of any WNV cases in the focal year (not restricted to first occurrence) |
| `first_year` | integer / NA | Year, or NA | Calendar year of first recorded WNV invasion in the region. NA for regions that remained disease-free throughout the study period. |

------------------------------------------------------------------------

### Spatial diffusion metrics

| Variable | Type | Units | Description |
|----|----|----|----|
| `nearest_invaded_km` | numeric | km | Euclidean distance from the NUTS3 region centroid to the centroid of the nearest region with documented WNV transmission prior to the focal year. Calculated in EPSG:3035 (Lambert Azimuthal Equal Area). Updated annually as new regions are invaded. |
| `distance100` | numeric | per 100 km | Geographic distance to nearest invaded region, scaled to units of 100 km (= nearest_invaded_km / 100). This is the form used in all Cox models; HRs are interpreted per 100 km increase. |
| `invasion_pressure` | numeric | dimensionless | Raw neighbourhood invasion pressure index (Φᵢₜ). Calculated as the distance-decayed sum of contributions from all regions with documented WNV transmission prior to year t: Φᵢₜ = Σⱼ exp(−dᵢⱼ / λ), where λ = 50 km (profile-likelihood optimised). Right-skewed distribution; standardised version (pressure_z) used in models. |
| `pressure_z` | numeric | SD units | Neighbourhood invasion pressure standardised to mean = 0, SD = 1. Used in all Cox models; HRs are interpreted per one standard deviation increase in the unstandardised index (Φ). |
| `cases` | numeric | count | Reported human WNV case count in the NUTS3 region in the focal year (used to compute log_cases_nearest for neighbouring regions in subsequent years). |
| `log_cases_nearest` | numeric | log-scale | Natural log of (1 + case count) in the nearest previously invaded region in year t−1. The log(1+x) transformation handles zero-case years. Used as the source-region transmission intensity metric (Sᵢₜ). |

------------------------------------------------------------------------

### Climatic covariates (ERA5, 2006–2024)

All climatic variables derived from the ERA5 global atmospheric reanalysis dataset (Hersbach et al. 2020, Q J Roy Meteor Soc) at 0.25° × 0.25° spatial resolution, spatially aggregated to NUTS3 using area-weighted averaging.

| Variable | Type | Units | Description |
|----|----|----|----|
| `tas` | numeric | °C | Annual mean air temperature (2 metres above surface). Primary indicator of climatic suitability for mosquito development and viral replication. Used in original units in Cox models. |
| `tasmax` | numeric | °C | Annual mean daily maximum air temperature. Not used in primary models; available for sensitivity analyses. |
| `tasmin` | numeric | °C | Annual mean daily minimum air temperature. Not used in primary models; available for sensitivity analyses. |
| `pr` | numeric | mm/s | Annual total precipitation rate. Standardised using scale(pr) in Cox models (mean = 0, SD = 1). |
| `hurs` | numeric | \% | Annual mean near-surface relative humidity. Used in original units (%) in Cox models. |

------------------------------------------------------------------------

### Landscape and geographic covariates

| Variable | Type | Values | Description |
|----|----|----|----|
| `URBN_TYPE` | integer | 1 = urban, 2 = intermediate, 3 = rural | Urbanisation category derived from the Eurostat urban-rural typology at NUTS3 scale. Time-invariant covariate. Entered as factor(URBN_TYPE) in models with urban (1) as the reference category. |
| `COAST_TYPE` | numeric | Coded | Coastal typology classification (Eurostat). Not used in primary models. |
| `MOUNT_TYPE` | numeric | Coded | Mountain/topographic classification (Eurostat). Not used in primary models. |

------------------------------------------------------------------------

### Analytical groupings

| Variable | Type | Values | Description |
|----|----|----|----|
| `invasion_phase` | factor | Early (2006–2010), Middle (2011–2015), Pre-COVID (2016–2019), Post-COVID (2020–2024) | Temporal phase of WNV expansion. Defined a priori based on distinct epidemiological periods in European WNV history. Used in phase-stratified PH sensitivity check (Analysis 2.11) and Supplementary Figure S4. |

------------------------------------------------------------------------

## Notes on data structure

**Counting-process format:** The dataset is structured in the (start, stop] counting-process format required by the `survival` package in R. Each region contributes one row per calendar year until either (i) first invasion or (ii) end of study period (2024). The `stop` column equals `start + 1` for all observations.

**Right-censoring:** Regions remaining disease-free through 2024 have `event = 0` for all rows. These are treated as right-censored at `stop = 18` (year 2024 = t₀ + 18).

**Temporal causality:** Both `distance100` (Dᵢₜ) and `invasion_pressure` / `pressure_z` (Φᵢₜ) are computed exclusively from regions with documented WNV transmission **before** the focal year t (using I\_{t-1}). This ensures that future invasion events cannot influence past exposure estimates.

**Subset for descriptive analyses:** `first_invaded <- analysis_cases %>% filter(event == 1)` gives the 269 first-invasion events used for Figures 4A, S3, S4, and the KW tests.

------------------------------------------------------------------------

## Missing data

No missing values are expected in the primary analytic variables (distance100, pressure_z, log_cases_nearest, tas, hurs, pr, URBN_TYPE) for the 16,318 region-years included in the analysis. Regions with incomplete covariate data were excluded during data preparation (see `00_data_preparation.R`).

------------------------------------------------------------------------

## Provenance

This derived dataset was constructed from the sources listed in the README.md. The full construction pipeline is documented in `code/00_data_preparation.R`. Key construction decisions:

-   **Invasion definition:** First calendar year with ≥ 1 confirmed locally acquired human WNV infection (ECDC definition)
-   **Spatial projection:** EPSG:3035 (Lambert Azimuthal Equal Area) for all distance calculations
-   **Decay parameter:** λ = 50 km, selected by profile-likelihood minimising AIC across 8 candidate values
-   **ERA5 aggregation:** Area-weighted averaging using exactextractr package

------------------------------------------------------------------------

*Codebook version 1.0 — July 2025*

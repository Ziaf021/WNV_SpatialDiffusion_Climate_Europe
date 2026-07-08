# Spatial diffusion and climatic suitability shape the West Nile Virus invasion risk across Europe

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R
version](https://img.shields.io/badge/R-%3E%3D4.5.0-blue.svg)](https://cran.r-project.org/)

**Author:** Zia Farooq\
**Affiliation:** Department of Epidemiology and Global Health, Umeå
University, Umeå, Sweden\
**Contact:** [zia.farooq\@umu.se](mailto:zia.farooq@umu.se){.email}\
**Manuscript status:** Under review, *Nature Communications*

------------------------------------------------------------------------

## Overview

This repository contains all analysis code and derived data necessary to
reproduce the primary findings reported in:

> Farooq Z. *Spatial diffusion and climatic suitability shape the West
> Nile Virus invasion risk across Europe.*

West Nile virus (WNV) has undergone one of the most remarkable
geographic expansions in modern European epidemiology, yet the
mechanisms governing its invasion of previously disease-free
regions remain poorly quantified. This study develops a spatial
diffusion framework, integrated with time-to-event survival modelling,
to quantify how geographic proximity to existing viral foci, cumulative
neighbourhood exposure to surrounding invaded regions, and climatic
suitability jointly govern regional WNV emergence across 1,011 NUTS3
administrative regions in 18 European countries between 2006 and 2024.

**Key findings:**

- 71% of first emergence events occurred within 150 km of an existing
  focus, with invasion risk declining sharply at greater distances
- Each 100 km of geographic isolation reduced the hazard of invasion by
  44% on average across the study period (HR = 0.56, 95% CI 0.50–0.62)
- Neighbourhood invasion pressure independently increased invasion
  hazard by 26% on average across the study period, per standard
  deviation (HR = 1.26, 95% CI 1.17–1.35)
- Warmer temperatures eroded the protective barrier of geographic
  isolation, enabling viral establishment at distances that would
  otherwise remain protected (interaction HR = 1.07, 95% CI 1.04–1.09)

------------------------------------------------------------------------

## Repository structure

```         
WNV_SpatialDiffusion_Europe/
│
├── README.md                          # This file
├── LICENSE                            # MIT licence
│
├── code/
│   ├── 00_data_preparation.R          # Event-history construction, spatial
│   │                                  # metrics, ERA5 climate aggregation
│   ├── 02_survival_analyses.R         # Univariable Cox, spline, KM curves,
│   │                                  # primary multivariable model
│   ├── 03_interaction_models.R        # Distance × temperature interaction,
│   │                                  # Figure 4A–B, Kruskal–Wallis tests
│   ├── 04_sensitivity_prespecified.R  # λ calibration, alternative kernels,
│   │                                  # shared frailty, climate anomaly, LOCO
│   ├── 05_ph_diagnostics.R            # Schoenfeld residuals, outlier check,
│   │                                  # time-varying coefficients, case-
│   │                                  # weighted kernels, phase stratification
│   └── 06_supplementary_figures.R     # Supplementary figures, verification summary
│
├── data/
│   ├── analysis_cases.rds             # Main longitudinal panel dataset
│   │                                  # (16,318 region-years, 269 events)
│   ├── dist_mat.rds                   # Pairwise NUTS3 centroid distance
│   │                                  # matrix (km, EPSG:3035)
│   └── codebook.md                    # Variable definitions and units
│
├── outputs/
│   ├── figures/                       # All main-text and supplementary figures
│   └── tables/                        # Supplementary tables S1–S10
│
└── renv.lock                          # Package version snapshot (renv)
```

------------------------------------------------------------------------

## Data

### Primary data sources

| Source | Description | Access |
|------------------------|------------------------|------------------------|
| [ECDC WNV Surveillance](https://www.ecdc.europa.eu/en/west-nile-virus-infection) | Annual human WNV case counts at NUTS3 level, 2006–2024 | Freely downloadable |
| [ERA5 Reanalysis](https://cds.climate.copernicus.eu/) | Annual mean temperature, precipitation, relative humidity at 0.25° × 0.25° | Freely downloadable (registration required) |
| [Eurostat NUTS3 boundaries](https://ec.europa.eu/eurostat/web/gisco/geodata) | NUTS3 administrative shapefiles, 2021 release, EPSG:3035 | Freely downloadable |
| [Eurostat urban-rural typology](https://ec.europa.eu/eurostat/statistics-explained/index.php/Urban-rural_typology) | NUTS3 urbanisation category (urban / intermediate / rural) | Freely downloadable |

### Derived analytical dataset (`data/analysis_cases.rds`)

The primary analytical dataset is a longitudinal panel of 16,318 NUTS3
region-year observations across 18 European countries (2006–2024),
constructed from the raw data sources above. It is provided as a
compressed R data object and can be loaded directly into R without
re-running the full data preparation pipeline.

**Raw ECDC surveillance data are not redistributed in this repository**
in accordance with ECDC terms of use. Raw data can be downloaded
directly from the ECDC website. The `00_data_preparation.R` script
documents all processing steps required to reproduce the derived dataset
from the raw inputs.

### Variable codebook

Full variable definitions are provided in `data/codebook.md`. Key
variables:

| Variable | Type | Description |
|------------------------|------------------------|------------------------|
| `NUTS_ID` | character | NUTS3 region identifier |
| `CNTR_CODE` | character | ISO country code (2-letter) |
| `year` | integer | Calendar year (2007–2024) |
| `event` | binary | First locally acquired human WNV infection (1 = yes) |
| `start` / `stop` | integer | Counting-process interval boundaries |
| `distance100` | numeric | Distance to nearest invaded region (units of 100 km) |
| `invasion_pressure` | numeric | Raw neighbourhood invasion pressure index (Φ, dimensionless) |
| `pressure_z` | numeric | Standardised invasion pressure (mean=0, SD=1) |
| `log_cases_nearest` | numeric | log(1 + cases) in nearest invaded source region |
| `tas` | numeric | Annual mean air temperature (°C, ERA5) |
| `pr` | numeric | Annual total precipitation (mm/s, ERA5) |
| `hurs` | numeric | Annual mean relative humidity (%, ERA5) |
| `URBN_TYPE` | integer | Urbanisation category (1=urban, 2=intermediate, 3=rural) |
| `invasion_phase` | factor | Study phase (Early/Middle/Pre-COVID/Post-COVID) |

------------------------------------------------------------------------

## Reproducing the analyses

### Prerequisites

**R version:** ≥ 4.5.0\
**Package management:** This repository uses
[`renv`](https://rstudio.github.io/renv/) for reproducible package
environments. All package versions are recorded in `renv.lock`.

**Required R packages:**

``` r
# Core analysis
library(survival)    # Cox PH models, survfit, cox.zph
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping (LOCO analysis)

# Spatial
library(sf)          # Centroid computation, distance calculation

# Climate data extraction
library(terra)       # ERA5 raster handling
library(exactextractr) # Area-weighted spatial aggregation to NUTS3

# Visualisation
library(ggplot2)
library(ggthemes)

# Reproducibility
library(renv)
library(broom)
```

**Install all packages from the lockfile:**

``` r
install.packages("renv")
renv::restore()
```

### Running order

Scripts must be run in numerical order. Each script saves intermediate R
objects (`.rds` files) that are loaded by subsequent scripts.

``` r
# Step 1: Build derived dataset from raw inputs
source("code/00_data_preparation.R")

# Step 2: Univariable Cox, spline, KM, primary multivariable model
source("code/02_survival_analyses.R")

# Step 3: Interaction models, Figure 4A–B, KW tests
source("code/03_interaction_models.R")

# Step 4: Five pre-specified sensitivity analyses
source("code/04_sensitivity_prespecified.R")

# Step 5: Post-hoc PH diagnostics
source("code/05_ph_diagnostics.R")

# Step 6: Supplementary figures + verification summary
source("code/06_supplementary_figures.R")
```

**If you already have `data/analysis_cases.rds`**, you can skip
`00_data_preparation.R` and begin at step 2. Each script includes
internal safety checks and prints expected values from the manuscript
alongside computed values for direct comparison.

### Expected run times (approximate, on a standard laptop)

| Script | Approximate time |
|------------------------------------|------------------------------------|
| 00_data_preparation.R | 10–30 min (ERA5 extraction is the bottleneck) |
| 02_survival_analyses.R | \< 2 min |
| 03_interaction_models.R | \< 2 min |
| 04_sensitivity_prespecified.R | 5–15 min (LOCO × 18 countries) |
| 05_ph_diagnostics.R | 3–5 min (time-varying model) |
| 06_supplementary_figures.R | \< 1 min |

------------------------------------------------------------------------

## Key results verification

The final script (`06_supplementary_figures.R`) prints a complete
verification summary comparing every key manuscript number against the
computed output. Expected values for the primary model are:

```         
PRIMARY MODEL (cox_source_full):
  Distance HR:    0.56 (0.50–0.62) → 44% reduction per 100 km, on average across the study period
  Pressure HR:    1.26 (1.17–1.35) → 26% increase per SD, on average across the study period
  Temperature HR: 1.25 (1.19–1.31) → 25% increase per 1°C
  Humidity HR:    0.96 (0.93–0.99) → p = 0.014
  Concordance:    0.875

INTERACTION:
  Distance × Temperature HR: 1.07 (1.04–1.09) → p < 0.001

SPATIAL SCALE:
  Within 75 km:  121/269 (45.0%)
  Within 150 km: 192/269 (71.4%)
  Within 300 km: 236/269 (87.7%)
```

------------------------------------------------------------------------

## Citation

If you use this code or data in your research, please cite:

**Paper:** Farooq Z. Spatial diffusion and climatic suitability shape
the West Nile Virus invasion risk across Europe. *Nature Communications*
(under review, 2026).

**Code and data:** Farooq Z. (2026). *WNV Spatial Diffusion Europe*
[Code and data]. Zenodo. <https://doi.org/10.5281/zenodo.XXXXXXX>

------------------------------------------------------------------------

## Related publications

This work extends and complements the following published studies:

-   Farooq Z, Segelmark L, Rocklöv J, et al. Impact of climate and
    *Aedes albopictus* establishment on dengue and chikungunya outbreaks
    in Europe: a time-to-event analysis. *Lancet Planetary Health* 2025;
    9:e374–e383.
-   Farooq Z, Rocklöv J, Wallin J, et al. Artificial intelligence to
    predict West Nile virus outbreaks with eco-climatic drivers. *Lancet
    Regional Health – Europe* 2022; 17:100370.
-   Farooq Z, Sjödin H, Semenza JC, et al. European projections of West
    Nile virus transmission under climate change scenarios. *One Health*
    2023; 16:100509.

------------------------------------------------------------------------

## Licence

This repository is licensed under the [MIT Licence](LICENSE). You are
free to use, modify, and distribute the code with appropriate
attribution. The derived analytical dataset (`data/analysis_cases.rds`)
is shared under [CC BY
4.0](https://creativecommons.org/licenses/by/4.0/).

------------------------------------------------------------------------

## Contributing and issues

This repository is associated with a submitted manuscript and is not
actively accepting contributions. If you identify an error or have a
question about the analysis, please open a [GitHub
Issue](https://github.com/Ziaf021/WNV_SpatialDiffusion_Climate_Europe/issues)
or contact the author directly at
[zia.farooq\@umu.se](mailto:zia.farooq@umu.se){.email}.

------------------------------------------------------------------------

## Acknowledgements

The author thanks the European Centre for Disease Prevention and Control
(ECDC) for providing publicly available WNV surveillance data, the
Copernicus Climate Change Service for ERA5 reanalysis data, and Eurostat
for NUTS3 administrative boundary shapefiles. This work was supported by
the Horizon Europe project IDAlert.

------------------------------------------------------------------------

*Last updated: July 2026*

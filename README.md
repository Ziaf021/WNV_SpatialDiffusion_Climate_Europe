# Spatial diffusion and climatic suitability shape West Nile virus invasion risk across Europe

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D4.5.0-blue.svg)](https://cran.r-project.org/)

**Authors:** Zia Farooq, Jan C. Semenza\
**Affiliations:** Department of Epidemiology and Global Health, Umeå University, Umeå, Sweden (ZF, JCS); Heidelberg Institute of Global Health, University of Heidelberg, Heidelberg, Germany (JCS)\
**Contact:** [zia.farooq@umu.se](mailto:zia.farooq@umu.se)\
**Manuscript status:** revised version under review at *The Lancet Regional Health – Europe* (September 2026)

------------------------------------------------------------------------

## Versions

| Tag | Content |
|---|---|
| `v2-revision` | **This version.** Code and derived data behind the revised manuscript: 21 countries, 1088 NUTS3 regions followed 2006 to 10 September 2026, model fitted 2007–2024, 2025 and 2026 held out. |
| `v1-original-submission` | The original submission: 18 countries, 1050 regions, 16 318 region-years, 269 first invasions, fitted 2006–2024. Superseded; kept for the record. |

------------------------------------------------------------------------

## Overview

We followed 1088 NUTS3 regions in the 21 European countries that have reported at least one locally acquired human West Nile virus (WNV) infection, from 2006 to 2026. Three diffusion metrics were derived for each region and year: distance to the nearest previously invaded region, neighbourhood invasion pressure (a distance-decayed sum over all invaded regions, λ = 50 km) and source-region transmission intensity. Time to first invasion was modelled with a counting-process Cox model over 2007–2024, adjusted for temperature, humidity, precipitation and urbanisation. The fitted model was then applied, without refitting, to the 2025 and 2026 seasons.

Of the 1088 regions, 1084 entered the analysis: three were invaded in 2006, the first year of follow-up, and contribute no time at risk, and one Portuguese region has no climatic covariates. The fitted cohort comprises 16 998 region-years and 323 first invasions.

**Main results (primary model, 2007–2024)**

| Quantity | Estimate |
|---|---|
| Distance to nearest invaded region, HR per 100 km | 0·58 (0·53–0·64) |
| Neighbourhood invasion pressure, HR per SD | 1·35 (1·26–1·45) |
| Source-region intensity, HR per unit log(cases + 1) | 1·39 (1·19–1·61) |
| Annual mean temperature, HR per °C | 1·20 (1·15–1·26) |
| Distance × temperature interaction, HR | 1·06 (1·04–1·08) |
| Harrell's C | 0·871 |
| Out-of-sample AUC, full model | 0·66 (2025), 0·82 (2026) |

Of the 323 first invasions in the fitted period, 239 (74·0%) occurred within 150 km of an existing focus; across 2007–2026, 285 of 381 first invasions with a source (74·8%) did so.

------------------------------------------------------------------------

## Repository structure

```
WNV_SpatialDiffusion_Climate_Europe/
├── README.md
├── LICENSE
├── run_revision_v3.sh              Full pipeline driver (steps 1–23)
├── run_fig.sh                      Builds one figure at a time
├── run_withBalkans_v2.R            Step 2: data preparation and model fit (00, 02, 03)
├── run_supplementary_v2.R          Step 5: sensitivity analyses and diagnostics (04, 05, 06)
├── merge_balkan_cases.py           Step 1: merges the Western Balkans series into the case file
├── code/
│   ├── 00_data_preparation_withBalkans.R   Event history, diffusion metrics, climate join
│   ├── 02_survival_analyses.R              Univariable models, distance spline, KM, primary model
│   ├── 03_interaction_models.R             Distance × temperature interaction, LR test for intensity
│   ├── 04_sensitivity_prespecified_02.R    λ profile, kernels, frailty, climatic anomaly, leave-one-country-out
│   ├── 05_ph_diagnostics_2.R               Schoenfeld tests, time-varying and phase-stratified models
│   ├── 06_supplementary_figures.R          Supplementary Figures S5 and S6
│   ├── 08_revision_analyses.R              Region area, boundary distance, pressure deciles,
│   │                                       mountain/coast, cross-border, direction of invasion
│   ├── 09_prospective_validation_v2.R      Out-of-sample test, 2025 and 2026
│   ├── 10_climate_nuts_crosswalk.R         Re-keys the climate file to NUTS 2021
│   ├── 12_regional_composition.R           Regional composition of first invasions
│   ├── 12b_direction_latitude_4periods.R   Latitude and regional composition, four periods
│   ├── 13_table1_univariable.R             Table 1
│   ├── _wnv_style.R, _wnv_figdata.R        Shared figure style and data loader
│   ├── fig01*.R … fig05*.R                 Main-text figures
│   ├── figS10_directionality_rose.R        Supplementary Figure S10
│   └── fig_urbanisation_*.R, urbanisation_phase_enrichment.R
│                                           Supplementary Figures S7–S9
└── data/
    ├── WNV_SurvivalDataset_withBalkans.csv   Derived analytical dataset (CSV)
    ├── analysis_cases_withBalkans.rds        The same dataset as an R object
    ├── spatial_objects_withBalkans.rds       Distance matrix, decay weights, first-invasion years
    ├── codebook.md                           Variable definitions
    └── provenance/
        ├── ecdc_patch_rows.csv               The 81 region-years added from ECDC tables (Appendix S3)
        ├── WNV_Serbia_ECDC_verification.xlsx
        ├── WNV_Balkans_ECDC_verification.xlsx
        └── sessionInfo_v3_refit.txt          R session of the model refit
```

------------------------------------------------------------------------

## Data

### Sources

| Source | Content | Access |
|---|---|---|
| [ECDC WNV surveillance](https://www.ecdc.europa.eu/en/west-nile-virus-infection) | Locally acquired human WNV infections by NUTS3 region and season; region-level web tables 2011–2026 | Public; subject to ECDC terms of reuse |
| [ISIMIP3a obsclim](https://doi.org/10.48364/ISIMIP.982724) (20crv3-era5, bias-adjusted) | Annual temperature, precipitation and relative humidity, 1991–2024, aggregated to NUTS3 | Public |
| [Eurostat GISCO](https://ec.europa.eu/eurostat/web/gisco/geodata) | NUTS 2021 boundaries (EPSG:3035), retrieved with the `eurostat` R package | Public |
| [Eurostat urban-rural typology](https://ec.europa.eu/eurostat/statistics-explained/index.php/Urban-rural_typology) | NUTS3 urbanisation category | Public |

**Raw ECDC case files are not redistributed.** The repository holds the analysis code and the derived analytical dataset only. The derived dataset contains the annual case count for each region-year (`wnf_sum`, `cases`), because the source-intensity metric is built from it.

### Derived analytical dataset

`data/WNV_SurvivalDataset_withBalkans.csv` (and the identical `analysis_cases_withBalkans.rds`) holds one row per region and year at risk: 16 998 rows, 27 columns, 1084 regions, 21 countries, 323 events, 2007–2024. Variables are defined in [`data/codebook.md`](data/codebook.md). The 2025 and 2026 seasons used for the out-of-sample test are not part of this file.

### Provenance of the case data

- **ECDC extract repair (Appendix, Supplementary Methods S3).** The NUTS3 extract lacked records that ECDC had received under pre-2021 NUTS codes (mainly Croatia, Budapest and Pest, France and Sardinia) and held no Turkish record before 2024. Comparing it season by season with ECDC's published tables for 2011–2024 added 81 region-years; no existing count was changed. Each added row, with its ECDC source, is listed in `data/provenance/ecdc_patch_rows.csv`.
- **Western Balkans.** Region-level series for Serbia, Albania, North Macedonia and Montenegro were compiled from ECDC publications. The two verification workbooks in `data/provenance/` list the source of every row and the NUTS3 or GAUL unit it was assigned to. ECDC published no region-level table for these countries in 2020, 2022 and 2023, and Serbian cases only at country level in 2017.

------------------------------------------------------------------------

## Quick check from the derived dataset

The primary model can be refitted from the CSV alone:

```r
library(survival)
d <- read.csv("data/WNV_SurvivalDataset_withBalkans.csv")

fit <- coxph(Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
               tas + hurs + pr_z + factor(URBN_TYPE),
             data = d, ties = "efron")
summary(fit)

fit_int <- coxph(Surv(start, stop, event) ~ distance100 * tas + pressure_z +
                   log_cases_nearest + hurs + pr_z + factor(URBN_TYPE),
                 data = d, ties = "efron")
summary(fit_int)$conf.int["distance100:tas", ]
```

Expected output:

```
distance100         0.580  (0.528–0.637)
pressure_z          1.353  (1.260–1.453)
log_cases_nearest   1.388  (1.193–1.615)
tas                 1.203  (1.148–1.259)
hurs                0.973  (0.948–0.998)
pr_z                0.914  (0.806–1.037)
URBN_TYPE 2         0.789  (0.589–1.056)
URBN_TYPE 3         0.724  (0.532–0.986)
Concordance         0.871
distance100:tas     1.060  (1.037–1.084)

Distance at first invasion (nearest_invaded_km, event == 1, n = 323):
  within  75 km   150/323  (46·4%)
  within 150 km   239/323  (74·0%)
  within 300 km   287/323  (88·9%)
  beyond 500 km     9/323   (2·8%)
  median           80·1 km
```

------------------------------------------------------------------------

## Reproducing the full analysis

### Requirements

- R ≥ 4.5.0 with `survival`, `dplyr`, `tidyr`, `tidyverse`, `sf`, `units`, `eurostat`, `broom`, `ggplot2`, `ggthemes`, `ggpubr`, `patchwork`, `cowplot`, `ggnewscale`, `scales`, `forcats`, `survminer`, `officer`, `flextable`; optionally `svglite`, `ragg` and `ggspatial`. The versions used are listed in `data/provenance/sessionInfo_v3_refit.txt`.
- Python 3 (standard library only) for step 1.
- Internet access on first run, to download NUTS 2021 geometry through `eurostat`.

### Inputs you must supply

Place these in the repository root. They are not included for licensing reasons.

| File | Content |
|---|---|
| `WNV_NUTS3cases_2006to2025.csv` | ECDC NUTS3 case extract, 2006–2025, with the rows in `data/provenance/ecdc_patch_rows.csv` added |
| `WNV_Serbia_NUTS3cases_2006to2025.csv`, `WNV_Balkans_NUTS3cases_2006to2025.csv` | Western Balkans series (sources in the verification workbooks) |
| `NUTS3ObsClimate_1991to2024.csv` | ISIMIP3a obsclim variables aggregated to NUTS3 |
| `WNV_2025_ECDC_regions.csv`, `WNV_2026_ECDC.csv` | ECDC region-level tables for the 2025 season and for 2026 to 10 September |

### Running

```bash
Rscript code/10_climate_nuts_crosswalk.R
bash run_revision_v3.sh
```

`10_climate_nuts_crosswalk.R` writes the NUTS 2021-keyed climate file to `run_withBalkans/output/revision/`. `run_revision_v3.sh` then runs, in order: the Balkans merge (step 1), data preparation and model fit (2), revision analyses (3), out-of-sample test (4), sensitivity analyses and diagnostics (5), Table 1 (6), regional composition (7), direction and latitude (8), each main and supplementary figure (9–22) and the urbanisation figures (23). All output goes to `run_withBalkans_v3/`. A failed run can be resumed with `START=<step> bash run_revision_v3.sh`.

A single figure is rebuilt with, for example, `WNV_RUN_DIR=run_withBalkans_v3 bash run_fig.sh 4`; `bash run_fig.sh list` shows the keys. The scripts read the environment variables `WNV_RUN_DIR`, `WNV_EXPECTED_EVENTS`, `WNV_CLIMATE_FILE`, `WNV_LAMBDA` and `WNV_FIG4A_VIOLIN`.

### Where each result comes from

| Manuscript item | Script | Output (under `run_withBalkans_v3/output/`) |
|---|---|---|
| Table 1 | `13_table1_univariable.R` | `revision/rev_table1_univariable.*` |
| Primary and interaction models | `02`, `03` via `run_withBalkans_v2.R` | `multivariable_model_*.csv`, `headline_numbers_*.csv` |
| Figure 1 | `fig01a_facet_periods.R`, `fig01b_country_year.R`, `fig01f_combined.R` | `revision/figures/Figure1_Facet_AB.*` |
| Figure 2 | `fig02a_*`, `fig02b_*`, `fig02_combined.R` | `revision/figures/Figure2_Main_AB.*` |
| Figure 3 | `fig03_forest_spline_km.R` | `revision/figures/Figure3_forest_spline_km.*` |
| Figure 4 | `fig04a_*`, `fig04b_*`, `fig04_combined.R` | `revision/figures/Figure4_Main_AB.*` |
| Figure 5, out-of-sample AUC | `09_prospective_validation_v2.R`, `fig05_out_of_sample.R` | `revision/rev_prospective_*.csv`, `revision/figures/Figure5_out_of_sample.*` |
| Supplementary Figures S2, S3; Tables S1–S5 | `04_sensitivity_prespecified_02.R` | `supplementary/FigureS1_lambda_profile.pdf`, `FigureS_LOCO_DotPlot.pdf`, `Supplementary_Tables.docx` |
| Supplementary Figure S4; Tables S6–S10 | `05_ph_diagnostics_2.R` | `supplementary/FigureS2_Schoenfeld.pdf`, `Supplementary_Tables_PH.docx` |
| Supplementary Figures S5, S6 | `06_supplementary_figures.R` | `supplementary/Supplementary_Figure_S3_TempByPressure.pdf`, `…_S4_DistanceByPhase.pdf` |
| Supplementary Figures S7, S8, S9 | `fig_urbanisation_firstinvaded.R`, `fig_urbanisation_by_phase.R`, `urbanisation_phase_enrichment.R` | `Fig_Urbanisation_*.pdf` |
| Supplementary Figure S10 | `figS10_directionality_rose.R` | `revision/figures/FigureS10_directionality_rose.*` |
| Supplementary Tables S11–S13, S15 | `08_revision_analyses.R` | `revision/rev_area_*.csv`, `rev_pressure_decile_rates.csv`, `rev_mountain_coast_models.csv`, `rev_boundary_distance_models.csv` |
| Supplementary Table S14 and its note | `08_revision_analyses.R`, `12_regional_composition.R`, `12b_direction_latitude_4periods.R` | `revision/rev_directionality*.csv`, `rev_latitude_change_2007_2026.csv`, `rev_regional_*_4p.csv` |

Output file names from `04`, `05` and `06` keep the analysis numbering used while the code was written; the table above gives the corresponding figure numbers in the revised appendix. `fig01a_chronology_map.R` and `fig01_combined.R` build an alternative, unfaceted layout of Figure 1 that was not submitted.

------------------------------------------------------------------------

## Citation

Farooq Z, Semenza JC. Spatial diffusion and climatic suitability shape West Nile virus invasion risk across Europe. Under review, 2026.

Code and data: Farooq Z, Semenza JC. WNV_SpatialDiffusion_Climate_Europe, release `v2-revision` (2026). <https://github.com/Ziaf021/WNV_SpatialDiffusion_Climate_Europe>

## Related publications

- Farooq Z, Segelmark L, Rocklöv J, et al. Impact of climate and *Aedes albopictus* establishment on dengue and chikungunya outbreaks in Europe: a time-to-event analysis. *Lancet Planet Health* 2025; 9: e374–83.
- Farooq Z, Rocklöv J, Wallin J, et al. Artificial intelligence to predict West Nile virus outbreaks with eco-climatic drivers. *Lancet Reg Health Eur* 2022; 17: 100370.
- Farooq Z, Sjödin H, Semenza JC, et al. European projections of West Nile virus transmission under climate change scenarios. *One Health* 2023; 16: 100509.

## Licence

Code: [MIT](LICENSE). Derived dataset and provenance files: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Case counts originate from ECDC and remain subject to ECDC's terms of reuse.

## Issues

Errors or questions about the analysis: open a [GitHub issue](https://github.com/Ziaf021/WNV_SpatialDiffusion_Climate_Europe/issues) or write to [zia.farooq@umu.se](mailto:zia.farooq@umu.se).

## Acknowledgements

We thank the European Centre for Disease Prevention and Control for its public WNV surveillance data, the ISIMIP project for the observational climate forcing, and Eurostat for the NUTS boundaries. This work was supported by the Horizon Europe project IDAlert.

------------------------------------------------------------------------

*Last updated: September 2026*

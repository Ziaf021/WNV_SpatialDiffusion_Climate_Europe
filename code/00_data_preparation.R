###############################################################################
# 01_data_prep.R
#
# Foundational data preparation for:
#   "Spatial diffusion and climatic suitability shape the West Nile Virus
#    invasion risk across Europe"
#   Author: Dr. Zia Farooq, Umea University
#
# Design:
#   Outcome  - first WNV outbreak in a NUTS3 region (time-to-event)
#   Exposure - lagged invasion pressure from previously invaded regions
#   Format   - counting process (start, stop, event) per region-year
#
# Objects exposed for downstream scripts (02-06, figures):
#   analysis_cases  - modelling data frame (see column contract below)
#   dist_mat        - region x region distance matrix (km); dimnames = NUTS_ID
#   region_lookup   - data.frame(NUTS_ID); row order matches dist_mat
#   region_index    - named vector NUTS_ID -> column index in dist_mat
#   first_outbreak  - data.frame(NUTS_ID, first_year)
#   W               - connectivity kernel used to build invasion_pressure
#
# analysis_cases column contract expected downstream:
#   NUTS_ID, CNTR_CODE, year, event, start, stop,
#   distance100, nearest_invaded_km, invasion_pressure, pressure_z,
#   log_cases_nearest, tas, hurs, pr, URBN_TYPE, invasion_phase
###############################################################################


## ---- 0. Configuration -------------------------------------------------------

data_dir   <- "data"
output_dir <- "output"

case_file    <- file.path(data_dir, "WNV_NUTS3cases_2006to2025.csv")
climate_file <- file.path(data_dir, "NUTS3ObsClimate_1991to2024.csv")

# PRIMARY invasion-pressure kernel (keep consistent with 04_sensitivity).
lambda_primary        <- 50    # exponential decay length (km)
distance_threshold_km <- 1000   # links beyond this truncated to 0

# Analysis window (inclusive).
analysis_year_min <- 2007
analysis_year_max <- 2024

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


## ---- 1. Packages ------------------------------------------------------------

library(tidyverse)   # dplyr, tidyr, readr, purrr, ggplot2, ...
library(sf)
library(eurostat)
# units::drop_units() is called via namespace; package must be installed.


## ---- 2. Load case data ------------------------------------------------------

dat <- read.csv(case_file, stringsAsFactors = FALSE) %>%
  select(-any_of("geometry"))


## ---- 3. NUTS3 geography (Eurostat 2021, level 3) ----------------------------

EU_NUTS3 <- get_eurostat_geospatial(
  resolution = 10,
  nuts_level = 3,
  year       = 2021
) %>%
  filter(!grepl("^FRY", NUTS_ID))

nuts_sf <- EU_NUTS3 %>%
  select(NUTS_ID, CNTR_CODE, URBN_TYPE, COAST_TYPE, MOUNT_TYPE, geometry)


## ---- 4. Restrict to countries that ever reported WNV ------------------------

dat <- dat %>%
  mutate(outbreak = ifelse(wnf_sum > 0, 1, 0))

wnv_countries <- dat %>%
  group_by(CNTR_CODE) %>%
  summarise(ever_outbreak = max(outbreak), .groups = "drop") %>%
  filter(ever_outbreak == 1) %>%
  pull(CNTR_CODE)

dat     <- dat     %>% filter(CNTR_CODE %in% wnv_countries)
nuts_sf <- nuts_sf %>% filter(CNTR_CODE %in% wnv_countries)


## ---- 5. First-outbreak year per region --------------------------------------

first_outbreak <- dat %>%
  group_by(NUTS_ID) %>%
  summarise(
    first_year = if (any(outbreak == 1)) min(year[outbreak == 1]) else NA_real_,
    .groups = "drop"
  )


## ---- 6. At-risk (survival) dataset + event indicator ------------------------

surv_dat <- dat %>%
  left_join(first_outbreak, by = "NUTS_ID") %>%
  filter(is.na(first_year) | year <= first_year) %>%
  mutate(event = ifelse(!is.na(first_year) & year == first_year, 1L, 0L))

stopifnot(!any(is.na(surv_dat$event)))


## ---- 7. Attach static region attributes -------------------------------------

surv_dat <- surv_dat %>%
  select(-any_of(c("URBN_TYPE", "COAST_TYPE", "MOUNT_TYPE"))) %>%
  left_join(
    nuts_sf %>%
      st_drop_geometry() %>%
      select(NUTS_ID, URBN_TYPE, COAST_TYPE, MOUNT_TYPE),
    by = "NUTS_ID"
  )


## ---- 8. Centroids, distance matrix, connectivity kernel ---------------------

nuts_proj <- nuts_sf %>%
  st_transform(3035) %>%
  st_make_valid()

centroids <- st_centroid(nuts_proj)

dist_mat <- st_distance(centroids)
dist_mat <- as.matrix(units::drop_units(dist_mat)) / 1000   # -> km

W <- exp(-dist_mat / lambda_primary)
W[dist_mat > distance_threshold_km] <- 0
diag(W) <- 0


## ---- 9. Region lookup / index (row-order contract for dist_mat) -------------

region_lookup <- nuts_sf %>%
  st_drop_geometry() %>%
  select(NUTS_ID)

region_index <- seq_len(nrow(region_lookup))
names(region_index) <- region_lookup$NUTS_ID

rownames(dist_mat) <- region_lookup$NUTS_ID
colnames(dist_mat) <- region_lookup$NUTS_ID
stopifnot(
  nrow(dist_mat) == nrow(region_lookup),
  all(rownames(dist_mat) == region_lookup$NUTS_ID)
)

first_outbreak <- first_outbreak %>%
  filter(NUTS_ID %in% region_lookup$NUTS_ID)


## ---- 10. Annual loop: invasion pressure, nearest distance, nearest source ---
# One pass yields all three, so nearest_invaded_km and the source region whose
# cases define log_cases_nearest are guaranteed to be the same region.

years <- sort(unique(surv_dat$year))
pressure_list <- vector("list", length(years))

for (i in seq_along(years)) {
  
  yy <- years[i]
  cat("Processing year:", yy, "\n")
  
  invaded_regions <- first_outbreak %>%
    filter(!is.na(first_year), first_year < yy) %>%
    pull(NUTS_ID)
  
  invaded_idx <- na.omit(region_index[invaded_regions])   # named positions
  
  if (length(invaded_idx) == 0) {
    pressure       <- rep(0,  nrow(W))
    nearest_dist   <- rep(NA_real_, nrow(W))
    nearest_source <- rep(NA_character_, nrow(W))
  } else {
    sub            <- dist_mat[, invaded_idx, drop = FALSE]
    pressure       <- rowSums(W[, invaded_idx, drop = FALSE])
    nearest_pos    <- max.col(-sub, ties.method = "first")   # argmin per row
    nearest_dist   <- sub[cbind(seq_len(nrow(sub)), nearest_pos)]
    nearest_source <- names(invaded_idx)[nearest_pos]
  }
  
  pressure_list[[i]] <- data.frame(
    NUTS_ID            = region_lookup$NUTS_ID,
    year               = yy,
    invasion_pressure  = pressure,
    nearest_invaded_km = nearest_dist,
    nearest_source     = nearest_source,
    stringsAsFactors   = FALSE
  )
}

pressure_df <- bind_rows(pressure_list)


## ---- 11. Merge pressure; counting-process intervals; standardise ------------

surv_dat <- surv_dat %>%
  left_join(pressure_df, by = c("NUTS_ID", "year")) %>%
  arrange(NUTS_ID, year) %>%
  group_by(NUTS_ID) %>%
  mutate(
    start = year - min(year),
    stop  = start + 1
  ) %>%
  ungroup()



## ---- 12. Nearest-source case intensity (log_cases_nearest) ------------------
# log_cases_nearest = log1p(cases in the nearest prior-invaded region, year-1).
#
# NOTE: case_series defaults to surv_dat, which is truncated at each region's
# first-outbreak year; a source therefore contributes cases only up to its own
# invasion year and 0 afterwards. This feeds the PRIMARY model, so it affects
# reported HRs. Set case_series <- dat to use the full observed case series.


## ---- 12. Nearest-source case intensity (log_cases_nearest) ------------------
# log_cases_nearest = log1p(cases in the nearest prior-invaded region, year-1).
#
# Source intensity uses the FULL observed case series (dat), not surv_dat.
# Rationale: surv_dat's truncation at each region's first-outbreak year is an
# OUTCOME (censoring) rule for the focal region; log_cases_nearest is an
# EXPOSURE measured on a different region (the nearest prior-invaded source).
# Using dat gives the source's actual prior-year burden; the one-year lag
# prevents temporal leakage and the source/focal separation prevents
# circularity. dat is already restricted to WNV countries (section 4) and is
# otherwise untruncated, which is exactly what this definition requires.

case_series <- dat

region_cases <- case_series %>%
  group_by(NUTS_ID, year) %>%
  summarise(cases = max(wnf_sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(cases = ifelse(is.finite(cases), cases, 0))   # guard -Inf (all-NA group)

# Lag one year and rename key to match nearest_source.
region_cases_lag <- region_cases %>%
  mutate(year = year + 1) %>%
  rename(nearest_source = NUTS_ID)

surv_dat <- surv_dat %>%
  left_join(region_cases_lag, by = c("nearest_source", "year")) %>%
  mutate(
    # NA here = the nearest source has no record at year-1 (genuinely no
    # observed activity), correctly treated as zero.
    cases             = replace_na(cases, 0),
    log_cases_nearest = log1p(cases)
  )


## ---- 13. Attach observed climate --------------------------------------------

climate_surv <- read.csv(climate_file, stringsAsFactors = FALSE) %>%
  filter(scenario == "obsclim", year >= 2006) %>%
  select(NUTS_ID, year, tas, tasmax, tasmin, pr, hurs)

surv_dat <- surv_dat %>%
  left_join(climate_surv, by = c("NUTS_ID", "year")) %>%
  mutate(distance100 = nearest_invaded_km / 100)


## ---- 14. Build analysis_cases -----------------------------------------------
# Analysis window + completeness. The distance100 filter conditions the
# diffusion model on an existing invaded source.

analysis_cases <- surv_dat %>%
  filter(year >= analysis_year_min, year <= analysis_year_max) %>%
  filter(
    !is.na(distance100),
    !is.na(invasion_pressure),
    !is.na(tas),
    !is.na(pr),
    !is.na(hurs)
  ) %>%
  mutate(
    invasion_phase = factor(
      case_when(
        year <= 2010 ~ "Early",
        year <= 2015 ~ "Middle",
        year <= 2019 ~ "Pre-COVID",
        TRUE         ~ "Post-COVID"
      ),
      levels = c("Early", "Middle", "Pre-COVID", "Post-COVID")
    )
  )

analysis_cases$pressure_z <- as.numeric(scale(analysis_cases$invasion_pressure))

urbn_levels <- levels(factor(analysis_cases$URBN_TYPE))
urbn_ref    <- urbn_levels[which.max(tabulate(match(analysis_cases$URBN_TYPE, urbn_levels)))]  # modal category

## ---- 15. Assertions (frozen-dataset anchors) --------------------------------

n_events <- sum(analysis_cases$event)
message("Events in analysis_cases: ", n_events)
stopifnot(n_events == 269L)

# URBN_TYPE NAs are dropped by coxph(factor(URBN_TYPE)), reducing model N below
# the 269 events; report so the discrepancy is not silent.
n_urbn_na <- sum(is.na(analysis_cases$URBN_TYPE))
if (n_urbn_na > 0) {
  message("URBN_TYPE missing on ", n_urbn_na,
          " analysis rows (these drop out of models using factor(URBN_TYPE)).")
}


## ---- 16. Persist outputs ----------------------------------------------------
# Delete this block if downstream scripts source() this file in a shared session.

write.csv(
  analysis_cases,
  file.path(output_dir, "WNV_SurvivalDataset.csv"),
  row.names = FALSE
)

saveRDS(analysis_cases, file.path(output_dir, "analysis_cases.rds"))
saveRDS(
  list(
    dist_mat       = dist_mat,
    region_lookup  = region_lookup,
    region_index   = region_index,
    first_outbreak = first_outbreak,
    W              = W
  ),
  file.path(output_dir, "spatial_objects.rds")
)

###############################################################################
# End of 01_data_prep.R
###############################################################################
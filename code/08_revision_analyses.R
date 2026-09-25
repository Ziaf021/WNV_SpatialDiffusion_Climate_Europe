###############################################################################
# 08_revision_analyses.R
#
# Analyses requested by the referees of TLRHEUROPE-D-26-01167 that were not in
# the submitted paper. Everything here reads the objects produced by
# run_withBalkans.R and writes to run_withBalkans/output/revision/. No existing
# file is modified.
#
# Sections
#   A  Region area, and the area distribution by country            (editor P2)
#   B  Mountain and coast terms, distance x mountain interaction    (Rev 1.7)
#   C  Area-adjusted and area-stratified models                     (Rev 1.1, P2)
#   D  Boundary-to-boundary distance, full refit of the primary model (Rev 1.1, P2)
#   E  Direction of invasion: bearings, Rayleigh test, rose table   (Rev 3.3)
#   F  Share of first invasions whose nearest source is abroad      (Rev 4.2)
#   G  Invasion rate by neighbourhood-pressure decile               (Rev 1.5)
#
# Run from the project folder:
#   Rscript code/08_revision_analyses.R
#
# Section D builds a 1049 x 1049 polygon-to-polygon distance matrix and is the
# slow step, roughly 10 to 40 minutes depending on the machine. It is cached in
# run_withBalkans/output/revision/boundary_dist_mat.rds and reused on later
# runs. To skip it entirely:
#   WNV_SKIP_BOUNDARY=1 Rscript code/08_revision_analyses.R
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(sf)
  library(eurostat)
  library(survival)
})

options(stringsAsFactors = FALSE)
set.seed(1)

ROOT <- normalizePath(getwd())
RUN  <- file.path(ROOT, Sys.getenv("WNV_RUN_DIR", "run_withBalkans"))
OUT  <- file.path(RUN, "output")
REV  <- file.path(OUT, "revision")
dir.create(REV, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(REV, paste0("run_log_08_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
con <- file(log_path, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

message("08_revision_analyses.R started ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

skip_boundary <- toupper(Sys.getenv("WNV_SKIP_BOUNDARY", "0")) %in% c("1", "TRUE")

## ---- 0. Inputs --------------------------------------------------------------

need_files <- c(
  file.path(OUT, "analysis_cases_withBalkans.rds"),
  file.path(OUT, "spatial_objects_withBalkans.rds"),
  file.path(RUN, "cox_source_full.rds"),
  file.path(ROOT, "WNV_NUTS3cases_2006to2025_withBalkans.csv")
)
missing <- need_files[!file.exists(need_files)]
if (length(missing)) {
  stop("Not found:\n  ", paste(missing, collapse = "\n  "),
       "\nRun run_withBalkans.sh first.", call. = FALSE)
}

ac  <- readRDS(file.path(OUT, "analysis_cases_withBalkans.rds"))
sp  <- readRDS(file.path(OUT, "spatial_objects_withBalkans.rds"))
fit_primary <- readRDS(file.path(RUN, "cox_source_full.rds"))

dist_mat       <- sp$dist_mat
W              <- sp$W
region_lookup  <- sp$region_lookup
region_index   <- sp$region_index
first_outbreak <- sp$first_outbreak

# Kernel settings, kept identical to 00_data_preparation.
lambda_primary        <- as.numeric(Sys.getenv("WNV_LAMBDA", "50"))
distance_threshold_km <- 1000
analysis_year_min     <- min(ac$year)
analysis_year_max     <- max(ac$year)

message("analysis_cases: ", nrow(ac), " region-years, ",
        sum(ac$event), " events, ",
        dplyr::n_distinct(ac$NUTS_ID), " regions, ",
        dplyr::n_distinct(ac$CNTR_CODE), " countries, years ",
        analysis_year_min, " to ", analysis_year_max)

dat_all <- read.csv(file.path(ROOT, "WNV_NUTS3cases_2006to2025_withBalkans.csv")) %>%
  select(any_of(c("NUTS_ID", "CNTR_CODE", "year", "wnf_sum")))

## ---- 0b. Geometry, reproduced exactly as in 00_data_preparation -------------

message("\nFetching NUTS3 2021 geometry from Eurostat ...")
EU_NUTS3 <- get_eurostat_geospatial(resolution = 10, nuts_level = 3, year = 2021) %>%
  filter(!grepl("^FRY", NUTS_ID))

nuts_sf <- EU_NUTS3 %>%
  select(NUTS_ID, CNTR_CODE, URBN_TYPE, COAST_TYPE, MOUNT_TYPE, geometry) %>%
  filter(NUTS_ID %in% region_lookup$NUTS_ID)

# Row order must match dist_mat.
nuts_sf <- nuts_sf[match(region_lookup$NUTS_ID, nuts_sf$NUTS_ID), ]
stopifnot(identical(nuts_sf$NUTS_ID, region_lookup$NUTS_ID))

nuts_proj <- nuts_sf %>% st_transform(3035) %>% st_make_valid()
cent_proj <- suppressWarnings(st_centroid(nuts_proj))
cent_ll   <- st_transform(cent_proj, 4326)
cent_xy   <- st_coordinates(cent_ll)

region_country <- nuts_sf %>% st_drop_geometry() %>% select(NUTS_ID, CNTR_CODE)

write_csv_q <- function(x, name) {
  write.csv(x, file.path(REV, name), row.names = FALSE)
  message("wrote ", name)
}

hr_table <- function(fit, label) {
  s <- summary(fit)
  data.frame(
    model   = label,
    term    = rownames(s$conf.int),
    HR      = round(s$conf.int[, "exp(coef)"], 3),
    lower95 = round(s$conf.int[, "lower .95"], 3),
    upper95 = round(s$conf.int[, "upper .95"], 3),
    p       = signif(s$coefficients[, "Pr(>|z|)"], 3),
    n       = s$n,
    events  = s$nevent,
    concordance = round(unname(s$concordance["C"]), 3),
    AIC     = round(AIC(fit), 1),
    row.names = NULL
  )
}

all_models <- list()

## =============================================================================
## A. Region area and its distribution by country
## =============================================================================

message("\n", strrep("-", 70), "\nA. Region area")

area_df <- data.frame(
  NUTS_ID   = nuts_proj$NUTS_ID,
  area_km2  = as.numeric(units::set_units(st_area(nuts_proj), "km^2"))
) %>%
  mutate(log_area = log(area_km2))

area_df$area_quartile <- cut(
  area_df$area_km2,
  breaks = quantile(area_df$area_km2, probs = seq(0, 1, 0.25), na.rm = TRUE),
  include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4")
)

# The descriptive table describes the analysis cohort, not the full geometry.
area_cohort <- area_df %>% filter(NUTS_ID %in% unique(ac$NUTS_ID))

area_by_country <- area_cohort %>%
  left_join(region_country, by = "NUTS_ID") %>%
  group_by(CNTR_CODE) %>%
  summarise(
    n_regions  = n(),
    median_km2 = round(median(area_km2), 0),
    q25_km2    = round(quantile(area_km2, 0.25), 0),
    q75_km2    = round(quantile(area_km2, 0.75), 0),
    min_km2    = round(min(area_km2), 0),
    max_km2    = round(max(area_km2), 0),
    .groups = "drop"
  ) %>%
  arrange(median_km2)

overall_area <- area_cohort %>%
  summarise(
    CNTR_CODE  = "ALL",
    n_regions  = n(),
    median_km2 = round(median(area_km2), 0),
    q25_km2    = round(quantile(area_km2, 0.25), 0),
    q75_km2    = round(quantile(area_km2, 0.75), 0),
    min_km2    = round(min(area_km2), 0),
    max_km2    = round(max(area_km2), 0)
  )

write_csv_q(bind_rows(area_by_country, overall_area), "rev_area_by_country.csv")
print(as.data.frame(overall_area))

ac <- ac %>% left_join(area_df, by = "NUTS_ID")

## =============================================================================
## B. Mountain and coast terms (Reviewer 1.7)
## =============================================================================

message("\n", strrep("-", 70), "\nB. Mountain and coast")

# Eurostat typologies. MOUNT_TYPE: 1 = more than 50 % of population in mountain
# areas, 2 = more than 50 % of surface, 3 = both, 4 = non-mountain.
# COAST_TYPE: 1 = coastal, 2 = coastal because more than half the population is
# within 50 km of the coast, 3 = non-coastal.
ac <- ac %>%
  mutate(
    mountain = case_when(
      MOUNT_TYPE %in% c(1, 2, 3) ~ "Mountain",
      MOUNT_TYPE %in% c(4)       ~ "Non-mountain",
      TRUE                       ~ NA_character_
    ),
    coastal = case_when(
      COAST_TYPE %in% c(1, 2) ~ "Coastal",
      COAST_TYPE %in% c(3)    ~ "Non-coastal",
      TRUE                    ~ NA_character_
    )
  ) %>%
  mutate(
    mountain = relevel(factor(mountain), ref = "Non-mountain"),
    coastal  = relevel(factor(coastal),  ref = "Non-coastal")
  )

message("Mountain classification of region-years:")
print(table(ac$mountain, useNA = "ifany"))
message("Coast classification of region-years:")
print(table(ac$coastal, useNA = "ifany"))

base_rhs <- "distance100 + pressure_z + log_cases_nearest + tas + hurs + pr_z + factor(URBN_TYPE)"

fit_cox <- function(rhs, data = ac, label) {
  f <- as.formula(paste("Surv(start, stop, event) ~", rhs))
  tryCatch({
    fit <- coxph(f, data = data, ties = "efron")
    message("  fitted: ", label)
    fit
  }, error = function(e) { message("  FAILED (", label, "): ", conditionMessage(e)); NULL })
}

m_base      <- fit_cox(base_rhs, label = "primary (refit here for comparability)")
m_mount     <- fit_cox(paste(base_rhs, "+ mountain"), label = "primary + mountain")
m_mountcst  <- fit_cox(paste(base_rhs, "+ mountain + coastal"), label = "primary + mountain + coast")
m_mountint  <- fit_cox(paste("distance100 * mountain + pressure_z + log_cases_nearest + tas + hurs + pr_z + factor(URBN_TYPE)"),
                       label = "distance x mountain interaction")

geo_tab <- bind_rows(
  if (!is.null(m_base))     hr_table(m_base, "Primary")                    else NULL,
  if (!is.null(m_mount))    hr_table(m_mount, "Primary + mountain")        else NULL,
  if (!is.null(m_mountcst)) hr_table(m_mountcst, "Primary + mountain + coast") else NULL,
  if (!is.null(m_mountint)) hr_table(m_mountint, "Distance x mountain")    else NULL
)
write_csv_q(geo_tab, "rev_mountain_coast_models.csv")
all_models$geography <- geo_tab

if (!is.null(m_base) && !is.null(m_mount)) {
  lrt <- anova(m_base, m_mount)
  message("LRT primary vs primary+mountain: p = ", signif(lrt$`Pr(>|Chi|)`[2], 3),
          "; dAIC = ", round(AIC(m_base) - AIC(m_mount), 2))
}
if (!is.null(m_mount) && !is.null(m_mountint)) {
  message("dAIC mountain main-effect vs interaction model: ",
          round(AIC(m_mount) - AIC(m_mountint), 2))
}

## =============================================================================
## C. Area-adjusted and area-stratified models (editor P2)
## =============================================================================

message("\n", strrep("-", 70), "\nC. Area adjustment")

m_area     <- fit_cox(paste(base_rhs, "+ log_area"), label = "primary + log(area)")
m_areastr  <- fit_cox(paste(base_rhs, "+ strata(area_quartile)"), label = "stratified by area quartile")

# Within-quartile fits, to show the distance effect is not an artefact of size.
by_quartile <- lapply(levels(ac$area_quartile), function(q) {
  sub <- ac %>% filter(area_quartile == q)
  f <- fit_cox(base_rhs, data = sub, label = paste("area quartile", q))
  if (is.null(f)) return(NULL)
  hr_table(f, paste0("Area quartile ", q)) %>% filter(term == "distance100")
})

area_tab <- bind_rows(
  if (!is.null(m_base))    hr_table(m_base, "Primary")                  else NULL,
  if (!is.null(m_area))    hr_table(m_area, "Primary + log(area)")      else NULL,
  if (!is.null(m_areastr)) hr_table(m_areastr, "Stratified by area quartile") else NULL,
  bind_rows(by_quartile)
)
write_csv_q(area_tab, "rev_area_models.csv")
all_models$area <- area_tab

## =============================================================================
## D. Boundary-to-boundary distance (Reviewer 1.1, editor P2)
## =============================================================================

message("\n", strrep("-", 70), "\nD. Boundary-to-boundary distance")

bmat_path <- file.path(REV, "boundary_dist_mat.rds")

# The boundary matrix is built over the full analysis geometry, which does not
# change when the covariate file is repaired, so a cache written by an earlier
# run directory is reused. The dimname check below rejects it if it does not
# match, so a stale cache cannot pass silently.
if (!file.exists(bmat_path)) {
  alt <- file.path(ROOT, "run_withBalkans", "output", "revision", "boundary_dist_mat.rds")
  if (file.exists(alt) && !identical(normalizePath(alt), normalizePath(bmat_path, mustWork = FALSE))) {
    message("Reusing the boundary matrix cached at ", alt)
    bmat_path <- alt
  }
}

if (skip_boundary) {
  message("WNV_SKIP_BOUNDARY set; section D skipped.")
  dist_b <- NULL
} else if (file.exists(bmat_path)) {
  message("Reusing cached ", basename(bmat_path))
  dist_b <- readRDS(bmat_path)
} else {
  message("Computing polygon-to-polygon minimum distances for ",
          nrow(nuts_proj), " regions. This is the slow step.")
  t0 <- Sys.time()
  dist_b <- matrix(NA_real_, nrow(nuts_proj), nrow(nuts_proj),
                   dimnames = list(nuts_proj$NUTS_ID, nuts_proj$NUTS_ID))
  chunk <- 25L
  starts <- seq(1L, nrow(nuts_proj), by = chunk)
  for (k in seq_along(starts)) {
    i0 <- starts[k]
    i1 <- min(i0 + chunk - 1L, nrow(nuts_proj))
    d  <- st_distance(nuts_proj[i0:i1, ], nuts_proj)
    dist_b[i0:i1, ] <- as.matrix(units::drop_units(d)) / 1000
    if (k %% 5 == 0 || k == length(starts)) {
      message("  rows ", i1, "/", nrow(nuts_proj), "  (",
              round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
    }
  }
  dist_b <- (dist_b + t(dist_b)) / 2   # enforce exact symmetry
  diag(dist_b) <- 0
  saveRDS(dist_b, bmat_path)
  message("Done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

if (!is.null(dist_b)) {

  stopifnot(identical(rownames(dist_b), region_lookup$NUTS_ID))

  W_b <- exp(-dist_b / lambda_primary)
  W_b[dist_b > distance_threshold_km] <- 0
  diag(W_b) <- 0

  # Annual loop, identical in structure to 00_data_preparation section 10.
  years <- sort(unique(ac$year))
  pres_list <- vector("list", length(years))
  for (i in seq_along(years)) {
    yy <- years[i]
    invaded_regions <- first_outbreak %>%
      filter(!is.na(first_year), first_year < yy) %>% pull(NUTS_ID)
    invaded_idx <- na.omit(region_index[invaded_regions])
    if (length(invaded_idx) == 0) {
      pres <- rep(0, nrow(W_b)); nd <- rep(NA_real_, nrow(W_b)); ns <- rep(NA_character_, nrow(W_b))
    } else {
      sub <- dist_b[, invaded_idx, drop = FALSE]
      pres <- rowSums(W_b[, invaded_idx, drop = FALSE])
      np <- max.col(-sub, ties.method = "first")
      nd <- sub[cbind(seq_len(nrow(sub)), np)]
      ns <- names(invaded_idx)[np]
    }
    pres_list[[i]] <- data.frame(NUTS_ID = region_lookup$NUTS_ID, year = yy,
                                 invasion_pressure_b = pres,
                                 nearest_invaded_km_b = nd,
                                 nearest_source_b = ns)
  }
  pres_b <- bind_rows(pres_list)

  # Source intensity on the new nearest source, same one-year lag and same
  # untruncated case series as the primary pipeline.
  region_cases_lag <- dat_all %>%
    group_by(NUTS_ID, year) %>%
    summarise(cases_b = max(wnf_sum, na.rm = TRUE), .groups = "drop") %>%
    mutate(cases_b = ifelse(is.finite(cases_b), cases_b, 0), year = year + 1) %>%
    rename(nearest_source_b = NUTS_ID)

  ac_b <- ac %>%
    select(-any_of(c("invasion_pressure_b", "nearest_invaded_km_b", "nearest_source_b"))) %>%
    left_join(pres_b, by = c("NUTS_ID", "year")) %>%
    left_join(region_cases_lag, by = c("nearest_source_b", "year")) %>%
    mutate(
      cases_b            = ifelse(is.na(cases_b), 0, cases_b),
      log_cases_nearest_b = log1p(cases_b),
      distance100_b      = nearest_invaded_km_b / 100,
      pressure_z_b       = as.numeric(scale(invasion_pressure_b))
    )

  message("Boundary distance among first invasions: median ",
          round(median(ac_b$nearest_invaded_km_b[ac_b$event == 1], na.rm = TRUE), 1), " km; ",
          "centroid distance median ",
          round(median(ac_b$nearest_invaded_km[ac_b$event == 1], na.rm = TRUE), 1), " km")

  m_bound <- fit_cox(
    "distance100_b + pressure_z_b + log_cases_nearest_b + tas + hurs + pr_z + factor(URBN_TYPE)",
    data = ac_b, label = "boundary-to-boundary distance")

  m_bound_area <- fit_cox(
    "distance100_b + pressure_z_b + log_cases_nearest_b + tas + hurs + pr_z + factor(URBN_TYPE) + log_area",
    data = ac_b, label = "boundary distance + log(area)")

  bound_tab <- bind_rows(
    if (!is.null(m_base))       hr_table(m_base, "Primary (centroid distance)") else NULL,
    if (!is.null(m_bound))      hr_table(m_bound, "Boundary-to-boundary distance") else NULL,
    if (!is.null(m_bound_area)) hr_table(m_bound_area, "Boundary distance + log(area)") else NULL
  )
  write_csv_q(bound_tab, "rev_boundary_distance_models.csv")
  all_models$boundary <- bound_tab

  write_csv_q(
    ac_b %>% filter(event == 1) %>%
      transmute(NUTS_ID, CNTR_CODE, year,
                centroid_km = round(nearest_invaded_km, 1),
                boundary_km = round(nearest_invaded_km_b, 1),
                source_centroid = nearest_source, source_boundary = nearest_source_b),
    "rev_event_distances_centroid_vs_boundary.csv")

  saveRDS(ac_b, file.path(REV, "analysis_cases_boundary.rds"))
}

## =============================================================================
## E. Direction of invasion (Reviewer 3.3)
## =============================================================================

message("\n", strrep("-", 70), "\nE. Directionality")

# Initial great-circle bearing from source centroid to invaded centroid,
# degrees clockwise from true north.
bearing_deg <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  dl <- (lon2 - lon1) * r
  y <- sin(dl) * cos(lat2 * r)
  x <- cos(lat1 * r) * sin(lat2 * r) - sin(lat1 * r) * cos(lat2 * r) * cos(dl)
  (atan2(y, x) / r + 360) %% 360
}

cent_df <- data.frame(NUTS_ID = nuts_sf$NUTS_ID, lon = cent_xy[, 1], lat = cent_xy[, 2])

events <- ac %>%
  filter(event == 1, !is.na(nearest_source)) %>%
  select(NUTS_ID, CNTR_CODE, year, invasion_phase, nearest_source, nearest_invaded_km) %>%
  left_join(cent_df, by = "NUTS_ID") %>%
  rename(lon_to = lon, lat_to = lat) %>%
  left_join(cent_df %>% rename(nearest_source = NUTS_ID, lon_from = lon, lat_from = lat),
            by = "nearest_source") %>%
  filter(!is.na(lon_from), !is.na(lon_to)) %>%
  mutate(bearing = bearing_deg(lon_from, lat_from, lon_to, lat_to))

message("Events with a computable bearing: ", nrow(events))

rayleigh <- function(theta_deg) {
  th <- theta_deg * pi / 180
  n  <- length(th)
  C  <- sum(cos(th)); S <- sum(sin(th))
  R  <- sqrt(C^2 + S^2)
  Rbar <- R / n
  mu <- (atan2(S, C) * 180 / pi + 360) %% 360
  z  <- R^2 / n
  p  <- exp(sqrt(1 + 4 * n + 4 * (n^2 - R^2)) - (1 + 2 * n))
  data.frame(n = n, mean_bearing_deg = round(mu, 1),
             resultant_length = round(Rbar, 3),
             rayleigh_z = round(z, 3), p_value = signif(p, 3))
}

dir_overall <- rayleigh(events$bearing) %>% mutate(stratum = "All events", .before = 1)

dir_phase <- bind_rows(lapply(levels(events$invasion_phase), function(ph) {
  e <- events %>% filter(invasion_phase == ph)
  if (nrow(e) < 5) return(NULL)
  rayleigh(e$bearing) %>% mutate(stratum = as.character(ph), .before = 1)
}))

# Short-range events only: long jumps are the ones the diffusion model does not
# claim to describe, so the anisotropy test is also reported within 150 km.
e_short <- events %>% filter(nearest_invaded_km <= 150)
dir_short <- if (nrow(e_short) >= 5) {
  rayleigh(e_short$bearing) %>% mutate(stratum = "Events within 150 km", .before = 1)
} else NULL

dir_tab <- bind_rows(dir_overall, dir_phase, dir_short)
write_csv_q(dir_tab, "rev_directionality.csv")
print(as.data.frame(dir_tab))

# 16-sector rose, and a simple north/south/east/west summary.
sector_labels <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                   "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
events <- events %>%
  mutate(
    sector = sector_labels[((floor((bearing + 11.25) / 22.5)) %% 16) + 1],
    quadrant = case_when(
      bearing >= 315 | bearing < 45  ~ "Northward",
      bearing >= 45  & bearing < 135 ~ "Eastward",
      bearing >= 135 & bearing < 225 ~ "Southward",
      TRUE                           ~ "Westward"
    )
  )

rose <- events %>%
  count(sector, name = "n") %>%
  right_join(data.frame(sector = sector_labels), by = "sector") %>%
  mutate(n = ifelse(is.na(n), 0L, n),
         sector = factor(sector, levels = sector_labels),
         pct = round(100 * n / sum(n), 1)) %>%
  arrange(sector)
write_csv_q(rose, "rev_directionality_rose.csv")

quad <- events %>% count(quadrant, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% arrange(desc(n))
write_csv_q(quad, "rev_directionality_quadrants.csv")
print(as.data.frame(quad))

# Binomial test of a north-versus-south asymmetry in the latitude change.
events <- events %>% mutate(dlat = lat_to - lat_from)
bt <- binom.test(sum(events$dlat > 0), sum(events$dlat != 0), p = 0.5)
message("Northward displacement: ", sum(events$dlat > 0), " of ", sum(events$dlat != 0),
        " (", round(100 * bt$estimate, 1), " %), binomial p = ", signif(bt$p.value, 3))
write_csv_q(data.frame(
  n_northward = sum(events$dlat > 0),
  n_directional = sum(events$dlat != 0),
  pct_northward = round(100 * unname(bt$estimate), 1),
  ci_lo = round(100 * bt$conf.int[1], 1),
  ci_hi = round(100 * bt$conf.int[2], 1),
  p_value = signif(bt$p.value, 3),
  median_dlat_deg = round(median(events$dlat), 3)
), "rev_directionality_latitude.csv")

write_csv_q(events %>%
              transmute(NUTS_ID, CNTR_CODE, year, invasion_phase, nearest_source,
                        km = round(nearest_invaded_km, 1),
                        bearing = round(bearing, 1), sector, quadrant,
                        dlat = round(dlat, 3)),
            "rev_event_bearings.csv")

## =============================================================================
## F. Cross-border sources (Reviewer 4.2)
## =============================================================================

message("\n", strrep("-", 70), "\nF. Cross-border sources")

src_country <- region_country %>% rename(nearest_source = NUTS_ID, source_CNTR = CNTR_CODE)

xb <- ac %>%
  filter(event == 1, !is.na(nearest_source)) %>%
  select(NUTS_ID, CNTR_CODE, year, invasion_phase, nearest_source, nearest_invaded_km) %>%
  left_join(src_country, by = "nearest_source") %>%
  mutate(cross_border = source_CNTR != CNTR_CODE)

xb_overall <- data.frame(
  stratum = "All first invasions",
  n_events = nrow(xb),
  n_cross_border = sum(xb$cross_border, na.rm = TRUE),
  pct_cross_border = round(100 * mean(xb$cross_border, na.rm = TRUE), 1)
)

xb_phase <- xb %>%
  group_by(stratum = as.character(invasion_phase)) %>%
  summarise(n_events = n(), n_cross_border = sum(cross_border),
            pct_cross_border = round(100 * mean(cross_border), 1), .groups = "drop")

xb_country <- xb %>%
  group_by(stratum = CNTR_CODE) %>%
  summarise(n_events = n(), n_cross_border = sum(cross_border),
            pct_cross_border = round(100 * mean(cross_border), 1), .groups = "drop") %>%
  arrange(desc(n_events))

write_csv_q(bind_rows(xb_overall, xb_phase), "rev_crossborder.csv")
write_csv_q(xb_country, "rev_crossborder_by_country.csv")
print(as.data.frame(xb_overall))
print(as.data.frame(xb_phase))

# Which national pairs the cross-border links run between.
pairs <- xb %>% filter(cross_border) %>%
  count(source_CNTR, CNTR_CODE, name = "n") %>% arrange(desc(n))
write_csv_q(pairs, "rev_crossborder_pairs.csv")

## =============================================================================
## G. Invasion rate by neighbourhood-pressure decile (Reviewer 1.5)
## =============================================================================

message("\n", strrep("-", 70), "\nG. Pressure deciles")

dec <- ac %>%
  filter(!is.na(invasion_pressure)) %>%
  mutate(pressure_decile = ntile(invasion_pressure, 10)) %>%
  group_by(pressure_decile) %>%
  summarise(
    region_years = n(),
    events = sum(event),
    rate_pct_per_region_year = round(100 * sum(event) / n(), 2),
    pressure_min = round(min(invasion_pressure), 4),
    pressure_max = round(max(invasion_pressure), 4),
    median_distance_km = round(median(nearest_invaded_km, na.rm = TRUE), 1),
    .groups = "drop"
  )
write_csv_q(dec, "rev_pressure_decile_rates.csv")
print(as.data.frame(dec))

lo <- dec %>% filter(pressure_decile <= 4)
hi <- dec %>% filter(pressure_decile == 10)
message("Deciles 1-4: ", sum(lo$events), " events in ", sum(lo$region_years),
        " region-years (", round(100 * sum(lo$events) / sum(lo$region_years), 2), " %)")
message("Decile 10: ", hi$events, " events in ", hi$region_years,
        " region-years (", hi$rate_pct_per_region_year, " %)")

kw <- kruskal.test(invasion_pressure ~ factor(event), data = ac)
message("Kruskal-Wallis, pressure by event status: chi2 = ", round(kw$statistic, 1),
        ", df = ", kw$parameter, ", p = ", signif(kw$p.value, 3))

## ---- close ------------------------------------------------------------------

write_csv_q(bind_rows(all_models), "rev_all_models.csv")
writeLines(capture.output(sessionInfo()), file.path(REV, "sessionInfo_08.txt"))

message("\n", strrep("=", 70))
message("Section outputs are in ", REV)
message("08_revision_analyses.R finished ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

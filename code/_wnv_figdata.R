## -----------------------------------------------------------------------------
## _wnv_figdata.R
##
## Everything the figure scripts share: the fitted objects, the invasion history
## through 2026, the NUTS3 geometry and the map frame. Each figure script sources
## this and then draws one figure, so any figure can be rebuilt on its own.
##
##     source(file.path("code", "_wnv_figdata.R"))
##
## THE 2025 TRAP. WNV_NUTS3cases_2006to2025_withBalkans.csv contains no 2025
## cases despite its name; the 2025 season lives only in
## WNV_2025_ECDC_regions.csv, and 2026 only in the prospective-validation output.
## first_outbreak therefore ends at 2024, and any figure that takes its invasion
## history from first_outbreak alone loses all 38 first invasions of 2025 without
## saying so. The invasion table built here stitches all three sources together
## and prints the count per year so the omission cannot recur unnoticed.
##
## Written 9 September 2026.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(sf)
  library(eurostat)
  library(survival)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(1)

if (!exists("ROOT")) ROOT <- normalizePath(getwd())
RUN <- file.path(ROOT, Sys.getenv("WNV_RUN_DIR", "run_withBalkans"))
OUT <- file.path(RUN, "output")
REV <- file.path(OUT, "revision")

source(file.path(ROOT, "code", "_wnv_style.R"))
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

OOS_LABEL <- "Data for 2026 are to 10 September 2026."

## ---- Fitted objects and observed seasons ------------------------------------

ac  <- readRDS(file.path(OUT, "analysis_cases_withBalkans.rds"))
sp  <- readRDS(file.path(OUT, "spatial_objects_withBalkans.rds"))
fit_primary <- readRDS(file.path(RUN, "cox_source_full.rds"))

region_lookup  <- sp$region_lookup
first_outbreak <- sp$first_outbreak

p25 <- read.csv(file.path(REV, "rev_prospective_2025_regions.csv"))
p26 <- read.csv(file.path(REV, "rev_prospective_2026_regions.csv"))
inv25 <- p25 %>% filter(invaded == 1)
inv26 <- p26 %>% filter(invaded == 1)

## ---- Invasion history, 2006 to 2026 -----------------------------------------

invasion <- bind_rows(
  first_outbreak %>% filter(!is.na(first_year)) %>%
    transmute(NUTS_ID, first_year = as.integer(first_year)),
  inv25 %>% transmute(NUTS_ID, first_year = 2025L),
  inv26 %>% transmute(NUTS_ID, first_year = 2026L)
) %>%
  group_by(NUTS_ID) %>%
  summarise(first_year = min(first_year), .groups = "drop") %>%
  mutate(period = period_of(first_year))

message("fitted events (2007-2024): ", sum(ac$event),
        " | 2025: ", nrow(inv25), " | 2026: ", nrow(inv26))
message("regions with a recorded first invasion by 2026: ", nrow(invasion))

yr_tab <- invasion %>% count(first_year, name = "regions") %>% arrange(first_year)
message("first invasions per year:")
print(as.data.frame(yr_tab), row.names = FALSE)

stopifnot(any(invasion$first_year == 2025), any(invasion$first_year == 2026))

## ---- Geometry ---------------------------------------------------------------

geo_cache <- file.path(REV, "figdata_geometry.rds")
refresh <- toupper(Sys.getenv("WNV_FIGDATA_REFRESH", "0")) %in% c("1", "TRUE")

if (!refresh && file.exists(geo_cache)) {
  message("reusing cached geometry: ", basename(geo_cache))
  nuts_all <- readRDS(geo_cache)
} else {
  message("fetching NUTS3 geometry from Eurostat ...")
  nuts_all <- get_eurostat_geospatial(resolution = 10, nuts_level = 3, year = 2021) %>%
    filter(!grepl("^FRY", NUTS_ID)) %>%
    select(NUTS_ID, CNTR_CODE, geometry)
  saveRDS(nuts_all, geo_cache)
}

# The cohort layer, for figures that report the modelled population.
nuts <- nuts_all %>% filter(NUTS_ID %in% region_lookup$NUTS_ID)

# Whole-of-Europe frame, matching the submitted Figure 1A: Iberia to eastern
# Anatolia, the Mediterranean to northern Scandinavia.
BBOX <- st_bbox(c(xmin = 2400000, xmax = 6950000, ymin = 1320000, ymax = 5460000),
                crs = st_crs(3035))

## ---- Small helpers ----------------------------------------------------------

fig_log <- function(slug) {
  path <- file.path(REV, paste0("run_log_", slug, "_",
                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  # Output is teed to the file; messages stay on stderr, which the shell wrapper
  # captures. A message sink cannot split, so sinking it would hide every error.
  con <- file(path, open = "wt")
  sink(con, split = TRUE)
  message("log: ", path)
  invisible(list(con = con, close = function() { sink(); close(con) }))
}

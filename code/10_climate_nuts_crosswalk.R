###############################################################################
# 10_climate_nuts_crosswalk.R
#
# The covariate file NUTS3ObsClimate_1991to2024.csv is not keyed on the same
# NUTS release as the case file and the analysis geometry. For the Netherlands,
# Portugal and part of Germany it carries the later codes (NL350, NL361, NL363,
# NL415, NL32B, PT1A0, PT1B0, PT1C1-4, PT1D1-3 and others), so those regions
# join to nothing and drop out of the modelled cohort without a message. The
# loss is 15 of 40 Dutch regions, including NL310 which is the Dutch focus,
# 14 of 25 Portuguese regions, and 7 German regions.
#
# This script repairs the join by geography rather than by string matching.
# For every NUTS 2021 analysis region with no climate series, it locates the
# region's centroid inside the polygon of whichever NUTS release the climate
# file is keyed on, and copies that unit's series. Running the match in this
# direction, from the target region outward, handles both splits and merges and
# gives each target region exactly one source.
#
# It also reports, without changing anything, the case-file regions that have
# no unit in NUTS 2021 at all and therefore never receive diffusion metrics.
# That is a separate defect with a separate fix and it is left to the author.
#
# Outputs, all under run_withBalkans/output/revision/:
#   NUTS3ObsClimate_1991to2024_NUTS2021.csv   the repaired covariate file
#   rev_climate_crosswalk.csv                 every mapping, with its method
#   rev_climate_gap_report.csv                what was recovered and what was not
#
# Run from the project folder:
#   Rscript code/10_climate_nuts_crosswalk.R
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(eurostat)
})

options(stringsAsFactors = FALSE)

ROOT <- normalizePath(getwd())
RUN  <- file.path(ROOT, "run_withBalkans")
OUT  <- file.path(RUN, "output")
REV  <- file.path(OUT, "revision")
dir.create(REV, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(REV, paste0("run_log_10_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
con <- file(log_path, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

message("10_climate_nuts_crosswalk.R started ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

clim_in  <- file.path(ROOT, "NUTS3ObsClimate_1991to2024.csv")
case_in  <- file.path(ROOT, "WNV_NUTS3cases_2006to2025_withBalkans.csv")
clim_out <- file.path(REV, "NUTS3ObsClimate_1991to2024_NUTS2021.csv")

for (f in c(clim_in, case_in)) if (!file.exists(f)) stop("Not found: ", f, call. = FALSE)

## ---- 1. Inputs --------------------------------------------------------------

clim <- read.csv(clim_in)
message("climate file: ", nrow(clim), " rows, ",
        dplyr::n_distinct(clim$NUTS_ID), " distinct codes")

cases <- read.csv(case_in) %>% select(any_of(c("NUTS_ID", "CNTR_CODE", "year", "wnf_sum")))

## ---- 2. Geometry, matching 00_data_preparation exactly ----------------------

message("\nFetching NUTS3 geometry ...")

nuts21 <- get_eurostat_geospatial(resolution = 10, nuts_level = 3, year = 2021) %>%
  filter(!grepl("^FRY", NUTS_ID)) %>%
  select(NUTS_ID, CNTR_CODE, geometry)

# Cohort definition, as in the preparation script: countries that ever reported.
wnv_countries <- cases %>%
  mutate(outbreak = ifelse(wnf_sum > 0, 1, 0)) %>%
  group_by(CNTR_CODE) %>%
  summarise(ever = max(outbreak), .groups = "drop") %>%
  filter(ever == 1) %>% pull(CNTR_CODE)

cohort21 <- nuts21 %>% filter(CNTR_CODE %in% wnv_countries)
message("cohort regions in NUTS 2021 geometry: ", nrow(cohort21),
        " across ", dplyr::n_distinct(cohort21$CNTR_CODE), " countries")

have_clim <- unique(clim$NUTS_ID)
missing21 <- cohort21 %>% filter(!(NUTS_ID %in% have_clim))
message("cohort regions with no climate series: ", nrow(missing21))
print(as.data.frame(table(missing21$CNTR_CODE)))

## ---- 3. Candidate source geometries -----------------------------------------
# Codes in the climate file that NUTS 2021 does not contain. Their polygons are
# looked for in the later release first, then the earlier one.

orphan_codes <- setdiff(have_clim, nuts21$NUTS_ID)
message("\nclimate codes absent from NUTS 2021: ", length(orphan_codes))

fetch_year <- function(y) {
  tryCatch(
    get_eurostat_geospatial(resolution = 10, nuts_level = 3, year = y) %>%
      select(NUTS_ID, geometry) %>% mutate(release = as.character(y)),
    error = function(e) { message("  could not fetch NUTS ", y, ": ", conditionMessage(e)); NULL })
}

src_list <- list()
for (y in c(2024, 2016, 2013)) {
  g <- fetch_year(y)
  if (is.null(g)) next
  hit <- g %>% filter(NUTS_ID %in% orphan_codes)
  message("  NUTS ", y, " supplies ", nrow(hit), " of the orphan codes")
  src_list[[as.character(y)]] <- hit
  orphan_codes <- setdiff(orphan_codes, hit$NUTS_ID)
  if (!length(orphan_codes)) break
}
if (length(orphan_codes))
  message("  still unlocated: ", paste(utils::head(orphan_codes, 20), collapse = ", "))

src <- if (length(src_list)) do.call(rbind, src_list) else NULL

## ---- 4. Match each missing 2021 region to a source unit ---------------------

crosswalk <- data.frame(target_NUTS_ID = character(), source_NUTS_ID = character(),
                        source_release = character(), method = character())

if (!is.null(src) && nrow(src) && nrow(missing21)) {

  tgt_cent <- suppressWarnings(st_centroid(st_transform(missing21, 3035) %>% st_make_valid()))
  src_proj <- st_transform(src, 3035) %>% st_make_valid()

  inside <- suppressMessages(st_within(tgt_cent, src_proj))
  pick <- vapply(seq_len(nrow(tgt_cent)), function(i) {
    h <- inside[[i]]
    if (length(h) >= 1) h[1] else NA_integer_
  }, integer(1))

  meth <- ifelse(is.na(pick), "no containing source unit", "centroid within")

  crosswalk <- data.frame(
    target_NUTS_ID = missing21$NUTS_ID,
    target_country = missing21$CNTR_CODE,
    source_NUTS_ID = ifelse(is.na(pick), NA_character_, src_proj$NUTS_ID[pick]),
    source_release = ifelse(is.na(pick), NA_character_, src_proj$release[pick]),
    method         = meth
  )
}

write.csv(crosswalk, file.path(REV, "rev_climate_crosswalk.csv"), row.names = FALSE)
message("\nClimate crosswalk:")
print(as.data.frame(crosswalk))

## ---- 5. Write the repaired covariate file -----------------------------------
# Rows are copied, not averaged: each target region takes the full series of the
# single source unit that contains its centroid. Where two target regions fall
# inside one source unit they receive the same series, which is the honest
# consequence of the source release being coarser there.

cw_ok <- crosswalk %>% filter(!is.na(source_NUTS_ID))

added <- clim %>%
  inner_join(cw_ok %>% select(source_NUTS_ID, target_NUTS_ID),
             by = c("NUTS_ID" = "source_NUTS_ID"), relationship = "many-to-many") %>%
  mutate(NUTS_ID = target_NUTS_ID) %>%
  select(-target_NUTS_ID)

clim_fixed <- bind_rows(clim, added) %>% distinct(NUTS_ID, scenario, model, year, .keep_all = TRUE)

write.csv(clim_fixed, clim_out, row.names = FALSE)
message("\nwrote ", clim_out)
message("rows: ", nrow(clim), " -> ", nrow(clim_fixed),
        "; distinct codes: ", dplyr::n_distinct(clim$NUTS_ID), " -> ",
        dplyr::n_distinct(clim_fixed$NUTS_ID))

still_missing <- cohort21 %>%
  st_drop_geometry() %>%
  filter(!(NUTS_ID %in% unique(clim_fixed$NUTS_ID)))

message("cohort regions still without climate: ", nrow(still_missing))
if (nrow(still_missing)) print(as.data.frame(table(still_missing$CNTR_CODE)))

## ---- 6. The other defect: case-file codes with no NUTS 2021 unit ------------
# These never enter region_lookup, so they never receive a distance or a
# pressure value and drop out of the cohort for a reason unrelated to climate.
# Reported, not repaired: reassigning their cases changes the outcome data and
# is the author's call.

case_codes <- cases %>% filter(CNTR_CODE %in% wnv_countries) %>% distinct(NUTS_ID, CNTR_CODE)
orphan_cases <- case_codes %>% filter(!(NUTS_ID %in% nuts21$NUTS_ID))

if (nrow(orphan_cases)) {
  burden <- cases %>%
    filter(NUTS_ID %in% orphan_cases$NUTS_ID) %>%
    group_by(NUTS_ID) %>%
    summarise(total_cases = sum(wnf_sum, na.rm = TRUE),
              first_case_year = suppressWarnings(min(year[wnf_sum > 0])),
              .groups = "drop") %>%
    mutate(first_case_year = ifelse(is.finite(first_case_year), first_case_year, NA))
  orphan_cases <- orphan_cases %>% left_join(burden, by = "NUTS_ID")
}

message("\ncase-file regions with no NUTS 2021 unit: ", nrow(orphan_cases))
if (nrow(orphan_cases)) print(as.data.frame(orphan_cases))

report <- bind_rows(
  data.frame(item = "cohort regions in NUTS 2021 geometry", n = nrow(cohort21)),
  data.frame(item = "cohort regions with no climate before repair", n = nrow(missing21)),
  data.frame(item = "cohort regions recovered by the crosswalk", n = nrow(cw_ok)),
  data.frame(item = "cohort regions still without climate", n = nrow(still_missing)),
  data.frame(item = "case-file regions with no NUTS 2021 unit", n = nrow(orphan_cases))
)
write.csv(report, file.path(REV, "rev_climate_gap_report.csv"), row.names = FALSE)
if (nrow(orphan_cases))
  write.csv(orphan_cases, file.path(REV, "rev_cases_without_nuts2021_unit.csv"), row.names = FALSE)

print(report)
writeLines(capture.output(sessionInfo()), file.path(REV, "sessionInfo_10.txt"))
message("\n10_climate_nuts_crosswalk.R finished ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

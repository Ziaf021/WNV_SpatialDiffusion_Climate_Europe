###############################################################################
# 09_prospective_validation_v2.R
#
# Prospective validation, with the two fixes agreed for the revision:
#
#   1. Codes reported by ECDC under a later NUTS release that have no match in
#      the NUTS 2021 analysis geometry are resolved by a geometry-based
#      crosswalk (the 2024 region's centroid is located inside the 2021 region
#      that contains it) instead of being dropped silently. This recovers the
#      five Dutch 2026 areas NL350, NL363, NL361, NL32B and NL415. Kosovo
#      (XK003, XK007) has no unit in NUTS 2021 at all and stays excluded, for
#      the reason already documented.
#
#   2. A cluster-aware sensitivity analysis. First invasions in the same season
#      whose centroids lie within 50 km of each other are treated as one
#      introduction, so the metrics do not reward or penalise a country simply
#      for having small NUTS3 units. The member closest to its nearest invaded
#      source is kept as the representative of the cluster; the other members
#      are removed from the evaluation set rather than counted as negatives.
#
# Nothing is refitted. The model is the one fitted through 2024.
# 07_prospective_validation.R is left as it is; this writes to
# run_withBalkans/output/revision/ under new names.
#
# Run from the project folder:
#   Rscript code/09_prospective_validation_v2.R
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(eurostat)
  library(survival)
})

options(stringsAsFactors = FALSE)

ROOT <- normalizePath(getwd())
RUN  <- file.path(ROOT, Sys.getenv("WNV_RUN_DIR", "run_withBalkans"))
OUT  <- file.path(RUN, "output")
REV  <- file.path(OUT, "revision")
dir.create(REV, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(REV, paste0("run_log_09_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
con <- file(log_path, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)

message("09_prospective_validation_v2.R started ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

CLUSTER_KM <- as.numeric(Sys.getenv("WNV_CLUSTER_KM", "50"))

sp  <- readRDS(file.path(OUT, "spatial_objects_withBalkans.rds"))
ac  <- readRDS(file.path(OUT, "analysis_cases_withBalkans.rds"))
fit <- readRDS(file.path(RUN, "cox_source_full.rds"))

dist_mat       <- sp$dist_mat
W              <- sp$W
region_lookup  <- sp$region_lookup
region_index   <- sp$region_index
first_outbreak <- sp$first_outbreak

## ---- 1. Observed seasons ----------------------------------------------------

cases_2025 <- read.csv(file.path(ROOT, "WNV_2025_ECDC_regions.csv"))

raw_2026 <- read.csv(file.path(ROOT, "WNV_2026_ECDC.csv"), check.names = FALSE)
code_col <- grep("NUTS|GAUL", names(raw_2026), value = TRUE)[1]
case_col <- grep("^Total human cases", names(raw_2026), value = TRUE)[1]
cases_2026 <- data.frame(NUTS_ID = trimws(raw_2026[[code_col]]),
                         wnf_sum = suppressWarnings(as.numeric(raw_2026[[case_col]])))

merged <- read.csv(file.path(ROOT, "WNV_NUTS3cases_2006to2025_withBalkans.csv"))[, c("NUTS_ID", "year", "wnf_sum")]

## ---- 2. Geometry and the NUTS 2024 to 2021 crosswalk ------------------------

message("\nFetching NUTS3 geometry (2021 and 2024) from Eurostat ...")

nuts21 <- get_eurostat_geospatial(resolution = 10, nuts_level = 3, year = 2021) %>%
  filter(!grepl("^FRY", NUTS_ID)) %>%
  select(NUTS_ID, CNTR_CODE, geometry)

nuts21_analysis <- nuts21 %>% filter(NUTS_ID %in% region_lookup$NUTS_ID)
nuts21_analysis <- nuts21_analysis[match(region_lookup$NUTS_ID, nuts21_analysis$NUTS_ID), ]
stopifnot(identical(nuts21_analysis$NUTS_ID, region_lookup$NUTS_ID))

cent21 <- suppressWarnings(st_centroid(st_transform(nuts21_analysis, 3035)))
cent21_xy <- st_coordinates(cent21)

nuts24 <- tryCatch(
  get_eurostat_geospatial(resolution = 10, nuts_level = 3, year = 2024) %>%
    select(NUTS_ID, CNTR_CODE, geometry),
  error = function(e) { message("NOTE: could not fetch NUTS 2024 geometry: ", conditionMessage(e)); NULL }
)

# Codes reported by ECDC that the analysis geometry does not contain.
unmatched <- union(
  setdiff(cases_2025$NUTS_ID[cases_2025$wnf_sum > 0], region_lookup$NUTS_ID),
  setdiff(cases_2026$NUTS_ID[!is.na(cases_2026$wnf_sum) & cases_2026$wnf_sum > 0],
          region_lookup$NUTS_ID)
)
unmatched <- setdiff(unmatched, NA_character_)
message("\nReported codes with no NUTS 2021 analysis unit: ",
        if (length(unmatched)) paste(unmatched, collapse = ", ") else "none")

crosswalk <- data.frame(reported_code = character(), mapped_NUTS_ID = character(),
                        method = character(), note = character())

if (length(unmatched) && !is.null(nuts24)) {
  cand <- nuts24 %>% filter(NUTS_ID %in% unmatched)
  if (nrow(cand)) {
    cand_proj <- st_transform(cand, 3035) %>% st_make_valid()
    cand_cent <- suppressWarnings(st_centroid(cand_proj))
    nuts21_proj <- st_transform(nuts21, 3035) %>% st_make_valid()

    # Containment only. A nearest-feature fallback would place a code from a
    # country with no NUTS 2021 unit at all (Kosovo) inside a neighbouring
    # country's region and manufacture an invasion event there, so codes whose
    # centroid falls in no NUTS 2021 region stay unresolved and excluded.
    inside <- suppressMessages(st_within(cand_cent, nuts21_proj))
    mapped <- vapply(seq_len(nrow(cand_cent)), function(i) {
      hit <- inside[[i]]
      if (length(hit) >= 1) nuts21_proj$NUTS_ID[hit[1]] else NA_character_
    }, character(1))

    crosswalk <- data.frame(
      reported_code  = cand$NUTS_ID,
      mapped_NUTS_ID = mapped,
      method         = ifelse(is.na(mapped), "unresolved", "centroid within"),
      note           = ifelse(is.na(mapped), "centroid in no NUTS 2021 region; excluded",
                       ifelse(mapped %in% region_lookup$NUTS_ID,
                              "in analysis cohort", "outside analysis cohort"))
    )
  }
  still <- setdiff(unmatched, crosswalk$reported_code)
  if (length(still)) {
    crosswalk <- bind_rows(crosswalk, data.frame(
      reported_code = still, mapped_NUTS_ID = NA_character_,
      method = "unresolved",
      note = "no unit in NUTS 2021 or NUTS 2024; excluded"))
  }
}

write.csv(crosswalk, file.path(REV, "rev_nuts_crosswalk.csv"), row.names = FALSE)
message("Crosswalk:")
print(as.data.frame(crosswalk))

apply_crosswalk <- function(df) {
  if (!nrow(crosswalk)) return(df)
  cw <- crosswalk %>%
    filter(method == "centroid within", !is.na(mapped_NUTS_ID),
           mapped_NUTS_ID %in% region_lookup$NUTS_ID)
  df %>%
    left_join(cw %>% select(reported_code, mapped_NUTS_ID),
              by = c("NUTS_ID" = "reported_code")) %>%
    mutate(NUTS_ID = ifelse(is.na(mapped_NUTS_ID), NUTS_ID, mapped_NUTS_ID)) %>%
    select(-mapped_NUTS_ID) %>%
    group_by(NUTS_ID) %>%
    summarise(wnf_sum = sum(wnf_sum, na.rm = TRUE), .groups = "drop")
}

cases_2025_x <- apply_crosswalk(cases_2025 %>% select(NUTS_ID, wnf_sum))
cases_2026_x <- apply_crosswalk(cases_2026 %>% select(NUTS_ID, wnf_sum))

## ---- 3. Scoring machinery (unchanged from 07) -------------------------------

pressure_mean <- mean(ac$invasion_pressure, na.rm = TRUE)
pressure_sd   <- sd(ac$invasion_pressure,   na.rm = TRUE)

# Precipitation is standardised against the fitting data too, for the same reason
# pressure is. Previously the model carried scale(pr) in its formula and
# predict() recomputed it on the scoring set, so out-of-sample precipitation was
# centred and scaled by the wrong quantities.
pr_mean <- mean(ac$pr, na.rm = TRUE)
pr_sd   <- sd(ac$pr,   na.rm = TRUE)

region_static <- ac %>%
  filter(year >= 2020) %>%
  group_by(NUTS_ID) %>%
  summarise(tas = mean(tas, na.rm = TRUE), hurs = mean(hurs, na.rm = TRUE),
            pr = mean(pr, na.rm = TRUE),
            URBN_TYPE = URBN_TYPE[which(!is.na(URBN_TYPE))[1]], .groups = "drop")

auc <- function(score, label) {
  pos <- score[label == 1]; neg <- score[label == 0]
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

build_season <- function(invaded_ids, observed_ids, prior_cases) {
  naive <- setdiff(region_lookup$NUTS_ID, invaded_ids)
  idx_naive   <- region_index[naive]
  idx_invaded <- na.omit(region_index[invaded_ids])
  sub  <- dist_mat[idx_naive, idx_invaded, drop = FALSE]
  near <- max.col(-sub, ties.method = "first")
  d    <- sub[cbind(seq_len(nrow(sub)), near)]
  src  <- names(idx_invaded)[near]
  pres <- rowSums(W[idx_naive, idx_invaded, drop = FALSE])

  df <- data.frame(NUTS_ID = naive, nearest_invaded_km = d, nearest_source = src,
                   invasion_pressure = pres) %>%
    mutate(distance100 = nearest_invaded_km / 100,
           pressure_z  = (invasion_pressure - pressure_mean) / pressure_sd) %>%
    left_join(prior_cases, by = c("nearest_source" = "NUTS_ID")) %>%
    mutate(cases = ifelse(is.na(cases), 0, cases),
           log_cases_nearest = log1p(cases)) %>%
    left_join(region_static, by = "NUTS_ID") %>%
    mutate(pr_z    = (pr - pr_mean) / pr_sd,
           invaded = as.integer(NUTS_ID %in% observed_ids))

  ok <- !is.na(df$tas) & !is.na(df$hurs) & !is.na(df$pr) & !is.na(df$URBN_TYPE)
  df$model_score <- NA_real_
  df$model_score[ok] <- predict(fit, newdata = df[ok, ], type = "lp")
  df
}

summarise_season <- function(df, season_label, set_label) {
  df <- df %>% mutate(pressure_decile = ntile(invasion_pressure, 10))
  dd <- df$nearest_invaded_km[df$invaded == 1]

  # The distance and pressure columns use every observed invasion, but the
  # full-model AUC can only use regions with complete covariates. If those two
  # sets differ, the row is not internally consistent and must say so rather
  # than let a reader assume one denominator throughout.
  unscored <- df$NUTS_ID[df$invaded == 1 & !is.finite(df$model_score)]
  if (length(unscored))
    message("  ", season_label, " (", set_label, "): ", length(unscored),
            " observed invasions have no model score and are excluded from ",
            "auc_full_model only: ", paste(unscored, collapse = ", "))

  data.frame(
    season = season_label,
    evaluation_set = set_label,
    regions_at_risk = nrow(df),
    new_invasions = sum(df$invaded),
    invasions_scored = sum(df$invaded == 1 & is.finite(df$model_score)),
    median_distance_km = round(median(dd, na.rm = TRUE), 1),
    pct_within_75km  = round(100 * mean(dd <= 75,  na.rm = TRUE), 1),
    pct_within_150km = round(100 * mean(dd <= 150, na.rm = TRUE), 1),
    pct_within_300km = round(100 * mean(dd <= 300, na.rm = TRUE), 1),
    pct_in_top_pressure_decile = round(100 * mean(df$pressure_decile[df$invaded == 1] == 10, na.rm = TRUE), 1),
    pct_in_top_two_deciles     = round(100 * mean(df$pressure_decile[df$invaded == 1] >= 9,  na.rm = TRUE), 1),
    auc_pressure   = round(auc(df$invasion_pressure, df$invaded), 3),
    auc_proximity  = round(auc(-df$nearest_invaded_km, df$invaded), 3),
    auc_full_model = round(auc(df$model_score, df$invaded), 3),
    naive_median_km = round(median(df$nearest_invaded_km, na.rm = TRUE), 1)
  )
}

## ---- 4. Cluster collapse ----------------------------------------------------
# Single-linkage grouping of the season's first invasions at CLUSTER_KM.
# Distances come from the same centroid matrix used everywhere else.

collapse_clusters <- function(df, threshold_km = CLUSTER_KM) {
  ev <- df %>% filter(invaded == 1)
  if (nrow(ev) <= 1) return(list(keep = ev$NUTS_ID, drop = character(0), table = NULL))
  ids <- ev$NUTS_ID
  d   <- dist_mat[ids, ids, drop = FALSE]
  hc  <- hclust(as.dist(d), method = "single")
  cl  <- cutree(hc, h = threshold_km)
  ev$cluster <- cl
  reps <- ev %>% group_by(cluster) %>%
    slice_min(nearest_invaded_km, n = 1, with_ties = FALSE) %>%
    ungroup()
  list(
    keep  = reps$NUTS_ID,
    drop  = setdiff(ids, reps$NUTS_ID),
    table = ev %>% arrange(cluster, nearest_invaded_km) %>%
      transmute(cluster, NUTS_ID, km = round(nearest_invaded_km, 1),
                nearest_source, representative = NUTS_ID %in% reps$NUTS_ID)
  )
}

## ---- 5. Run both seasons ----------------------------------------------------

invaded_2024 <- first_outbreak %>% filter(!is.na(first_year), first_year <= 2024) %>% pull(NUTS_ID)
prior_2024 <- merged %>% filter(year == 2024) %>%
  group_by(NUTS_ID) %>% summarise(cases = max(wnf_sum, na.rm = TRUE), .groups = "drop")

obs_2025_all <- cases_2025_x %>% filter(wnf_sum > 0) %>% pull(NUTS_ID)
obs_2025 <- setdiff(intersect(obs_2025_all, region_lookup$NUTS_ID), invaded_2024)
dropped_2025 <- setdiff(obs_2025_all, region_lookup$NUTS_ID)

d25 <- build_season(invaded_2024, obs_2025, prior_2024)

invaded_2025 <- union(invaded_2024, obs_2025)
prior_2025 <- cases_2025_x %>% filter(wnf_sum > 0) %>%
  group_by(NUTS_ID) %>% summarise(cases = max(wnf_sum), .groups = "drop")

obs_2026_all <- cases_2026_x %>% filter(!is.na(wnf_sum), wnf_sum > 0) %>% pull(NUTS_ID)
obs_2026 <- setdiff(intersect(obs_2026_all, region_lookup$NUTS_ID), invaded_2025)
dropped_2026 <- setdiff(obs_2026_all, region_lookup$NUTS_ID)

d26 <- build_season(invaded_2025, obs_2026, prior_2025)

message("\n2025: ", length(obs_2025), " new invasions; 2026 to date: ", length(obs_2026))
if (length(dropped_2025)) message("2025 codes still unresolved: ", paste(dropped_2025, collapse = ", "))
if (length(dropped_2026)) message("2026 codes still unresolved: ", paste(dropped_2026, collapse = ", "))

cl25 <- collapse_clusters(d25)
cl26 <- collapse_clusters(d26)

message("\nCluster collapse at ", CLUSTER_KM, " km:")
message("  2025: ", length(obs_2025), " invaded regions -> ", length(cl25$keep), " introductions")
message("  2026: ", length(obs_2026), " invaded regions -> ", length(cl26$keep), " introductions")

if (!is.null(cl25$table)) write.csv(cl25$table, file.path(REV, "rev_clusters_2025.csv"), row.names = FALSE)
if (!is.null(cl26$table)) write.csv(cl26$table, file.path(REV, "rev_clusters_2026.csv"), row.names = FALSE)

d25_cl <- d25 %>% filter(!(NUTS_ID %in% cl25$drop))
d26_cl <- d26 %>% filter(!(NUTS_ID %in% cl26$drop))

res <- bind_rows(
  summarise_season(d25,    2025, "all regions"),
  summarise_season(d25_cl, 2025, paste0("cluster-collapsed at ", CLUSTER_KM, " km")),
  summarise_season(d26,    2026, "all regions"),
  summarise_season(d26_cl, 2026, paste0("cluster-collapsed at ", CLUSTER_KM, " km"))
)

cat("\n============= PROSPECTIVE VALIDATION, REVISED =============\n")
cat("Model fitted through 2024. Nothing refitted.\n")
cat("NUTS crosswalk applied to reported codes from later NUTS releases.\n\n")
print(as.data.frame(res), row.names = FALSE)

write.csv(res, file.path(REV, "rev_prospective_validation_summary.csv"), row.names = FALSE)
write.csv(d25 %>% mutate(pressure_decile = ntile(invasion_pressure, 10)),
          file.path(REV, "rev_prospective_2025_regions.csv"), row.names = FALSE)
write.csv(d26 %>% mutate(pressure_decile = ntile(invasion_pressure, 10)),
          file.path(REV, "rev_prospective_2026_regions.csv"), row.names = FALSE)

show_events <- function(df, label) {
  cat("\nNew invasions, ", label, ":\n", sep = "")
  tab <- df %>% filter(invaded == 1) %>%
    mutate(pressure_decile = ntile(invasion_pressure, 10)) %>%
    transmute(NUTS_ID, km = round(nearest_invaded_km, 1), nearest_source, pressure_decile) %>%
    arrange(km) %>% as.data.frame()
  print(tab, row.names = FALSE)
}
show_events(d25, "2025")
show_events(d26, "2026 to date")

writeLines(capture.output(sessionInfo()), file.path(REV, "sessionInfo_09.txt"))
message("\nWrote outputs to ", REV)
message("09_prospective_validation_v2.R finished ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

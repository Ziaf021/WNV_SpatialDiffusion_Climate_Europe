###############################################################################
# 12b_direction_latitude_4periods.R
#
# Reviewer 3.2, reconciled to one basis: every first invasion to 2026 and four
# display periods, 2006-2010, 2011-2015, 2016-2020 and 2021-2026. This matches
# Supplementary Figure S10 (figS10_directionality_rose.R).
#
# Nothing here is modelled, so the 2025 and 2026 seasons are included.
#
#   1. Latitude change relative to the event's own source: binomial test of the
#      share of first invasions lying north of their nearest source. Needs a
#      source, so 362 events (the three regions invaded in 2006 have none).
#      Centroids and source as in 08 section E.
#   2. Absolute latitude of newly invaded regions: median by period, Spearman
#      correlation with year and linear trend. Needs no source, so all 365
#      regions. Centroid latitude as in 12_regional_composition.R.
#   3. Composition between the European regional groupings by period: chi-
#      squared and logistic trend in the West share, as in 12, on 365 regions.
#
# Writes to <run>/output/revision:
#   rev_latitude_change_2007_2026.csv
#   rev_latitude_by_period_4p.csv
#   rev_regional_composition_4p.csv
#   rev_regional_tests_4p.csv
#
# Run from the project folder, one line:
#   WNV_RUN_DIR=run_withBalkans_v2 Rscript code/12b_direction_latitude_4periods.R
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("direction_latitude_4p")
on.exit(lg$close(), add = TRUE)

PERIODS <- c("2006-2010", "2011-2015", "2016-2020", "2021-2026")
period4 <- function(y) factor(as.character(cut(y, c(2005, 2010, 2015, 2020, 2026),
                                               labels = PERIODS)), levels = PERIODS)

## ---- 1. Latitude change relative to the source (362 events) -------------------
cent_ll <- suppressWarnings(nuts %>% st_transform(3035) %>% st_make_valid() %>%
                              st_centroid() %>% st_transform(4326))
cent_df <- data.frame(NUTS_ID = cent_ll$NUTS_ID, st_coordinates(cent_ll)) %>%
  rename(lon = X, lat = Y)

src_ev <- bind_rows(
  ac %>% filter(event == 1) %>% transmute(NUTS_ID, year, nearest_source),
  inv25 %>% transmute(NUTS_ID, year = 2025L, nearest_source),
  inv26 %>% transmute(NUTS_ID, year = 2026L, nearest_source)
) %>%
  left_join(cent_df %>% select(NUTS_ID, lat_to = lat), by = "NUTS_ID") %>%
  left_join(cent_df %>% select(nearest_source = NUTS_ID, lat_from = lat), by = "nearest_source") %>%
  mutate(dlat = lat_to - lat_from, period = period4(year))

stopifnot(!anyNA(src_ev$dlat), !anyNA(src_ev$period))
message("events with a source: ", nrow(src_ev))

# The fitted-period part must reproduce 08.
ref <- read.csv(file.path(REV, "rev_event_bearings.csv"), stringsAsFactors = FALSE)
chk <- inner_join(src_ev %>% filter(year <= 2024) %>% select(NUTS_ID, dlat),
                  ref %>% select(NUTS_ID, dlat_ref = dlat), by = "NUTS_ID")
stopifnot(nrow(chk) == nrow(ref), max(abs(chk$dlat - chk$dlat_ref)) < 0.001)
message("fitted-period latitude changes match 08 for all ", nrow(chk), " events")

binom_row <- function(d, label) {
  k <- sum(d$dlat > 0); n <- sum(d$dlat != 0)
  bt <- binom.test(k, n, p = 0.5)
  data.frame(stratum = label, n_events = nrow(d), n_northward = k, n_directional = n,
             pct_northward = round(100 * k / n, 1),
             ci_lo = round(100 * bt$conf.int[1], 1), ci_hi = round(100 * bt$conf.int[2], 1),
             p_value = signif(bt$p.value, 3),
             median_dlat_deg = round(median(d$dlat), 3))
}
lat_change <- bind_rows(
  binom_row(src_ev, "All first invasions 2007-2026"),
  bind_rows(lapply(PERIODS, function(p) binom_row(src_ev[src_ev$period == p, ], p))))
message("\nLatitude change relative to the nearest source:")
print(lat_change, row.names = FALSE)
write.csv(lat_change, file.path(REV, "rev_latitude_change_2007_2026.csv"), row.names = FALSE)

## ---- 2. Absolute latitude of newly invaded regions (365 regions) --------------
cent_abs <- nuts_all %>% st_transform(4326) %>% { suppressWarnings(st_centroid(.)) }
lat_abs <- data.frame(NUTS_ID = cent_abs$NUTS_ID, lat = st_coordinates(cent_abs)[, 2])

ev <- invasion %>%
  mutate(period = period4(first_year)) %>%
  left_join(lat_abs, by = "NUTS_ID")
stopifnot(!anyNA(ev$lat), !anyNA(ev$period))
message("\ninvaded regions with a latitude: ", nrow(ev))

lat_tab <- ev %>% group_by(period) %>%
  summarise(n = n(), median_lat = round(median(lat), 2),
            q25 = round(quantile(lat, 0.25), 2), q75 = round(quantile(lat, 0.75), 2),
            .groups = "drop")
print(as.data.frame(lat_tab), row.names = FALSE)
write.csv(lat_tab, file.path(REV, "rev_latitude_by_period_4p.csv"), row.names = FALSE)

sp <- suppressWarnings(cor.test(ev$first_year, ev$lat, method = "spearman"))
sl <- summary(lm(lat ~ first_year, data = ev))$coefficients["first_year", ]
message("Spearman rho = ", round(unname(sp$estimate), 3), ", p = ", signif(sp$p.value, 3))
message("Linear trend ", round(sl["Estimate"], 4), " degrees per year (SE ",
        round(sl["Std. Error"], 4), "), p = ", signif(sl["Pr(>|t|)"], 3))

## ---- 3. Composition between regional groupings (365 regions) ------------------
cls <- read.csv(file.path(ROOT, "Europe_region_classification.csv")) %>%
  transmute(NUTS2_ID = trimws(NUTS2_ID), CNTR_CODE = trimws(CNTR_CODE),
            europe_region = trimws(europe_region)) %>%
  filter(nzchar(europe_region)) %>%
  distinct(NUTS2_ID, .keep_all = TRUE)
modal <- cls %>% count(CNTR_CODE, europe_region, name = "k") %>%
  group_by(CNTR_CODE) %>% slice_max(k, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(CNTR_CODE, fallback = europe_region)

ev <- ev %>%
  mutate(NUTS2_ID = substr(NUTS_ID, 1, 4), CNTR_CODE = substr(NUTS_ID, 1, 2)) %>%
  left_join(cls %>% select(NUTS2_ID, europe_region), by = "NUTS2_ID") %>%
  left_join(modal, by = "CNTR_CODE") %>%
  mutate(europe_region = ifelse(is.na(europe_region), fallback, europe_region))
stopifnot(!anyNA(ev$europe_region))

comp <- ev %>% count(period, europe_region, name = "n") %>%
  group_by(period) %>% mutate(total = sum(n), pct = round(100 * n / total, 1)) %>%
  ungroup() %>% arrange(period, europe_region)
print(as.data.frame(comp), row.names = FALSE)
write.csv(comp, file.path(REV, "rev_regional_composition_4p.csv"), row.names = FALSE)

tab <- table(ev$period, ev$europe_region)
tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
cs <- suppressWarnings(chisq.test(tab))
mc_p <- if (any(cs$expected < 5)) chisq.test(tab, simulate.p.value = TRUE, B = 20000)$p.value else NA_real_
message("Chi-squared = ", round(unname(cs$statistic), 2), ", df = ", unname(cs$parameter),
        ", p = ", signif(cs$p.value, 3), if (!is.na(mc_p)) paste0("; Monte Carlo p = ", signif(mc_p, 3)))

ev <- ev %>% mutate(is_west = as.integer(europe_region == "West"),
                    period_idx = as.integer(period))
share <- ev %>% group_by(period) %>%
  summarise(n = n(), n_west = sum(is_west), pct_west = round(100 * mean(is_west), 1),
            .groups = "drop")
print(as.data.frame(share), row.names = FALSE)
ct <- summary(glm(is_west ~ period_idx, family = binomial, data = ev))$coefficients
or <- exp(ct["period_idx", "Estimate"])
ci <- exp(ct["period_idx", "Estimate"] + c(-1.96, 1.96) * ct["period_idx", "Std. Error"])
message("West share trend: odds ratio per period ", round(or, 2), " (", round(ci[1], 2),
        " to ", round(ci[2], 2), "), p = ", signif(ct["period_idx", "Pr(>|z|)"], 3))

write.csv(data.frame(
  test = c("binomial, share north of source (381 events)",
           "chi-squared, period by grouping", "chi-squared, Monte Carlo",
           "logistic trend, share in the West grouping",
           "Spearman, year against latitude", "linear trend, latitude on year"),
  estimate = c(lat_change$pct_northward[1], round(unname(cs$statistic), 2), NA, round(or, 3),
               round(unname(sp$estimate), 3), round(unname(sl["Estimate"]), 4)),
  interval = c(paste0(lat_change$ci_lo[1], " to ", lat_change$ci_hi[1]),
               paste("df", unname(cs$parameter)), NA,
               paste0(round(ci[1], 3), " to ", round(ci[2], 3)), NA,
               paste0("SE ", round(unname(sl["Std. Error"]), 4))),
  p_value = c(lat_change$p_value[1], signif(cs$p.value, 3), signif(mc_p, 3),
              signif(ct["period_idx", "Pr(>|z|)"], 3), signif(sp$p.value, 3),
              signif(unname(sl["Pr(>|t|)"]), 3))
), file.path(REV, "rev_regional_tests_4p.csv"), row.names = FALSE)
write.csv(share, file.path(REV, "rev_regional_trend_4p.csv"), row.names = FALSE)

message("\n12b done")

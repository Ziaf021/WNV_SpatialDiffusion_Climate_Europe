###############################################################################
# 12_regional_composition.R
#
# Where in Europe the newly invaded regions were, period by period, and whether
# they moved north.
#
# Three questions, deliberately kept apart, because the manuscript has been
# conflating them:
#
#   1. Does invasion have a preferred compass bearing relative to its OWN
#      nearest source? Answered in section E of 08_revision_analyses.R: no.
#      Bearings are near-isotropic and the north-south split is 50/50.
#
#   2. Did the composition of newly invaded regions shift between the European
#      regional groupings? That is what the submitted caption meant by
#      "expanded westward". Answered here.
#
#   3. Did newly invaded regions move to higher latitudes over time? That is
#      what the caption meant by "northward", and it is not the same as
#      question 1, which is relative to each event's own source rather than
#      absolute. Also answered here.
#
# The groupings come from Europe_region_classification.csv, extracted from the
# author's Europe_regions shapefile and used in the earlier WNV projection work.
# In this cohort every country falls wholly within one grouping, so the
# composition analysis is in practice a country-composition analysis; the script
# checks and reports that rather than assuming it.
#
# Run from the project folder:
#   WNV_RUN_DIR=run_withBalkans_v2 Rscript code/12_regional_composition.R
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("regional_composition")
on.exit(lg$close(), add = TRUE)

cls_file <- file.path(ROOT, "Europe_region_classification.csv")
if (!file.exists(cls_file)) stop("Not found: ", cls_file, call. = FALSE)

cls <- read.csv(cls_file) %>%
  transmute(NUTS2_ID = trimws(NUTS2_ID), CNTR_CODE = trimws(CNTR_CODE),
            europe_region = trimws(europe_region)) %>%
  filter(nzchar(europe_region)) %>%
  distinct(NUTS2_ID, .keep_all = TRUE)

message("classification: ", nrow(cls), " NUTS2 regions across ",
        dplyr::n_distinct(cls$CNTR_CODE), " countries")
print(as.data.frame(cls %>% count(europe_region, name = "nuts2_regions")), row.names = FALSE)

## ---- Assign every invaded region to a grouping ------------------------------

ev <- invasion %>%
  mutate(NUTS2_ID = substr(NUTS_ID, 1, 4),
         CNTR_CODE = substr(NUTS_ID, 1, 2)) %>%
  left_join(cls %>% select(NUTS2_ID, europe_region), by = "NUTS2_ID")

unmatched <- ev %>% filter(is.na(europe_region))
message("\nfirst invasions: ", nrow(ev), "; unmatched on NUTS2 prefix: ", nrow(unmatched))

if (nrow(unmatched)) {
  print(as.data.frame(unmatched %>% count(CNTR_CODE, NUTS2_ID, name = "regions")),
        row.names = FALSE)
  modal <- cls %>%
    count(CNTR_CODE, europe_region, name = "n") %>%
    group_by(CNTR_CODE) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(CNTR_CODE, fallback = europe_region)
  ev <- ev %>%
    left_join(modal, by = "CNTR_CODE") %>%
    mutate(assigned_by = ifelse(is.na(europe_region), "country fallback", "NUTS2 code"),
           europe_region = ifelse(is.na(europe_region), fallback, europe_region)) %>%
    select(-fallback)
  message("after fallback, still unassigned: ", sum(is.na(ev$europe_region)))
} else {
  ev$assigned_by <- "NUTS2 code"
}

ev <- ev %>%
  filter(!is.na(europe_region)) %>%
  mutate(europe_region = factor(europe_region,
                                levels = c("South", "Central", "West", "North")),
         period = period_of(first_year))

print(as.data.frame(ev %>% count(assigned_by, name = "regions")), row.names = FALSE)

## ---- Is the grouping country-uniform in this cohort? ------------------------

per_country <- ev %>% distinct(CNTR_CODE, europe_region) %>%
  count(CNTR_CODE, name = "groupings")
split_countries <- per_country %>% filter(groupings > 1)
message("\ncohort countries spanning more than one grouping: ", nrow(split_countries))
if (nrow(split_countries)) print(as.data.frame(split_countries), row.names = FALSE)

message("grouping membership of the cohort:")
print(as.data.frame(ev %>% distinct(CNTR_CODE, europe_region) %>%
                      arrange(europe_region, CNTR_CODE) %>%
                      group_by(europe_region) %>%
                      summarise(countries = paste(CNTR_CODE, collapse = " "),
                                .groups = "drop")), row.names = FALSE)

# Groupings that can never receive an event, because no cohort country belongs
# to them. Any claim about movement into such a grouping is untestable here and
# must not be made.
empty_groups <- setdiff(levels(ev$europe_region), unique(as.character(ev$europe_region)))
if (length(empty_groups)) {
  message("\nGROUPINGS WITH NO COHORT COUNTRY, and therefore no possible event: ",
          paste(empty_groups, collapse = ", "))
  message("A claim of movement into these groupings cannot be supported by these data.")
}

## ---- Composition by period --------------------------------------------------

comp <- ev %>%
  count(period, europe_region, name = "n", .drop = FALSE) %>%
  group_by(period) %>%
  mutate(total = sum(n), pct = round(100 * n / total, 1)) %>%
  ungroup() %>%
  arrange(period, europe_region)

wide <- comp %>%
  select(period, europe_region, n) %>%
  tidyr::pivot_wider(names_from = europe_region, values_from = n, values_fill = 0) %>%
  arrange(period)

message("\nFirst invasions by period and grouping:")
print(as.data.frame(comp), row.names = FALSE)
message("\nCounts:")
print(as.data.frame(wide), row.names = FALSE)

write.csv(comp, file.path(REV, "rev_regional_composition.csv"), row.names = FALSE)
write.csv(wide, file.path(REV, "rev_regional_composition_wide.csv"), row.names = FALSE)

## ---- Did the composition change? --------------------------------------------
# Empty groupings are dropped before testing; a column of structural zeros gives
# an expected count of zero and returns NaN rather than a test.

ev_test <- ev %>% mutate(europe_region = droplevels(europe_region))
tab <- table(ev_test$period, ev_test$europe_region)
tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]

cs <- suppressWarnings(chisq.test(tab))
message("\nChi-squared, period by grouping (empty groupings dropped): chi2 = ",
        round(unname(cs$statistic), 2), ", df = ", unname(cs$parameter),
        ", p = ", signif(cs$p.value, 3))
mc_p <- NA_real_
if (any(cs$expected < 5)) {
  sim <- chisq.test(tab, simulate.p.value = TRUE, B = 20000)
  mc_p <- sim$p.value
  message("  some expected counts below 5; Monte Carlo p = ", signif(mc_p, 3))
}

## ---- The westward claim ------------------------------------------------------
# Tested as a trend in the share falling in the West grouping across the ordered
# periods, which is the specific claim, rather than as a general association.

ev <- ev %>% mutate(is_west = as.integer(europe_region == "West"),
                    period_idx = as.integer(period))

share <- ev %>% group_by(period) %>%
  summarise(n = n(), n_west = sum(is_west),
            pct_west = round(100 * sum(is_west) / n(), 1), .groups = "drop")

message("\nShare of first invasions in the West grouping, by period:")
print(as.data.frame(share), row.names = FALSE)

trend <- glm(is_west ~ period_idx, family = binomial, data = ev)
ct <- summary(trend)$coefficients
or <- exp(ct["period_idx", "Estimate"])
ci <- exp(ct["period_idx", "Estimate"] + c(-1.96, 1.96) * ct["period_idx", "Std. Error"])
message("Trend in the West share: odds ratio per period ", round(or, 2),
        " (", round(ci[1], 2), " to ", round(ci[2], 2), "), p = ",
        signif(ct["period_idx", "Pr(>|z|)"], 3))

write.csv(share, file.path(REV, "rev_regional_trend.csv"), row.names = FALSE)

## ---- The northward claim -----------------------------------------------------
# Absolute latitude of the newly invaded region, which is what "northward" means
# in plain reading, and is a different question from the bearing of each event
# relative to its own source.

cent <- nuts_all %>%
  st_transform(4326) %>%
  { suppressWarnings(st_centroid(.)) }
lat <- data.frame(NUTS_ID = cent$NUTS_ID, lat = st_coordinates(cent)[, 2])

ev <- ev %>% left_join(lat, by = "NUTS_ID")
message("\nevents with a centroid latitude: ", sum(!is.na(ev$lat)), " of ", nrow(ev))

lat_tab <- ev %>% filter(!is.na(lat)) %>%
  group_by(period) %>%
  summarise(n = n(), median_lat = round(median(lat), 2),
            q25 = round(quantile(lat, 0.25), 2),
            q75 = round(quantile(lat, 0.75), 2), .groups = "drop")

message("\nLatitude of newly invaded regions, by period:")
print(as.data.frame(lat_tab), row.names = FALSE)

sp <- suppressWarnings(cor.test(ev$first_year, ev$lat, method = "spearman"))
lm_lat <- lm(lat ~ first_year, data = ev)
sl <- summary(lm_lat)$coefficients["first_year", ]

message("Spearman correlation between year of first invasion and latitude: rho = ",
        round(unname(sp$estimate), 3), ", p = ", signif(sp$p.value, 3))
message("Linear trend: ", round(sl["Estimate"], 4), " degrees of latitude per year",
        " (SE ", round(sl["Std. Error"], 4), "), p = ", signif(sl["Pr(>|t|)"], 3),
        "; over the 20 years that is ", round(sl["Estimate"] * 20, 2), " degrees.")

write.csv(lat_tab, file.path(REV, "rev_latitude_by_period.csv"), row.names = FALSE)

## ---- One tidy table of every test -------------------------------------------

write.csv(data.frame(
  test = c("chi-squared, period by grouping",
           "chi-squared, Monte Carlo",
           "logistic trend, share in the West grouping",
           "Spearman, year against latitude",
           "linear trend, latitude on year"),
  estimate = c(round(unname(cs$statistic), 2), NA,
               round(or, 3), round(unname(sp$estimate), 3),
               round(unname(sl["Estimate"]), 4)),
  interval = c(paste("df", unname(cs$parameter)), NA,
               paste0(round(ci[1], 3), " to ", round(ci[2], 3)), NA,
               paste0("SE ", round(unname(sl["Std. Error"]), 4))),
  p_value = c(signif(cs$p.value, 3), signif(mc_p, 3),
              signif(ct["period_idx", "Pr(>|z|)"], 3),
              signif(sp$p.value, 3), signif(unname(sl["Pr(>|t|)"]), 3))
), file.path(REV, "rev_regional_tests.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(REV, "sessionInfo_12.txt"))
message("\n12_regional_composition.R done")

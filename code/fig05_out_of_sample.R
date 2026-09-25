###############################################################################
# fig05_out_of_sample.R
#
# Figure 5. Out-of-sample test of the fitted model on the 2025 and 2026 seasons.
# The model was fitted to the end of 2024 and nothing was refitted.
#
#   A  Europe, showing where the model placed its risk before the 2025 season and
#      which regions were then invaded
#   B  observed invasion rate by decile of the model's predicted risk, fitted
#      period against out of sample: calibration
#   C  receiver operating characteristic, in sample against out of sample:
#      discrimination
#
#   bash run_fig.sh 5
#
# Rewritten 10 September 2026. What changed and why:
#
#   - Panel A now follows the Figure 1A map rules: a country-level grey base for
#     everything outside the cohort, the cohort drawn at region level with no
#     internal boundaries, ETRS89 Lambert azimuthal equal-area, a scale bar, and
#     the legend anchored to a piece of ground rather than to a fraction of the
#     panel so it cannot drift onto the data.
#   - The shading is the FITTED MODEL's predicted hazard, not neighbourhood
#     pressure. The previous version cut its tiers from ntile(invasion_pressure)
#     while the panel was captioned as the model's out-of-sample test, so it was
#     showing one thing and claiming another. model_score in the prospective
#     output is the model's linear predictor for each still-naive region,
#     computed from the invasion state at the end of 2024.
#   - Tier boundaries are stated as percentiles of that score and the legend says
#     so, rather than leaving "risk tier" undefined.
#   - Regions already invaded, and regions invaded in 2025 or 2026, are marked
#     with small symbols inside the region rather than with outlines around it.
#     Previously the latter were unexplained dark boundaries.
#   - Panel B contrasts the fitted period with the out-of-sample seasons pooled,
#     rather than splitting 2025 from 2026. The question the panel answers is
#     whether the gradient the model learned held up afterwards, which needs the
#     fitted period on the same axes as its reference.
#   - Panel C was a density of distances with a rug of events beneath it. It is
#     now the ROC for the fitted period against the two out-of-sample seasons
#     pooled, which is the comparison that says whether discrimination held up.
#     WNV_FIG5C=density restores the old panel.
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig05")
on.exit(lg$close(), add = TRUE)

check_greyscale()

MAP_CRS <- as.integer(Sys.getenv("WNV_MAP_CRS", "3035"))
message("map CRS: EPSG:", MAP_CRS)

## =============================================================================
## Panel A. Where the model placed its risk, and what happened
## =============================================================================

invaded_2024 <- first_outbreak %>%
  filter(!is.na(first_year), first_year <= 2024) %>% pull(NUTS_ID)
message("regions invaded by the end of 2024: ", length(invaded_2024))

TIERS <- c("Lowest 50%", "50-75%", "75-90%", "Highest 10%")

# Regions already invaded carry no predicted risk: the model ranks only the
# regions still naive, which is what makes this a prospective test. They are
# filled with the project teal lightened towards white, which sits well clear of
# both the red scale and the country base and needs no symbol of its own. Only
# the regions invaded in the two following seasons are marked, with a small
# circle drawn inside the region rather than an outline around it.
PRIOR_KEY  <- "Invaded, 2006-2024"
OOS_KEY    <- "Invaded, 2025-2026"
OUTSIDE    <- "Not included in the study"
PRIOR_FILL <- "#5E8791"      # COL_P1 at 72% over white
MARK_SIZE  <- 0.75

LEVELS <- c(TIERS, PRIOR_KEY, OUTSIDE)

# Country base: a neutral grey with no boundaries drawn at all, inside or out.
# Nothing on this map is outlined; every edge is the edge of a fill. Kept in step
# with Figure 1A, where it has to stay clearly darker than the cohort's "no
# recorded invasion" tone.
BASE_FILL <- "#D4D8DB"

# Red rather than blue: the quantity is risk, and red is what a reader reads as
# danger without being told. The already-invaded class is a level of the same
# scale rather than a missing value, so that it carries a legend key.
FILLS <- stats::setNames(
  c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", PRIOR_FILL, BASE_FILL), LEVELS)


# Percentile of the fitted model's predicted hazard among the regions still
# naive at the end of 2024. p25 is the prospective scoring set for 2025.
tiers <- p25 %>%
  mutate(pct  = dplyr::percent_rank(model_score),
         tier = cut(pct, breaks = c(-Inf, 0.5, 0.75, 0.9, Inf), labels = TIERS)) %>%
  select(NUTS_ID, tier)

message("regions per predicted-risk tier:")
print(as.data.frame(table(tiers$tier)), row.names = FALSE)

base0 <- tryCatch(
  get_eurostat_geospatial(resolution = 10, nuts_level = 0, year = 2021) %>%
    select(CNTR_CODE, geometry),
  error = function(e) { message("could not fetch NUTS0: ", conditionMessage(e)); NULL })
stopifnot(!is.null(base0))

region_layer <- nuts %>%
  left_join(tiers, by = "NUTS_ID") %>%
  mutate(cls = ifelse(NUTS_ID %in% invaded_2024, PRIOR_KEY, as.character(tier)),
         cls = factor(cls, levels = LEVELS))

# The seam colour is precomputed rather than derived with after_scale(fill).
# PT170 has no covariates and therefore no model score, so its class is NA;
# after_scale returned an NA colour for it and ggplot2 rejected the entire layer,
# which is what produced empty output files. A plain vector cannot fail that way.
# The region falls back to the country base colour, so it reads as unmodelled
# rather than as a risk class.
region_layer <- region_layer %>%
  mutate(seam = ifelse(is.na(cls), BASE_FILL, unname(FILLS[as.character(cls)])))
if (any(is.na(region_layer$cls))) {
  message("regions with no model score, drawn in the base colour: ",
          paste(region_layer$NUTS_ID[is.na(region_layer$cls)], collapse = " "))
}

print(as.data.frame(region_layer %>% st_drop_geometry() %>%
                      count(cls, name = "regions")), row.names = FALSE)

prior_sf <- region_layer %>% filter(NUTS_ID %in% invaded_2024)
obs_sf   <- region_layer %>% filter(NUTS_ID %in% c(inv25$NUTS_ID, inv26$NUTS_ID))
message("regions already invaded: ", nrow(prior_sf))
message("regions invaded in 2025 or 2026: ", nrow(obs_sf))

# The Figure 1A window, carried into whichever CRS is in use.
win <- st_as_sfc(st_bbox(c(xmin = -10, xmax = 45, ymin = 32, ymax = 70), crs = 4326))
if (MAP_CRS != 4326) win <- st_transform(win, MAP_CRS)
bb <- st_bbox(win)

# The key for this layer is a COLUMN of base0, not a loose object referenced from
# aes(). saveRDS stores a plot with a reference to the environment it was built
# in, so a loose name resolves here and is gone when the combiner reopens the
# panel in a new session; that is what emptied the combined Figure 1.
base0$key <- factor(OUTSIDE, levels = OUTSIDE)

# The legend is built as three groups, each from its own scale, so that the map
# reads as three separate statements: what the model predicts, what has in fact
# been invaded, and what lies outside the cohort. One scale can only produce one
# legend box, so the risk tiers keep the fill scale, the two invaded keys share
# an alpha scale, and the country base takes the colour scale.
INV_KEYS <- c(PRIOR_KEY, OOS_KEY)
tier_layer  <- region_layer %>% filter(as.character(cls) %in% TIERS)
prior_layer <- region_layer %>%
  filter(as.character(cls) == PRIOR_KEY) %>%
  mutate(key = factor(PRIOR_KEY, levels = INV_KEYS))
none_layer  <- region_layer %>% filter(is.na(cls))

base0        <- st_transform(base0, MAP_CRS)
tier_layer   <- st_transform(tier_layer, MAP_CRS)
prior_layer  <- st_transform(prior_layer, MAP_CRS)
none_layer   <- st_transform(none_layer, MAP_CRS)
region_layer <- st_transform(region_layer, MAP_CRS)
prior_sf     <- st_transform(prior_sf, MAP_CRS)
obs_sf       <- st_transform(obs_sf, MAP_CRS)

# Only the out-of-sample invasions are marked. The already-invaded regions are
# identified by their fill, so a symbol on each of the 300-odd of them added
# stipple without adding information.
marks_sf <- suppressWarnings(st_point_on_surface(obs_sf)) %>%
  transmute(key = factor(OOS_KEY, levels = INV_KEYS)) %>%
  st_as_sf()
message("markers drawn: ", nrow(marks_sf))

# Anchored to ground, left edge, as in Figure 1A. Nudge without editing:
#   WNV_LEG_LON=38 WNV_LEG_LAT=57 bash run_fig.sh 5
LEG_LON <- as.numeric(Sys.getenv("WNV_LEG_LON", "31"))
LEG_LAT <- as.numeric(Sys.getenv("WNV_LEG_LAT", "56"))
LEG <- legend_at_lonlat(lon = LEG_LON, lat = LEG_LAT, bbox = bb, crs = MAP_CRS)
message("legend left edge anchored at lon ", LEG_LON, ", lat ", LEG_LAT,
        " -> panel position x ", round(LEG["x"], 3), ", y ", round(LEG["y"], 3))

fA <- ggplot() +
  geom_sf(data = base0, aes(colour = key), fill = BASE_FILL, linewidth = 0.06) +
  geom_sf(data = none_layer, fill = BASE_FILL, colour = BASE_FILL,
          linewidth = 0.06) +
  geom_sf(data = tier_layer, aes(fill = cls), colour = tier_layer$seam,
          linewidth = 0.06) +
  geom_sf(data = prior_layer, aes(alpha = key), fill = PRIOR_FILL,
          colour = PRIOR_FILL, linewidth = 0.06) +
  geom_sf(data = marks_sf, aes(size = key), shape = 16, colour = COL_OOS,
          stroke = 0.3) +
  scale_fill_manual(name = "Predicted West Nile virus invasion risk",
                    values = FILLS[TIERS], limits = TIERS, drop = FALSE,
                    na.value = BASE_FILL) +
  # The two invaded keys share one legend box because the alpha and size scales
  # carry the same name and the same breaks, and ggplot2 merges legends that
  # agree on both. Each key then shows only the mark that belongs to it: the
  # polygon is drawn at alpha 0 on the out-of-sample key, and the point at size 0
  # on the already-invaded key, so neither leaves a ghost on the other.
  scale_alpha_manual(name = NULL,
                     values = stats::setNames(c(1, 0), INV_KEYS),
                     limits = INV_KEYS, drop = FALSE) +
  scale_size_manual(name = NULL,
                    values = stats::setNames(c(0, MARK_SIZE), INV_KEYS),
                    limits = INV_KEYS, drop = FALSE) +
  scale_colour_manual(name = NULL,
                      values = stats::setNames(BASE_FILL, OUTSIDE)) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  add_scale_bar() +
  # The fill keys would otherwise be drawn with the geom's default black border,
  # since the per-region seam colour is a vector and cannot be carried into a key.
  guides(fill = guide_legend(byrow = TRUE, order = 1,
                             override.aes = list(colour = NA)),
         alpha = guide_legend(byrow = TRUE, order = 2,
                              override.aes = list(size = c(0, 1.6))),
         size  = guide_legend(byrow = TRUE, order = 2,
                              override.aes = list(size = c(0, 1.6))),
         colour = guide_legend(order = 3,
                               override.aes = list(fill = BASE_FILL,
                                                   colour = NA))) +
  labs(tag = "A") +
  theme_wnv_map(base = 9) +
  theme(legend.position = LEG,
        legend.justification = c(0, 0.5),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.margin = margin(2, 4, 2, 2),
        legend.spacing.y = unit(2.6, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4.5, "mm"),
        legend.title = element_text(size = 8.5, colour = INK),
        legend.text = element_text(size = 8, colour = INK),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 10, colour = INK),
        plot.tag.position = c(0.01, 0.99),
        plot.margin = margin(2, 2, 2, 2))

## =============================================================================
## The model score, shared by panels B and C
## =============================================================================
## Both the score and the outcome for the fitted period are taken from the fitted
## object rather than from ac. predict() drops any row with a missing covariate,
## so scoring the model's own rows against ac$event could misalign the two
## vectors without raising an error; the model's linear predictors and the status
## column of its own response are aligned by construction.

lp_fit <- as.numeric(fit_primary$linear.predictors)
y_fit  <- as.numeric(fit_primary$y[, ncol(fit_primary$y)])
stopifnot(length(lp_fit) == length(y_fit))
message("fitted-period rows scored: ", length(lp_fit), ", events ", sum(y_fit))

oos_all <- bind_rows(p25, p26)
message("out-of-sample rows scored: ", nrow(oos_all),
        ", invasions ", sum(oos_all$invaded))

## =============================================================================
## Panel B. Invasion rate by decile of the model's predicted risk
## =============================================================================
## Calibration. Panel C asks whether the model ranks regions correctly; this panel
## asks whether the regions it ranked high were in fact invaded at a higher rate,
## and whether that held out of sample. Each series is cut into deciles of its own
## predicted score, so decile 10 is the tenth of that series the model considered
## most at risk.
##
## This panel used to plot deciles of neighbourhood invasion pressure, which made
## it a coarsened redraw of Figure 2B: same quantity, same covariate, and the ten
## out-of-sample values were the identical ten numbers. Figure 2B keeps the
## empirical Phi relationship; this panel now carries the full model, so the two
## no longer overlap and the figure set gains a calibration check it did not have.

# Bare year ranges, matching Figures 2A, 2B and 4. Both series in panel B are
# OBSERVED invasion rates: regions invaded over regions at risk. Nothing there is
# fitted or predicted, and the model's only role is to define the covariate on the
# horizontal axis. Labelling one series "model fit" invited the reader to take the
# teal bars for fitted values and the gap between the two for model error. Which
# period the model was fitted on is stated in the caption, where it belongs.
# 2007, not 2006: the fitted period begins in 2007 because distance to the
# nearest invaded region is undefined in the first year of follow-up.
SER <- c(fit = "2007-2024", oos = "2025-2026")

# Panel C labels its series separately. In C the two curves ARE the model, scoring
# its own data against seasons it never saw, so the qualifier is true there; in B
# both series are observed counts and it is not. One shared vector for both panels
# was the reason C lost these when B was corrected.
SER_C <- c(fit = "2007-2024, model fit", oos = "2025-2026, predicted")
SER_COLS <- stats::setNames(c(TEAL, COL_OOS), unname(SER))

decile_rate <- function(score, invaded, label) {
  data.frame(score = score, invaded = invaded) %>%
    filter(is.finite(score)) %>%
    mutate(dec = ntile(score, 10)) %>%
    group_by(dec) %>%
    summarise(rate = 100 * mean(invaded), n = n(), ev = sum(invaded),
              .groups = "drop") %>%
    mutate(series = label)
}

dec_fit <- decile_rate(lp_fit, y_fit, unname(SER["fit"]))
dec_oos <- decile_rate(oos_all$model_score, oos_all$invaded, unname(SER["oos"]))

dec <- bind_rows(dec_fit, dec_oos) %>%
  mutate(series = factor(series, levels = unname(SER)))

print(as.data.frame(dec), row.names = FALSE)

# The two series sit on different baselines: 306 invasions in 17 305 region-years
# against 56 in 1520, so 2025 and 2026 were more active seasons than the fitted
# average. The rates are printed to the log and stated in the caption rather than
# carried in the key.
ov <- dec %>% group_by(series) %>%
  summarise(ov = 100 * sum(ev) / sum(n), .groups = "drop")
message("overall invasion rate by series:")
print(as.data.frame(ov), row.names = FALSE)

fB <- ggplot(dec, aes(factor(dec), rate, fill = series)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.74) +
  scale_fill_manual(values = SER_COLS, limits = unname(SER), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Decile of predicted invasion risk (1 = lowest, 10 = highest)",
       y = "Observed annual NUTS3 region invasion rate (%)",
       tag = "B") +
  theme_wnv_plot(base = 7) +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(b = 0), legend.box.spacing = unit(1.5, "mm"),
        legend.key.height = unit(3, "mm"), legend.key.width = unit(3.4, "mm"),
        plot.tag = element_text(size = 10, colour = INK),
        plot.tag.position = c(0.01, 0.99))

## =============================================================================
## Panel C. Discrimination out of sample
## =============================================================================

as_density <- identical(tolower(Sys.getenv("WNV_FIG5C", "")), "density")

## An ROC sweeps every possible cut point on a single ranking score; no threshold
## is chosen. Regions are sorted by score from highest to lowest, and at each
## position the true positive rate is the share of regions that were invaded
## lying above the cut and the false positive rate is the share of regions that
## were not. The area under the curve is then the probability that a randomly
## chosen invaded region scored above a randomly chosen uninvaded one, with 0.5
## the diagonal.
##
## The score is the fitted model's linear predictor throughout, not any one
## covariate. For the fitted period it is taken from the model on its own data,
## one value per region-year, against whether that region-year was a first
## invasion. For the out-of-sample seasons it is model_score in the prospective
## output, one value per still-naive region per season, against whether that
## region was invaded in the season. The two seasons are pooled, since the same
## model produced both scores.

roc_of <- function(score, y) {
  keep <- is.finite(score) & !is.na(y)
  score <- score[keep]; y <- y[keep]
  o <- order(score, decreasing = TRUE)
  y <- y[o]
  data.frame(fpr = c(0, cumsum(1 - y) / sum(1 - y)),
             tpr = c(0, cumsum(y) / sum(y)))
}
auc_of <- function(r) sum(diff(r$fpr) * (head(r$tpr, -1) + tail(r$tpr, -1)) / 2)

if (as_density) {
  bg <- bind_rows(p25 %>% mutate(season = "2025"), p26 %>% mutate(season = "2026"))
  fC <- ggplot(bg, aes(nearest_invaded_km)) +
    geom_density(fill = COL_NEVER, colour = SOFT, linewidth = 0.3, adjust = 1.1) +
    geom_rug(data = bg %>% filter(invaded == 1), sides = "b",
             colour = COL_OOS, alpha = 0.8, length = unit(3, "pt")) +
    coord_cartesian(xlim = c(0, 800)) +
    labs(x = "Distance to the nearest region already invaded by West Nile virus (kilometres)",
         y = "Density, all regions at risk", tag = "C") +
    theme_wnv_plot(base = 7) +
    theme(plot.tag = element_text(size = 10, colour = INK),
          plot.tag.position = c(0.01, 0.99))
} else {
  oos <- oos_all

  roc_dat <- bind_rows(
    roc_of(lp_fit, y_fit)                %>% mutate(series = unname(SER_C["fit"])),
    roc_of(oos$model_score, oos$invaded) %>% mutate(series = unname(SER_C["oos"]))) %>%
    mutate(series = factor(series, levels = unname(SER_C)))

  aucs <- roc_dat %>% group_by(series) %>%
    group_modify(~ data.frame(auc = auc_of(.x))) %>% ungroup()
  message("area under the ROC curve:")
  print(as.data.frame(aucs), row.names = FALSE)
  message("Harrell's C for the fitted model, for comparison: ",
          round(summary(fit_primary)$concordance[1], 3))

  lab <- aucs %>% mutate(lab = sprintf("%s (AUC %.2f)", series, auc))

  fC <- ggplot(roc_dat, aes(fpr, tpr, colour = series)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = RULE,
                linewidth = 0.3) +
    geom_step(linewidth = 0.6) +
    scale_colour_manual(values = stats::setNames(unname(SER_COLS), unname(SER_C)),
                        limits = unname(SER_C),
                        labels = stats::setNames(lab$lab, lab$series),
                        name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = c(0.01, 0)) +
    coord_fixed() +
    labs(x = "False positive rate", y = "True positive rate", tag = "C",
         caption = "AUC, area under the receiver operating characteristic curve") +
    theme_wnv_plot(base = 7) +
    theme(legend.position = c(0.98, 0.05), legend.justification = c(1, 0),
          legend.background = element_blank(),
          legend.key.height = unit(3, "mm"), legend.key.width = unit(4, "mm"),
          plot.caption = element_text(size = 5.6, colour = MID, hjust = 0,
                                      margin = margin(t = 2)),
          plot.tag = element_text(size = 10, colour = INK),
          plot.tag.position = c(0.01, 0.99))
}

## =============================================================================
## Assembly
## =============================================================================
## The map takes the full width, because it carries a legend and a scale bar and
## because it is the panel a reader looks at first. B and C share the row
## beneath, each wrapped so that patchwork aligns panels within a row rather than
## across the whole figure.

f5 <- wrap_elements(fA) / wrap_elements(fB | fC) +
  plot_layout(heights = c(1.45, 1))

save_fig(f5, "Figure5_out_of_sample", WIDTH_2COL, 215)

v <- read.csv(file.path(REV, "rev_prospective_validation_summary.csv")) %>%
  filter(evaluation_set == "all regions")

write_caption("Figure5_out_of_sample", paste(
  "Figure 5. Out-of-sample test of the diffusion model for the 2025 and 2026",
  "seasons. The model was fitted on data to the end of 2024 and nothing was",
  "refitted; every region was scored using only the invasion history available",
  "before the season began. (A) Regions still naive at the end of 2024 are shaded by",
  "the percentile of the model's predicted invasion hazard for the 2025 season, so",
  "the deepest red is the tenth of naive regions the model ranked highest. The score",
  "uses only information available to the end of 2024. Regions already invaded",
  "during 2006 to 2024 carry no predicted risk, since the model ranks only regions",
  "still at risk; they take the pale teal fill and their own legend entry. Circles",
  "mark the regions invaded in 2025 or 2026, across both seasons. Countries outside the study",
  "cohort are shown in light grey at country level. No boundaries are drawn anywhere",
  "on the map, so every edge is the edge of a filled area, and the projection is",
  "ETRS89 Lambert azimuthal equal-area (EPSG:3035). (B) Observed invasion rate by",
  "decile of the model's predicted invasion risk, decile 1 being the tenth of",
  "regions the model ranked lowest and decile 10 the tenth it ranked highest. Each",
  "series is cut into deciles of its own scores: for the fitted period the model's",
  "linear predictor across the 16 998 region-years of 2007 to 2024, and for the",
  "out-of-sample seasons the score assigned to each still-naive region before the",
  "season began, with 2025 and 2026 pooled. Both series are observed rates, not",
  "predictions; the panel asks whether the regions the model ranked high were in",
  "fact invaded more often, which is calibration, where panel C asks whether the",
  "ranking is correct at all, which is discrimination. The out-of-sample bars sit",
  "higher than the fitted bars across the range because 2025 and 2026 were more",
  "active seasons than the fitted average, 58 first invasions in 1487 region-seasons",
  "against 323 in 16 998 region-years, and not because the model underpredicts.",
  if (as_density)
    "(C) Distance to the nearest invaded region across all regions at risk, with the observed first invasions marked beneath."
  else
    paste("(C) Receiver operating characteristic curves for the model's ranking of",
          "regions, in the fitted period and in the two out-of-sample seasons pooled.",
          "The score is the model's linear predictor: for the fitted period one value",
          "per region-year against whether that year was a first invasion, and for the",
          "out-of-sample seasons the score assigned to each still-naive region before",
          "the season began against whether it was then invaded. The curve traces every",
          "possible cut point on that ranking rather than any chosen threshold, and the",
          "area beneath it is the probability that an invaded region scored above an",
          "uninvaded one. Taken season by season the area under the curve was",
          sprintf("%.2f in 2025 and %.2f in 2026, against %.2f and %.2f for neighbourhood pressure alone and %.2f and %.2f for proximity alone.",
                  v$auc_full_model[v$season == 2025], v$auc_full_model[v$season == 2026],
                  v$auc_pressure[v$season == 2025], v$auc_pressure[v$season == 2026],
                  v$auc_proximity[v$season == 2025], v$auc_proximity[v$season == 2026])),
  paste("The out-of-sample curve falls below the diagonal across its upper range.",
        "That is not noise: it means a small group of regions the model ranked low",
        "were nonetheless invaded. These are long-distance introductions, which a",
        "model built on diffusion from neighbouring regions cannot anticipate by",
        "construction, and they are the same events behind the raised invasion rate",
        "in the lowest pressure decile in panel B. The weaker 2025 result is driven",
        "by one such introduction into Ile-de-France counted across eight adjacent",
        "regions; the cluster-collapsed figures are given in the appendix."),
  OOS_LABEL))

message("fig05 done")

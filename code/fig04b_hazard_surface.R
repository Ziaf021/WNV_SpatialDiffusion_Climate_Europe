###############################################################################
# fig04b_hazard_surface.R
#
# Figure 4B. Predicted distance by temperature invasion-hazard surface, from the
# interaction model, with the observed first invasions drawn over it.
#
#   bash run_fig.sh 4b
#
# Written 10 September 2026 against the submitted panel in
# Fig4_Main_TempDistInteractions.R. This panel had not been regenerated on the
# refitted cohort at all.
#
# Two things are done differently from the submitted code, both deliberate.
#
# 1. The surface is computed from the fitted coefficients rather than through
#    predict(). The interaction model carries scale(pr) in its formula, and
#    scale() has no makepredictcall method, so it is recomputed on whatever
#    newdata is supplied. The submitted grid holds pr at a single constant, which
#    gives a zero standard deviation and therefore NaN for every row. Every
#    covariate other than distance and temperature is held fixed in this panel,
#    so those terms contribute a constant to the linear predictor which cancels
#    in both the percentile rank and the relative hazard. Computing the linear
#    predictor directly is exact here and cannot silently return NaN.
#
# 2. The grid is drawn with geom_raster rather than re-binned with
#    stat_summary_2d(bins = 60). The prediction grid is already regular, so
#    binning it a second time only blurs it.
#
# Scale. The submitted panel showed the within-grid percentile rank of predicted
# risk. That is kept as the default because it is what the manuscript describes,
# but note that a percentile rank depends on the extent of the grid: widen the
# distance axis and every percentile changes. WNV_FIG4B=hazard draws the relative
# hazard against a stated reference instead, which does not have that property.
#
# Colour. The submitted ramp was RdYlBu reversed. The default here keeps the same
# cool-to-warm structure in the project hues, which is smoother in the middle and
# safer for deuteranopia. WNV_FIG4B_RAMP=submitted restores the original five
# stops exactly.
#
# Points. The submitted panel carried a shape legend reading "Observed WNV
# invasion". That legend is kept, split into the two seasons, in the teal and
# reserved plum of Figure 2A and at the size and opacity used in panel A.
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig04b")
on.exit(lg$close(), add = TRUE)

as_hazard <- identical(tolower(Sys.getenv("WNV_FIG4B", "")), "hazard")
use_sub   <- identical(tolower(Sys.getenv("WNV_FIG4B_RAMP", "")), "submitted")

SURF_RAMP <- if (use_sub) {
  c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
} else {
  c("#2E6FA7", "#9DC3DE", "#FCE79B", "#F2A48F", "#CB181D")
}
message("surface ramp: ", if (use_sub) "submitted RdYlBu" else "project cool-to-warm")
message("surface quantity: ", if (as_hazard) "relative hazard" else "percentile rank")

## ---- The fitted interaction model -------------------------------------------

int_path <- file.path(RUN, "cox_interaction.rds")
if (!file.exists(int_path)) {
  stop("Not found: ", int_path, "\nRun run_revision_v2.sh first.", call. = FALSE)
}
fit_int <- readRDS(int_path)

cf <- coef(fit_int)
need <- c("distance100", "tas", "distance100:tas")
absent <- setdiff(need, names(cf))
if (length(absent)) {
  stop("Interaction model is missing: ", paste(absent, collapse = ", "),
       call. = FALSE)
}
b_d <- cf[["distance100"]]; b_t <- cf[["tas"]]; b_dt <- cf[["distance100:tas"]]

message(sprintf("distance %.4f | tas %.4f | distance x tas %.4f", b_d, b_t, b_dt))
message(sprintf("HR for the interaction term: %.3f", exp(b_dt)))

## ---- Prediction grid ---------------------------------------------------------

D_MAX <- 750                     # km, the submitted horizontal extent
t_rng <- range(ac$tas, na.rm = TRUE)

grid <- expand.grid(
  km  = seq(0, D_MAX, length.out = 240),
  tas = seq(t_rng[1], t_rng[2], length.out = 240)) %>%
  mutate(distance100 = km / 100,
         lp = b_d * distance100 + b_t * tas + b_dt * distance100 * tas)

# Reference for the relative-hazard version: the median distance and median
# temperature among the observed first invasions.
ev_fit <- ac %>% filter(event == 1)
d_ref <- median(ev_fit$nearest_invaded_km, na.rm = TRUE) / 100
t_ref <- median(ev_fit$tas, na.rm = TRUE)
lp_ref <- b_d * d_ref + b_t * t_ref + b_dt * d_ref * t_ref

grid <- grid %>%
  mutate(pct = dplyr::percent_rank(lp) * 100,
         rel = exp(lp - lp_ref),
         z   = if (as_hazard) rel else pct)

stopifnot(all(is.finite(grid$z)))

message(sprintf("reference: %.0f km, %.1f C | relative hazard range %.2f to %.2f",
                d_ref * 100, t_ref, min(grid$rel), max(grid$rel)))

## ---- Observed events ---------------------------------------------------------
## The fitted and out-of-sample pair, matching Figure 2A and panel A exactly.

# 2007, not 2006: the fitted period begins in 2007 because distance to the
# nearest invaded region is undefined in the first year of follow-up.
SETS <- c(fitted = "2007-2024", oos = "2025-2026")
SET_COLS <- stats::setNames(c(TEAL, COL_OOS), unname(SETS))
PT_SIZE  <- 1.1
PT_ALPHA <- 0.95

obs <- bind_rows(
  ev_fit %>% transmute(km = nearest_invaded_km, tas, set = unname(SETS["fitted"])),
  inv25 %>% transmute(km = nearest_invaded_km, tas, set = unname(SETS["oos"])),
  inv26 %>% transmute(km = nearest_invaded_km, tas, set = unname(SETS["oos"]))) %>%
  filter(!is.na(km), !is.na(tas), km <= D_MAX) %>%
  mutate(set = factor(set, levels = unname(SETS)))

message("observed points drawn: ")
print(as.data.frame(table(obs$set)), row.names = FALSE)

## ---- Panel -------------------------------------------------------------------

## The white marks across the colourbar are its tick marks. guide_colourbar(ticks
## = FALSE) stopped governing them in ggplot2 3.5, where they became the
## legend.ticks theme element, so both are set. The tryCatch keeps the script
## running on an older ggplot2, where naming that element is an error.
tick_theme <- tryCatch(theme(legend.ticks = element_blank()),
                       error = function(e) theme())

fill_scale <- if (as_hazard) {
  scale_fill_gradientn(colours = SURF_RAMP, trans = "log10",
                       name = sprintf("Hazard relative to\n%.0f km and %.1f °C",
                                      d_ref * 100, t_ref),
                       guide = guide_colourbar(order = 1, ticks = FALSE,
                                               ticks.colour = NA,
                                               frame.colour = NA))
} else {
  scale_fill_gradientn(colours = SURF_RAMP, limits = c(0, 100),
                       breaks = c(0, 25, 50, 75, 100),
                       name = "Predicted invasion\nhazard percentile",
                       guide = guide_colourbar(order = 1, ticks = FALSE,
                                               ticks.colour = NA,
                                               frame.colour = NA))
}

## Contours are off by default. On a continuous surface they impose steps that
## the quantity does not have, and the eye follows the line rather than the
## gradient. WNV_FIG4B_CONTOUR=1 draws them.
show_contours  <- identical(Sys.getenv("WNV_FIG4B_CONTOUR", "0"), "1")
contour_breaks <- if (as_hazard) c(0.25, 0.5, 1, 2, 4) else c(25, 50, 75, 90)

f4b <- ggplot(grid, aes(km, tas)) +
  geom_raster(aes(fill = z), interpolate = TRUE) +
  {if (show_contours)
    geom_contour(aes(z = z), breaks = contour_breaks, colour = INK,
                 linewidth = 0.2, alpha = 0.35)} +
  geom_point(data = obs, aes(km, tas, colour = set), inherit.aes = FALSE,
             shape = 16, size = PT_SIZE, alpha = PT_ALPHA) +
  scale_colour_manual(values = SET_COLS, limits = unname(SETS),
                      name = "Reported West Nile virus\ninvasion in humans",
                      guide = guide_legend(order = 2,
                        override.aes = list(size = 2, alpha = 1))) +
  fill_scale +
  scale_x_continuous(limits = c(0, D_MAX), breaks = seq(0, D_MAX, 150),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Distance to the nearest region already invaded by West Nile virus (kilometres)",
       y = "Annual mean temperature (°C)") +
  theme_wnv_plot(base = 7) +
  theme(legend.position = "top", legend.box = "horizontal",
        legend.key.height = unit(2.6, "mm"), legend.key.width = unit(9, "mm"),
        legend.title = element_text(size = 6.2, colour = MID),
        legend.text = element_text(size = 5.8, colour = MID)) +
  tick_theme

saveRDS(f4b, file.path(DIR_FIG, "panel_4B.rds"))
save_fig(f4b, "Figure4B_hazard_surface", 120, 108)

write_caption("Figure4B_hazard_surface", paste(
  "Figure 4B. Predicted hazard of first West Nile virus invasion across the joint",
  "range of distance to the nearest invaded region and annual mean temperature,",
  "from the Cox model containing their interaction, with neighbourhood invasion",
  "pressure, source-region transmission intensity, relative humidity and",
  "precipitation held at their means and urbanization at its reference level.",
  if (as_hazard)
    sprintf(paste("Shading is the hazard relative to a region %0.f km from the",
                  "nearest invaded region at %.1f degrees C, the medians among",
                  "observed first invasions, on a log colour scale."),
            d_ref * 100, t_ref)
  else
    paste("Shading is the percentile rank of the predicted hazard within the",
          "plotted grid, so it orders risk across the surface rather than giving",
          "an absolute rate."),
  if (show_contours)
    paste("Contours mark", paste(contour_breaks, collapse = ", "),
          if (as_hazard) "fold." else "percent.")
  else "",
  "Points are observed first invasions, in teal for the fitted period and in the",
  "reserved colour for the 2025 and 2026 seasons, which were not used to fit the",
  "model. Warmth and proximity act together: at a given distance a warmer region",
  "carries the higher hazard, and the temperature gradient steepens as distance",
  "grows, which is the interaction the model estimates.", OOS_LABEL))

message("fig04b done")

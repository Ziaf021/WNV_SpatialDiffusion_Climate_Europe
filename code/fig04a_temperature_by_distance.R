###############################################################################
# fig04a_temperature_by_distance.R
#
# Figure 4A. Annual mean temperature at first invasion, by distance class.
#
#   bash run_fig.sh 4a
#
# Rewritten 10 September 2026 against the submitted panel in
# Fig4_Main_TempDistInteractions.R. Carried back from the submitted version: the
# per-class sample sizes above each violin, the Kruskal-Wallis annotation on the
# panel, and the white interior boxplot.
#
# Colour. The submitted pale-to-deep red ramp is kept, running light at the
# shortest distances to dark at the longest, so the fill deepens with distance.
# WNV_FIG4A_RAMP=amber draws the project amber distance ramp instead, which is
# what Figure 3C uses for the same five classes.
#
# Points. Both seasons are drawn at the same size, shape and opacity and carry a
# legend, in the teal and reserved plum of Figure 2A.
#
# The violins have no outline. An outline on a filled shape that already carries
# the variable in its fill adds a second contour for the eye to follow and, at
# 0.3 pt, prints as a grey fringe rather than a line.
#
# NOTE ON saveRDS. Anything referenced inside aes() must be a column of the
# layer's own data. A loose object in the calling environment resolves while this
# script runs and is gone when fig04_combined.R reopens the saved panel in a new
# session, which is what broke the combined figure on 10 September. Values inside
# annotate() and coord_cartesian() are evaluated at build time and are safe.
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig04a")
on.exit(lg$close(), add = TRUE)

use_amber <- identical(tolower(Sys.getenv("WNV_FIG4A_RAMP", "")), "amber")
FILL_RAMP <- if (use_amber) {
  DIST_COLOURS
} else {
  stats::setNames(colorRampPalette(c("#FEE5D9", "#CB181D"))(length(DIST_LEVELS)),
                  DIST_LEVELS)
}
message("panel A fill ramp: ", if (use_amber) "project amber" else "reds, light near to dark far")

# The fitted and out-of-sample pair, matching Figure 2A exactly.
# 2007, not 2006: the fitted period begins in 2007 because distance to the
# nearest invaded region is undefined in the first year of follow-up.
SETS <- c(fitted = "2007-2024", oos = "2025-2026")
SET_COLS <- stats::setNames(c(TEAL, COL_OOS), unname(SETS))
SET_EDGE <- stats::setNames(c("#1C4E59", COL_OOS_EDGE), unname(SETS))
PT_SIZE  <- 1.1
PT_ALPHA <- 0.55

## Violin scaling. "area", the ggplot default and the scaling of the submitted
## panel, gives every violin equal area, so a class with more events is drawn
## wider. WNV_FIG4A_VIOLIN=width gives every violin the same maximum width
## instead, which isolates the shape of the density from the sample size. Either
## is defensible here because the n labels state the sample sizes outright.
## 22 Sep 2026: default changed to "width". With "area" every violin has the same
## area, so the 9-event >500 km class, whose temperatures are tightly clustered, is
## drawn widest and the 186-event <=75 km class is squeezed to a sliver. "width"
## gives every class the same maximum width; n is printed above each violin.
## WNV_FIG4A_VIOLIN=area restores the previous behaviour.
VIOLIN_SCALE <- if (identical(tolower(Sys.getenv("WNV_FIG4A_VIOLIN", "")), "area"))
  "area" else "width"
message("violin scaling: ", VIOLIN_SCALE)

## The interior boxplot: a very light grey body with a single slate outline. The
## median is the box's own median line, thickened through fatten, rather than a
## mark drawn over it; overplotting a second line left a short dash sitting on
## top of the first.
BOX_FILL <- "#F2F1EF"
BOX_LINE <- "#4D5A66"

## Horizontal spread of the points. Wide enough that they sit across the body of
## the violin rather than piling onto the box in the middle of it.
JIT_W <- 0.2

## ---- Data --------------------------------------------------------------------

fitted_ev <- ac %>% filter(event == 1) %>%
  transmute(km = nearest_invaded_km, tas, set = unname(SETS["fitted"]))

# The bias-adjusted covariate product ends in 2024, so the out-of-sample events
# are placed at each region's own 2020 to 2024 mean temperature. Only the
# temperature axis is substituted; the distance is fully observed.
oos_ev <- bind_rows(
  inv25 %>% transmute(km = nearest_invaded_km, tas),
  inv26 %>% transmute(km = nearest_invaded_km, tas)) %>%
  filter(!is.na(tas)) %>%
  mutate(set = unname(SETS["oos"]))

ev <- bind_rows(fitted_ev, oos_ev) %>%
  mutate(cls = distance_class(km),
         set = factor(set, levels = unname(SETS)))

fit <- ev %>% filter(set == unname(SETS["fitted"]))
oos <- ev %>% filter(set == unname(SETS["oos"]))

## One kernel bandwidth for all five classes. With a per-class bandwidth the
## >500 km class (9 events, 5 of them near 16 degC) got a very small bandwidth
## and was drawn as separate beads. The pooled Silverman bandwidth of the
## fitted-period temperatures smooths every class alike.
VIOLIN_BW <- bw.nrd0(fit$tas[is.finite(fit$tas)])
message("violin bandwidth (pooled, degC): ", round(VIOLIN_BW, 3))

kw <- kruskal.test(tas ~ cls, data = fit)
message("Kruskal-Wallis, fitted period only: chi2 = ", round(kw$statistic, 3),
        ", df = ", kw$parameter, ", p = ", signif(kw$p.value, 3))

pw <- pairwise.wilcox.test(fit$tas, fit$cls, p.adjust.method = "BH")
message("pairwise Wilcoxon, BH adjusted:")
print(round(pw$p.value, 4))

message("out-of-sample events plotted: ", nrow(oos), " of ",
        nrow(inv25) + nrow(inv26), " (the remainder have no covariates)")

y_lo   <- min(ev$tas, na.rm = TRUE)
y_hi   <- max(ev$tas, na.rm = TRUE)
y_n    <- y_hi + 0.55
y_test <- y_hi + 1.75

# y is a column, not a loose scalar, so the panel survives saveRDS.
#
# Both series are plotted, so the count above each class has to describe both.
# Earlier it counted the fitted period alone while out-of-sample points were
# drawn in the same column, which understated every class. The split is drawn as
# two pieces, each in the colour of the series it counts, matching Figure 2A.
n_lab <- ev %>%
  group_by(cls, .drop = FALSE) %>%
  summarise(n_fitted = sum(set == unname(SETS["fitted"])),
            n_oos    = sum(set == unname(SETS["oos"])),
            .groups = "drop") %>%
  mutate(y = y_n,
         fit_lab = as.character(n_fitted),
         oos_lab = ifelse(n_oos > 0, paste0("+ ", n_oos), ""),
         h_fit = ifelse(n_oos > 0, 1.15, 0.5))

# df is reported in the caption; on the panel it only lengthens the line.
kw_lab <- sprintf("Kruskal-Wallis chi-squared = %.1f, %s",
                  kw$statistic,
                  ifelse(kw$p.value < 0.0001, "p < 0.0001",
                         paste0("p = ", signif(kw$p.value, 2))))

## ---- Panel -------------------------------------------------------------------

f4a <- ggplot(ev, aes(cls, tas)) +
  geom_violin(data = fit, aes(fill = cls), colour = NA, alpha = 0.85,
              scale = VIOLIN_SCALE, bw = VIOLIN_BW, trim = TRUE) +
  geom_boxplot(data = fit, width = 0.11, outlier.shape = NA, colour = BOX_LINE,
               fill = BOX_FILL, linewidth = 0.3, fatten = 2.4) +
  geom_point(data = fit, aes(fill = NULL, colour = set),
             position = position_jitter(width = JIT_W, height = 0, seed = 1),
             size = PT_SIZE, alpha = PT_ALPHA, shape = 16) +
  geom_point(data = oos, aes(fill = NULL, colour = set),
             position = position_jitter(width = JIT_W, height = 0, seed = 2),
             size = PT_SIZE, alpha = PT_ALPHA, shape = 16) +
  geom_text(data = n_lab, aes(cls, y, label = fit_lab, hjust = h_fit),
            inherit.aes = FALSE, size = 2.1, colour = TEAL,
            fontface = "bold") +
  geom_text(data = n_lab %>% filter(nzchar(oos_lab)),
            aes(cls, y, label = oos_lab), inherit.aes = FALSE,
            hjust = -0.20, size = 2.1, colour = COL_OOS, fontface = "bold") +
  annotate("text", x = 0.45, y = y_test, label = kw_lab, hjust = 0, size = 2.1,
           colour = MID) +
  scale_fill_manual(values = FILL_RAMP, limits = DIST_LEVELS, guide = "none") +
  scale_colour_manual(values = SET_COLS, limits = unname(SETS),
                      name = "Reported West Nile virus\ninvasion in humans",
                      guide = guide_legend(
                        override.aes = list(size = 2, alpha = 1))) +
  scale_x_discrete(labels = c("≤75", "76-150", "151-300", "301-500", ">500")) +
  coord_cartesian(ylim = c(y_lo - 0.5, y_test + 0.9), clip = "off") +
  labs(x = "Distance to the nearest region already invaded by West Nile virus (kilometres)",
       y = "Annual mean temperature (°C)") +
  theme_wnv_plot(base = 7) +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 5.8, colour = MID, lineheight = 1.1),
        legend.text = element_text(size = 5.8, colour = MID),
        legend.key.height = unit(3, "mm"), legend.key.width = unit(3.4, "mm"),
        legend.margin = margin(0, 0, 0, 0))

saveRDS(f4a, file.path(DIR_FIG, "panel_4A.rds"))
save_fig(f4a, "Figure4A_temperature_by_distance", 120, 102)

write_caption("Figure4A_temperature_by_distance", paste(
  "Figure 4A. Annual mean temperature at the time of first human West Nile virus",
  "invasion, by the distance from the newly invaded region to the nearest region",
  "already invaded at the start of that season. Violins show the density within each",
  "distance class and deepen in colour with distance, the interior box gives the",
  "median and interquartile range, and points are the individual first invasions.",
  "Violins and boxes summarise the fitted period alone, on the same events as the",
  "test;",
  "above each class, n is the number of first invasions in the fitted period, followed in the out-of-sample colour by the number added in 2025 and 2026. Regions",
  "invaded at greater distances were warmer than those invaded through short-range",
  "diffusion",
  sprintf("(Kruskal-Wallis chi-squared = %.1f, df = %d, p = %s;",
          kw$statistic, kw$parameter, format.pval(kw$p.value, digits = 2)),
  "test computed on the fitted period only). Invasions in 2025 and 2026 are shown in",
  "the reserved out-of-sample colour and are placed at each region's own 2020 to 2024",
  "mean temperature, because the bias-adjusted covariate product ends in 2024; their",
  "distances are fully observed and only the temperature axis is substituted. These",
  "points are excluded from the test.", OOS_LABEL))

message("fig04a done")

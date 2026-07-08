############################################################
# Figure 4 (main text):
#   Panel A - temperature by distance-to-invaded violin
#   Panel B - predicted distance x temperature hazard surface
#
# Panel B is a model-predicted percentile surface: a distance x
# temperature grid is scored with cox_interaction (other covariates
# held at their mean, urbanisation at its reference level) and the
# predicted risk is shown as a within-grid percentile rank.
#
# Bug-fix / update notes:
#   - analysis_dat -> analysis_cases throughout (incl. the Panel B
#     geom_point layer, which had reverted to analysis_dat)
#   - added library(ggpubr) for theme_classic2() and ggarrange()
#   - Panel A: removed duplicate n_labels and dead 'reds' definitions
#     (kept the colorRampPalette version), na.rm = TRUE on tas ranges
#   - Panel B: consolidated the fill-legend title to a single source
#   - export bug fix retained (Fig_TempDistance saved under its own name)
#
# Expects in the environment:
#   analysis_cases   counting-process data
#   cox_interaction  interaction Cox model used for the Panel B surface
#                    (newdat must supply every covariate in its formula)
############################################################

## ---- Libraries ----
library(dplyr)
library(ggplot2)
library(ggthemes)
library(scales)
library(ggpubr)     # theme_classic2(), ggarrange()

## ---- Distance categories on first-invasion events ----
first_invaded <- analysis_cases %>%
  filter(event == 1)

first_invaded$distance_cat <- cut(
  first_invaded$nearest_invaded_km,
  breaks = c(0, 75, 150, 300, 500, Inf),
  include.lowest = TRUE
)

# Not used by the figures below; retained in case it is referenced elsewhere.
first_invaded$distance_class <- cut(
  first_invaded$nearest_invaded_km,
  breaks = c(0, 50, 100, 250, 500, Inf),
  include.lowest = TRUE
)

## ---- Group tests (reported in text) ----
kruskal.test(tas ~ distance_cat, data = first_invaded)

pairwise.wilcox.test(
  first_invaded$tas,
  first_invaded$distance_cat,
  p.adjust.method = "BH"
)

## ---- Sample-size labels and fill palette ----
n_labels <- first_invaded %>%
  group_by(distance_cat) %>%
  summarise(n = n(), .groups = "drop")

n_labels$y <- max(first_invaded$tas, na.rm = TRUE) + 0.5

reds <- colorRampPalette(c("#FEE5D9", "#CB181D"))(
  nlevels(first_invaded$distance_cat)
)

## ---- Panel A: temperature by distance category ----
Fig_TempDistance <-
  ggplot(first_invaded, aes(x = distance_cat, y = tas, fill = distance_cat)) +
  geom_violin(alpha = 0.8, linewidth = 0.5, colour = "grey30") +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1.5, colour = "black") +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_text(
    data = n_labels,
    aes(x = distance_cat, y = y, label = paste0("n=", n)),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_fill_manual(values = reds, guide = "none") +
  annotate(
    "text",
    x = 1.2,
    y = max(first_invaded$tas, na.rm = TRUE) + 1.8,
    label = "Kruskal-Wallis p < 0.001",
    hjust = 0,
    size = 4
  ) +
  labs(
    x = "Distance from nearest invaded region (km)",
    y = "Annual mean temperature (°C)"
  ) +
  coord_cartesian(
    ylim = c(
      min(first_invaded$tas, na.rm = TRUE) - 0.5,
      max(first_invaded$tas, na.rm = TRUE) + 2.5
    )
  ) +
  theme_classic2(base_size = 13) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

Fig_TempDistance

## ---- Panel B: predicted distance x temperature hazard surface ----

# Prediction grid: vary distance and temperature, hold everything else fixed.
newdat <- expand.grid(
  distance100 = seq(0, 10, length.out = 100),
  tas = seq(
    min(analysis_cases$tas, na.rm = TRUE),
    max(analysis_cases$tas, na.rm = TRUE),
    length.out = 100
  )
)

newdat$pressure_z        <- mean(analysis_cases$pressure_z, na.rm = TRUE)
newdat$log_cases_nearest <- mean(analysis_cases$log_cases_nearest, na.rm = TRUE)
newdat$hurs              <- mean(analysis_cases$hurs, na.rm = TRUE)
newdat$pr                <- mean(analysis_cases$pr, na.rm = TRUE)
newdat$URBN_TYPE <- factor(
  levels(factor(analysis_cases$URBN_TYPE))[1],
  levels = levels(factor(analysis_cases$URBN_TYPE))
)

# newdat must supply every covariate in cox_interaction's formula.
newdat$pred_hazard <- predict(cox_interaction, newdata = newdat, type = "risk")
newdat$hazard_pct  <- percent_rank(newdat$pred_hazard) * 100

Fig_DistanceTempHazard <-
  ggplot(newdat, aes(x = distance100 * 100, y = tas, z = hazard_pct)) +
  stat_summary_2d(bins = 60, aes(fill = after_stat(value))) +
  scale_fill_gradientn(
    colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
    limits = c(0, 100),
    name = "Predicted invasion\nhazard percentile"
  ) +
  labs(
    x = "Distance to nearest invaded region (km)",
    y = "Annual mean temperature (°C)"
  ) +
  scale_x_continuous(limits = c(0, 1000 + 10), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_point(
    data = analysis_cases %>% filter(event == 1),
    aes(x = nearest_invaded_km, y = tas, shape = "Observed WNV invasion"),
    inherit.aes = FALSE,
    colour = "black",   # shape 19 uses colour, not fill
    size = 2,
    alpha = 0.8
  ) +
  scale_shape_manual(name = "", values = c("Observed WNV invasion" = 19)) +
  guides(
    fill  = guide_colorbar(order = 1),
    shape = guide_legend(order = 2, override.aes = list(size = 3, alpha = 1))
  ) +
  theme_classic2(base_size = 13) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )

Fig_DistanceTempHazard

## ---- Export individual panels ----
ggsave("Fig_DistanceTempHazard.svg", Fig_DistanceTempHazard, width = 8, height = 6)
ggsave("Fig_DistanceTempHazard.pdf", Fig_DistanceTempHazard, width = 8, height = 6)

# Saved under its own name (originally clobbered the surface plot's file).
ggsave("Fig_TempDistance.svg", Fig_TempDistance, width = 8, height = 6)
ggsave("Fig_TempDistance.pdf", Fig_TempDistance, width = 8, height = 6)

## ---- Combined A/B panel ----
Fig4_Main_TempDist_interaction_AB <-
  ggarrange(
    Fig_TempDistance,
    Fig_DistanceTempHazard,
    labels = c("A", "B"),
    ncol = 1
  )

Fig4_Main_TempDist_interaction_AB

ggsave("Fig4_Main_TempDist_interaction_AB.svg",
       Fig4_Main_TempDist_interaction_AB, width = 7, height = 9)
ggsave("Fig4_Main_TempDist_interaction_AB.pdf",
       Fig4_Main_TempDist_interaction_AB, width = 7, height = 9)
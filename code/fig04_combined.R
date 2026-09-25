###############################################################################
# fig04_combined.R
#
# Figure 4, panels A and B stacked, at the proportions of the submitted figure
# (9 by 7 inches, A above B).
#
#   bash run_fig.sh 4a     then
#   bash run_fig.sh 4b     then
#   bash run_fig.sh 4
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig04")
on.exit(lg$close(), add = TRUE)

pa <- file.path(DIR_FIG, "panel_4A.rds")
pb <- file.path(DIR_FIG, "panel_4B.rds")

for (f in c(pa, pb)) {
  if (!file.exists(f)) {
    stop("Not found: ", f,
         "\nBuild the panels first:  bash run_fig.sh 4a   and   bash run_fig.sh 4b",
         call. = FALSE)
  }
}

fig4 <- (readRDS(pa) / readRDS(pb)) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 12, face = "bold", colour = INK),
                  plot.margin = margin(4, 4, 4, 4)))

# Narrower and taller than the two-column maximum. The panels are stacked, so
# width beyond what the axis labels need only flattens them; 168 mm still holds
# the longest label in either panel, which is the shared x-axis title.
save_fig(fig4, "Figure4_Main_AB", 168, 210)

write_caption("Figure4_Main_AB", paste(
  "Figure 4. Temperature and distance act together on West Nile virus invasion",
  "risk. (A) Annual mean temperature at the time of first invasion, by the distance",
  "from the newly invaded region to the nearest region already invaded at the start",
  "of that season. Violins show the density within each distance class, the interior",
  "box the median and interquartile range, and points the individual first",
  "invasions; above each class, n is the number of first invasions in the fitted period, followed in the out-of-sample colour by the number added in 2025 and 2026.",
  "Regions invaded at greater distances were warmer than those invaded through",
  "short-range diffusion. (B) Predicted hazard of first invasion across the joint",
  "range of distance and annual mean temperature, from the Cox model containing",
  "their interaction, with the remaining covariates held at their means and",
  "urbanization at its reference level; white contours are drawn on the shaded",
  "quantity and points are observed first invasions. At a given distance a warmer",
  "region carries the higher hazard, and the temperature gradient steepens as",
  "distance grows. Invasions in the 2025 and 2026 seasons are shown throughout in",
  "the reserved out-of-sample colour and were not used to fit the model; in panel A",
  "they are placed at each region's own 2020 to 2024 mean temperature, because the",
  "bias-adjusted covariate product ends in 2024, and they are excluded from the",
  "Kruskal-Wallis test.", OOS_LABEL))

message("fig04 combined done")

###############################################################################
# fig01_combined.R
#
# Figure 1, panels A and B side by side on one canvas, at the proportions of the
# submitted figure (15 by 6 inches, two columns).
#
# Reads the panel objects saved by fig01a and fig01b, so each panel can be
# adjusted on its own and only this script re-run to see them together.
#
#   bash run_fig.sh 1a     then
#   bash run_fig.sh 1b     then
#   bash run_fig.sh 1
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig01")
on.exit(lg$close(), add = TRUE)

pa <- file.path(DIR_FIG, "panel_1A.rds")
pb <- file.path(DIR_FIG, "panel_1B.rds")

for (f in c(pa, pb)) {
  if (!file.exists(f)) {
    stop("Not found: ", f,
         "\nBuild the panels first:  bash run_fig.sh 1a   and   bash run_fig.sh 1b",
         call. = FALSE)
  }
}

f1a <- readRDS(pa)
f1b <- readRDS(pb)

# The map carries more white space than the heatmap, so it takes slightly more
# width; the two panels then read as the same size on the page.
fig1 <- (f1a | f1b) +
  plot_layout(widths = c(1.05, 1)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "NUTS3, Nomenclature of Territorial Units for Statistics (Eurostat, 2021",
      "release), level 3. GAUL, Global Administrative Unit Layers (for Serbia and",
      "Albania). *Data to 10 September 2026."),
    theme = theme(
      plot.tag = element_text(size = 12, face = "bold", colour = INK),
      plot.caption = element_text(size = 7, colour = MID, hjust = 0,
                                  margin = margin(t = 4)),
      plot.margin = margin(4, 4, 4, 4)))

# 15 by 6 inches, as submitted.
save_fig(fig1, "Figure1_Main_AB", 381, 152)

## ---- The whole-figure caption ----------------------------------------------
# The submitted caption is revised rather than replaced. Five things in it are no
# longer true: the period runs to 2026, the bands changed, the palette is no
# longer a single light-to-dark ramp, and absence is now described as no recorded
# HUMAN infection rather than shown in blue.
#
# The fifth is a deliberate removal. The submitted caption said emergence
# "expanded westward and northward". Direction is not something this figure
# shows, and the two analyses that do bear on it point different ways: movement
# into the western grouping is strong, movement into the northern grouping cannot
# occur at all because no cohort country belongs to it, and invasion relative to
# each event's own source is near-isotropic. All of that belongs in the response
# to Reviewer 3, not in a caption. The claim is therefore dropped here rather
# than restated in a weaker form.

write_caption("Figure1_Main_AB", paste(
  "Figure 1. Chronology of West Nile virus invasion across Europe, 2006 to 2026.",
  "(A) Geographic distribution of first recorded locally acquired human WNV",
  "infections. Countries outside the study cohort are shown in grey at country",
  "level. Within the 21 cohort countries, each region is shaded by the period of",
  "first documented emergence, from pale yellow for 2006 to 2010 through gold",
  "(2011 to 2015) and ochre (2016 to 2020) to dark bronze for 2021 to 2024, and plum for 2025 to 2026,",
  "the colour reserved throughout for the two seasons that lie outside the period",
  "used to fit the model. Internal region boundaries are not drawn, so the invaded",
  "area reads as a continuous surface. Regions in the pale neutral tone have no",
  "recorded locally acquired human WNV infection over the study period, which does",
  "not exclude avian, equine or entomological detection, since human surveillance",
  "was the only outcome available at this resolution across all countries.",
  "Projection ETRS89 Lambert azimuthal equal-area (EPSG:3035).",
  "(B) Annual count of regions recording their first locally acquired human WNV",
  "infection, by country and year. Countries run from the earliest invaded at the",
  "bottom to the most recently invaded at the top, with ties broken by the number",
  "of regions eventually invaded. Tile shading and the number within each tile give",
  "the count of regions newly invaded in that country and year, with darker tiles",
  "representing larger annual counts; blank tiles are years in which no first",
  "invasion was recorded in that country.",
  "NUTS3, Nomenclature of Territorial Units for Statistics (Eurostat, 2021",
  "release), level 3. GAUL, Global Administrative Unit Layers (for Serbia and",
  "Albania). *Data to 10 September 2026."))

message("fig01 combined done")

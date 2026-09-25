###############################################################################
# fig01f_combined.R
#
# The experimental Figure 1: the period-faceted maps as panel A, the country by
# year heatmap as panel B, stacked in one column rather than side by side.
#
# Vertical because panel A is now a grid of five maps and has an aspect ratio of
# its own; putting it beside the heatmap would either squeeze the maps to
# illegibility or stretch the heatmap across a width its rows do not need.
#
# This does not touch Figure1_Main_AB. It reads panel_1A_facet.rds and
# panel_1B.rds and writes its own file, so the submitted layout stands.
#
#   bash run_fig.sh 1af
#   bash run_fig.sh 1b
#   bash run_fig.sh 1f
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig01f")
on.exit(lg$close(), add = TRUE)

pa <- file.path(DIR_FIG, "panel_1A_facet.rds")
pb <- file.path(DIR_FIG, "panel_1B.rds")

for (f in c(pa, pb)) {
  if (!file.exists(f)) {
    stop("Not found: ", f,
         "\nBuild the panels first:  bash run_fig.sh 1af   and   bash run_fig.sh 1b",
         call. = FALSE)
  }
}

f1a <- readRDS(pa)
f1b <- readRDS(pb)

# wrap_elements around each row, so patchwork aligns nothing across the two.
# Without it the heatmap's y axis labels pull the map grid rightward to share an
# axis that the maps do not have, which opens a gap down the left of panel A.
## Panel A's height is not a taste question. The maps are drawn in an equal-area
## projection with coord_sf, so each holds a fixed ratio of height to width; give
## the row more height than that ratio asks for and the surplus comes back as
## white space above and below the grid, which is what opened the gap between the
## two panels. So the split is computed from the window rather than guessed.

W <- as.numeric(Sys.getenv("WNV_FIG1F_W", "215"))

MAP_CRS <- as.integer(Sys.getenv("WNV_MAP_CRS", "3035"))
win <- st_as_sfc(st_bbox(c(xmin = -10, xmax = 45, ymin = 32, ymax = 70),
                         crs = 4326))
if (MAP_CRS != 4326) win <- st_transform(win, MAP_CRS)
bbw <- st_bbox(win)
asp <- unname((bbw["ymax"] - bbw["ymin"]) / (bbw["xmax"] - bbw["xmin"]))

NCOL <- as.integer(Sys.getenv("WNV_FACET_NCOL", "3"))
NROW <- ceiling(5 / NCOL)
cell_w <- (W - 8) / NCOL
H_A <- NROW * cell_w * asp + 6 * NROW
## The split is 60/40 by default. Panel A keeps the height its map aspect asks
## for, since any surplus there comes straight back as white space, so it is the
## canvas that follows the split rather than the panel: the total is whatever
## makes A's natural height 60% of it, and B takes the remaining 40%.
SPLIT <- as.numeric(strsplit(Sys.getenv("WNV_FIG1F_SPLIT", "0.55,0.45"),
                             ",", fixed = TRUE)[[1]])
SPLIT <- SPLIT / sum(SPLIT)

H <- as.numeric(Sys.getenv("WNV_FIG1F_H", as.character(round(H_A / SPLIT[1] + 6))))

HEIGHTS <- if (nzchar(Sys.getenv("WNV_FIG1F_HEIGHTS", ""))) {
  as.numeric(strsplit(Sys.getenv("WNV_FIG1F_HEIGHTS"), ",", fixed = TRUE)[[1]])
} else SPLIT

message("map aspect (h/w): ", round(asp, 3),
        " | canvas ", round(W), " by ", round(H), " mm",
        " | split ", paste(round(HEIGHTS * 100), collapse = "/"),
        " | panel A natural height ", round(H_A, 1), " mm")

# Panel A is itself a patchwork once the inset is added, so a tag cannot be put
# on it with labs(): on a patchwork that reaches the last nested plot, which is
# the inset. Each panel is wrapped instead and the tags are placed on the two
# wrappers, which are the only elements plot_annotation() counts here.
fig1f <- wrap_elements(full = f1a) / wrap_elements(full = f1b) +
  plot_layout(heights = HEIGHTS) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "NUTS3, Nomenclature of Territorial Units for Statistics (Eurostat, 2021",
      "release), level 3. GAUL, Global Administrative Unit Layers (for Serbia and",
      "Albania). *Data to 10 September 2026."),
    theme = theme(
      plot.tag = element_text(size = 12, face = "bold", colour = INK),
      plot.caption = element_text(size = 7, colour = MID, hjust = 0,
                                  margin = margin(t = 3)),
      plot.margin = margin(2, 3, 2, 3)))

save_fig(fig1f, "Figure1_Facet_AB", W, H)

write_caption("Figure1_Facet_AB", paste(
  "Figure 1, experimental layout. (A) First recorded locally acquired human West",
  "Nile virus infection by NUTS3 region, one map per period. Within each map,",
  "regions invaded during that period are shaded by their position in the period,",
  "first year darkest; regions",
  "invaded in an earlier period are shown in a muted dark gold, so each map reads as",
  "an increment on the accumulated history. Regions with no recorded invasion by",
  "the end of the period take the pale neutral, and countries outside the study",
  "cohort the base grey. The 2006-2010",
  "panel includes the three regions invaded in 2006, which seed the diffusion",
  "process and contribute no time at risk to the model; the 2025-2026 panel is out",
  "of sample and is not used in fitting. The scale bar applies to every panel and",
  "is drawn once, in the 2016-2020 map. Equal-area projection, EPSG:3035.",
  "(B) First invasions by country and year."))

###############################################################################
# fig01a_chronology_map.R
#
# Figure 1A. Chronology of West Nile virus invasion across Europe, 2006 to 2026.
#
# Built to the design of the submitted fig1_chronology_map.R: a country-level
# base layer in dark grey with black boundaries, and the study countries drawn
# over it at region level with no internal boundary colour, so the invaded area
# reads as one continuous surface rather than as a mosaic of administrative
# units. Two departures from the submitted version, both deliberate:
#
#   1. The region layer covers only the countries in the study cohort. In the
#      submitted figure every European region was drawn and those outside the
#      cohort fell into the "no invasion" class, which conflated "not followed"
#      with "followed and never invaded". The grey base now carries the former.
#
#   2. The map is drawn in ETRS89 Lambert azimuthal equal-area rather than in
#      degrees. On an unprojected map the scale varies by roughly a factor of two
#      between southern Greece and northern Scandinavia, so a scale bar on it
#      would be wrong almost everywhere, and this paper is about distance.
#      Set WNV_MAP_CRS=4326 to revert to the submitted plate carree.
#
#   bash run_fig.sh 1a
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig01a")
on.exit(lg$close(), add = TRUE)

check_greyscale()

MAP_CRS <- as.integer(Sys.getenv("WNV_MAP_CRS", "3035"))
message("map CRS: EPSG:", MAP_CRS)

## ---- Layers -----------------------------------------------------------------

# Country-level base for the whole of Europe.
base0 <- tryCatch(
  get_eurostat_geospatial(resolution = 10, nuts_level = 0, year = 2021) %>%
    select(CNTR_CODE, geometry),
  error = function(e) { message("could not fetch NUTS0: ", conditionMessage(e)); NULL })
stopifnot(!is.null(base0))

cohort_countries <- sort(unique(nuts$CNTR_CODE))
message("cohort countries drawn at region level: ", length(cohort_countries), " (",
        paste(cohort_countries, collapse = " "), ")")

# The outcome is a first locally acquired HUMAN infection. A region in this class
# may well have had avian, equine or entomological detections; we simply have no
# human case on record for it. The label has to say so, or the map claims more
# than the data support.
# "No recorded human invasion" reads as an invasion carried out by humans. The
# outcome is a first locally acquired infection IN humans, so the label says so.
ABSENT  <- "No recorded invasion in humans"
OUTSIDE <- "Not included in the study"
LEVELS  <- c(PERIOD_LEVELS, ABSENT, OUTSIDE)

region_layer <- nuts %>%
  left_join(invasion, by = "NUTS_ID") %>%
  mutate(period = factor(ifelse(is.na(as.character(period)), ABSENT,
                                as.character(period)), levels = LEVELS))

print(as.data.frame(region_layer %>% st_drop_geometry() %>%
                      count(period, name = "regions")), row.names = FALSE)

# Country base and seams, matched to Figure 5A so the two maps read as one pair.
# No boundaries are drawn anywhere; each region takes a hairline stroke in its
# own fill colour, which closes the slivers between adjacent polygons at this
# generalisation without introducing a boundary of a different colour.
# Country base. Deliberately darker than the cohort's "no recorded invasion"
# tone (#EAE8E3, a warm off-white). At #E6E8EA the two were within one step of
# each other in lightness and a reader could not tell a country outside the
# cohort from a cohort region that was followed and never invaded, which is the
# one distinction this map exists to make.
BASE_FILL <- "#D4D8DB"

## Period colours. The four fitted periods are ordered bands of time, so they
## take a single-hue sequential ramp, darkest for the earliest, rather than the
## teal-then-amber pair the project palette carries; two hue families across an
## ordered scale make the reader ask whether the change of family means
## something. Luminance rises monotonically, so the order survives greyscale and
## the check_greyscale() test above. 2025-2026 keeps the reserved plum, which is
## the point: it is out of sample and should not sit on the ramp at all.
## Four candidate ramps, switched without editing:
##   WNV_PERIOD_RAMP=amber   bronze to saturated lemon (default)
##   WNV_PERIOD_RAMP=lemon   olive-gold to pale lemon, single hue throughout
##   WNV_PERIOD_RAMP=warm    lemon earliest, deep amber latest (reversed)
##   WNV_PERIOD_RAMP=style   the shared project palette
RAMPS <- list(
  amber = c("#6B4410", "#B07A18", "#E0B02E", "#F4DC63"),
  lemon = c("#7A6A14", "#B49F26", "#DCCB55", "#F1E58C"),
  warm  = c("#F4DC63", "#E0B02E", "#B07A18", "#6B4410"))

ramp_key <- tolower(Sys.getenv("WNV_PERIOD_RAMP", "amber"))
PERIOD_RAMP <- if (identical(ramp_key, "style")) {
  PERIOD_COLOURS
} else {
  if (!ramp_key %in% names(RAMPS)) {
    message("unknown ramp '", ramp_key, "'; using amber")
    ramp_key <- "amber"
  }
  stats::setNames(c(RAMPS[[ramp_key]], COL_OOS), PERIOD_LEVELS)
}
message("ramp key: ", ramp_key)
message("period ramp: ", paste(PERIOD_RAMP, collapse = " "))

fills <- c(PERIOD_RAMP,
           stats::setNames(COL_NEVER, ABSENT),
           stats::setNames(BASE_FILL, OUTSIDE))

region_layer <- region_layer %>%
  mutate(seam = unname(fills[as.character(period)]))
stopifnot(!any(is.na(region_layer$seam)))

## ---- Frame ------------------------------------------------------------------
# The submitted window, longitude -10 to 45 and latitude 32 to 70, carried into
# whichever CRS is in use so the two versions frame the same ground.

win <- st_as_sfc(st_bbox(c(xmin = -10, xmax = 45, ymin = 32, ymax = 70),
                         crs = 4326))
if (MAP_CRS != 4326) win <- st_transform(win, MAP_CRS)
bb <- st_bbox(win)

# The key for this layer is a COLUMN of base0, not a loose object referenced from
# aes(). saveRDS stores a plot with a reference to the environment it was built
# in, so a loose name resolves here and is gone when the combiner reopens the
# panel in a new session; that is what emptied the combined Figure 1.
base0$key <- OUTSIDE
base0        <- st_transform(base0, MAP_CRS)
region_layer <- st_transform(region_layer, MAP_CRS)

# The legend is anchored to a piece of ground rather than to a fraction of the
# panel, so it stays put if the projection or the window changes. The anchor is
# the legend's LEFT edge, not its centre: centring it at a point over Ukraine
# pushed the box back over Romania and Bulgaria, which is what was overlapping
# the data. Anchoring the left edge lets the box grow rightward into the blank
# area east and north of the Black Sea instead.
#
# Nudge it without editing this file:
#   WNV_LEG_LON=38 WNV_LEG_LAT=57 bash run_fig.sh 1a
LEG_LON <- as.numeric(Sys.getenv("WNV_LEG_LON", "31"))
LEG_LAT <- as.numeric(Sys.getenv("WNV_LEG_LAT", "56"))
LEG <- legend_at_lonlat(lon = LEG_LON, lat = LEG_LAT, bbox = bb, crs = MAP_CRS)
message("legend left edge anchored at lon ", LEG_LON, ", lat ", LEG_LAT,
        " -> panel position x ", round(LEG["x"], 3), ", y ", round(LEG["y"], 3))

## ---- Figure -----------------------------------------------------------------

f1a <- ggplot() +
  geom_sf(data = base0, aes(fill = key), colour = NA) +
  geom_sf(data = region_layer, aes(fill = period), colour = region_layer$seam,
          linewidth = 0.06) +
  scale_fill_manual(name = "First WNV invasion", values = fills,
                    limits = LEVELS, drop = FALSE) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  add_scale_bar() +
  # Without this the fill keys carry the geom's default black border, since the
  # per-region seam colour is a vector and cannot be carried into a key.
  guides(fill = guide_legend(byrow = TRUE, override.aes = list(colour = NA))) +
  theme_wnv_map(base = 9) +
  theme(
    legend.position = LEG,
    legend.justification = c(0, 0.5),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.margin = margin(2, 4, 2, 2),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(4.5, "mm"),
    legend.title = element_text(size = 8.5, colour = INK),
    legend.text = element_text(size = 8, colour = INK),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(2, 2, 2, 2))

saveRDS(f1a, file.path(DIR_FIG, "panel_1A.rds"))
save_fig(f1a, "Figure1A_chronology_map", 190, 152)

write_caption("Figure1A_chronology_map", paste(
  "Figure 1A. Chronology of West Nile virus invasion across Europe, 2006 to 2026.",
  "Countries outside the study cohort are shown in grey at country level, with no",
  "boundaries drawn anywhere on the map. Within",
  "the 21 cohort countries, each region is shaded by the period in which its first",
  "locally acquired human infection was recorded; internal region boundaries are",
  "not drawn, so the invaded area reads as a continuous surface. Regions in the",
  "pale neutral tone have no recorded locally acquired human infection over the",
  "study period, which does not exclude avian, equine or entomological detection,",
  "since only human surveillance data were available at this resolution across all",
  "countries. The 2025 and 2026 seasons lie outside the period used to fit the",
  "model.", OOS_LABEL,
  "Regions are NUTS3 units (Nomenclature of Territorial Units for Statistics,",
  "level 3, 2021 release) and the equivalent Eurostat statistical regions for",
  "candidate countries. Projection ETRS89 Lambert azimuthal equal-area",
  "(EPSG:3035)."))

message("fig01a done")

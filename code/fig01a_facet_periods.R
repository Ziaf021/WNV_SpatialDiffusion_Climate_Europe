###############################################################################
# fig01a_facet_periods.R
#
# An experiment, not a manuscript panel. Figure 1A shows the whole chronology on
# one map, which forces a single colour to stand for a five-year band and hides
# the order of arrival inside it. This version gives each period its own map and
# resolves the individual years within it, so the reader can follow the sequence
# rather than the aggregate.
#
# Five facets: the four fitted periods and the two out-of-sample seasons. Within
# each facet a region invaded during that period is filled by its position in
# the period (first year darkest) and carries the year of invasion written over
# its geometry. Regions invaded in an earlier period are shown in a flat neutral
# so each panel reads as an increment on the history rather than as a fresh map,
# and regions not yet invaded by the end of the period take the cohort's
# never-invaded tone.
#
#   bash run_fig.sh 1af
#
# Year text was tried over each invaded region and dropped. At NUTS3 size the
# labels either overlapped or had to be set so small that they were unreadable,
# and suppressing the overlaps meant a region invaded in a crowded panel carried
# no number at all, which is worse than carrying none anywhere.
#
# Switches, none of which require editing this file:
#   WNV_FACET_NCOL=2       facet columns, default 3
#   WNV_MAP_CRS=4326       revert to plate carree
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig01a_facet")
on.exit(lg$close(), add = TRUE)

MAP_CRS   <- as.integer(Sys.getenv("WNV_MAP_CRS", "3035"))
NCOL      <- as.integer(Sys.getenv("WNV_FACET_NCOL", "3"))
message("map CRS: EPSG:", MAP_CRS, " | facet columns: ", NCOL)

## ---- Period bounds ----------------------------------------------------------
## Taken from the level names rather than restated, so this file cannot drift
## from the banding the rest of the figures use.

bounds <- do.call(rbind, lapply(PERIOD_LEVELS, function(p) {
  yr <- as.integer(strsplit(p, "-", fixed = TRUE)[[1]])
  data.frame(period = p, y0 = yr[1], y1 = yr[2], stringsAsFactors = FALSE)
}))
print(bounds, row.names = FALSE)

## ---- Classes and colours ----------------------------------------------------
## The ordinal ramp runs over five steps because the longest band is five years.
## The 2025-2026 facet uses the first two of them, which is correct: a reader
## comparing panels should see the same colour mean the same position in the
## band, not a rescaled version of it.

## Yellows, with age reading as depth of ink: the accumulated history is darkest,
## and within the facet's own period the first year is darker than the last. One
## hue throughout, so the only thing a reader has to decode is lightness. The
## saturation is held high at every step; a desaturated yellow ramp has almost no
## range left to separate the bands and turns to mud.
ORD_LABS <- c("First year of the period", "Second", "Third", "Fourth", "Fifth")
ORD_COLS <- c("#B8860B", "#DBA80F", "#EFC93A", "#F7E066", "#FCF2A8")

EARLIER <- "Invaded in an earlier period"
NOT_YET <- "Not yet invaded"
OUTSIDE <- "Not included in the study"

EARLIER_FILL <- "#7A6B47"   # greyed dark gold: context, not subject
BASE_FILL    <- "#AEB6BB"   # clearly darker than the never-invaded tone

LEVELS_F <- c(ORD_LABS, EARLIER, NOT_YET, OUTSIDE)
FILLS_F  <- stats::setNames(
  c(ORD_COLS, EARLIER_FILL, COL_NEVER, BASE_FILL), LEVELS_F)

## ---- Layers -----------------------------------------------------------------

base0 <- tryCatch(
  get_eurostat_geospatial(resolution = 10, nuts_level = 0, year = 2021) %>%
    select(CNTR_CODE, geometry),
  error = function(e) { message("could not fetch NUTS0: ", conditionMessage(e)); NULL })
stopifnot(!is.null(base0))

regions <- nuts %>%
  left_join(invasion %>% select(NUTS_ID, first_year), by = "NUTS_ID") %>%
  st_transform(MAP_CRS)
base0 <- st_transform(base0, MAP_CRS)

# One copy of every layer per facet. The facet variable is a COLUMN of each
# layer's own data, never a loose object referenced from aes(); a loose name
# resolves while the script runs and is gone the moment the plot is reopened
# elsewhere, which is what emptied earlier combined panels.
facet_of <- function(df, p) { df$period <- factor(p, levels = PERIOD_LEVELS); df }

region_facets <- do.call(rbind, lapply(seq_len(nrow(bounds)), function(i) {
  y0 <- bounds$y0[i]; y1 <- bounds$y1[i]
  d <- regions
  ord <- d$first_year - y0 + 1L
  cls <- ifelse(is.na(d$first_year) | d$first_year > y1, NOT_YET,
         ifelse(d$first_year < y0, EARLIER, ORD_LABS[pmin(pmax(ord, 1L), 5L)]))
  d$cls <- factor(cls, levels = LEVELS_F)
  facet_of(d, bounds$period[i])
}))

base_facets <- do.call(rbind, lapply(bounds$period, function(p) {
  d <- base0
  d$key <- factor(OUTSIDE, levels = OUTSIDE)
  facet_of(d, p)
}))

message("regions per facet and class:")
print(as.data.frame(region_facets %>% st_drop_geometry() %>%
                      count(period, cls, name = "regions")), row.names = FALSE)

# The seam is a precomputed vector, not after_scale(fill). A region with no
# class would give after_scale an NA colour and ggplot2 rejects the whole layer
# on that, which produces an empty output file rather than an error you can read.
region_facets$seam <- unname(FILLS_F[as.character(region_facets$cls)])
stopifnot(!any(is.na(region_facets$seam)))

## ---- Frame ------------------------------------------------------------------

win <- st_as_sfc(st_bbox(c(xmin = -10, xmax = 45, ymin = 32, ymax = 70),
                         crs = 4326))
if (MAP_CRS != 4326) win <- st_transform(win, MAP_CRS)
bb <- st_bbox(win)

## ---- Legend placement -------------------------------------------------------
## The empty cell of the facet grid, in panel coordinates. With three columns and
## five maps it is bottom right; with two columns the grid is full, so the key
## falls back under the block. Nudge without editing:
##   WNV_FACET_LEG=0.86,0.20 bash run_fig.sh 1af

LEG_WHERE <- tolower(Sys.getenv("WNV_FACET_LEGPOS", "bottom"))
LEG_NROW  <- as.integer(Sys.getenv("WNV_FACET_LEGROWS", "2"))

n_free <- (NCOL * ceiling(nrow(bounds) / NCOL)) - nrow(bounds)
LEG_POS <- if (LEG_WHERE %in% c("top", "right", "bottom")) LEG_WHERE else if (n_free > 0) {
  env <- Sys.getenv("WNV_FACET_LEG", "")
  if (nzchar(env)) {
    as.numeric(strsplit(env, ",", fixed = TRUE)[[1]])
  } else {
    col_i <- nrow(bounds) %% NCOL + 1
    n_row <- ceiling(nrow(bounds) / NCOL)
    c((col_i - 0.5) / NCOL, (0.5) / n_row)
  }
} else "bottom"

top_leg  <- LEG_WHERE %in% c("top", "right", "bottom")
LEG_JUST <- if (identical(LEG_WHERE, "right")) {
  "top"
} else if (LEG_WHERE %in% c("top", "bottom")) {
  "left"
} else if (is.character(LEG_POS)) {
  "left"
} else {
  c(0.5, 0.5)
}
# Left-justified, then nudged in from the edge so the key sits under the body of
# the grid rather than hard against its left margin.
LEG_INDENT <- as.numeric(Sys.getenv("WNV_FACET_LEGINDENT", "22"))
LEG_BOX_MAR <- if (LEG_WHERE %in% c("top", "bottom")) {
  ggplot2::margin(l = LEG_INDENT)
} else {
  ggplot2::margin()
}

row_leg  <- LEG_WHERE %in% c("top", "bottom")
LEG_NCOL <- if (row_leg) NULL else 1L
if (!row_leg) LEG_NROW <- NULL
message("legend position: ",
        if (is.character(LEG_POS)) LEG_POS
        else paste(round(LEG_POS, 3), collapse = ", "))

## ---- Layer split for the two keys -------------------------------------------
## The five within-period shades are one ordered quantity and belong in a single
## continuous strip, not in five separate keys with five separate labels. The
## other three classes are not on that scale at all and keep ordinary keys. One
## ggplot2 scale produces one legend box, so the ordered shades keep the fill
## scale and the three context classes move onto an alpha scale whose keys take
## their colours from override.aes.

CTX <- c(EARLIER, NOT_YET, OUTSIDE)
CTX_FILL <- c(EARLIER_FILL, COL_NEVER, BASE_FILL)

ord_layer     <- region_facets %>% filter(as.character(cls) %in% ORD_LABS)
earlier_layer <- region_facets %>% filter(as.character(cls) == EARLIER) %>%
  mutate(ctx = factor(EARLIER, levels = CTX))
notyet_layer  <- region_facets %>% filter(as.character(cls) == NOT_YET) %>%
  mutate(ctx = factor(NOT_YET, levels = CTX))
base_facets$key <- factor(OUTSIDE, levels = CTX)

## ggplot2 gained legend.key.spacing in 3.5; on an older build the argument is
## rejected and the whole theme call fails, so it is added only if it exists.
tight_keys <- tryCatch(theme(legend.key.spacing.x = unit(0, "pt"),
                             legend.key.spacing.y = unit(0, "pt")),
                       error = function(e) theme())

## ---- Figure -----------------------------------------------------------------

f <- ggplot() +
  # Every edge on this map is the edge of a fill. The base takes a hairline in
  # its own colour and each region takes a hairline in its own, which closes the
  # slivers of white page that otherwise show through between adjacent polygons
  # at this generalisation and read as boundaries that were never drawn.
  geom_sf(data = base_facets, aes(alpha = key), fill = BASE_FILL,
          colour = BASE_FILL, linewidth = 0.12) +
  geom_sf(data = notyet_layer, aes(alpha = ctx), fill = COL_NEVER,
          colour = COL_NEVER, linewidth = 0.16) +
  geom_sf(data = earlier_layer, aes(alpha = ctx), fill = EARLIER_FILL,
          colour = EARLIER_FILL, linewidth = 0.16) +
  geom_sf(data = ord_layer, aes(fill = cls), colour = ord_layer$seam,
          linewidth = 0.16)

## A scale bar in one panel, not five. ggspatial::annotation_scale draws itself
## in every facet, which repeats the same statement four times over; a bar built
## from a segment and a label carries a facet column of its own and so appears
## exactly once. Metres are the projection's units, so this is skipped under
## WNV_MAP_CRS=4326, where a single bar would be wrong almost everywhere anyway.

SCALE_PANEL <- Sys.getenv("WNV_FACET_SCALE_PANEL", "2016-2020")
SCALE_KM    <- as.numeric(Sys.getenv("WNV_FACET_SCALE_KM", "500"))
xr <- unname(bb["xmax"] - bb["xmin"])
yr <- unname(bb["ymax"] - bb["ymin"])

if (MAP_CRS != 4326 && SCALE_PANEL %in% PERIOD_LEVELS) {
  # Down in the corner, off the Moroccan coast, where nothing is drawn. Nudge:
  #   WNV_FACET_SCALE_XY=0.06,0.05 bash run_fig.sh 1af
  xy  <- Sys.getenv("WNV_FACET_SCALE_XY", "0.035,0.025")
  fxy <- as.numeric(strsplit(xy, ",", fixed = TRUE)[[1]])
  x0 <- unname(bb["xmin"]) + fxy[1] * xr
  y0 <- unname(bb["ymin"]) + fxy[2] * yr
  x1 <- x0 + SCALE_KM * 1000
  tick <- 0.012 * yr
  fp <- factor(SCALE_PANEL, levels = PERIOD_LEVELS)

  bar <- data.frame(x = c(x0, x0, x1), xend = c(x1, x0, x1),
                    y = c(y0, y0, y0), yend = c(y0, y0 + tick, y0 + tick),
                    period = fp)
  bar_lab <- data.frame(x = (x0 + x1) / 2, y = y0 + tick * 1.5,
                        lab = paste0(SCALE_KM, " km"), period = fp)

  f <- f +
    geom_segment(data = bar, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = MID, linewidth = 0.35) +
    geom_text(data = bar_lab, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, colour = MID, size = 2.1, vjust = 0)
  message("scale bar: ", SCALE_KM, " km in the ", SCALE_PANEL, " panel")
}

ORD_TICKS <- strsplit(Sys.getenv("WNV_FACET_ORDLAB", "1st,2nd,3rd,4th,5th"),
                      ",", fixed = TRUE)[[1]]

f <- f +
  scale_fill_manual(name = "Year of the period",
                    values = stats::setNames(ORD_COLS, ORD_LABS),
                    limits = ORD_LABS, labels = ORD_TICKS, drop = FALSE) +
  scale_alpha_manual(name = NULL,
                     values = stats::setNames(rep(1, length(CTX)), CTX),
                     limits = CTX, drop = FALSE) +
  facet_wrap(~ period, ncol = NCOL) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  # The ordered shades are drawn as one continuous strip with the labels beneath
  # it, so it reads as a scale rather than as five unrelated categories.
  guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE,
                             title.position = "top", title.hjust = 0,
                             label.position = "bottom",
                             keywidth = unit(7, "mm"),
                             keyheight = unit(3.2, "mm"),
                             override.aes = list(colour = NA)),
         alpha = guide_legend(order = 2, byrow = TRUE, ncol = 1,
                              override.aes = list(fill = CTX_FILL,
                                                  colour = NA, alpha = 1))) +
  theme_wnv_map(base = 8) +
  # The key runs under the maps in two rows, headed by the figure's title, clear
  # of both the panel tag and the inset. WNV_FACET_LEGPOS=top moves it above the
  # maps, =right down the side, =cell back into the empty sixth cell, which then
  # leaves no room for the inset. WNV_FACET_LEGROWS sets the number of rows.
  theme(legend.position = LEG_POS,
        legend.justification = LEG_JUST,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(3.4, "mm"),
        legend.key.width = unit(4, "mm"),
        legend.text = element_text(size = 7, colour = INK),
        legend.margin = margin(1, 2, 1, 2),
        legend.box.margin = LEG_BOX_MAR,
        legend.box = "horizontal",
        legend.box.just = "top",
        legend.spacing.x = unit(6, "mm"),
        legend.title = element_text(size = 7.4, colour = INK, hjust = 0,
                                    margin = margin(b = 2)),
        strip.text = element_text(size = 8.5, colour = INK,
                                  margin = margin(b = 2)),
        panel.spacing = unit(1, "mm"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(0, 2, 0, 2)) +
  tight_keys

## ---- The empty cell: cumulative regions invaded ------------------------------
## With the key moved out of the grid, the sixth cell is free. A small curve there
## answers the question the five maps cannot: how fast the total grew. It is the
## same data as the maps, summed, so it adds no claim; it restates one.

cum_dat <- invasion %>%
  count(first_year, name = "new") %>%
  arrange(first_year) %>%
  mutate(total = cumsum(new),
         period = period_of(first_year))

message("cumulative invaded regions by ", max(cum_dat$first_year), ": ",
        max(cum_dat$total))

## Two quantities, two axes. Plotted on one axis the yearly bars are a few pixels
## tall beside a total that reaches the high hundreds, so the bars carry no
## information at all. The cumulative line is rescaled onto the bars' axis and
## read off a right-hand axis of its own; the bars keep a neutral grey, since
## they are the same data the maps already colour by period.

# The rescaling constant is applied to the data, not written into aes(), and the
# right-hand axis transform carries its own copy in a local environment. Anything
# inside aes() or a formula resolves against the global environment, which is a
# reference rather than a value: it is present while this script runs and gone
# when the combiner reopens the saved panel, and the layer then fails on a
# missing object. That is the same trap that emptied Figures 4 and 5 earlier.
k <- max(cum_dat$new) / max(cum_dat$total)
cum_dat$total_scaled <- cum_dat$total * k
cum_inv <- local({ kk <- k; function(x) x / kk })

mini <- ggplot(cum_dat, aes(first_year)) +
  geom_col(aes(y = new), fill = COL_OBS, width = 0.72) +
  geom_line(aes(y = total_scaled), colour = EARLIER_FILL, linewidth = 0.5) +
  geom_point(aes(y = total_scaled), colour = EARLIER_FILL, size = 0.5) +
  scale_x_continuous(breaks = seq(2006, 2026, 5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(
    name = "Newly invaded regions",
    expand = expansion(mult = c(0, 0.10)),
    sec.axis = sec_axis(cum_inv, name = "Cumulative invasions")) +
  labs(x = NULL) +
  theme_wnv_plot(base = 6) +
  theme(axis.title.y = element_text(size = 6, colour = MID),
        axis.title.y.right = element_text(size = 6, colour = EARLIER_FILL),
        axis.text = element_text(size = 5.6, colour = MID),
        axis.text.y.right = element_text(size = 5.6, colour = EARLIER_FILL),
        plot.margin = margin(2, 2, 2, 2))

## Placed in the free cell rather than beside the block. inset_element measures
## from the panel area, so the cell's edges are simple fractions of the grid.
INSET <- as.numeric(strsplit(
  Sys.getenv("WNV_FACET_INSET", "0.69,0.05,0.99,0.44"), ",", fixed = TRUE)[[1]])

if (top_leg && n_free > 0) {
  f <- f + inset_element(mini, left = INSET[1], bottom = INSET[2],
                         right = INSET[3], top = INSET[4],
                         align_to = "panel", on_top = TRUE)
  message("inset placed at ", paste(INSET, collapse = ", "))
} else {
  message("no free cell or key is in it; inset not drawn")
}

saveRDS(f, file.path(DIR_FIG, "panel_1A_facet.rds"))

# The canvas follows the map window rather than a round number. coord_sf holds
# each map to a fixed height-to-width ratio, so any height beyond what that ratio
# asks for is returned as white space above and below the grid.
asp    <- unname((bb["ymax"] - bb["ymin"]) / (bb["xmax"] - bb["xmin"]))
W_1AF  <- as.numeric(Sys.getenv("WNV_FACET_W", "215"))
NROW   <- ceiling(nrow(bounds) / NCOL)
H_1AF  <- NROW * ((W_1AF - 8) / NCOL) * asp + 6 * NROW
message("canvas: ", round(W_1AF), " by ", round(H_1AF), " mm (aspect ",
        round(asp, 3), ")")
save_fig(f, "Figure1A_facet_periods", W_1AF, H_1AF)

write_caption("Figure1A_facet_periods", paste(
  "Figure 1A, faceted variant. First recorded locally acquired human West Nile",
  "virus infection by NUTS3 region, one map per period. Within each map, regions",
  "invaded during that period are shaded by their position in the period, first",
  "year darkest, and carry the calendar year of first invasion. Regions invaded",
  "in an earlier period are shown in a muted dark gold, so each map reads as an",
  "increment on the accumulated history; regions with no recorded invasion by the",
  "end of the period take the cohort's never-invaded tone, and countries outside",
  "the study cohort the darker base grey. The 2006-2010 panel includes the three",
  "regions invaded in 2006, which seed the diffusion process and contribute no",
  "time at risk to the model. The 2025-2026 panel is out of sample: those seasons",
  "are not used in fitting. The inset gives the number of regions recording a",
  "first invasion in each year (bars, left axis) and the running total (line,",
  "right axis). The scale bar applies to every panel",
  "and is drawn",
  "once, in the 2016-2020 map. Equal-area projection, EPSG:3035."))

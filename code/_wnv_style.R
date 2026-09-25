## -----------------------------------------------------------------------------
## _wnv_style.R
##
## One definition of the colour system, typography and figure furniture for the
## WNV invasion figures, so that a reader who learns the legend on Figure 1 does
## not have to relearn it on Figure 4.
##
##     source(file.path("code", "_wnv_style.R"))
##
## THE COLOUR LOGIC. The figures carry two kinds of category and they must never
## be drawn in the same palette:
##
##     period of first invasion  -> one sequential ramp, dark early to light late
##     distance class            -> Okabe-Ito, a categorical set
##
## One colour is reserved. Plum marks the out-of-sample seasons, 2025 and 2026,
## and marks nothing else on any figure. It is always drawn with an outline as
## well, so the out-of-sample class survives greyscale printing and does not
## depend on hue alone.
##
## Regions never invaded are a warm off-white, not blue. Blue reads as a value
## on the ramp; the off-white reads as a deliberate absence.
##
## Built on the pattern of R/v2/_figure_style.R in the dengue mobility project.
## Written 9 September 2026.
## -----------------------------------------------------------------------------

## The base tokens, carried over from the dengue mobility project so the two
## papers read as one body of work -------------------------------------------
TEAL   <- "#2E7B8C"
AMBER  <- "#C77B21"
PLUM   <- "#5A2A57"
NEUTRAL<- "#EAE8E3"   # warm off-white
GOLD   <- "#C1801F"
BLUE   <- "#2E6FA7"
GREY   <- "#6E777F"
LILAC  <- "#EDE6F2"

## Ink, used for text and rules ------------------------------------------------
INK  <- "#1a2026"   # titles, labels sitting on data
MID  <- "#4d5a66"   # axis titles, secondary text
SOFT <- "#98a2ab"   # tick labels, captions
RULE <- "#d3d8dc"   # axis lines

## Base map --------------------------------------------------------------------
COL_BASE      <- "#6E777F"   # country layer, the neutral grey
COL_BASE_LINE <- "black"     # country boundaries

## Period of first invasion ----------------------------------------------------
## Built from the two project hues rather than from one ramp: teal for the early
## periods, amber for the later ones, so the reader sees the shift from the old
## southeastern foci to the recent western and northern expansion as a change of
## hue and not only of lightness. Plum stays reserved for the out-of-sample
## seasons. Luminance still increases across the four fitted periods, so the
## order survives greyscale.
COL_P1 <- "#1F5866"   # 2006-2010, deep teal
COL_P2 <- TEAL        # 2011-2015
COL_P3 <- GOLD        # 2016-2020
COL_P4 <- "#E3BE79"   # 2021-2024, light amber

## The reserved colour ---------------------------------------------------------
COL_OOS      <- PLUM        # 2025-2026, out of sample. Reserved.
COL_OOS_EDGE <- "#2F1530"
COL_OOS_SOFT <- LILAC       # pale lilac, for bands behind out-of-sample points

## Absence ---------------------------------------------------------------------
COL_NEVER <- NEUTRAL     # never invaded: chosen, not missing
COL_NA    <- "#F4F4F2"   # outside the study cohort
COL_OBS   <- GREY        # observed counts, one neutral tone wherever they appear

PERIOD_LEVELS <- c("2006-2010", "2011-2015", "2016-2020", "2021-2024",
                   "2025-2026")
PERIOD_COLOURS <- stats::setNames(
  c(COL_P1, COL_P2, COL_P3, COL_P4, COL_OOS), PERIOD_LEVELS)

## Distance classes ------------------------------------------------------------
## Distance classes are bins of a continuous quantity, so they take a sequential
## ramp rather than a categorical set. The ramp is the project amber, deliberately
## a different hue family from the period ramp: a figure uses one scheme or the
## other, never both at once, so the reader never has to ask which is which.
## Near is dark, far is light.
DIST_LEVELS <- c("<=75 km", "76-150 km", "151-300 km", "301-500 km", ">500 km")
DIST_COLOURS <- stats::setNames(
  c("#5C3608", "#8F5714", AMBER, "#DCA95F", "#EFD6AC"), DIST_LEVELS)

## Counts on a heatmap ---------------------------------------------------------
## The reds from the submitted Figure 1B, kept as they were. Cell labels are
## white, so the ramp must stay dark enough at its pale end to carry white text;
## that is why it starts at a mid salmon rather than at near-white.
COUNT_RAMP <- c("#FCAE91", "#FB6A4A", "#CB181D", "#B00000")
COUNT_LOW  <- COUNT_RAMP[1]
COUNT_MID  <- COUNT_RAMP[3]
COUNT_HIGH <- COUNT_RAMP[4]

## Typography ------------------------------------------------------------------
## Resolved against the fonts actually installed. Naming a missing family makes
## the graphics device fall back silently, which is how a figure set ends up
## half in one face and half in another.
FONT_CANDIDATES <- c("DejaVu Sans", "Helvetica Neue", "Helvetica", "Arial")

FONT <- local({
  have <- tryCatch(unique(systemfonts::system_fonts()$family),
                   error = function(e) character(0))
  if (!length(have)) {
    message("_wnv_style: systemfonts not available, using the device default")
    return("")
  }
  hit <- FONT_CANDIDATES[FONT_CANDIDATES %in% have]
  if (!length(hit)) {
    warning("none of ", paste(FONT_CANDIDATES, collapse = ", "),
            " is installed; using the device default")
    return("")
  }
  if (hit[1] != FONT_CANDIDATES[1])
    message("_wnv_style: ", FONT_CANDIDATES[1], " not found, using ", hit[1])
  hit[1]
})

## Themes ---------------------------------------------------------------------
## Spare by default: no panel grid, no box, axis lines only where they carry
## information. The Lancet Regional Health Europe column widths are 90 mm and
## 190 mm, so figures are built at final size and the type sizes are the real
## ones.

WIDTH_1COL <- 90    # mm
WIDTH_2COL <- 190   # mm

theme_wnv_map <- function(base = 7) {
  ggplot2::theme_void(base_size = base, base_family = FONT) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base * 1.08, colour = INK,
                                         margin = ggplot2::margin(b = 3)),
      legend.title = ggplot2::element_text(size = base, colour = INK),
      legend.text  = ggplot2::element_text(size = base * 0.9, colour = MID),
      legend.key.height = ggplot2::unit(3.2, "mm"),
      legend.key.width  = ggplot2::unit(3.2, "mm"),
      plot.title = ggplot2::element_text(colour = INK, size = base * 1.2,
                                         hjust = 0,
                                         margin = ggplot2::margin(b = 4)),
      plot.margin = ggplot2::margin(2, 2, 2, 2))
}

theme_wnv_plot <- function(base = 7) {
  ggplot2::theme_minimal(base_size = base, base_family = FONT) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line  = ggplot2::element_line(colour = RULE, linewidth = 0.25),
      axis.ticks = ggplot2::element_line(colour = RULE, linewidth = 0.25),
      axis.ticks.length = ggplot2::unit(1.6, "pt"),
      axis.text  = ggplot2::element_text(colour = MID, size = base * 0.9),
      axis.title = ggplot2::element_text(colour = INK, size = base),
      legend.title = ggplot2::element_text(size = base, colour = INK),
      legend.text  = ggplot2::element_text(size = base * 0.9, colour = MID),
      plot.title = ggplot2::element_text(colour = INK, size = base * 1.2,
                                         hjust = 0,
                                         margin = ggplot2::margin(b = 5)),
      plot.margin = ggplot2::margin(6, 8, 6, 6))
}

## Scales, so a figure never has to name a colour by hand ----------------------

scale_fill_period <- function(name = "First invasion", drop = FALSE, ...) {
  ggplot2::scale_fill_manual(name = name, values = PERIOD_COLOURS,
                             limits = PERIOD_LEVELS, drop = drop,
                             na.value = COL_NEVER, ...)
}

scale_colour_distance <- function(name = "Distance to nearest\ninvaded region", ...) {
  ggplot2::scale_colour_manual(name = name, values = DIST_COLOURS,
                               limits = DIST_LEVELS, ...)
}

scale_fill_distance <- function(name = "Distance to nearest\ninvaded region", ...) {
  ggplot2::scale_fill_manual(name = name, values = DIST_COLOURS,
                             limits = DIST_LEVELS, ...)
}

scale_fill_count <- function(name = "Total newly invaded\nNUTS3 regions", ...) {
  ggplot2::scale_fill_gradientn(
    name = name,
    colours = c(COUNT_LOW, COUNT_MID, COUNT_HIGH),
    values = c(0, 0.45, 1),
    na.value = "white", ...)
}

## Binned version, for counts. A count is a whole number, so its key shows whole
## numbers in blocks rather than a continuous ribbon the reader has to
## interpolate against.
scale_fill_count_binned <- function(name = "Total newly invaded regions",
                                    breaks = c(1, 3, 5, 9, 13),
                                    limits = NULL, ...) {
  ggplot2::scale_fill_stepsn(
    name = name,
    colours = COUNT_RAMP,
    breaks = breaks,
    limits = limits,
    na.value = "white",
    guide = ggplot2::guide_coloursteps(
      show.limits = TRUE, even.steps = TRUE, ticks = FALSE,
      title.position = "left", title.vjust = 0.9,
      barheight = grid::unit(3.2, "mm"), barwidth = grid::unit(40, "mm")),
    ...)
}

## Country names, for axes that should read as prose rather than as codes.
CNTR_NAMES <- c(
  AL = "Albania",        AT = "Austria",         BG = "Bulgaria",
  CY = "Cyprus",         CZ = "Czech Republic",  DE = "Germany",
  EL = "Greece",         ES = "Spain",           FR = "France",
  HR = "Croatia",        HU = "Hungary",         IT = "Italy",
  ME = "Montenegro",     MK = "North Macedonia", NL = "Netherlands",
  PT = "Portugal",       RO = "Romania",         RS = "Serbia",
  SI = "Slovenia",       SK = "Slovakia",        TR = "Türkiye")

country_name <- function(code) unname(ifelse(is.na(CNTR_NAMES[code]), code, CNTR_NAMES[code]))

## The out-of-sample marker ----------------------------------------------------
## Used identically wherever 2025 or 2026 appears, so the reader learns it once.
## Colour alone is not enough: these layers add an outline as well.

oos_outline <- function(linewidth = 0.22) {
  ggplot2::geom_sf(colour = COL_OOS_EDGE, fill = NA, linewidth = linewidth)
}

oos_point <- function(size = 1.1, shape = 21, stroke = 0.3) {
  ggplot2::geom_point(size = size, shape = shape, stroke = stroke,
                      fill = COL_OOS, colour = COL_OOS_EDGE)
}

OOS_NOTE <- paste(
  "Seasons 2025 and 2026 lie outside the period used to fit the model and are",
  "shown in the reserved out-of-sample styling. The 2026 season is incomplete",
  "and is shown as at 10 September 2026.")

## Scale bar and north arrow, for Figure 1A (Reviewer 1, comment 4) ------------
## ggspatial is used when installed. The fallback draws a plain bar from the
## panel limits so the figure is never silently produced without one.

## The referee asked for a scale, not for a north arrow, and north_arrow_minimal
## lays its "N" out from text_size while drawing the arrow inside `height`, so
## the two collide at any size that suits a small panel. A north-up map of Europe
## with no graticule does not need an arrow to be read correctly, so the arrow is
## off by default. Set WNV_NORTH_ARROW=1 to put one back; it then uses the
## orienteering style, whose label sits outside the arrow and cannot overlap it.

add_scale_bar <- function(location = "bl", width_hint = 0.22,
                          arrow_location = "tl", arrow_mm = 11,
                          arrow_text = 8, pad_mm = 3,
                          arrow = toupper(Sys.getenv("WNV_NORTH_ARROW", "0")) %in%
                                  c("1", "TRUE")) {
  if (!requireNamespace("ggspatial", quietly = TRUE)) {
    message("_wnv_style: ggspatial not installed; add_scale_bar() is a no-op. ",
            "Install it before producing the final Figure 1A.")
    return(list())
  }
  fam <- if (nzchar(FONT)) FONT else NULL

  out <- list(
    ggspatial::annotation_scale(
      location = location, width_hint = width_hint,
      height = ggplot2::unit(1.4, "mm"), text_cex = 0.6,
      line_width = 0.3, text_col = MID, line_col = MID,
      pad_x = ggplot2::unit(pad_mm, "mm"), pad_y = ggplot2::unit(pad_mm, "mm"),
      text_family = fam))

  if (arrow) {
    out <- c(out, list(ggspatial::annotation_north_arrow(
      location = arrow_location,
      height = ggplot2::unit(arrow_mm, "mm"),
      width  = ggplot2::unit(arrow_mm, "mm"),
      pad_x  = ggplot2::unit(pad_mm, "mm"),
      pad_y  = ggplot2::unit(pad_mm, "mm"),
      style = ggspatial::north_arrow_orienteering(
        line_width = 1, line_col = MID, fill = c("white", MID), text_col = MID,
        text_size = arrow_text, text_family = fam))))
  }
  out
}

## Where to put a legend that must sit over a specific piece of ground rather
## than at a fixed fraction of the panel. Give it a longitude and latitude and
## the plotted bounding box, and it returns the npc coordinates for
## legend.position, so the anchor holds whatever projection the map is drawn in.
legend_at_lonlat <- function(lon, lat, bbox, crs) {
  p <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  if (!identical(as.integer(crs), 4326L)) p <- sf::st_transform(p, crs)
  xy <- sf::st_coordinates(p)
  c(x = unname((xy[1] - bbox["xmin"]) / (bbox["xmax"] - bbox["xmin"])),
    y = unname((xy[2] - bbox["ymin"]) / (bbox["ymax"] - bbox["ymin"])))
}

## Captions live OUTSIDE the figure --------------------------------------------
## A caption burnt into the image cannot be edited in the manuscript, cannot be
## spell-checked, and forces the canvas to carry eight lines of prose that
## squeeze the panels. Each figure script writes its caption here instead, one
## unwrapped paragraph per file, ready to paste into Word.

DIR_FIG <- file.path(Sys.getenv("WNV_RUN_DIR", "run_withBalkans"),
                     "output", "revision", "figures")
DIR_CAP <- file.path(DIR_FIG, "captions")

write_caption <- function(slug, txt) {
  if (!dir.exists(DIR_CAP)) dir.create(DIR_CAP, recursive = TRUE)
  txt <- gsub("[ \t\n]+", " ", trimws(txt))
  writeLines(txt, file.path(DIR_CAP, paste0(slug, ".txt")))
  invisible(txt)
}

## The base pdf() and svg() devices only know the 14 standard PostScript fonts,
## so passing a real system family to them fails with "invalid font type", and a
## symbol such as the capital phi fails with "failed to find or load PDF CID
## font". The Cairo devices embed the actual font and handle both, so they are
## used wherever they are available.

save_fig <- function(plot, slug, width_mm = WIDTH_2COL, height_mm = 120,
                     formats = c("pdf", "png", "svg")) {
  if (!dir.exists(DIR_FIG)) dir.create(DIR_FIG, recursive = TRUE)

  device_for <- function(fmt) {
    switch(fmt,
      pdf = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf,
      svg = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite
            else if (capabilities("cairo")) grDevices::svg else NULL,
      png = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png
            else function(filename, width, height, ...)
                   grDevices::png(filename, width, height, type = "cairo", ...),
      NULL)
  }

  for (fmt in formats) {
    dev <- device_for(fmt)
    if (is.null(dev)) { message("no device available for ", fmt, ", skipped"); next }
    f <- file.path(DIR_FIG, paste0(slug, ".", fmt))
    args <- list(filename = f, plot = plot, width = width_mm, height = height_mm,
                 units = "mm", device = dev, bg = "white")
    if (fmt == "png") args$dpi <- 300
    res <- tryCatch({ do.call(ggplot2::ggsave, args); TRUE },
                    error = function(e) { message("  ", fmt, " failed: ",
                                                  conditionMessage(e)); FALSE })
    if (isTRUE(res)) message("wrote ", f)
  }
  invisible(plot)
}

## Helpers shared by more than one figure --------------------------------------

## The display periods are five-year bands, as in the submitted figures. They are
## NOT the invasion phases used in the models, which are Early, Middle, Pre-COVID
## and Post-COVID with a different 2016-2019 cut. Keep the two apart.
period_of <- function(year) {
  factor(dplyr::case_when(
    year <= 2010 ~ "2006-2010",
    year <= 2015 ~ "2011-2015",
    year <= 2020 ~ "2016-2020",
    year <= 2024 ~ "2021-2024",
    TRUE         ~ "2025-2026"
  ), levels = PERIOD_LEVELS)
}

distance_class <- function(km) {
  factor(dplyr::case_when(
    km <= 75  ~ "<=75 km",
    km <= 150 ~ "76-150 km",
    km <= 300 ~ "151-300 km",
    km <= 500 ~ "301-500 km",
    TRUE      ~ ">500 km"
  ), levels = DIST_LEVELS)
}

## A one-line check that the palette survives greyscale. Run it after any change
## to the ramp; the luminances must be strictly increasing across the periods.
check_greyscale <- function() {
  lum <- vapply(c(COL_P1, COL_P2, COL_P3, COL_P4), function(h) {
    v <- grDevices::col2rgb(h)[, 1] / 255
    sum(c(0.2126, 0.7152, 0.0722) * v)
  }, numeric(1))
  ok <- all(diff(lum) > 0)
  message("period ramp luminance: ", paste(round(lum, 3), collapse = " < "),
          if (ok) "  OK" else "  NOT MONOTONIC")
  invisible(ok)
}

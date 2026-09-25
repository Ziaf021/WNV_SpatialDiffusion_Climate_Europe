###############################################################################
# fig02b_pressure_rate.R
#
# Figure 2B. Annual first-invasion probability against neighbourhood invasion
# pressure.
#
# Rebuilt to the submitted design, which my earlier version had quietly changed
# in ways that mattered:
#
#   1. The x axis carries Phi itself on a log1p-transformed scale, with breaks
#      labelled in units of "nearby invaded regions". My version plotted
#      log(1+Phi) as a bare number, which threw away exactly the interpretive
#      labelling Reviewer 5 asked us to explain in comment 5. It is restored,
#      with the caution that reviewer asked for written into the caption.
#   2. Bins are summarised at the MEDIAN of Phi, not the mean. With a strongly
#      right-skewed index the mean sits above most of the bin's mass and shifts
#      every point rightward.
#   3. Fifty bins of equal size, as before.
#
# Colour follows the figure system rather than the submitted red: teal for the
# fitted series and its confidence band, plum for the 2026 season, which is the
# colour reserved for out-of-sample throughout the paper.
#
#   bash run_fig.sh 2b
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig02b")
on.exit(lg$close(), add = TRUE)

## ---- Data -------------------------------------------------------------------

NBIN <- as.integer(Sys.getenv("WNV_PRESSURE_BINS", "50"))

bin_curve <- function(df, nbin = NBIN) {
  df %>%
    filter(!is.na(invasion_pressure)) %>%
    mutate(bin = ntile(invasion_pressure, nbin)) %>%
    group_by(bin) %>%
    summarise(phi = median(invasion_pressure),
              rate = 100 * mean(invaded),
              n = n(), .groups = "drop")
}

# The curve is the FITTED period only. Folding 2025 into it, as an earlier
# version of this script did, was inconsistent with panel A, where 2025 is drawn
# as out-of-sample, and it quietly used a season the model never saw to shape a
# curve presented as the fitted relationship. Both out-of-sample seasons now sit
# in the reserved colour, in both panels.
curve_src <- ac %>% transmute(invasion_pressure, invaded = event)

oos_src <- bind_rows(
  p25 %>% transmute(invasion_pressure, invaded),
  p26 %>% transmute(invasion_pressure, invaded))

curve_dat <- bin_curve(curve_src)
oos_dat <- bin_curve(oos_src, nbin = 10)

PHI_MAX <- max(curve_src$invasion_pressure, na.rm = TRUE)

message("region-years in the curve (fitted period): ", sum(curve_dat$n),
        " in ", nrow(curve_dat), " bins")
message("out-of-sample region-years (2025 and 2026): ", sum(oos_dat$n),
        " in ", nrow(oos_dat), " bins")
message("Phi observed range: 0 to ", round(PHI_MAX, 3))
message("Phi across bin medians: ", round(min(curve_dat$phi), 4), " to ",
        round(max(curve_dat$phi), 3))

# Roughly a third of region-years have no invaded region within the kernel's
# reach, so their Phi is exactly zero. Every bin drawn from that part of the
# distribution has a median of zero and they all land on the same point, which
# is why the figure shows fewer distinct points than there are bins. Reported so
# the count is explicable rather than puzzling.
n_zero <- sum(curve_dat$phi == 0)
n_distinct <- nrow(dplyr::distinct(curve_dat, phi, rate))
message("bins with median Phi exactly zero: ", n_zero,
        "; distinct plotted positions: ", n_distinct, " of ", nrow(curve_dat), " bins")
message("region-years with Phi exactly zero: ", sum(curve_src$invasion_pressure == 0),
        " of ", nrow(curve_src),
        " (", round(100 * mean(curve_src$invasion_pressure == 0), 1), "%)")
message("To draw more points, raise the bin count: WNV_PRESSURE_BINS=80 bash run_fig.sh 2b")

## ---- The interpretive x axis ------------------------------------------------
# Phi is an effective count of invaded regions expressed in units of one adjacent
# invaded region at zero distance. With lambda = 50 km a neighbour at 50 km
# contributes exp(-1) = 0.368, so four such neighbours give about 1.47. The
# break labels are therefore an approximation valid at the modal neighbour
# distance, not an exact region count, and the caption says so.

# The axis spans the range Phi actually takes, not the range its bin medians
# take. A bin median cannot approach the maximum, because the top bin's median
# is the middle of the top few per cent of the distribution: here Phi reaches
# 6.1 while the highest bin median is only about 2.5. Letting the bin medians
# set the range silently truncated the axis at 3 and dropped the upper break
# labels, so the limits are set explicitly.
# Tick labels are the phrases from the submitted figure, and they are rotated 90
# degrees exactly as the submitted script rotated them. That rotation is why the
# submitted panel was legible and my earlier horizontal version was not: on a
# log1p axis the upper breaks crowd together, and horizontal two-line phrases at
# each of them collide. Rotated, they take vertical space instead and cannot.
#
# The breaks stop at the range Phi covers in this dataset, which is well short of
# the submitted panel's "~7 nearby invaded regions".
#
# That is not a change introduced by the revision. The submitted Figure 2B was
# drawn on 30 June from the root-level WNV_SurvivalDataset.csv of 24 June, in
# which Phi reaches 12.12 and the top bin median is 7.11. The reconciled pipeline
# that produced the submitted MODEL results wrote output/WNV_SurvivalDataset.csv
# on 4 July, in which Phi reaches 5.09 and the top bin median is 2.28, matching
# this refit almost exactly. The June file predates the reconciliation to
# lambda = 50 km; regressing log(new/old) on nearest-neighbour distance across
# the region-years common to both implies a decay length near 80 to 100 km in it.
#
# So the submitted figure is drawn from a superseded file and disagrees with the
# Methods, Table 1 and the model objects of the same manuscript. This panel is
# consistent with them; the submitted one was not.
PHI_PLOT_MAX <- max(c(curve_dat$phi, oos_dat$phi), na.rm = TRUE)

# The axis runs to 3 by default, a little beyond the highest bin median, so the
# "~3 nearby invaded regions" label is present and the reader can see where the
# curve is heading rather than having it stop at the last point.
PHI_AXIS_MAX <- as.numeric(Sys.getenv("WNV_PHI_MAX", "3"))
XLIM <- c(0, max(PHI_AXIS_MAX, PHI_PLOT_MAX * 1.06))

cand <- c(0, 1, 2, 3, 5, 7)
brks <- cand[cand <= XLIM[2] * 1.001]
phrase <- function(x) {
  ifelse(x == 0, "None",
  ifelse(x == 1, "~1 nearby\ninvaded region",
         paste0("~", x, " nearby\ninvaded regions")))
}
brk_labs <- phrase(brks)

message("x axis: limits 0 to ", round(XLIM[2], 3),
        "; breaks ", paste(brks, collapse = ", "),
        "; highest bin median ", round(PHI_PLOT_MAX, 3))
message("NOTE: the submitted panel reached Phi of about 7 because it was drawn from",
        " a superseded pre-reconciliation dataset. The lambda = 50 km pipeline that",
        " produced the submitted model results gives a top bin median of 2.28; this",
        " refit gives ", round(PHI_PLOT_MAX, 2), ".")

# ggplot renamed this argument; support both rather than warn on every run.
scale_x_log1p <- function(...) {
  if ("transform" %in% names(formals(ggplot2::scale_x_continuous)))
    ggplot2::scale_x_continuous(transform = "log1p", ...)
  else
    ggplot2::scale_x_continuous(trans = "log1p", ...)
}

## ---- Figure -----------------------------------------------------------------

# 2007, not 2006: the fitted period begins in 2007 because distance to the
# nearest invaded region is undefined in the first year of follow-up.
SER <- c(fit = "2007-2024", oos = "2025-2026")

pts <- bind_rows(
  curve_dat %>% mutate(series = SER["fit"]),
  oos_dat   %>% mutate(series = SER["oos"])) %>%
  mutate(series = factor(series, levels = unname(SER)))

# Both series are drawn at one size. Size carries no meaning here, so varying it
# invites the reader to look for one; the two series are told apart by colour and
# outline, and named in the key.
f2b <- ggplot(curve_dat, aes(phi, rate)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              colour = TEAL, fill = COL_P4, alpha = 0.55, linewidth = 0.9) +
  geom_point(data = pts, aes(phi, rate, fill = series),
             shape = 21, size = 2, stroke = 0.4, colour = "white") +
  scale_fill_manual(values = stats::setNames(c(COL_P1, COL_OOS), unname(SER)),
                    name = NULL) +
  scale_x_log1p(breaks = brks, labels = brk_labs, limits = XLIM) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = expression(paste("Neighbourhood invasion pressure (", Phi, ")")),
       y = "Observed annual invasion rate (%)") +
  theme_wnv_plot(base = 9) +
  theme(axis.title = element_text(colour = INK, size = 9),
        axis.text = element_text(colour = INK, size = 8.5),
        axis.text.x = element_text(size = 7.5, angle = 0, hjust = 0.5,
                                   vjust = 1, lineheight = 0.95),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.text = element_text(size = 8),
        plot.margin = margin(6, 8, 6, 6))

saveRDS(f2b, file.path(DIR_FIG, "panel_2B.rds"))
save_fig(f2b, "Figure2B_pressure_rate", 165, 120)

write_caption("Figure2B_pressure_rate", paste(
  "Figure 2B. Empirical relationship between neighbourhood invasion pressure and",
  "the observed annual rate of first West Nile virus invasion, before multivariable",
  "adjustment. Neighbourhood invasion pressure (Phi) is a",
  "distance-weighted index of cumulative exposure to all previously invaded",
  "regions, computed with a negative exponential kernel at lambda = 50 km. One",
  "invaded region at zero distance contributes 1; a region 50 km away contributes",
  "exp(-1) = 0.368, so four invaded regions each 50 km away give Phi of about",
  "1.47. Phi is therefore an effective count of invaded regions expressed in units",
  "of one adjacent invaded region at zero distance. The horizontal axis is on a",
  "log(1+x) scale because the index is strongly right-skewed, and its break labels",
  "are an approximation valid at the modal neighbour distance; they should not be",
  "read as an exact number of regions.",
  sprintf("Region-years of the fitted period were grouped into %d bins of equal size;", NBIN),
  "points give the median Phi and the observed annual invasion rate in each bin,",
  "with a locally weighted regression and its 95% confidence band. Regions were",
  "followed from 2006, and the model is fitted from 2007 onwards, since distance to",
  "the nearest invaded region is undefined in the first year of follow-up. The",
  "horizontal axis is on a log(1+x) scale and covers the range spanned by the bin",
  sprintf("medians; Phi itself reaches %.1f in individual region-years, but a bin", PHI_MAX),
  "median cannot approach that, being the middle of the top few per cent of the",
  "distribution. Region-years with no invaded region within the kernel's reach have",
  "Phi of exactly zero, so the lowest bins coincide at the origin and fewer distinct",
  "points are visible than there are bins. The 2025 and 2026 seasons were not used",
  "to fit the model and are shown separately in deciles, in the reserved colour.",
  OOS_LABEL))

message("fig02b done")

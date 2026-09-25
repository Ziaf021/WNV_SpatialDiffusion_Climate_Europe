###############################################################################
# fig02_combined.R
#
# Figure 2, panels A and B side by side, at the proportions of the submitted
# figure (10 by 5 inches, two columns).
#
#   bash run_fig.sh 2a     then
#   bash run_fig.sh 2b     then
#   bash run_fig.sh 2
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig02")
on.exit(lg$close(), add = TRUE)

pa <- file.path(DIR_FIG, "panel_2A.rds")
pb <- file.path(DIR_FIG, "panel_2B.rds")

for (f in c(pa, pb)) {
  if (!file.exists(f)) {
    stop("Not found: ", f,
         "\nBuild the panels first:  bash run_fig.sh 2a   and   bash run_fig.sh 2b",
         call. = FALSE)
  }
}

fig2 <- (readRDS(pa) | readRDS(pb)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 12, face = "bold", colour = INK),
                  plot.margin = margin(4, 4, 4, 4)))

# 10 by 5 inches, as submitted.
save_fig(fig2, "Figure2_Main_AB", 254, 127)

write_caption("Figure2_Main_AB", paste(
  "Figure 2. Spatial scale of West Nile virus invasion risk across Europe.",
  "(A) Distribution of first invasions by the distance between the newly invaded",
  "region and the nearest region already invaded at the start of that season. Bars",
  "give counts; the label within each bar is the total for that class with its",
  "percentage of all 381 first invasions with a measurable distance, and above each",
  "bar the total is split into the fitted period, 2007 to 2024, and, in the reserved",
  "colour, the 2025 and 2026 seasons, which were not used to fit the model.",
  "The distribution is strongly concentrated within 150 km of an existing focus.",
  "(B) Empirical relationship between neighbourhood invasion pressure (Phi) and the",
  "observed annual rate of first invasion, before multivariable adjustment. Phi is a",
  "distance-weighted index of cumulative exposure to all previously invaded regions,",
  "computed with a negative exponential kernel at lambda = 50 km: one invaded region",
  "at zero distance contributes 1, a region 50 km away contributes exp(-1) = 0.368,",
  "and four invaded regions each 50 km away give Phi of about 1.47. Phi is therefore",
  "an effective count of invaded regions expressed in units of one adjacent invaded",
  "region at zero distance. The horizontal axis is on a log(1+x) scale because the",
  "index is strongly right-skewed, and its break labels are an approximation valid",
  "at the modal neighbour distance rather than an exact region count. Region-years",
  "from 2007 to 2025 were grouped into 50 bins of equal size; points give the median",
  "Phi and the observed annual invasion rate in each bin, with a locally weighted",
  "regression and its 95% confidence band. The 2026 season is shown separately in",
  "deciles, in the colour reserved for out-of-sample data, because a partial season",
  "biases an annual rate downward. Invasion probability stays low across the lower",
  "range of Phi and rises sharply once exposure accumulates beyond roughly one to",
  "three nearby invaded regions.", OOS_LABEL))

message("fig02 combined done")

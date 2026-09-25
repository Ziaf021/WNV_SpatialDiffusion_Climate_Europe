###############################################################################
# fig03_forest_spline_km.R
#
# Figure 3. Three panels.
#
#   A  multivariable forest plot with the hazard-ratio and p columns beside it
#   B  distance spline
#   C  Kaplan-Meier by distance class, with the numbers at risk beneath it
#
# The numbers at risk are part of panel C. They are drawn as a separate ggplot
# because a risk table has to align to the survival curve's x axis and cannot
# live inside it, but they carry no tag and the caption describes three panels.
#
#   bash run_fig.sh 3
#
# Revised 10 September 2026, third pass. What changed against the second pass:
#
#   - panel C keeps the time-varying stratification of the submitted figure, by
#     the author's decision: a region is classified each year by its distance to
#     the nearest region already invaded, so the panel measures exposure to the
#     moving front, which is what the paper is about. The consequence is stated
#     rather than hidden. A region contributes follow-up to more than one class,
#     the curves summarise class-specific annual invasion rates rather than the
#     experience of fixed groups, and the counts beneath rise as well as fall.
#     WNV_KM_CLASS=baseline fixes each region on its 2007 class instead, which
#     gives proper Kaplan-Meier curves of fixed groups at the cost of a nearest
#     class holding only 7 regions
#   - the axis is drawn in calendar years. The follow-up clock runs start = 1,
#     stop = 2 for 2007 through start = 18, stop = 19 for 2024, so an axis of 0
#     to 18 cut off the whole of 2024, the largest season at 56 first invasions
#   - panel A no longer uses amber. Amber is the middle step of the distance
#     ramp in panel C, so the climate domain and the 151-300 km class were the
#     same colour in one figure. Spatial diffusion now takes the red of the
#     panel B spline, which also ties the distance row in A to the curve in B
#   - the year row is back under the risk table, with the axis title there
#   - axis labels in C reworded, see the note above panel C
#   - urbanization spelled to match the manuscript text
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))
suppressPackageStartupMessages(library(splines))

lg <- fig_log("fig03")
on.exit(lg$close(), add = TRUE)

## ---- Colour ------------------------------------------------------------------
## DIST_LINE. DIST_COLOURS ends at #EFD6AC, chosen for fills. A 0.65 pt line in
## that colour on white is close to invisible in print, so the ramp is compressed
## at the pale end and stays monotone in lightness. WNV_DIST_RAMP=style falls
## back to the shared ramp exactly.

DIST_LINE <- if (identical(toupper(Sys.getenv("WNV_DIST_RAMP", "")), "STYLE")) {
  DIST_COLOURS
} else {
  stats::setNames(c("#4A2B06", "#8A5413", "#C77B21", "#D2994A", "#C7AE82"),
                  DIST_LEVELS)
}

darken <- function(hex, k = 0.72) {
  m <- grDevices::col2rgb(hex)
  grDevices::rgb(t(m * k), maxColorValue = 255)
}
DIST_TEXT <- stats::setNames(darken(DIST_LINE), DIST_LEVELS)

## Panel A domains. Three hues, none of them from the amber family, because the
## amber family is the distance ramp in panel C and a reader should never have to
## ask whether two identical colours in one figure mean the same thing. Red is
## the colour of the spline in panel B, so the distance row in A and the curve in
## B now read as the same quantity. Plum stays unused, reserved for the
## out-of-sample seasons.
DOM_RED   <- "#CB181D"
DOM_BLUE  <- "#2171B5"
DOM_GREEN <- "#238B45"
dom_cols <- c("Spatial diffusion" = DOM_RED,
              "Climate"           = DOM_BLUE,
              "Urbanization"      = DOM_GREEN)

## Panel B, the colours of the submitted figure.
SPLINE_LINE   <- DOM_RED
SPLINE_RIBBON <- "#FCAE91"

DIST_PRETTY <- stats::setNames(
  c("≤75 km", "76-150 km", "151-300 km", "301-500 km", ">500 km"),
  DIST_LEVELS)

## ---- Panel A, forest ---------------------------------------------------------

s <- summary(fit_primary)
labs_map <- c(
  distance100          = "Distance to nearest invaded region\n(per 100 km)",
  pressure_z           = "Neighbourhood invasion pressure\n(per SD)",
  log_cases_nearest    = "Source-region transmission intensity\n(per unit log(cases+1))",
  tas                  = "Annual mean temperature\n(per 1 °C)",
  hurs                 = "Relative humidity\n(per 1%)",
  pr_z                 = "Annual precipitation\n(per SD)",
  `factor(URBN_TYPE)2` = "Intermediate urbanization\n(vs urban)",
  `factor(URBN_TYPE)3` = "Rural\n(vs urban)")
dom_map <- c(distance100 = "Spatial diffusion", pressure_z = "Spatial diffusion",
             log_cases_nearest = "Spatial diffusion", tas = "Climate",
             hurs = "Climate", pr_z = "Climate",
             `factor(URBN_TYPE)2` = "Urbanization", `factor(URBN_TYPE)3` = "Urbanization")

# A term the label map does not know is dropped by the filter below without a
# word, which is how the precipitation row vanished from this panel on 10
# September when the figure was built against a model fitted before pr_z
# replaced scale(pr). Fail loudly instead.
model_terms <- rownames(summary(fit_primary)$conf.int)
missing_lab <- setdiff(model_terms, names(labs_map))
if (length(missing_lab)) {
  stop("The fitted model has terms with no label: ",
       paste(missing_lab, collapse = ", "),
       "\nThe model in ", RUN, " is older than this script. Refit before drawing.",
       call. = FALSE)
}

fd <- data.frame(term = rownames(s$conf.int),
                 HR = s$conf.int[, "exp(coef)"],
                 lo = s$conf.int[, "lower .95"],
                 hi = s$conf.int[, "upper .95"],
                 p  = s$coefficients[, "Pr(>|z|)"], row.names = NULL) %>%
  filter(term %in% names(labs_map)) %>%
  mutate(Variable = factor(labs_map[term], levels = rev(labs_map[names(labs_map)])),
         Domain = factor(dom_map[term],
                         levels = c("Spatial diffusion", "Climate", "Urbanization")),
         lab = sprintf("%.2f (%.2f–%.2f)", HR, lo, hi),
         plab = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

print(fd %>% select(term, HR, lo, hi, plab), row.names = FALSE)

# Both plots share this y scale exactly, which is what makes the rows of the
# text column line up with the rows of the forest.
Y_EXPAND <- expansion(add = c(0.55, 0.95))

pA <- ggplot(fd, aes(HR, Variable, colour = Domain)) +
  geom_vline(xintercept = 1, linetype = 2, colour = RULE, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.6) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = dom_cols, name = NULL) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2),
                labels = c("0.5", "0.75", "1.0", "1.5", "2.0")) +
  scale_y_discrete(expand = Y_EXPAND) +
  coord_cartesian(xlim = c(0.45, 2.05)) +
  labs(x = "Hazard ratio", y = NULL, tag = "A",
       caption = "HR, hazard ratio; CI, confidence interval; SD, standard deviation") +
  theme_wnv_plot(base = 7) +
  theme(axis.text.y = element_text(colour = MID, size = 6.3, lineheight = 1.05),
        legend.position = "bottom",
        legend.margin = margin(t = 0, b = 0),
        legend.box.spacing = unit(2, "mm"),
        legend.key.width = unit(4, "mm"),
        plot.caption = element_text(size = 5.6, colour = MID, hjust = 0,
                                    margin = margin(t = 2)),
        plot.tag = element_text(size = 8, colour = INK),
        plot.tag.position = c(0.002, 0.99),
        plot.margin = margin(6, 2, 4, 4))

# The text column. A plot rather than annotation outside the panel, so nothing
# depends on a margin being wide enough.
pAtab <- ggplot(fd, aes(y = Variable)) +
  geom_text(aes(x = 0.02, label = lab), hjust = 0, size = 2.2, colour = INK) +
  geom_text(aes(x = 0.78, label = plab), hjust = 0, size = 2.2, colour = INK) +
  annotate("text", x = 0.02, y = nrow(fd) + 0.62, label = "HR (95% CI)",
           hjust = 0, size = 2.2, colour = MID) +
  annotate("text", x = 0.78, y = nrow(fd) + 0.62, label = "p value",
           hjust = 0, size = 2.2, colour = MID) +
  scale_x_continuous(limits = c(0, 1.25), expand = c(0, 0)) +
  scale_y_discrete(expand = Y_EXPAND) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 7) +
  theme(plot.margin = margin(6, 4, 4, 2))

rowA <- pA + pAtab + plot_layout(widths = c(1, 0.46))

## ---- Panel B, distance spline ------------------------------------------------

cox_spline <- coxph(Surv(start, stop, event) ~ ns(distance100, df = 4), data = ac)
first_inv <- ac %>% filter(event == 1)
upper <- as.numeric(quantile(first_inv$nearest_invaded_km, 0.99, na.rm = TRUE))
nd <- data.frame(distance100 = seq(min(ac$distance100, na.rm = TRUE), upper / 100,
                                   length.out = 800))
refl <- predict(cox_spline, newdata = data.frame(distance100 = 1), type = "lp")
pr <- predict(cox_spline, newdata = nd, type = "lp", se.fit = TRUE)
nd <- nd %>% mutate(km = distance100 * 100,
                    HR = exp(pr$fit - refl),
                    lo = exp(pr$fit - refl - 1.96 * pr$se.fit),
                    hi = exp(pr$fit - refl + 1.96 * pr$se.fit))

pB <- ggplot(nd, aes(km, HR)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = SPLINE_RIBBON, alpha = 0.35) +
  geom_line(colour = SPLINE_LINE, linewidth = 0.75) +
  geom_hline(yintercept = 1, linetype = 2, colour = RULE, linewidth = 0.3) +
  geom_rug(data = first_inv, aes(x = nearest_invaded_km), inherit.aes = FALSE,
           sides = "b", alpha = 0.3, colour = MID, length = unit(2.2, "pt")) +
  scale_x_continuous(breaks = seq(0, 600, 100),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(breaks = scales::pretty_breaks(5),
                     expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = c(0, 600)) +
  annotate("text", x = 235, y = 1.5, label = "Reference 100 km",
           size = 2.2, colour = MID, hjust = 0) +
  labs(x = "Distance to nearest invaded region (km)",
       y = "Relative invasion hazard", tag = "B") +
  theme_wnv_plot(base = 7) +
  theme(plot.tag = element_text(size = 8, colour = INK),
        plot.tag.position = c(0.002, 0.99),
        plot.margin = margin(6, 6, 4, 4))

## ---- Panel C, Kaplan-Meier with the numbers at risk beneath ------------------
##
## Stratification. A region is classified in every year at risk by its distance
## to the nearest region already invaded at the start of that year, as in the
## submitted figure. A region therefore moves between classes as the front
## approaches it, contributes follow-up to more than one class, and the counts
## beneath the curves rise as well as fall. These are not Kaplan-Meier estimates
## for fixed groups; they summarise the class-specific annual invasion rates, and
## the caption says so. WNV_KM_CLASS=baseline fixes each region on the class it
## held in 2007 instead, which gives fixed groups and a monotone risk table at
## the cost of a nearest class of 7 regions.
##
## Axes. The follow-up clock in analysis_cases runs start = 1, stop = 2 for 2007
## through start = 18, stop = 19 for 2024, so model time is the calendar year
## minus 2005. An axis running 0 to 18 therefore cuts off the whole of 2024, the
## largest season in the series at 56 first invasions, and shifts every label by
## a year. The axis is drawn in calendar years, which also settles what "year"
## means without a clause in the caption.

base_class <- identical(tolower(Sys.getenv("WNV_KM_CLASS", "")), "baseline")

if (!base_class) {
  message("panel C: distance class recomputed each region-year (time-varying)")
  acc <- ac %>% mutate(cls0 = distance_class(nearest_invaded_km))
  message("region-years per distance class:")
  print(as.data.frame(table(acc$cls0)), row.names = FALSE)
} else {
  message("panel C: distance class fixed at each region's first year at risk")
  base_cls <- ac %>%
    arrange(NUTS_ID, start) %>%
    group_by(NUTS_ID) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(NUTS_ID, cls0 = distance_class(nearest_invaded_km))
  message("regions per baseline distance class:")
  print(as.data.frame(table(base_cls$cls0)), row.names = FALSE)
  acc <- ac %>% left_join(base_cls, by = "NUTS_ID")
}

stopifnot(!any(is.na(acc$cls0)))

km <- survfit(Surv(start, stop, event) ~ cls0, data = acc)

T_ORIGIN  <- 2005            # model time t corresponds to calendar year t + 2005
T_TIMES   <- c(1, 4, 7, 10, 13, 16, 19)
YR_BREAKS <- T_TIMES + T_ORIGIN
YR_RANGE  <- c(2005.7, 2024.35)

kmd <- data.frame(t = km$time, surv = km$surv,
                  cls = rep(sub("^cls0=", "", names(km$strata)), km$strata)) %>%
  mutate(cls = factor(cls, levels = DIST_LEVELS), yr = t + T_ORIGIN)

# survfit returns event and censoring times only, so each curve is anchored at
# the start of follow-up to make it begin at 1 rather than at its first event.
kmd <- bind_rows(
  data.frame(t = 1, surv = 1, cls = factor(DIST_LEVELS, levels = DIST_LEVELS),
             yr = 1 + T_ORIGIN),
  kmd) %>% arrange(cls, yr)

stopifnot(max(kmd$yr) >= 2024)

pC <- ggplot(kmd, aes(yr, surv, colour = cls)) +
  geom_step(linewidth = 0.65) +
  scale_colour_manual(values = DIST_LINE, limits = DIST_LEVELS, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.03))) +
  scale_x_continuous(breaks = YR_BREAKS, limits = YR_RANGE, expand = c(0.01, 0)) +
  labs(x = "Year", y = "Probability a region\nremains uninvaded",
       subtitle = paste("Grouped by distance to the nearest region",
                        "already\ninvaded in that year"),
       tag = "C") +
  theme_wnv_plot(base = 7) +
  theme(plot.subtitle = element_text(size = 5.9, colour = MID,
                                     lineheight = 1.1,
                                     margin = margin(b = 3)),
        plot.tag = element_text(size = 8, colour = INK),
        plot.tag.position = c(0.002, 0.99),
        plot.margin = margin(6, 6, 2, 4))

sm <- summary(km, times = T_TIMES, extend = TRUE)
risk <- data.frame(cls = sub("^cls0=", "", as.character(sm$strata)),
                   t = sm$time, n = sm$n.risk) %>%
  mutate(cls = factor(cls, levels = rev(DIST_LEVELS)), yr = t + T_ORIGIN)

message("regions at risk by year, summed over classes:")
print(risk %>% group_by(yr) %>% summarise(n = sum(n), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

lev_bottom_up <- levels(risk$cls)

pR <- ggplot(risk, aes(yr, cls, label = n)) +
  geom_text(size = 2.0, colour = INK) +
  scale_x_continuous(breaks = YR_BREAKS, limits = YR_RANGE, expand = c(0.01, 0)) +
  scale_y_discrete(labels = DIST_PRETTY[lev_bottom_up],
                   expand = expansion(add = 0.6)) +
  labs(x = NULL, y = NULL,
       title = "Regions at risk in each distance class, recomputed each year") +
  theme_wnv_plot(base = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = DIST_TEXT[lev_bottom_up], size = 5.9,
                                   hjust = 1),
        plot.title = element_text(size = 6.2, colour = MID,
                                  margin = margin(t = 1, b = 2)),
        plot.margin = margin(1, 6, 4, 4))

panelC <- pC / pR + plot_layout(heights = c(1, 0.5))

rowB <- pB + panelC + plot_layout(widths = c(1, 1.12))

## ---- Assembly ----------------------------------------------------------------
## Each row is wrapped so that patchwork aligns panels within a row and not
## across the whole figure. Without this, panel A's label column sets the left
## offset for panel B as well.

f3 <- wrap_elements(rowA) / wrap_elements(rowB) + plot_layout(heights = c(1, 0.98))

save_fig(f3, "Figure3_forest_spline_km", WIDTH_2COL, 188)

write_caption("Figure3_forest_spline_km", paste(
  "Figure 3. Determinants of West Nile virus invasion risk and the distance decay of",
  "invasion hazard, NUTS3 regions, 2007 to 2024. (A) Hazard ratios and 95% confidence",
  "intervals from the multivariable Cox proportional hazards model, grouped by domain;",
  "points left of 1 indicate a lower hazard of first invasion. HR, hazard ratio; CI,",
  "confidence interval; SD, standard deviation. The source-region term is",
  "the log of the case count in the nearest previously invaded region in the preceding",
  "season, so its hazard ratio is per unit increase on that log scale, which corresponds",
  "to roughly a 32% higher hazard for each doubling of that count. (B) Association between",
  "distance to the nearest previously invaded region and relative invasion hazard, from a",
  "model with a natural cubic spline on distance (4 degrees of freedom), expressed relative",
  "to 100 km; the band is the pointwise 95% confidence interval and the rug marks observed",
  "first invasions. (C) Survival curves for remaining free of a first invasion across the",
  "fitted period 2007 to 2024. A region is assigned to a distance class in each year at risk",
  "according to how far it lay from the nearest region already invaded at the start of that",
  "year, so a region moves between classes as the invasion front approaches it. The curves",
  "therefore summarise the class-specific annual invasion rates and are not survival",
  "estimates for fixed groups of regions. The number of regions in each class in each year is",
  "given beneath the curves; these counts rise as well as fall, because class membership is",
  "recomputed annually. Row labels carry the colour of the corresponding curve. The 2025 and",
  "2026 seasons are out of sample and are shown in Figure 5."))

message("fig03 done")

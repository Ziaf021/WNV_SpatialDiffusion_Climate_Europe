###############################################################################
# fig02a_distance_distribution.R
#
# Figure 2A. Distance from the nearest invaded region at first invasion.
#
# Built to the design of the submitted fig2_barplot_pressure.R panel A: five
# distance bins, counts on the y axis, and each bar labelled with its count and
# its percentage. Two changes:
#
#   1. The 2025 and 2026 seasons are stacked on top in the reserved colour, so
#      the out-of-sample events can be seen to fall in the same near bins.
#   2. The count and percentage label refers to the FITTED period only, so the
#      numbers on the figure are the ones quoted in the Results. The
#      out-of-sample addition is labelled separately above the bar as "+n" in
#      the reserved colour. Mixing the two into one percentage would have put a
#      number on the figure that appears nowhere in the text.
#
#   bash run_fig.sh 2a
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig02a")
on.exit(lg$close(), add = TRUE)

## ---- Data -------------------------------------------------------------------

fitted_ev <- ac %>% filter(event == 1) %>%
  transmute(km = nearest_invaded_km, set = "fitted")

oos_ev <- bind_rows(
  inv25 %>% transmute(km = nearest_invaded_km),
  inv26 %>% transmute(km = nearest_invaded_km)) %>%
  mutate(set = "oos")

ev <- bind_rows(fitted_ev, oos_ev) %>% mutate(cls = distance_class(km))

n_fit <- sum(ev$set == "fitted")
n_oos <- sum(ev$set == "oos")
message("fitted events: ", n_fit, "; out-of-sample events: ", n_oos)

# Denominator fixed here, outside the summarise: inside it, the per-class column
# n_oos masks the global n_oos, which gave total / (306 + class n_oos).
n_all <- n_fit + n_oos

lab <- ev %>%
  group_by(cls, .drop = FALSE) %>%
  summarise(n_fitted = sum(set == "fitted"),
            n_oos = sum(set == "oos"),
            total = n(), .groups = "drop") %>%
  # The label describes the bar. Earlier it gave the fitted-period count and its
  # percentage of 306 while the bar stood at the combined height, so a bar of 178
  # carried the label 47.7%, which is 146/306. Either the label describes what is
  # drawn or it does not belong on it.
  #
  # The denominator is every first invasion with a measurable distance, 362, and
  # not the 365 regions invaded by the end of 2026. The three regions first
  # invaded in 2006 have no distance to a previously invaded region, since nothing
  # had been invaded when they were, so they cannot be placed in any class here.
  mutate(pct = 100 * total / n_all,
         label = paste0(total, "\n(", sprintf("%.1f", pct), "%)"),
         fit_label = ifelse(n_oos > 0, as.character(n_fitted), ""),
         oos_label = ifelse(n_oos > 0, paste0("+ ", n_oos), ""))

print(as.data.frame(lab), row.names = FALSE)

ymax <- max(lab$total)

# A label sits inside its bar when the bar is tall enough to hold it, and above
# it otherwise. The >500 km bin holds a handful of events and would otherwise
# get a label wider than the bar.
lab <- lab %>%
  mutate(inside = total >= 0.16 * ymax,
         y_lab  = ifelse(inside, total / 2, total),
         v_lab  = ifelse(inside, 0.5, -0.45),
         # A bar with no out-of-sample events carries no coloured split above it,
         # so its total is the fitted-period count and is drawn in the fitted
         # colour and weight, like every other number on the panel.
         c_lab  = ifelse(inside, "white",
                         ifelse(n_oos == 0, TEAL, INK)),
         f_lab  = ifelse(inside | n_oos > 0, "plain", "bold"))

# 2007, not 2006. The fitted period begins in 2007 because distance to the
# nearest invaded region is undefined in the first year of follow-up; 2006 is the
# seed year, and the three regions invaded in it are sources rather than events.
SETS <- c(fitted = "2007-2024",
          oos    = "2025-2026")

ev <- ev %>% mutate(set_lab = factor(SETS[set], levels = unname(SETS)))

## ---- Figure -----------------------------------------------------------------

f2a <- ggplot(ev, aes(x = cls, fill = set_lab)) +
  geom_bar(width = 0.78) +
  geom_text(data = lab, aes(x = cls, y = y_lab, label = label,
                            vjust = v_lab, colour = c_lab,
                            fontface = f_lab),
            inherit.aes = FALSE, size = 3, lineheight = 0.95,
            show.legend = FALSE) +
  # The split above each bar is drawn as two pieces, each in the colour of the
  # series it counts, so the reader does not have to consult the key to know
  # which number is which.
  geom_text(data = lab %>% filter(nzchar(fit_label)),
            aes(x = cls, y = total, label = fit_label),
            inherit.aes = FALSE, vjust = -0.45, hjust = 1.35,
            size = 2.7, colour = TEAL, fontface = "bold") +
  geom_text(data = lab %>% filter(nzchar(oos_label)),
            aes(x = cls, y = total, label = oos_label),
            inherit.aes = FALSE, vjust = -0.45, hjust = -0.35,
            size = 2.7, colour = COL_OOS, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_manual(values = stats::setNames(c(TEAL, COL_OOS), unname(SETS)),
                    name = NULL) +
  scale_x_discrete(expand = c(0, 0.55)) +
  scale_y_continuous(limits = c(0, ymax * 1.12), expand = c(0, 0)) +
  labs(x = "Distance to nearest invaded region (km)",
       y = "Number of first invasion events") +
  theme_wnv_plot(base = 9) +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.text = element_text(size = 8),
        axis.title = element_text(colour = INK, size = 9),
        axis.text = element_text(colour = INK, size = 8.5))

saveRDS(f2a, file.path(DIR_FIG, "panel_2A.rds"))
save_fig(f2a, "Figure2A_distance_distribution", 127, 127)

write_caption("Figure2A_distance_distribution", paste(
  "Figure 2A. Distribution of first West Nile virus invasions by the distance",
  "between the newly invaded region and the nearest region already invaded at the",
  "start of that season. Bars give counts, and the label within each bar is the",
  "total for that class with its percentage of all",
  sprintf("%d first invasions with a measurable distance;", n_fit + n_oos),
  "above each bar the total is split into the fitted period and the 2025 and 2026",
  "seasons, which were not used to fit the model.",
  "Regions were followed from 2006 and the fitted period begins in 2007, because",
  "distance to the nearest invaded region is undefined in the first year of",
  "follow-up: the three regions first invaded in 2006 had no previously invaded",
  "region to be measured from, so they act as sources of exposure for other regions",
  "rather than as events and appear in no distance class here.", OOS_LABEL))

message("fig02a done")

###############################################################################
# fig01b_country_year.R
#
# Figure 1B. First invasions by country and year, 2006 to 2026.
#
# Built to the design of the submitted fig1_chronology_map.R panel B: red ramp,
# white cell labels, no tile borders, horizontal key above the panel, countries
# ordered by year of first invasion with the earliest at the bottom.
#
#   bash run_fig.sh 1b
###############################################################################

ROOT <- normalizePath(getwd())
source(file.path(ROOT, "code", "_wnv_figdata.R"))

lg <- fig_log("fig01b")
on.exit(lg$close(), add = TRUE)

## ---- Data -------------------------------------------------------------------

ctry <- nuts_all %>% st_drop_geometry() %>% select(NUTS_ID, CNTR_CODE)

tiles <- invasion %>%
  left_join(ctry, by = "NUTS_ID") %>%
  filter(!is.na(CNTR_CODE), first_year >= 2006, first_year <= 2026) %>%
  count(CNTR_CODE, first_year, name = "n")

message("regions represented: ", sum(tiles$n), " in ", nrow(tiles), " country-years")
print(as.data.frame(tiles %>% group_by(first_year) %>%
                      summarise(regions = sum(n), .groups = "drop")), row.names = FALSE)

# Earliest invaded at the bottom, most recent at the top. Ties on the first year
# are broken by the number of regions eventually invaded, which puts Romania
# below Hungary in 2006 and Türkiye at the top.
order_c <- tiles %>%
  group_by(CNTR_CODE) %>%
  summarise(first = min(first_year), tot = sum(n), .groups = "drop") %>%
  arrange(first, desc(tot)) %>%
  pull(CNTR_CODE)

grid <- expand.grid(CNTR_CODE = order_c, first_year = 2006:2026,
                    stringsAsFactors = FALSE) %>%
  left_join(tiles, by = c("CNTR_CODE", "first_year")) %>%
  mutate(country = factor(country_name(CNTR_CODE), levels = country_name(order_c)))

yrs <- 2006:2026
ylabs <- ifelse(yrs == 2026, "2026*", as.character(yrs))

max_n <- max(grid$n, na.rm = TRUE)
brks <- c(1, 3, 5, 9, 13)
brks <- brks[brks < max_n]
message("maximum regions in one country-year: ", max_n,
        "; legend breaks: ", paste(brks, collapse = ", "))

## ---- Figure -----------------------------------------------------------------
# No tile borders and no outline on the out-of-sample columns: every cell is
# drawn the same way, and the two seasons outside the fitted period are marked
# by the asterisk on the axis and by the caption instead.

f1b <- ggplot(grid, aes(x = factor(first_year), y = country)) +
  geom_tile(aes(fill = n), colour = NA) +
  geom_text(aes(label = ifelse(is.na(n), "", n)), colour = "white", size = 3) +
  scale_fill_count_binned(name = "Total newly invaded\nregions (NUTS3/GAUL)",
                          breaks = brks, limits = c(0, max_n)) +
  scale_x_discrete(labels = ylabs, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Year of first human West Nile virus invasion", y = "Invaded country") +
  theme_wnv_plot(base = 9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(colour = INK, size = 8.5),
        axis.text.x = element_text(colour = INK, size = 8.5, angle = 90,
                                   hjust = 1, vjust = 0.5),
        axis.title = element_text(colour = INK, size = 9),
        legend.position = c(0.02, 1.04), legend.justification = c(0, 1),
        legend.direction = "horizontal",
        legend.title = element_text(size = 8, colour = INK),
        legend.text = element_text(size = 7.5, colour = MID),
        plot.margin = margin(22, 6, 4, 4))

saveRDS(f1b, file.path(DIR_FIG, "panel_1B.rds"))
save_fig(f1b, "Figure1B_country_year", 190, 152)

write_caption("Figure1B_country_year", paste(
  "Figure 1B. Number of regions recording their first locally acquired human West",
  "Nile virus infection, by country and year, 2006 to 2026. Countries run from the",
  "earliest invaded at the bottom to the most recently invaded at the top, with",
  "ties broken by the number of regions eventually invaded. Shading and the number",
  "in each cell give the count of regions newly invaded in that country and year;",
  "blank cells are years with no first invasion in that country. The 2025 and 2026",
  "seasons lie outside the period used to fit the model.", OOS_LABEL,
  "Regions are NUTS3 units (Nomenclature of Territorial Units for Statistics,",
  "level 3) and the equivalent Eurostat statistical regions for candidate",
  "countries; for Serbia and Albania, ECDC reports cases under GAUL codes (Global",
  "Administrative Unit Layers), mapped to the 25 Serbian okruzi and the 12",
  "Albanian prefectures."))

message("fig01b done")

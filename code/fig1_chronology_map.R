############################################################
# Figure 1 (main text): chronology of WNV invasion across Europe
#   Panel A: NUTS3 map coloured by period of first human WNV emergence
#   Panel B: country x year tile heatmap of newly invaded NUTS3 regions
#
# Expects analysis_cases (counting-process data; see 00_data_preparation.R)
# in the environment. Downloads NUTS boundary shapefiles from Eurostat via
# the `eurostat` package on first run (cached locally by the package after
# that).
#
# Cleaned for the public repository from the author's original working
# script. Changes made during cleanup, noted here for transparency:
#   - `first_outbreak` and `nuts_sf` (NUTS3 boundaries) are now built
#     explicitly below; the original script assumed both were already
#     present in the session.
#   - `library(ggpubr)` added (theme_classic2(), ggarrange() both come
#     from it, not from ggthemes).
#   - Removed a dead code block that joined a `country_first` /
#     `country_names` lookup that was never used downstream; the actual
#     heatmap is built from `df_country`, which already carries country
#     names via `country_lookup_df`.
#   - The variable name `country_lookup` was reused for two different
#     objects in the original (a data.frame, then a named character
#     vector, overwriting the first). Split into `country_lookup_df`
#     (used for the join) and `country_lookup_labels` (used for the axis
#     labels) so both survive independently.
#   - Output filenames prefixed with Fig1_ for consistency with the
#     fig2_/fig3_/fig4_ scripts elsewhere in code/.
############################################################

## ---- Libraries ----
library(eurostat)
library(sf)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggnewscale)
library(ggthemes)
library(ggpubr)     # theme_classic2(), ggarrange()
library(svglite)

## ---- Data: first WNV emergence per NUTS3 region ----
# first_outbreak: one row per NUTS3 region that ever recorded a first
# locally acquired human WNV infection (event == 1), giving the calendar
# year of that first invasion. Deliberately NOT carrying CNTR_CODE here:
# nuts_sf (below) already has its own authoritative CNTR_CODE column, and
# joining two CNTR_CODE columns on NUTS_ID alone made dplyr silently
# rename them to CNTR_CODE.x/CNTR_CODE.y, breaking every downstream
# reference to plain CNTR_CODE (group_by(), left_join(), etc.).
first_outbreak <- analysis_cases %>%
  filter(event == 1) %>%
  group_by(NUTS_ID) %>%
  summarise(
    first_year = min(year),
    .groups = "drop"
  )

## ---- NUTS boundary shapefiles (Eurostat) ----
EU_NUTS0 <- eurostat::get_eurostat_geospatial(
  resolution = 10, nuts_level = 0, year = 2021
)

nuts_sf <- eurostat::get_eurostat_geospatial(
  resolution = 10, nuts_level = 3, year = 2021
)

## ---- Panel A: invasion chronology map ----
nuts_map <-
  nuts_sf %>%
  left_join(first_outbreak, by = "NUTS_ID") %>%
  mutate(
    `First WNV emergence` = case_when(
      is.na(first_year)   ~ "No invasion (in humans)",
      first_year <= 2010  ~ "2006-2010",
      first_year <= 2015  ~ "2011-2015",
      first_year <= 2020  ~ "2016-2020",
      TRUE                ~ "2021-2024"
    )
  )

Fig_Invasion_ChronologyMap <-
  ggplot(EU_NUTS0) +
  geom_sf(
    aes(geometry = geometry, fill = LEVL_CODE),
    color = "black", size = 0.1, alpha = 0.5, show.legend = FALSE
  ) +
  new_scale("fill") +
  geom_sf(
    data = nuts_map,
    aes(geometry = geometry, fill = `First WNV emergence`),
    color = NA
  ) +
  scale_fill_brewer(palette = "Spectral", direction = 1) +
  theme_map() +
  scale_y_continuous(limits = c(32, 70), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-10, 45), expand = c(0, 0)) +
  theme(
    legend.direction = "vertical",
    legend.position = c(0.7, 0.4),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.box = "vertical"
  )

Fig_Invasion_ChronologyMap

ggsave("Fig1_Invasion_ChronologyMap.svg", Fig_Invasion_ChronologyMap, height = 6, width = 8)
ggsave("Fig1_Invasion_ChronologyMap.pdf", Fig_Invasion_ChronologyMap, height = 6, width = 8)

## ---- Panel B: country x year tile heatmap ----

# Country code -> full name, used both to order countries by first
# appearance and to label the y-axis.
country_lookup_df <- data.frame(
  CNTR_CODE = c(
    "AT", "BG", "CY", "CZ", "DE", "EL", "ES", "FR",
    "HR", "HU", "IT", "MK", "NL", "PT", "RO", "SI",
    "SK", "TR"
  ),
  Country = c(
    "Austria", "Bulgaria", "Cyprus", "Czech Republic",
    "Germany", "Greece", "Spain", "France",
    "Croatia", "Hungary", "Italy",
    "North Macedonia", "Netherlands", "Portugal",
    "Romania", "Slovenia", "Slovakia", "Turkey"
  )
)

country_lookup_labels <- c(
  AT = "Austria", BG = "Bulgaria", HR = "Croatia", CY = "Cyprus",
  CZ = "Czech Republic", FR = "France", DE = "Germany", EL = "Greece",
  HU = "Hungary", IT = "Italy", NL = "Netherlands", MK = "North Macedonia",
  PT = "Portugal", RO = "Romania", SK = "Slovakia", SI = "Slovenia",
  ES = "Spain", TR = "Turkey"
)

df_country <- nuts_map %>%
  st_drop_geometry() %>%
  # Restrict to rows with a recorded invasion; using na.omit() here would
  # drop rows for any NA across nuts_sf's ~20 metadata columns (CAPT,
  # NAME_FREN, ISO3_CODE, etc., which are sparsely populated and
  # irrelevant to this figure), silently emptying the whole table.
  filter(!is.na(first_year), !is.na(CNTR_CODE)) %>%
  group_by(first_year, CNTR_CODE) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(country_lookup_df, by = "CNTR_CODE")

# Order countries on the y-axis by the year of their first invasion.
country_order <- df_country %>%
  group_by(CNTR_CODE) %>%
  summarise(first_country_year = min(first_year), .groups = "drop") %>%
  arrange(first_country_year)

df_country$CNTR_CODE <- factor(df_country$CNTR_CODE, levels = country_order$CNTR_CODE)

Fig_CountryHeatmap <-
  ggplot(
    df_country %>% mutate(Country = factor(Country)),
    aes(x = factor(first_year), y = CNTR_CODE, fill = n)
  ) +
  scale_y_discrete(labels = country_lookup_labels) +
  geom_tile(colour = NA, linewidth = 0.3) +
  scale_fill_gradientn(
    name = "Total newly invaded\nNUTS3 regions",
    colors = c("#FCAE91", "#FB6A4A", "#CB181D", "red"),
    values = scales::rescale(c(0, 5, 10, 15)),
    limits = c(1, 13),
    breaks = c(1, 5, 10, 13),
    labels = c("1", "5", "9", "13")
  ) +
  geom_text(aes(label = n), color = "white", size = 3) +
  labs(
    x = "Year of first human West Nile virus invasion",
    y = "Invaded country"
  ) +
  theme_classic2(base_size = 13) +
  theme(
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.35, 0.95),
    legend.direction = "horizontal",
    legend.text = element_text(angle = 0),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")
  )

Fig_CountryHeatmap

ggsave("Fig1_CountryHeatmap.svg", Fig_CountryHeatmap, height = 6, width = 9)
ggsave("Fig1_CountryHeatmap.pdf", Fig_CountryHeatmap, height = 6, width = 9)

## ---- Combine and export ----
Fig1_Main_AB <-
  ggarrange(
    Fig_Invasion_ChronologyMap,
    Fig_CountryHeatmap,
    labels = c("A", "B"),
    ncol = 2
  )

Fig1_Main_AB

ggsave("Fig1_Main_AB.svg", Fig1_Main_AB, height = 6, width = 15)
ggsave("Fig1_Main_AB.pdf", Fig1_Main_AB, height = 6, width = 15)

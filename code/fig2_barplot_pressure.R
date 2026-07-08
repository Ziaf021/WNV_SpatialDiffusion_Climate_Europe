############################################################
# Figure 2 (main text): distance barplot + pressure threshold
# Expects analysis_cases (counting-process data) in the environment.
############################################################

## ---- Libraries ----
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)     # theme_classic2(), ggarrange()

## ---- Data ----
#dat <- analysis_cases %>%
#  filter(year > 2005, year < 2025)

first_invaded <- analysis_cases%>%
  #dat %>%
  filter(event == 1)

## ---- Panel A: observed invasion distances ----
distance_summary <- first_invaded %>%
  mutate(
    distance_bin = cut(
      nearest_invaded_km,
      breaks = c(0, 75, 150, 300, 500, Inf),
      include.lowest = TRUE
    )
  ) %>%
  count(distance_bin) %>%
  mutate(
    percent = 100 * n / sum(n),
    label   = paste0(n, "\n(", round(percent, 1), "%)")
  )

Fig2_Main_DistBarplot <-
  ggplot(distance_summary, aes(distance_bin, n)) +
  geom_col(fill = "#CB181D", alpha = 0.85) +
  geom_text(aes(label = label), vjust = 0.5, size = 3.5) +
  labs(
    x = "Distance to nearest invaded region (km)",
    y = "Number of first invasion events"
  ) +
  theme_classic2(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(expand = c(0, 0)) +
  # Hardcoded upper limit tied to the current data (tallest bar ~121).
  scale_y_continuous(limits = c(0, 121 + 3), expand = c(0, 0))

Fig2_Main_DistBarplot

## ---- Panel B: empirical invasion rate vs neighbourhood pressure ----
pressure_curve <- analysis_cases %>%
  filter(!is.na(invasion_pressure)) %>%
  mutate(pressure_bin = ntile(invasion_pressure, 50)) %>%
  group_by(pressure_bin) %>%
  summarise(
    pressure_mid  = median(invasion_pressure),
    invasion_rate = mean(event),
    n             = n(),
    .groups = "drop"
  )

Fig_PressureThreshold <-
  ggplot(pressure_curve, aes(x = pressure_mid, y = invasion_rate * 100)) +
  geom_point(size = 2.5, colour = "#CB181D") +
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = "#CB181D",
    fill = "#FCAE91",
    linewidth = 1.2
  ) +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 1, 3, 7),
    labels = c(
      "None",
      "~1 nearby\ninvaded region",
      "~3 nearby\ninvaded regions",
      "~7 nearby\ninvaded regions"
    )
  ) +
  labs(
    x = expression(paste("Neighbourhood invasion pressure (", Phi, ")")),
    y = "Annual first WNV invasion probability (%)"
  ) +
  theme_classic2(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90))

Fig_PressureThreshold

## ---- Combine and export ----
Fig2_Main_AB_Barplot_PressureInvasion <-
  ggarrange(
    Fig2_Main_DistBarplot,
    Fig_PressureThreshold,
    labels = c("A", "B"),
    ncol = 2
  )

Fig2_Main_AB_Barplot_PressureInvasion

# Added: the combined panel was never written to disk in the original.
ggsave("Fig2_Main_AB_Barplot_PressureInvasion.svg",
       Fig2_Main_AB_Barplot_PressureInvasion, width = 10, height = 5)
ggsave("Fig2_Main_AB_Barplot_PressureInvasion.pdf",
       Fig2_Main_AB_Barplot_PressureInvasion, width = 10, height = 5)

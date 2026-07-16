## ============================================================
## fig_urbanisation_rate_by_phase.R
## Supplementary Figure S9: Phase-specific first-invasion rate by
## Urbanization category (per 1,000 at-risk region-years).
##
## This script was missing from the repository entirely; only its
## output PDF had been committed previously, from a run that predates
## this script's addition. Reconstructed here to match Supplementary
## Figure S9's caption exactly: first-invasion rate per 1,000 at-risk
## region-years (not compositional percentage, which is Figure S8),
## by Urbanization category and invasion phase, with a Poisson
## likelihood-ratio test for the Urbanization x phase interaction,
## repeated after excluding the 2020-2024 period as a sensitivity
## check (per the caption's reported p = 0.006 overall vs p = 0.54
## with Post-COVID excluded).
## ============================================================

library(dplyr)
library(ggplot2)
library(ggpubr)     # theme_classic2()

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!exists("analysis_cases"))
  analysis_cases <- readRDS("analysis_cases_with_cats.rds")

urbn_labels  <- c("1" = "Urban", "2" = "Intermediate", "3" = "Rural")
phase_levels <- c("Early", "Middle", "Pre-COVID", "Post-COVID")

## One row per region-year already represents one at-risk region-year
## (counting-process format from 00_data_preparation.R); event == 1
## marks that region-year's first invasion.
rate_dat <- analysis_cases %>%
  filter(!is.na(URBN_TYPE), !is.na(invasion_phase)) %>%
  mutate(
    URBN  = factor(urbn_labels[as.character(URBN_TYPE)],
                   levels = c("Urban", "Intermediate", "Rural")),
    phase = factor(invasion_phase, levels = phase_levels)
  )

## ---- Phase x Urbanization: at-risk region-years and events ----
rate_summary <- rate_dat %>%
  group_by(phase, URBN) %>%
  summarise(
    region_years = n(),
    events       = sum(event == 1),
    .groups = "drop"
  ) %>%
  mutate(
    rate_per_1000 = 1000 * events / region_years,
    # Exact Poisson 95% CI for the event count, scaled to the same rate units.
    ci_low  = 1000 * qgamma(0.025, shape = events) / region_years,
    ci_high = 1000 * qgamma(0.975, shape = events + 1) / region_years
  )

cat("Phase x Urbanization: region-years, events, rate per 1,000 region-years\n")
print(rate_summary %>% mutate(rate_per_1000 = round(rate_per_1000, 2)))

## ---- Poisson interaction test: Urbanization x phase, all phases ----
pois_full <- glm(
  event ~ URBN * phase, family = poisson(), data = rate_dat
)
pois_main <- glm(
  event ~ URBN + phase, family = poisson(), data = rate_dat
)
lrt_full <- anova(pois_main, pois_full, test = "LRT")
p_full <- lrt_full$`Pr(>Chi)`[2]
cat(sprintf(
  "\nUrbanization x phase interaction (all phases), Poisson LRT: p = %.3g\n",
  p_full
))
cat("Expected: p = 0.006\n")

## ---- Sensitivity: exclude 2020-2024 (Post-COVID) ----
rate_dat_excl <- rate_dat %>% filter(phase != "Post-COVID")
pois_full_excl <- glm(
  event ~ URBN * phase, family = poisson(),
  data = droplevels(rate_dat_excl)
)
pois_main_excl <- glm(
  event ~ URBN + phase, family = poisson(),
  data = droplevels(rate_dat_excl)
)
lrt_excl <- anova(pois_main_excl, pois_full_excl, test = "LRT")
p_excl <- lrt_excl$`Pr(>Chi)`[2]
cat(sprintf(
  "Urbanization x phase interaction (Post-COVID excluded), Poisson LRT: p = %.3g\n",
  p_excl
))
cat("Expected: p = 0.54\n\n")

## ---- Figure S9 ----
green3 <- c("Urban" = "#006D2C", "Intermediate" = "#41AB5D", "Rural" = "#A1D99B")

Fig_Urban_Rate_Phase <-
  ggplot(rate_summary, aes(x = phase, y = rate_per_1000, colour = URBN, group = URBN)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.6, alpha = 0.6) +
  scale_colour_manual(values = green3, name = "Urbanisation") +
  scale_y_log10() +
  labs(
    x = "Invasion phase",
    y = "First-invasion rate\n(per 1,000 at-risk region-years, log scale)"
  ) +
  theme_classic2(base_size = 13) +
  theme(legend.position = "top")

Fig_Urban_Rate_Phase
ggsave(file.path(output_dir, "FigureS9_Urbanisation_rate_by_phase.pdf"),
       Fig_Urban_Rate_Phase, width = 8, height = 5)
cat("Saved: FigureS9_Urbanisation_rate_by_phase.pdf\n")

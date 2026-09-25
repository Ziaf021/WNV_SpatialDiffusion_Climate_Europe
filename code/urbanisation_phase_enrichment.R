## ============================================================
## urbanisation_phase_enrichment.R
## Phase-specific urbanisation invasion RATE and the decisive
## urbanisation x phase interaction (Poisson, offset = at-risk
## region-years). Green palette; year-range x-axis.
## Includes a Post-COVID-excluded robustness refit.
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
green3 <- c("Urban" = "#006D2C", "Intermediate" = "#41AB5D", "Rural" = "#A1D99B")

## Person-time (at-risk region-years) and events per phase x urbanisation.
strata <- analysis_cases %>%
  filter(!is.na(URBN_TYPE)) %>%
  mutate(URBN  = factor(urbn_labels[as.character(URBN_TYPE)],
                        levels = c("Urban", "Intermediate", "Rural")),
         phase = factor(invasion_phase, levels = phase_levels)) %>%
  group_by(phase, URBN) %>%
  summarise(events = sum(event), pyears = dplyr::n(), .groups = "drop") %>%
  mutate(rate1000 = 1000 * events / pyears,
         lo = 1000 * qgamma(0.025, events)     / pyears,
         hi = 1000 * qgamma(0.975, events + 1) / pyears)

cat("Events, at-risk region-years, and rate per 1,000:\n")
print(as.data.frame(strata %>% mutate(across(c(rate1000, lo, hi), \(x) round(x, 2)))))

## ---- Decisive test: urbanisation x phase interaction --------------------
m_add  <- glm(events ~ URBN + phase, offset = log(pyears), family = poisson, data = strata)
m_full <- glm(events ~ URBN * phase, offset = log(pyears), family = poisson, data = strata)
lr <- anova(m_add, m_full, test = "LRT")
int_p <- lr$`Pr(>Chi)`[2]
cat(sprintf("\nUrbanisation x phase interaction: LR chi-square = %.2f, df = %d, p = %.4f\n",
            lr$Deviance[2], lr$Df[2], int_p))
p_lab <- ifelse(int_p < 0.001, "p < 0.001", sprintf("p = %.3f", int_p))

## ---- Robustness: exclude the Post-COVID phase ---------------------------
## Tests whether the urbanisation gradient is present across Early/Middle/
## Pre-COVID alone, i.e. before the surveillance-disrupted period.
strata_pre <- strata %>% filter(phase != "Post-COVID") %>% mutate(phase = droplevels(phase))
m_add_e  <- glm(events ~ URBN + phase, offset = log(pyears), family = poisson, data = strata_pre)
m_full_e <- glm(events ~ URBN * phase, offset = log(pyears), family = poisson, data = strata_pre)
lr_e <- anova(m_add_e, m_full_e, test = "LRT")
cat(sprintf("\nInteraction EXCLUDING Post-COVID (Early/Middle/Pre-COVID only):\n  LR chi-square = %.2f, df = %d, p = %.4f\n",
            lr_e$Deviance[2], lr_e$Df[2], lr_e$`Pr(>Chi)`[2]))

## ---- Per-phase enrichment (event share / at-risk share) -----------------
enrich <- strata %>%
  group_by(phase) %>%
  mutate(enrichment = (events / sum(events)) / (pyears / sum(pyears))) %>%
  ungroup()
cat("\nPer-phase enrichment (>1 = over-represented vs at-risk pool):\n")
print(as.data.frame(enrich %>% transmute(phase, URBN, enrichment = round(enrichment, 2))))

## ---- Year-range x-axis labels -------------------------------------------
phase_years <- analysis_cases %>%
  filter(!is.na(invasion_phase)) %>%
  transmute(phase = factor(invasion_phase, levels = phase_levels), year) %>%
  group_by(phase) %>% summarise(y0 = min(year), y1 = max(year), .groups = "drop") %>%
  # Early is labelled by its definitional lower bound (study start, 2006);
  # the modelling data begin in 2007, so no 2006 events enter this analysis.
  mutate(y0 = ifelse(as.character(phase) == "Early", 2006, y0),
         lab = paste0(as.character(phase), "\n(", y0, "\u2013", y1, ")"))
phase_labels <- setNames(phase_years$lab, as.character(phase_years$phase))

## ---- Figure: invasion rate by phase and urbanisation --------------------
Fig_rate <- ggplot(strata, aes(x = phase, y = rate1000, colour = URBN, group = URBN)) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18,
                position = position_dodge(width = 0.3)) +
  scale_colour_manual(values = green3, name = "Urbanisation") +
  scale_x_discrete(labels = phase_labels) +
  scale_y_log10() +
  annotate("text", x = 0.6, y = max(strata$hi),
           label = paste0("Urbanisation \u00D7 phase interaction: ", p_lab),
           hjust = 0, size = 3.6) +
  labs(x = "Invasion phase",
       y = "First-invasion rate (per 1,000 region-years, log scale)") +
  theme_classic2(base_size = 13) +
  theme(legend.position = "top")

Fig_rate
ggsave(file.path(output_dir, "Fig_Urbanisation_rate_by_phase.pdf"),
       Fig_rate, width = 8, height = 5)

## ============================================================
## 06_supplementary_figures.R
## Analyses 2.12 – 2.13
## Figure S3: Temperature × invasion pressure quintiles (violin)
## Figure S4: Invasion distance × temporal phase (violin)
## WNV Spatial Diffusion Manuscript — Farooq, Umeå University
## ============================================================

library(dplyr)
library(ggplot2)
library(ggthemes)

if (!exists("analysis_cases")) {
  analysis_cases <- readRDS("analysis_cases_with_cats.rds")
}

first_invaded <- analysis_cases %>% filter(event == 1)
cat("First invasion events:", nrow(first_invaded), "(expected: 269)\n\n")

## ============================================================
## ANALYSIS 2.12 — Supplementary Figure S3
##   Temperature at first invasion by invasion pressure quintile
## ============================================================

cat("── Analysis 2.12: Temperature × invasion pressure (Figure S3) ──\n")

first_invaded <- first_invaded %>%
  mutate(
    pressure_cat = cut(
      invasion_pressure,
      breaks = quantile(invasion_pressure, probs = seq(0, 1, 0.2), na.rm = TRUE),
      labels = c("Very low", "Low", "Moderate", "High", "Very high"),
      include.lowest = TRUE
    )
  )

# Check group sizes
n_s3 <- first_invaded %>%
  group_by(pressure_cat) %>%
  summarise(n = n(), .groups = "drop")
cat("Group sizes (expected ~53-54 each):\n")
print(n_s3)

# Kruskal-Wallis test
kw_s3 <- kruskal.test(tas ~ pressure_cat, data = first_invaded)
cat(sprintf("\nKruskal-Wallis: chi-sq = %.1f, df = %d, p = %.2e\n",
            kw_s3$statistic, kw_s3$parameter, kw_s3$p.value))
cat("Expected: chi-sq = 43.1, df = 4, p = 9.9e-09\n")

# Pairwise tests
pairwise_s3 <- pairwise.wilcox.test(
  x               = first_invaded$tas,
  g               = first_invaded$pressure_cat,
  p.adjust.method = "BH"
)
cat("\nPairwise Wilcoxon (BH-adjusted) — Supplementary Figure S3:\n")
print(pairwise_s3)
cat("Expected: Very Low differs from all higher categories (all adj p < 2.1e-05)\n\n")

# Figure S3
reds_s3 <- c("#FEE5D9","#FCAE91","#FB6A4A","#CB181D","#7F0000")

Fig_S3 <- ggplot(
  first_invaded,
  aes(x = pressure_cat, y = tas, fill = pressure_cat)
) +
  geom_violin(alpha = 0.8, colour = "grey30") +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1.4) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_text(
    data = n_s3,
    aes(x = pressure_cat, label = paste0("n=", n)),
    y = max(first_invaded$tas, na.rm = TRUE) + 0.5,
    inherit.aes = FALSE, size = 3.8
  ) +
  annotate(
    "text", x = 3, y = min(first_invaded$tas, na.rm = TRUE) - 0.5,
    label = sprintf("Kruskal-Wallis chi-sq=%.1f, p=%.2e",
                    kw_s3$statistic, kw_s3$p.value),
    size = 3.2, colour = "grey30"
  ) +
  scale_fill_manual(values = reds_s3, guide = "none") +
  labs(
    x = "Neighbourhood invasion pressure category",
    y = "Annual mean temperature (°C)"
  ) +
  theme_classic2(base_size = 13)

ggsave("Supplementary_Figure_S3_TempByPressure.pdf", Fig_S3,
       width = 8, height = 5)
cat("Saved: Supplementary_Figure_S3_TempByPressure.pdf\n\n")

## ============================================================
## ANALYSIS 2.13 — Supplementary Figure S4
##   Invasion distance × temporal phase (violin)
## ============================================================

cat("── Analysis 2.13: Invasion distance × phase (Figure S4) ──\n")

first_invaded <- first_invaded %>%
  mutate(
    invasion_phase = factor(
      invasion_phase,
      levels = c("Early", "Middle", "Pre-COVID", "Post-COVID")
    )
  )

# Check group sizes
n_s4 <- first_invaded %>%
  group_by(invasion_phase) %>%
  summarise(n = n(), .groups = "drop")
cat("Events by invasion phase:\n")
print(n_s4)

# Kruskal-Wallis test
kw_s4 <- kruskal.test(nearest_invaded_km ~ invasion_phase, data = first_invaded)
cat(sprintf("\nKruskal-Wallis: chi-sq = %.1f, df = %d, p = %.2e\n",
            kw_s4$statistic, kw_s4$parameter, kw_s4$p.value))
cat("Expected: chi-sq = 25.1, df = 3, p = 1.5e-05\n")

# Pairwise tests
pairwise_s4 <- pairwise.wilcox.test(
  x               = first_invaded$nearest_invaded_km,
  g               = first_invaded$invasion_phase,
  p.adjust.method = "BH"
)
cat("\nPairwise Wilcoxon (BH-adjusted) — Supplementary Figure S4:\n")
print(pairwise_s4)
cat("Expected: Early differs from Middle (p=7.1e-04), Pre-COVID (p=9.2e-06),\n")
cat("          Post-COVID (p=1.8e-04); Later phases indistinguishable (p>0.16)\n\n")

# Figure S4
phase_cols <- c("#FEE5D9","#FC9272","#FB6A4A","#CB181D")

Fig_S4 <- ggplot(
  first_invaded,
  aes(x = invasion_phase, y = nearest_invaded_km, fill = invasion_phase)
) +
  geom_violin(alpha = 0.85, colour = "grey30") +
  geom_jitter(width = 0.10, alpha = 0.30, size = 1.5) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_text(
    data = n_s4,
    aes(x = invasion_phase, label = paste0("n=", n)),
    y = max(first_invaded$nearest_invaded_km, na.rm = TRUE) * 1.15,
    inherit.aes = FALSE, size = 4
  ) +
  annotate(
    "text", x = 2.5,
    y = max(first_invaded$nearest_invaded_km, na.rm = TRUE) * 1.25,
    label = sprintf("Kruskal-Wallis chi-sq=%.1f, p=%.2e",
                    kw_s4$statistic, kw_s4$p.value),
    size = 3.2, colour = "grey30"
  ) +
  scale_y_log10(
    breaks = c(10, 25, 50, 100, 250, 500, 1000)
  ) +
  scale_fill_manual(values = phase_cols, guide = "none") +
  labs(
    x = "Invasion phase",
    y = "Distance to nearest invaded region (km)"
  ) +
  theme_classic2(base_size = 13)

ggsave("Supplementary_Figure_S4_DistanceByPhase.pdf", Fig_S4,
       width = 8, height = 5)
cat("Saved: Supplementary_Figure_S4_DistanceByPhase.pdf\n\n")

## ============================================================
## FINAL VERIFICATION SUMMARY
## ============================================================

cat("=" , strrep("=", 60), "\n", sep="")
cat("FINAL VERIFICATION SUMMARY\n")
cat("=" , strrep("=", 60), "\n", sep="")
cat("Run all scripts in order 02 → 03 → 04 → 05 → 06\n\n")
cat("KEY NUMBERS TO VERIFY AGAINST MANUSCRIPT:\n\n")
cat("RESULTS SECTION:\n")
cat("  Total first invasions: 269\n")
cat("  Within 75 km:  121 (45.0%)\n")
cat("  Within 150 km: 192 (71.4%)\n")
cat("  Within 300 km: 236 (87.7%)\n")
cat("  Beyond 500 km:   6  (2.2%)\n\n")
cat("PRIMARY MODEL (cox_source_full):\n")
cat("  Distance HR:     0.56 (0.51-0.63), 44% reduction\n")
cat("  Pressure HR:     1.37 (1.23-1.52), 37% increase\n")
cat("  Intensity HR:    1.32 (1.09-1.61), p=0.005\n")
cat("  Temperature HR:  1.24 (1.18-1.31), 24% increase\n")
cat("  Humidity HR:     0.95 (0.92-0.99), p=0.004\n")
cat("  Precipitation:   0.97 (0.84-1.12), p=0.688 (NOT SIGNIFICANT)\n")
cat("  Concordance:     0.872\n\n")
cat("INTERACTION MODEL:\n")
cat("  Distance x Temperature HR: 1.07 (1.04-1.09), p<0.001\n\n")
cat("SENSITIVITY ANALYSES:\n")
cat("  Optimal lambda: 50 km (lowest AIC)\n")
cat("  Climate anomaly DELTA AIC: ~120 (absolute values better)\n")
cat("  LOCO: stable across all 18 exclusions\n\n")
cat("PH DIAGNOSTICS:\n")
cat("  distance100 Schoenfeld p: 0.011\n")
cat("  pressure_z  Schoenfeld p: <0.001\n")
cat("  Outlier exclusion: violation persists (p unchanged)\n")
cat("  Time-varying model: pressure tt() p=0.035 (significant)\n")
cat("  Case-weighted kernel: DELTA AIC >= 21 vs unweighted\n")
cat("  Phase stratification: HRs unchanged (violation not due to baseline)\n\n")
cat("SUPPLEMENTARY FIGURES:\n")
cat("  Fig S3 KW: chi-sq=43.1, p=9.9e-09\n")
cat("  Fig S4 KW: chi-sq=25.1, p=1.5e-05\n")

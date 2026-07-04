## ============================================================
## 03_interaction_models.R
## Analyses 1.13 – 1.15
## Transmission intensity LR test | Distance × temperature interaction |
## Kruskal–Wallis temperature by distance category
## WNV Spatial Diffusion Manuscript — Farooq, Umeå University
## ============================================================
## REQUIRES: analysis_cases (with distance_cat column from 02_survival_analyses.R)
##           cox_source_full (saved from 02_survival_analyses.R)
##           first_invaded   (subset where event == 1)
## ============================================================

library(survival)
library(dplyr)
library(ggplot2)
library(ggthemes)

# ── Load objects ─────────────────────────────────────────────
if (!exists("cox_source_full")) {
  cox_source_full <- readRDS("cox_source_full.rds")
}
if (!exists("analysis_cases")) {
  analysis_cases <- readRDS("analysis_cases_with_cats.rds")
}

first_invaded <- analysis_cases %>% filter(event == 1)

## ============================================================
## ANALYSIS 1.13 — LR test: does transmission intensity
##                 improve primary model fit?
## ============================================================

cat("── Analysis 1.13: LR test for source-region transmission intensity ──\n")

# Model WITHOUT log_cases_nearest
cox_without_intensity <- coxph(
  Surv(start, stop, event) ~
    distance100 + pressure_z +
    tas + hurs + scale(pr) + factor(URBN_TYPE),
  data  = analysis_cases,
  ties  = "efron"
)

# LR test
lr_result <- anova(cox_without_intensity, cox_source_full)
cat("\nLikelihood-ratio test — with vs without log_cases_nearest:\n")
print(lr_result)

# AIC comparison
aic_without <- AIC(cox_without_intensity)
aic_with    <- AIC(cox_source_full)
delta_aic   <- aic_without - aic_with
cat(sprintf("\nAIC without intensity: %.1f\n", aic_without))
cat(sprintf("AIC with intensity:    %.1f\n", aic_with))
cat(sprintf("DELTA AIC:             %.1f\n", delta_aic))
cat("Expected: DELTA AIC = 4.7, LR p = 0.01\n\n")

## ============================================================
## ANALYSIS 1.14 — Distance × temperature interaction model
##                 (Figure 4B: predicted hazard surface)
## ============================================================

cat("── Analysis 1.14: Distance × temperature interaction ──\n")

cox_interaction <- coxph(
  Surv(start, stop, event) ~
    distance100 * tas +
    pressure_z + log_cases_nearest +
    hurs + scale(pr) + factor(URBN_TYPE),
  data  = analysis_cases,
  ties  = "efron"
)

summary_int <- summary(cox_interaction)
cat("\nInteraction model summary:\n")
print(summary_int$coefficients)

# Extract interaction term
int_row <- grep("distance100:tas", rownames(summary_int$coefficients))
int_HR  <- summary_int$conf.int[int_row, "exp(coef)"]
int_lo  <- summary_int$conf.int[int_row, "lower .95"]
int_hi  <- summary_int$conf.int[int_row, "upper .95"]
int_p   <- summary_int$coefficients[int_row, "Pr(>|z|)"]

cat(sprintf(
  "\nInteraction HR: %.2f (95%% CI %.2f-%.2f), p = %.4f\n",
  int_HR, int_lo, int_hi, int_p
))
cat("Expected: HR=1.07 (1.04-1.09), p<0.001\n\n")

# ── Predicted hazard surface for Figure 4B ──

newdat <- expand.grid(
   distance100 =seq(0,10,length.out = 100),
   tas =seq( min(analysis_cases$tas, na.rm = TRUE), max(analysis_cases$tas, na.rm = TRUE),
      length.out = 100))

newdat$pressure_z <- mean(analysis_cases$pressure_z,na.rm = TRUE)

newdat$log_cases_nearest <- mean( analysis_cases$log_cases_nearest,na.rm = TRUE)

newdat$hurs <- mean( analysis_cases$hurs,na.rm = TRUE)
newdat$pr <- mean( analysis_cases$pr,na.rm = TRUE)
newdat$URBN_TYPE <- factor(
  levels(factor(analysis_cases$URBN_TYPE))[1],
  levels = levels(factor(analysis_cases$URBN_TYPE)))

newdat$pred_hazard <-predict( cox_interaction,newdata = newdat, type = "risk"  )

newdat$hazard_pct <-percent_rank( newdat$pred_hazard) * 100

cat("Hazard surface computed:", nrow(newdat), "combinations\n")
cat("Max relative hazard at:",
    round(newdat$distance100[which.max(newdat$pred_hazard)] * 100, 0),
    "km,",
    round(newdat$tas[which.max(newdat$pred_hazard)], 1),
    "deg C\n\n")

# Merge back for overlaying observed points in Figure 4B
first_invaded_fig4b <- first_invaded %>%
  select(nearest_invaded_km, tas) %>%
  mutate(distance100 = nearest_invaded_km / 100)

## ============================================================
## ANALYSIS 1.15 — Kruskal–Wallis and pairwise tests
##                 Temperature by distance category at first invasion
##                 (Figure 4A)
## ============================================================

cat("── Analysis 1.15: KW test — temperature by distance category ──\n")

first_invaded <- first_invaded %>%
  mutate(
    distance_cat = cut(
      nearest_invaded_km,
      breaks = c(0, 75, 150, 300, 500, Inf),
      labels = c("<=75 km", "76-150 km", "151-300 km", "301-500 km", ">500 km"),
      right  = TRUE
    )
  )

kw_result <- kruskal.test(tas ~ distance_cat, data = first_invaded)
cat("\nKruskal-Wallis test:\n")
print(kw_result)
cat("Expected: chi-sq = 29.4, df = 4, p < 0.001\n\n")

pairwise_result <- pairwise.wilcox.test(
  x                = first_invaded$tas,
  g                = first_invaded$distance_cat,
  p.adjust.method  = "BH"
)
cat("Pairwise Wilcoxon tests (BH-adjusted):\n")
print(pairwise_result)
cat("Expected: <=75 km vs all categories beyond 150 km: all adjusted p < 0.01\n\n")

# ── Figure 4A: Violin plot — temperature by distance category ──
n_labels_4a <- first_invaded %>%
  group_by(distance_cat) %>%
  summarise(n = n(), .groups = "drop")

Fig_4A <- ggplot(
  first_invaded,
  aes(x = distance_cat, y = tas, fill = distance_cat)
) +
  geom_violin(alpha = 0.8, colour = "grey30") +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1.4) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_text(
    data = n_labels_4a,
    aes(x = distance_cat, label = paste0("n=", n)),
    y = max(first_invaded$tas, na.rm = TRUE) + 0.5,
    inherit.aes = FALSE, size = 3.5
  ) +
  scale_fill_manual(
    values = c("#FEE5D9","#FCAE91","#FB6A4A","#CB181D","#7F0000"),
    guide  = "none"
  ) +
  labs(
    x = "Distance to nearest invaded region",
    y = "Annual mean temperature (°C)"
  ) +
  annotate("text", x = 3, y = 2,
           label = paste0("Kruskal-Wallis p < 0.001"),
           size = 3.5, colour = "grey30") +
  theme_classic2(base_size = 13)

ggsave("Figure_4A_TempByDistance.pdf", Fig_4A, width = 8, height = 5)
cat("Saved: Figure_4A_TempByDistance.pdf\n")

## ============================================================
## SAVE INTERACTION MODEL AND GRID FOR FIGURE 4B
## ============================================================

saveRDS(cox_interaction, file = "cox_interaction.rds")
saveRDS(pred_grid,       file = "pred_grid_fig4b.rds")
cat("Saved: cox_interaction.rds, pred_grid_fig4b.rds\n")

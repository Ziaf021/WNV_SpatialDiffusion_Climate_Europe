## ============================================================
## 02_survival_analyses.R
## Analyses 1.10 – 1.12
## Univariable Cox models | Spline distance-hazard | Primary multivariable model
## WNV Spatial Diffusion Manuscript — Farooq, Umeå University
## ============================================================
## REQUIRES: analysis_cases (loaded from 00_data_preparation.R or .rds)
##           first_invaded  (subset of analysis_cases where event == 1)
## ============================================================

library(survival)
library(dplyr)
library(broom)

# ── Safety checks ────────────────────────────────────────────
stopifnot(
  exists("analysis_cases"),
  all(c("start","stop","event","distance100","pressure_z",
        "log_cases_nearest","tas","hurs","pr",
        "URBN_TYPE","CNTR_CODE") %in% names(analysis_cases))
)

first_invaded <- analysis_cases %>% filter(event == 1)
cat("Total first invasion events:", nrow(first_invaded), "\n")
cat("Expected: 269\n\n")

## ============================================================
## ANALYSIS 1.10 — Univariable Cox proportional hazards models
## ============================================================

cat("── Analysis 1.10: Univariable Cox models ──\n")

fit_univar <- function(formula_str, label, data = analysis_cases) {
  fit  <- coxph(as.formula(formula_str), data = data, ties = "efron")
  s    <- summary(fit)
  coef <- s$coefficients
  ci   <- s$conf.int
  data.frame(
    label    = label,
    HR       = round(coef[1, "exp(coef)"], 3),
    lower95  = round(ci[1, "lower .95"],   3),
    upper95  = round(ci[1, "upper .95"],   3),
    p        = signif(coef[1, "Pr(>|z|)"], 3),
    stringsAsFactors = FALSE
  )
}

univar_results <- bind_rows(
  fit_univar("Surv(start,stop,event) ~ distance100",
             "Distance per 100 km"),
  fit_univar("Surv(start,stop,event) ~ pressure_z",
             "Invasion pressure per SD"),
  fit_univar("Surv(start,stop,event) ~ log_cases_nearest",
             "Source-region intensity (log scale)"),
  fit_univar("Surv(start,stop,event) ~ tas",
             "Temperature per 1 deg C"),
  fit_univar("Surv(start,stop,event) ~ scale(pr)",
             "Precipitation per SD"),
  fit_univar("Surv(start,stop,event) ~ hurs",
             "Relative humidity per 1 pct"),
  fit_univar("Surv(start,stop,event) ~ factor(URBN_TYPE)",
             "Urbanisation (vs urban)")
)

print(univar_results)

cat("\nEXPECTED VALUES (Table 1):\n")
cat("  Distance:    HR=0.49 (0.44-0.54),  p<0.001\n")
cat("  Pressure:    HR=1.89 (1.74-2.05),  p<0.001\n")
cat("  Intensity:   HR=1.10 (0.99-1.23),  p=0.088\n")
cat("  Temperature: HR=1.22 (1.17-1.26),  p<0.001\n")
cat("  Precip:      HR=0.83 (0.73-0.94),  p=0.004\n")
cat("  Humidity:    HR=0.95 (0.93-0.96),  p<0.001\n\n")

## ============================================================
## ANALYSIS 1.11 — Spline-based Cox model (Figure 3B)
## ============================================================

cat("── Analysis 1.11: Spline distance-hazard (Figure 3B) ──\n")

cox_spline <- coxph(
  Surv(start, stop, event) ~
    pspline(distance100, df = 4) +
    pressure_z + log_cases_nearest +
    tas + hurs + scale(pr) + factor(URBN_TYPE),
  data  = analysis_cases,
  ties  = "efron"
)

# Likelihood-ratio test vs linear distance specification
cox_linear_dist <- coxph(
  Surv(start, stop, event) ~
    distance100 + pressure_z + log_cases_nearest +
    tas + hurs + scale(pr) + factor(URBN_TYPE),
  data  = analysis_cases,
  ties  = "efron"
)

lr_test_spline <- anova(cox_linear_dist, cox_spline)
cat("LR test: linear vs spline distance specification:\n")
print(lr_test_spline)

# Predicted spline curve for Figure 3B
# Vary distance100 across observed range; hold other covariates at means
dist_grid <- data.frame(
  distance100      = seq(0.1, 10, length.out = 200),
  pressure_z       = mean(analysis_cases$pressure_z,       na.rm = TRUE),
  log_cases_nearest= mean(analysis_cases$log_cases_nearest, na.rm = TRUE),
  tas              = mean(analysis_cases$tas,               na.rm = TRUE),
  hurs             = mean(analysis_cases$hurs,              na.rm = TRUE),
  pr               = mean(analysis_cases$pr,                na.rm = TRUE),
  URBN_TYPE        = 3L    # reference: rural (most common)
)

# Reference value for HR = 1 at distance100 = 1 (100 km)
ref_row <- dist_grid[1, ]
ref_row$distance100 <- 1.0

pred_lp  <- predict(cox_spline, newdata = dist_grid,  type = "lp")
pred_ref <- predict(cox_spline, newdata = ref_row,     type = "lp")
dist_grid$log_HR <- pred_lp - as.numeric(pred_ref)
dist_grid$HR     <- exp(dist_grid$log_HR)

cat("\nSpline curve: % invasions by distance:\n")
cat("  Hazard at 75km relative to 100km ref:  ",
    round(exp(dist_grid$log_HR[which.min(abs(dist_grid$distance100 - 0.75))]), 3), "\n")
cat("  Hazard at 300km relative to 100km ref: ",
    round(exp(dist_grid$log_HR[which.min(abs(dist_grid$distance100 - 3.0))]), 3), "\n\n")

## ============================================================
## ANALYSIS 1.9 — Kaplan–Meier survival curves (Figure 3B inset)
## ============================================================

cat("── Analysis 1.9: Kaplan-Meier by distance category ──\n")

analysis_cases <- analysis_cases %>%
  mutate(
    distance_cat = cut(
      nearest_invaded_km,
      breaks = c(0, 75, 150, 300, 500, Inf),
      labels = c("<=75 km", "76-150 km", "151-300 km", "301-500 km", ">500 km"),
      right  = TRUE
    )
  )

km_fit <- survfit(
  Surv(start, stop, event) ~ distance_cat,
  data = analysis_cases
)

cat("Kaplan-Meier fit summary:\n")
print(km_fit)

## ============================================================
## ANALYSIS 1.12 — Primary multivariable Cox model (cox_source_full)
##                 (Figure 3A forest plot)
## ============================================================

cat("── Analysis 1.12: Primary multivariable Cox model ──\n")

cox_source_full <- coxph(
  Surv(start, stop, event) ~
    distance100 + pressure_z + log_cases_nearest +
    tas + hurs + scale(pr) + factor(URBN_TYPE),
  data  = analysis_cases,
  ties  = "efron"
)

summary_full <- summary(cox_source_full)
cat("\nPrimary model summary:\n")
print(summary_full)

cat("\n── Extracted HRs and 95% CIs ──\n")
coef_table <- data.frame(
  term    = rownames(summary_full$coefficients),
  HR      = round(summary_full$conf.int[, "exp(coef)"],  3),
  lower95 = round(summary_full$conf.int[, "lower .95"],  3),
  upper95 = round(summary_full$conf.int[, "upper .95"],  3),
  p       = signif(summary_full$coefficients[, "Pr(>|z|)"], 3)
)
print(coef_table)

cat("\nConcordance:", round(summary_full$concordance["C"], 3), "\n")

cat("\nEXPECTED VALUES (verified R output):\n")
expected <- data.frame(
  term    = c("distance100","pressure_z","log_cases_nearest",
              "tas","hurs","scale(pr)","URBN2","URBN3"),
  HR      = c(0.564, 1.368, 1.323, 1.244, 0.955, 0.971, 0.792, 0.737),
  lower95 = c(0.508, 1.233, 1.087, 1.184, 0.925, 0.839, 0.579, 0.527),
  upper95 = c(0.625, 1.518, 1.611, 1.307, 0.985, 1.123, 1.083, 1.031),
  p_label = c("<0.001","<0.001","0.005","<0.001","0.004","0.688","0.143","0.075")
)
print(expected)
cat("Concordance expected: 0.872\n\n")

# ── Percentage reductions / increases for manuscript text ──
cat("── Percentage effects for manuscript text ──\n")
hr_dist  <- summary_full$conf.int["distance100", "exp(coef)"]
hr_pres  <- summary_full$conf.int["pressure_z",  "exp(coef)"]
hr_temp  <- summary_full$conf.int["tas",          "exp(coef)"]
cat(sprintf("Distance:    %.1f%% reduction per 100km\n", (1 - hr_dist) * 100))
cat(sprintf("Pressure:    %.1f%% increase per SD\n",     (hr_pres - 1) * 100))
cat(sprintf("Temperature: %.1f%% increase per 1C\n",     (hr_temp - 1) * 100))

## ============================================================
## SAVE MODEL OBJECT FOR DOWNSTREAM ANALYSES
## ============================================================

saveRDS(cox_source_full, file = "cox_source_full.rds")
saveRDS(analysis_cases,  file = "analysis_cases_with_cats.rds")
cat("\nSaved: cox_source_full.rds, analysis_cases_with_cats.rds\n")

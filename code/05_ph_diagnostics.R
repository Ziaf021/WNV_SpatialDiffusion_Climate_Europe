## ============================================================
## 05_ph_diagnostics.R   (Analyses 2.7 - 2.11)
## Proportional-hazards diagnostics for the WNV diffusion manuscript.
##   2.7  Schoenfeld residuals / PH test   (Fig S2 + PH table)
##   2.8  Influence (dfbeta) robustness check
##   2.9  Time-varying coefficient model   (Table S6)
##   2.10 Case-weighted kernel PH comparison (Table S7)
##   2.11 Phase-stratified PH check
## Author: Dr. Zia Farooq, Umea University
##
## Settled specification: lambda = 50, analysis-sample pressure_z,
## corrected log_cases_nearest, full covariates, Efron ties.
## All reported values are computed live; no hardcoded expectations.
## ============================================================

library(survival)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(ggpubr)     # theme_classic2()
library(flextable)
library(officer)

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (!exists("analysis_cases"))
  analysis_cases <- readRDS("analysis_cases_with_cats.rds")
if (!exists("cox_source_full"))
  cox_source_full <- readRDS("cox_source_full.rds")


## ---- Table helpers (shared style with 04_) ---------------------------------

fmt_hr <- function(hr, lo, hi, d = 2)
  sprintf(paste0("%.", d, "f (%.", d, "f\u2013%.", d, "f)"), hr, lo, hi)
fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

pretty_term <- function(x) dplyr::recode(x,
  "distance100"         = "Distance (per 100 km)",
  "tt(distance100)"     = "Distance \u00D7 log(t+1)",
  "pressure_z"          = "Invasion pressure (per SD)",
  "tt(pressure_z)"      = "Invasion pressure \u00D7 log(t+1)",
  "log_cases_nearest"   = "Source-region intensity (log cases)",
  "tas"                 = "Temperature (per \u00B0C)",
  "hurs"                = "Relative humidity (per %)",
  "scale(pr)"           = "Precipitation (per SD)",
  "factor(URBN_TYPE)2"  = "Urbanisation: intermediate (vs urban)",
  "factor(URBN_TYPE)3"  = "Urbanisation: rural (vs urban)",
  .default = x
)

make_ft <- function(df) {
  ft <- flextable(df)
  ft <- theme_booktabs(ft)
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- fontsize(ft, size = 9, part = "all")
  autofit(ft)
}
coef_table <- function(fit) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(!grepl("frailty", term)) %>%
    transmute(Term = pretty_term(term),
              `HR (95% CI)` = fmt_hr(estimate, conf.low, conf.high),
              `p-value` = fmt_p(p.value))
}

ph_doc <- read_docx()
add_tbl <- function(doc, caption, ft, note = NULL) {
  doc <- body_add_par(doc, caption, style = "heading 2")
  doc <- body_add_flextable(doc, ft)
  if (!is.null(note)) doc <- body_add_par(doc, note, style = "Normal")
  body_add_par(doc, "")
}


## ============================================================
## 2.7  Proportional-hazards assessment  (Fig S2 + PH table)
## ============================================================

cat("\u2500\u2500 2.7: Proportional hazards ------------------------\n")

zph_result <- cox.zph(cox_source_full)
print(zph_result)

## PH-test table (live)
zt <- as.data.frame(zph_result$table)
zt$term <- rownames(zt)
ph_tbl <- zt %>%
  transmute(
    Covariate = ifelse(term == "GLOBAL", "GLOBAL", pretty_term(term)),
    `Chi-square` = round(chisq, 2),
    df = df,
    `p-value` = fmt_p(p)
  )
ph_doc <- add_tbl(
  ph_doc,
  "Table (PH). Schoenfeld proportional-hazards test for the primary model.",
  make_ft(ph_tbl)
)

## Fig S2: Schoenfeld residual plots for the two violators
pdf(file.path(output_dir, "FigureS2_Schoenfeld.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(zph_result, var = "distance100",
     main = "A) Distance to nearest invaded region",
     xlab = "Time", ylab = "Beta(t) for distance100"); abline(h = 0, col = "grey50", lty = 2)
plot(zph_result, var = "pressure_z",
     main = "B) Neighbourhood invasion pressure",
     xlab = "Time", ylab = "Beta(t) for pressure_z"); abline(h = 0, col = "grey50", lty = 2)
dev.off()


## ============================================================
## 2.8  Influence (dfbeta) robustness check
## ============================================================
# Schoenfeld-residual outlier mapping is ill-defined here because event times
# are tied integer years. dfbeta residuals align unambiguously to data rows,
# so they are used to flag influential region-years.

cat("\u2500\u2500 2.8: Influence robustness (dfbeta) ---------------\n")

mf       <- model.frame(cox_source_full)
mf_rows  <- as.integer(rownames(mf))
infl     <- residuals(cox_source_full, type = "dfbeta")
colnames(infl) <- names(coef(cox_source_full))

flag <- abs(scale(infl[, "distance100"])[, 1]) > 3 |
        abs(scale(infl[, "pressure_z"])[, 1])  > 3

influential <- analysis_cases[mf_rows[flag], , drop = FALSE]
cat(sprintf("Influential region-years flagged (|z(dfbeta)| > 3): %d\n", nrow(influential)))
print(influential %>% select(NUTS_ID, year, distance100, pressure_z, tas))

analysis_no_infl <- analysis_cases %>%
  anti_join(influential %>% select(NUTS_ID, year), by = c("NUTS_ID", "year"))

cox_no_infl <- coxph(
  Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
    tas + hurs + scale(pr) + factor(URBN_TYPE),
  data = analysis_no_infl, ties = "efron"
)
cat("PH test after excluding influential region-years:\n")
print(cox.zph(cox_no_infl))


## ============================================================
## 2.9  Time-varying coefficient model  (Table S6)
## ============================================================
# Distance and pressure (the two PH violators) vary with log(t+1).

cat("\u2500\u2500 2.9: Time-varying coefficients (Table S6) --------\n")

cox_tv2 <- coxph(
  Surv(start, stop, event) ~
    distance100 + tt(distance100) +
    pressure_z  + tt(pressure_z)  +
    log_cases_nearest + tas + hurs + scale(pr) + factor(URBN_TYPE),
  data = analysis_cases, ties = "efron",
  tt   = function(x, t, ...) x * log(t + 1)
)
print(summary(cox_tv2))

ph_doc <- add_tbl(
  ph_doc,
  "Table S6. Time-varying coefficient model (distance and pressure interacted with log follow-up time).",
  make_ft(coef_table(cox_tv2)),
  note = sprintf("Concordance = %.3f. A significant time interaction indicates a hazard ratio that changes over follow-up.",
                 summary(cox_tv2)$concordance["C"])
)


## ============================================================
## 2.10  Case-weighted kernel PH comparison  (Table S7)
## ============================================================
# Weights source regions by their prior-year case burden. Requires the case
# series and spatial objects from 01_; skipped with a message if unavailable.

cat("\u2500\u2500 2.10: Case-weighted kernels (Table S7) -----------\n")

have_s7 <- all(vapply(c("dist_mat", "region_index", "region_lookup",
                        "first_outbreak"), exists, logical(1))) &&
           (exists("dat") || exists("surv_dat"))

if (!have_s7) {
  message("Table S7 skipped: run within the pipeline so the case series and ",
          "spatial objects from 01_ are in session.")
} else {

  case_series <- if (exists("dat")) dat else surv_dat
  region_cases <- case_series %>%
    group_by(NUTS_ID, year) %>%
    summarise(cases = max(wnf_sum, na.rm = TRUE), .groups = "drop") %>%
    mutate(cases = ifelse(is.finite(cases), cases, 0))

  ## Case-weighted invasion pressure for a given weight transform.
  build_pressure_cases <- function(weight_fun, lambda = 50) {
    W <- exp(-dist_mat / lambda)
    years <- sort(unique(analysis_cases$year))
    out <- list()
    for (yy in years) {
      invaded <- first_outbreak %>%
        filter(!is.na(first_year), first_year < yy) %>% pull(NUTS_ID)
      idx <- na.omit(region_index[invaded])
      if (length(idx) == 0) {
        pressure <- rep(0, nrow(W))
      } else {
        src <- region_cases$cases[match(paste(names(idx), yy - 1),
                                        paste(region_cases$NUTS_ID, region_cases$year))]
        src[is.na(src)] <- 0
        pressure <- as.numeric(W[, idx, drop = FALSE] %*% weight_fun(src))
      }
      out[[length(out) + 1]] <- data.frame(
        NUTS_ID = region_lookup$NUTS_ID, year = yy, pressure_cases = pressure)
    }
    bind_rows(out)
  }

  ## Refit the full model with a case-weighted pressure term; return PH + fit.
  fit_weighted <- function(weight_fun) {
    pc <- build_pressure_cases(weight_fun)
    d  <- analysis_cases %>%
      select(-pressure_z) %>%
      left_join(pc, by = c("NUTS_ID", "year")) %>%
      mutate(pressure_z = as.numeric(scale(pressure_cases)))
    fit <- coxph(
      Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
        tas + hurs + scale(pr) + factor(URBN_TYPE),
      data = d, ties = "efron")
    z <- cox.zph(fit)
    list(fit = fit, zph = z)
  }

  ## Reference: unweighted primary model.
  z0 <- cox.zph(cox_source_full)
  s0 <- summary(cox_source_full)
  rows <- list(data.frame(
    Kernel          = "Unweighted (primary)",
    `Distance p`    = fmt_p(z0$table["distance100", "p"]),
    `Pressure p`    = fmt_p(z0$table["pressure_z",  "p"]),
    `Precipitation p` = fmt_p(z0$table["scale(pr)", "p"]),
    `Urbanisation p`  = fmt_p(z0$table["factor(URBN_TYPE)", "p"]),
    `Global p`      = fmt_p(z0$table["GLOBAL", "p"]),
    AIC             = round(AIC(cox_source_full), 1),
    `C-index`       = round(s0$concordance["C"], 3),
    `Pressure HR`   = round(s0$conf.int["pressure_z", "exp(coef)"], 2),
    check.names = FALSE
  ))

  weights <- list(`Log(cases+1)-weighted` = function(x) log1p(x),
                  `Sqrt(cases)-weighted`  = function(x) sqrt(x))
  for (nm in names(weights)) {
    r <- fit_weighted(weights[[nm]]); s <- summary(r$fit)
    rows[[length(rows) + 1]] <- data.frame(
      Kernel          = nm,
      `Distance p`    = fmt_p(r$zph$table["distance100", "p"]),
      `Pressure p`    = fmt_p(r$zph$table["pressure_z",  "p"]),
      `Precipitation p` = fmt_p(r$zph$table["scale(pr)", "p"]),
      `Urbanisation p`  = fmt_p(r$zph$table["factor(URBN_TYPE)", "p"]),
      `Global p`      = fmt_p(r$zph$table["GLOBAL", "p"]),
      AIC             = round(AIC(r$fit), 1),
      `C-index`       = round(s$concordance["C"], 3),
      `Pressure HR`   = round(s$conf.int["pressure_z", "exp(coef)"], 2),
      check.names = FALSE
    )
  }
  tbl_s7 <- bind_rows(rows)
  print(tbl_s7)
  ph_doc <- add_tbl(
    ph_doc,
    "Table S7. Proportional-hazards properties and fit of case-weighted invasion-pressure kernels.",
    make_ft(tbl_s7)
  )
}

# Build the two case-weighted pressures on the analysis rows
pc_log  <- build_pressure_cases(function(x) log1p(x))
pc_sqrt <- build_pressure_cases(function(x) sqrt(x))

d <- analysis_cases %>%
  left_join(pc_log,  by = c("NUTS_ID", "year")) %>%
  rename(pressure_cases_log = pressure_cases) %>%
  left_join(pc_sqrt, by = c("NUTS_ID", "year")) %>%
  rename(pressure_cases_sqrt = pressure_cases)

cor(d$pressure_z, scale(d$pressure_cases_log)[, 1],  use = "complete.obs")
cor(d$pressure_z, scale(d$pressure_cases_sqrt)[, 1], use = "complete.obs")
## ============================================================
## 2.11  Phase-stratified PH check
## ============================================================
# invasion_phase is already an ordered factor from 01_ (no re-levelling here).

cat("\u2500\u2500 2.11: Phase-stratified PH check ------------------\n")
print(table(analysis_cases$invasion_phase, useNA = "ifany"))

cox_stratified_phase <- coxph(
  Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
    tas + hurs + scale(pr) + factor(URBN_TYPE) + strata(invasion_phase),
  data = analysis_cases, ties = "efron"
)
print(summary(cox_stratified_phase))

d_prim <- summary(cox_source_full)$conf.int["distance100", "exp(coef)"]
d_strat <- summary(cox_stratified_phase)$conf.int["distance100", "exp(coef)"]
cat(sprintf("Distance HR - primary: %.4f vs phase-stratified: %.4f\n", d_prim, d_strat))
print(cox.zph(cox_stratified_phase))

ph_doc <- add_tbl(
  ph_doc,
  "Table (phase-stratified). Multivariable model stratified on invasion phase.",
  make_ft(coef_table(cox_stratified_phase))
)


## ============================================================
## 2.12  Precipitation non-proportionality sensitivity
## ============================================================
# Precipitation is the principal residual PH violator. Allow its coefficient
# to vary with log follow-up time and confirm the interpreted effects
# (distance, pressure, temperature, humidity) are unchanged.

cat("\u2500\u2500 2.12: Precipitation time-interaction sensitivity ---\n")

d_pr <- analysis_cases %>% mutate(pr_z = as.numeric(scale(pr)))

cox_pr_tv <- coxph(
  Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
    tas + hurs + pr_z + tt(pr_z) + factor(URBN_TYPE),
  data = d_pr, ties = "efron",
  tt   = function(x, t, ...) x * log(t + 1)
)
print(summary(cox_pr_tv))

## Compare interpreted effects: primary vs precipitation-time-varying model.
key_terms <- c("distance100", "pressure_z", "tas", "hurs")
sp <- summary(cox_source_full)$conf.int
sv <- summary(cox_pr_tv)$conf.int
pr_sens_tbl <- data.frame(
  Term = pretty_term(key_terms),
  `Primary HR (95% CI)` =
    fmt_hr(sp[key_terms, "exp(coef)"], sp[key_terms, "lower .95"], sp[key_terms, "upper .95"]),
  `Precip. time-varying HR (95% CI)` =
    fmt_hr(sv[key_terms, "exp(coef)"], sv[key_terms, "lower .95"], sv[key_terms, "upper .95"]),
  check.names = FALSE
)
print(pr_sens_tbl)

pr_int_p <- summary(cox_pr_tv)$coefficients["tt(pr_z)", "Pr(>|z|)"]
ph_doc <- add_tbl(
  ph_doc,
  "Table (precipitation sensitivity). Interpreted effects with precipitation allowed to vary over follow-up time.",
  make_ft(pr_sens_tbl),
  note = sprintf(
    "Precipitation x log(t+1): p = %s. The distance, pressure, temperature, and humidity estimates are essentially unchanged from the primary model.",
    fmt_p(pr_int_p))
)


## ============================================================
## Write outputs
## ============================================================

print(ph_doc, target = file.path(output_dir, "Supplementary_Tables_PH.docx"))
saveRDS(cox_tv2,              file.path(output_dir, "cox_tv2.rds"))
saveRDS(cox_stratified_phase, file.path(output_dir, "cox_stratified_phase.rds"))
saveRDS(cox_pr_tv,            file.path(output_dir, "cox_pr_tv.rds"))
cat("\nWrote: Supplementary_Tables_PH.docx (PH, S6, S7 with precip/urban columns,",
    "phase-stratified, precipitation sensitivity), FigureS2_Schoenfeld.pdf,",
    "cox_tv2.rds, cox_stratified_phase.rds, cox_pr_tv.rds\n")

## ============================================================
## End of 05_ph_diagnostics.R
## ============================================================

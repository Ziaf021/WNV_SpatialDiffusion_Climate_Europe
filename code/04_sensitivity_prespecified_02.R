# =============================================================================
# 04_sensitivity_prespecified.R
#
# Pre-specified sensitivity analyses for:
#   "Spatial diffusion and climatic suitability shape the West Nile Virus
#    invasion risk across Europe"
#
# Produces MS-Word-ready tables (one consolidated Supplementary_Tables.docx)
# and the associated supplementary figures:
#   Table S1 / Figure S1  Decay-parameter (lambda) sensitivity     [Section 1]
#   Table S2              Alternative spatial-kernel sensitivity    [Section 2]
#   Table S3              Country-level shared-frailty model        [Section 3]
#   Table S4              Climate-anomaly model + AIC comparison    [Section 4]
#   Table S5 / Figure     Leave-one-country-out (LOCO) analysis     [Section 5]
#   Table 1  (optional)   Primary multivariable model (main text)
#
# Settled specification (see 01_/02_ and the lambda profile below):
#   lambda_primary = 50 (AIC-minimising); pressure standardised on the
#   analysis sample; log_cases_nearest from the full case series.
#
# Upstream objects expected in the session (from 01_/02_; run those first):
#   analysis_cases, dist_mat, region_lookup, region_index, first_outbreak
#   cox_source_full  (from 02_; read from disk if not in session)
#   NUTS3ObsClimate_1991to2024  (read from climate_file if not in session)
# =============================================================================


## ---- Libraries --------------------------------------------------------------

library(survival)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(ggpubr)      # theme_classic2()
library(patchwork)
library(flextable)
library(officer)


## ---- Configuration ----------------------------------------------------------

data_dir     <- "data"
output_dir   <- "output"
climate_file <- file.path(data_dir, "NUTS3ObsClimate_1991to2024.csv")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Fixed decay length for the kernel-form comparison (Section 2). This equals
# the AIC-minimising primary value, so the comparison is on the primary scale.
lambda_kernel <- 50

# Primary model: use the in-session object if present, else read from disk.
cox_source_full <- if (exists("cox_source_full")) {
  cox_source_full
} else {
  readRDS("cox_source_full.rds")
}


## ---- Shared table helpers ---------------------------------------------------
# One formatter for every table, so all supplementary tables are consistent.

# Hazard ratio with 95% CI. En-dash (range), not em-dash; change if preferred.
fmt_hr <- function(hr, lo, hi, d = 2)
  sprintf(paste0("%.", d, "f (%.", d, "f\u2013%.", d, "f)"), hr, lo, hi)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# Human-readable covariate labels for coefficient tables.
pretty_term <- function(x) dplyr::recode(x,
  "distance100"        = "Distance (per 100 km)",
  "pressure_z"         = "Invasion pressure (per SD)",
  "log_cases_nearest"  = "Source-region intensity (log cases)",
  "tas"                = "Temperature (per \u00B0C)",
  "hurs"               = "Relative humidity (per %)",
  "pr_z"          = "Precipitation (per SD)",
  "tas_anom"           = "Temperature anomaly (per \u00B0C)",
  "hurs_anom"          = "Humidity anomaly (per \u00B0C)",
  "scale(pr_anom)"     = "Precipitation anomaly (per SD)",
  "factor(URBN_TYPE)2" = "Urbanisation: intermediate (vs urban)",
  "factor(URBN_TYPE)3" = "Urbanisation: rural (vs urban)",
  .default = x
)

# Data frame -> formatted flextable.
make_ft <- function(df) {
  ft <- flextable(df)
  ft <- theme_booktabs(ft)
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- padding(ft, padding = 3, part = "all")
  autofit(ft)
}

# Coefficient table (HR, 95% CI, p) from a fitted coxph model.
coef_table <- function(fit, drop_frailty = TRUE) {
  tb <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  if (drop_frailty) tb <- dplyr::filter(tb, !grepl("frailty", term))
  tb %>%
    transmute(
      Term          = pretty_term(term),
      `HR (95% CI)` = fmt_hr(estimate, conf.low, conf.high),
      `p-value`     = fmt_p(p.value)
    )
}

# Accumulate tables into one Supplementary document.
supp_doc <- read_docx()
add_supp <- function(doc, caption, ft, note = NULL) {
  doc <- body_add_par(doc, caption, style = "heading 2")
  doc <- body_add_flextable(doc, ft)
  if (!is.null(note)) doc <- body_add_par(doc, note, style = "Normal")
  body_add_par(doc, "")
}


# =============================================================================
# Section 1. Decay-parameter (lambda) sensitivity        -> Table S1 / Figure S1
# =============================================================================

lambda_values <- c(25, 50, 75, 100, 150, 200, 250, 500)

## Year-by-region invasion pressure for a given lambda (exponential kernel).
## NOTE: this kernel is untruncated, whereas the PRIMARY pressure in 01_ is
## truncated at 500 km. The difference is immaterial at lambda = 50
## (exp(-500/50) ~ 5e-5) but affects the high-lambda rows below. Add
## `W[dist_mat > 500] <- 0; diag(W) <- 0` here to match 01_ exactly.
build_pressure <- function(lambda) {
  W <- exp(-dist_mat / lambda)
  years <- sort(unique(analysis_cases$year))
  pressure_list <- list()
  for (yy in years) {
    invaded_regions <- first_outbreak %>%
      filter(!is.na(first_year), first_year < yy) %>%
      pull(NUTS_ID)
    invaded_idx <- na.omit(region_index[invaded_regions])
    if (length(invaded_idx) == 0) {
      pressure <- rep(0, nrow(W))
    } else {
      pressure <- rowSums(W[, invaded_idx, drop = FALSE])
    }
    pressure_list[[length(pressure_list) + 1]] <- data.frame(
      NUTS_ID = region_lookup$NUTS_ID, year = yy, invasion_pressure = pressure
    )
  }
  bind_rows(pressure_list)
}

lambda_results <- list()
for (lambda in lambda_values) {
  cat("Running lambda =", lambda, "\n")
  pressure_df <- build_pressure(lambda)
  tmp_dat <- analysis_cases %>%
    select(-pressure_z, -invasion_pressure) %>%
    left_join(pressure_df, by = c("NUTS_ID", "year"))
  tmp_dat$pressure_z <- as.numeric(scale(tmp_dat$invasion_pressure))  # analysis-sample

  fit <- coxph(
    Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
      tas + hurs + pr_z + factor(URBN_TYPE),
    data = tmp_dat, ties = "efron"
  )
  s <- summary(fit)
  lambda_results[[as.character(lambda)]] <- data.frame(
    lambda      = lambda,
    HR          = s$conf.int["pressure_z", "exp(coef)"],
    lower       = s$conf.int["pressure_z", "lower .95"],
    upper       = s$conf.int["pressure_z", "upper .95"],
    p           = s$coefficients["pressure_z", "Pr(>|z|)"],
    AIC         = AIC(fit),
    concordance = s$concordance[1]
  )
}

lambda_df <- bind_rows(lambda_results)   # keep the data frame intact for figures

## Word-ready Table S1
tbl_s1 <- lambda_df %>%
  transmute(
    `\u03BB (km)`   = lambda,
    `HR (95% CI)`   = fmt_hr(HR, lower, upper),
    `p-value`       = fmt_p(p),
    AIC             = round(AIC, 1),
    `C-index`       = round(concordance, 3)
  )
ft_s1 <- make_ft(tbl_s1)
supp_doc <- add_supp(
  supp_doc,
  "Table S1. Sensitivity of the invasion-pressure hazard ratio to the distance-decay length (\u03BB).",
  ft_s1
)

## ---- Figure S1: three-panel lambda profile ---------------------------------
P1 <- ggplot(lambda_df, aes(lambda, HR)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#FB6A4A", alpha = 0.25) +
  geom_line(linewidth = 1, colour = "#CB181D") +
  geom_point(size = 3, colour = "#CB181D") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_vline(xintercept = 50, linetype = 2, colour = "black") +
  annotate("text", x = 50, y = max(lambda_df$upper) * 0.98,
           label = expression(lambda == 50), size = 3) +
  theme_classic2(base_size = 12) +
  labs(x = expression(lambda ~ "(km)"), y = "Hazard ratio for\ninvasion pressure")

best_lambda <- lambda_df$lambda[which.min(lambda_df$AIC)]
best_AIC    <- min(lambda_df$AIC)
P2 <- ggplot(lambda_df, aes(lambda, AIC)) +
  geom_line(linewidth = 1, colour = "#2171B5") +
  geom_point(size = 3, colour = "#2171B5") +
  geom_vline(xintercept = 50, linetype = 2, colour = "black") +
  geom_point(aes(best_lambda, best_AIC), colour = "black", size = 4) +
  annotate("text", x = best_lambda + 30, y = best_AIC + 8,
           label = paste0("Minimum AIC\n(", best_lambda, " km)"), size = 3) +
  theme_classic2(base_size = 12) +
  labs(x = expression(lambda ~ "(km)"), y = "AIC")

P3 <- ggplot(lambda_df, aes(lambda, concordance)) +
  geom_line(linewidth = 1, colour = "#238B45") +
  geom_point(size = 3, colour = "#238B45") +
  geom_vline(xintercept = 50, linetype = 2, colour = "black") +
  annotate("text", x = 50 + 18, y = max(lambda_df$concordance) - 0.001,
           label = expression(lambda == 50), size = 3) +
  scale_y_continuous(
    limits = c(floor(min(lambda_df$concordance) * 100) / 100,
               ceiling(max(lambda_df$concordance) * 100) / 100),
    expand = c(0.02, 0)) +
  theme_classic2(base_size = 12) +
  labs(x = expression(lambda ~ "(km)"), y = "Concordance")

Fig_lambda_Sensitivity <- (P1 | P2 | P3) +
  plot_layout(widths = c(1.2, 1, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(output_dir, "FigureS1_lambda_profile.pdf"),
       Fig_lambda_Sensitivity, width = 11, height = 4)


# =============================================================================
# Section 2. Alternative spatial kernels                          -> Table S2
# =============================================================================

build_pressure_kernel <- function(kernel, lambda) {
  if      (kernel == "exp")      W <- exp(-dist_mat / lambda)
  else if (kernel == "gaussian") W <- exp(-(dist_mat^2) / (2 * lambda^2))
  else if (kernel == "inverse")  W <- 1 / (dist_mat + 1)
  else stop("Unknown kernel: ", kernel)

  years <- sort(unique(analysis_cases$year))
  pressure_list <- list()
  for (yy in years) {
    invaded_regions <- first_outbreak %>%
      filter(!is.na(first_year), first_year < yy) %>%
      pull(NUTS_ID)
    invaded_idx <- na.omit(region_index[invaded_regions])
    if (length(invaded_idx) == 0) {
      pressure <- rep(0, nrow(W))
    } else {
      pressure <- rowSums(W[, invaded_idx, drop = FALSE])
    }
    pressure_list[[length(pressure_list) + 1]] <- data.frame(
      NUTS_ID = region_lookup$NUTS_ID, year = yy, invasion_pressure = pressure
    )
  }
  bind_rows(pressure_list)
}

kernel_results <- list()
for (k in c("exp", "gaussian", "inverse")) {
  pressure_df <- build_pressure_kernel(kernel = k, lambda = lambda_kernel)
  tmp_dat <- analysis_cases %>%
    select(-pressure_z, -invasion_pressure) %>%
    left_join(pressure_df, by = c("NUTS_ID", "year"))
  tmp_dat$pressure_z <- as.numeric(scale(tmp_dat$invasion_pressure))

  fit <- coxph(
    Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
      tas + hurs + pr_z + factor(URBN_TYPE),
    data = tmp_dat, ties = "efron"
  )
  s <- summary(fit)
  kernel_results[[k]] <- data.frame(
    kernel = k,
    HR     = s$conf.int["pressure_z", "exp(coef)"],
    lower  = s$conf.int["pressure_z", "lower .95"],
    upper  = s$conf.int["pressure_z", "upper .95"],
    p      = s$coefficients["pressure_z", "Pr(>|z|)"],
    AIC    = AIC(fit),
    C      = s$concordance[1]
  )
}

kernel_df <- bind_rows(kernel_results)

tbl_s2 <- kernel_df %>%
  transmute(
    Kernel        = recode(kernel, exp = "Exponential",
                           gaussian = "Gaussian", inverse = "Inverse-distance"),
    `HR (95% CI)` = fmt_hr(HR, lower, upper),
    `p-value`     = fmt_p(p),
    AIC           = round(AIC, 1),
    `C-index`     = round(C, 3)
  )
ft_s2 <- make_ft(tbl_s2)
supp_doc <- add_supp(
  supp_doc,
  paste0("Table S2. Sensitivity to the functional form of the spatial-",
         "connectivity kernel (\u03BB = ", lambda_kernel, " km)."),
  ft_s2
)

Fig_kernel_SI <- kernel_df %>%
  select(kernel, HR, AIC, C) %>%
  pivot_longer(-kernel, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(kernel, Value, fill = kernel)) +
  geom_col(width = 0.7) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  theme_classic2() + theme(legend.position = "none")
ggsave(file.path(output_dir, "FigureS2_kernel_comparison.pdf"),
       Fig_kernel_SI, width = 9, height = 3.5)


# =============================================================================
# Section 3. Country-level shared-frailty model                   -> Table S3
# =============================================================================

cox_frailty <- coxph(
  Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
    tas + hurs + pr_z + factor(URBN_TYPE) + frailty(CNTR_CODE),
  data = analysis_cases
)
print(summary(cox_frailty))

frailty_var <- tryCatch(cox_frailty$history[[1]]$theta, error = function(e) NA_real_)
n_countries <- dplyr::n_distinct(analysis_cases$CNTR_CODE)

ft_s3 <- make_ft(coef_table(cox_frailty))
supp_doc <- add_supp(
  supp_doc,
  "Table S3. Multivariable model with a gamma shared frailty at country level.",
  ft_s3,
  note = sprintf(
    "Shared frailty over %d countries%s. Fixed-effect hazard ratios shown.",
    n_countries,
    if (is.na(frailty_var)) "" else sprintf("; estimated frailty variance = %.3f", frailty_var)
  )
)

AIC(cox_source_full, cox_frailty)

# =============================================================================
# Section 4. Climate-anomaly analysis                             -> Table S4
# =============================================================================

if (!exists("NUTS3ObsClimate_1991to2024"))
  NUTS3ObsClimate_1991to2024 <- read.csv(climate_file, stringsAsFactors = FALSE)

## Baseline climatology (observed climate, 2000-2006)
climatology <- NUTS3ObsClimate_1991to2024 %>%
  filter(scenario == "obsclim", year >= 2000, year <= 2006) %>%
  group_by(NUTS_ID) %>%
  summarise(tas_mean = mean(tas, na.rm = TRUE),
            hurs_mean = mean(hurs, na.rm = TRUE),
            pr_mean = mean(pr, na.rm = TRUE), .groups = "drop")

climate_anom <- NUTS3ObsClimate_1991to2024 %>%
  left_join(climatology, by = "NUTS_ID") %>%
  mutate(tas_anom = tas - tas_mean,
         hurs_anom = hurs - hurs_mean,
         pr_anom = pr - pr_mean)

analysis_anom <- analysis_cases %>%
  select(-tas, -hurs, -pr) %>%
  left_join(climate_anom %>% select(NUTS_ID, year, tas_anom, hurs_anom, pr_anom),
            by = c("NUTS_ID", "year"))

cox_anom <- coxph(
  Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
    tas_anom + hurs_anom + scale(pr_anom) + factor(URBN_TYPE),
  data = analysis_anom, ties = "efron"
)
print(summary(cox_anom))

# AIC comparison is valid only on identical observations.
if (cox_source_full$n != cox_anom$n)
  warning("Primary and anomaly models fitted on different n; AIC not comparable.")
aic_primary <- AIC(cox_source_full)
aic_anom    <- AIC(cox_anom)
delta_aic   <- aic_anom - aic_primary
cat(sprintf("\nAIC primary = %.1f; AIC anomaly = %.1f; delta = %.1f\n",
            aic_primary, aic_anom, delta_aic))

ft_s4a <- make_ft(coef_table(cox_anom))
tbl_s4b <- data.frame(
  Model             = c("Primary (absolute climate)", "Climate-anomaly"),
  AIC               = round(c(aic_primary, aic_anom), 1),
  `\u0394AIC vs primary` = round(c(0, delta_aic), 1),
  check.names = FALSE
)
ft_s4b <- make_ft(tbl_s4b)

supp_doc <- add_supp(
  supp_doc,
  "Table S4. Climate-anomaly specification (deviation from 2000-2006 baseline).",
  ft_s4a
)
supp_doc <- add_supp(
  supp_doc,
  "Table S4 (continued). Model-fit comparison, absolute climate versus anomalies.",
  ft_s4b
)


# =============================================================================
# Section 5. Leave-one-country-out (LOCO)                -> Table S5 / Figure
# =============================================================================
# Single full-model implementation, matching the primary specification.

countries    <- sort(unique(analysis_cases$CNTR_CODE))
loco_results <- data.frame()

for (ctry in countries) {
  cat(sprintf("  Excluding %s...\n", ctry))
  dat_excl <- analysis_cases %>% filter(CNTR_CODE != ctry)
  fit_loco <- tryCatch(
    coxph(
      Surv(start, stop, event) ~ distance100 + pressure_z + log_cases_nearest +
        tas + hurs + pr_z + factor(URBN_TYPE),
      data = dat_excl, ties = "efron"
    ), error = function(e) NULL)
  if (is.null(fit_loco)) next
  s <- summary(fit_loco)
  loco_results <- bind_rows(loco_results, data.frame(
    country_excluded = ctry,
    HR_distance = s$conf.int["distance100", "exp(coef)"],
    lo_distance = s$conf.int["distance100", "lower .95"],
    hi_distance = s$conf.int["distance100", "upper .95"],
    HR_pressure = s$conf.int["pressure_z",  "exp(coef)"],
    lo_pressure = s$conf.int["pressure_z",  "lower .95"],
    hi_pressure = s$conf.int["pressure_z",  "upper .95"],
    HR_temp     = s$conf.int["tas",          "exp(coef)"],
    lo_temp     = s$conf.int["tas",          "lower .95"],
    hi_temp     = s$conf.int["tas",          "upper .95"]
  ))
}

# Primary (full-dataset) reference row.
s_full <- summary(cox_source_full)
loco_results <- bind_rows(loco_results, data.frame(
  country_excluded = "ALL (primary)",
  HR_distance = s_full$conf.int["distance100", "exp(coef)"],
  lo_distance = s_full$conf.int["distance100", "lower .95"],
  hi_distance = s_full$conf.int["distance100", "upper .95"],
  HR_pressure = s_full$conf.int["pressure_z",  "exp(coef)"],
  lo_pressure = s_full$conf.int["pressure_z",  "lower .95"],
  hi_pressure = s_full$conf.int["pressure_z",  "upper .95"],
  HR_temp     = s_full$conf.int["tas",          "exp(coef)"],
  lo_temp     = s_full$conf.int["tas",          "lower .95"],
  hi_temp     = s_full$conf.int["tas",          "upper .95"]
))

## Word-ready Table S5
tbl_s5 <- loco_results %>%
  transmute(
    `Country excluded`        = country_excluded,
    `Distance HR (95% CI)`    = fmt_hr(HR_distance, lo_distance, hi_distance),
    `Pressure HR (95% CI)`    = fmt_hr(HR_pressure, lo_pressure, hi_pressure),
    `Temperature HR (95% CI)` = fmt_hr(HR_temp, lo_temp, hi_temp)
  )
ft_s5 <- make_ft(tbl_s5)
supp_doc <- add_supp(
  supp_doc,
  "Table S5. Leave-one-country-out sensitivity of the principal hazard ratios.",
  ft_s5,
  note = "Final row: primary full-dataset estimates for reference."
)

## LOCO dot plot
loco_long <- loco_results %>%
  pivot_longer(c(HR_distance, HR_pressure, HR_temp),
               names_to = "metric", values_to = "HR") %>%
  mutate(name_label = factor(recode(metric,
           HR_distance = "Distance", HR_pressure = "Invasion pressure",
           HR_temp = "Temperature"),
           levels = c("Distance", "Invasion pressure", "Temperature")))

mean_hrs <- loco_long %>%
  filter(country_excluded != "ALL (primary)") %>%
  group_by(name_label) %>%
  summarise(mean_HR = mean(HR, na.rm = TRUE), .groups = "drop")
primary_hrs <- loco_long %>%
  filter(country_excluded == "ALL (primary)") %>%
  transmute(name_label, primary_HR = HR)

Fig_LOCO <- ggplot(
  loco_long %>% filter(country_excluded != "ALL (primary)"),
  aes(HR, reorder(country_excluded, HR), colour = name_label)
) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_vline(data = mean_hrs, aes(xintercept = mean_HR, colour = name_label),
             linetype = "solid", linewidth = 1.0, alpha = 0.7) +
  geom_vline(data = primary_hrs, aes(xintercept = primary_HR, colour = name_label),
             linetype = "dashed", linewidth = 0.7, alpha = 0.5) +
  facet_wrap(~ name_label, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c("Distance" = "#7B3294",
    "Invasion pressure" = "#1B7837", "Temperature" = "#D95F0E"), guide = "none") +
  labs(x = "Hazard ratio", y = "Country excluded",
       caption = "Solid: mean HR across exclusions  |  Dashed: primary full-dataset HR") +
  theme_classic2(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 11),
        plot.caption = element_text(size = 8, colour = "grey50", hjust = 0))
ggsave(file.path(output_dir, "FigureS_LOCO_DotPlot.pdf"),
       Fig_LOCO, width = 10, height = 6)


# =============================================================================
# Optional: main-text Table 1 (primary multivariable model)
# =============================================================================

ft_main <- make_ft(coef_table(cox_source_full))
main_doc <- read_docx()
main_doc <- body_add_par(main_doc,
  "Table 1. Multivariable Cox model of first regional WNV invasion.",
  style = "heading 2")
main_doc <- body_add_flextable(main_doc, ft_main)
main_doc <- body_add_par(main_doc,
  sprintf("n = %d region-years; %d events; concordance = %.3f.",
          cox_source_full$n, cox_source_full$nevent, s_full$concordance[1]),
  style = "Normal")
print(main_doc, target = file.path(output_dir, "Table1_PrimaryModel.docx"))


# =============================================================================
# Write outputs
# =============================================================================

print(supp_doc, target = file.path(output_dir, "Supplementary_Tables.docx"))

saveRDS(loco_results, file.path(output_dir, "loco_results.rds"))
saveRDS(cox_anom,     file.path(output_dir, "cox_anom.rds"))
saveRDS(cox_frailty,  file.path(output_dir, "cox_frailty.rds"))

cat("\nWrote:\n",
    " ", file.path(output_dir, "Supplementary_Tables.docx"), "\n",
    " ", file.path(output_dir, "Table1_PrimaryModel.docx"), "\n",
    " supplementary figures (S1, S2, LOCO)\n")

# =============================================================================
# End of 04_sensitivity_prespecified.R
# =============================================================================

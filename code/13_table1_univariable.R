###############################################################################
# 13_table1_univariable.R
#
# Regenerates Table 1 of TLRHEUROPE-D-26-01167 on the refitted cohort.
#
# Table 1 in the submitted manuscript is a set of UNIVARIABLE Cox models, one
# covariate at a time, not the primary multivariable model. The v2 run already
# writes univariable_models_full_follow_up_climatefix.csv, but that file keeps
# only the first coefficient of each fit, so the urbanisation term appears as a
# single row and the rural-versus-urban contrast is lost. This script refits the
# same seven models and keeps every coefficient, then formats the result in the
# layout and the numeral style of the submitted table.
#
# Model specification is identical to 02_survival_analyses.R: counting-process
# Surv(start, stop, event), Efron ties, no strata, no robust variance.
#
# Run from the project folder:
#   Rscript code/13_table1_univariable.R
#
# Writes to run_withBalkans_v2/output/revision/
#   rev_table1_univariable.csv   machine readable, plain decimal points
#   rev_table1_univariable.md    the table as it should appear in the paper
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

options(stringsAsFactors = FALSE)

ROOT <- normalizePath(getwd())
RUN  <- file.path(ROOT, Sys.getenv("WNV_RUN_DIR", "run_withBalkans_v2"))
OUT  <- file.path(RUN, "output")
REV  <- file.path(OUT, "revision")
dir.create(REV, recursive = TRUE, showWarnings = FALSE)

ac_path <- file.path(OUT, "analysis_cases_withBalkans.rds")
if (!file.exists(ac_path)) {
  stop("Not found: ", ac_path, "\nRun run_revision_v2.sh first.", call. = FALSE)
}

analysis_cases <- readRDS(ac_path)

need <- c("start", "stop", "event", "distance100", "pressure_z",
          "log_cases_nearest", "tas", "hurs", "pr", "URBN_TYPE")
absent <- setdiff(need, names(analysis_cases))
if (length(absent)) {
  stop("Missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
}

n_ry     <- nrow(analysis_cases)
n_events <- sum(analysis_cases$event)
n_reg    <- dplyr::n_distinct(analysis_cases$NUTS_ID)
n_cnt    <- dplyr::n_distinct(analysis_cases$CNTR_CODE)
yr_min   <- min(analysis_cases$year)
yr_max   <- max(analysis_cases$year)

message("analysis_cases: ", n_ry, " region-years, ", n_events, " events, ",
        n_reg, " regions, ", n_cnt, " countries, years ", yr_min, " to ", yr_max)

if (n_events != as.integer(Sys.getenv("WNV_EXPECTED_EVENTS", "306"))) {
  message("NOTE: events = ", n_events,
          ", the 10 September handoff expects 306. Check WNV_RUN_DIR.")
}

## ---- fit ---------------------------------------------------------------------

fit_univar <- function(formula_str, terms_wanted) {
  fit <- coxph(as.formula(formula_str), data = analysis_cases, ties = "efron")
  s   <- summary(fit)
  cf  <- s$coefficients
  ci  <- s$conf.int
  out <- data.frame(
    term    = rownames(cf),
    HR      = unname(cf[, "exp(coef)"]),
    lower95 = unname(ci[, "lower .95"]),
    upper95 = unname(ci[, "upper .95"]),
    p       = unname(cf[, "Pr(>|z|)"])
  )
  out[match(terms_wanted, out$term), , drop = FALSE]
}

spec <- list(
  list(f = "Surv(start,stop,event) ~ distance100",
       t = "distance100",
       covariate = "Distance to nearest invaded region",
       unit = "per 100 km increase"),
  list(f = "Surv(start,stop,event) ~ pressure_z",
       t = "pressure_z",
       covariate = "Neighbourhood invasion pressure",
       unit = "per SD increase"),
  list(f = "Surv(start,stop,event) ~ log_cases_nearest",
       t = "log_cases_nearest",
       covariate = "Source-region transmission intensity",
       unit = "per unit increase in log(cases+1)"),
  list(f = "Surv(start,stop,event) ~ tas",
       t = "tas",
       covariate = "Annual mean temperature",
       unit = "per 1 degree C increase"),
  list(f = "Surv(start,stop,event) ~ pr_z",
       t = "pr_z",
       covariate = "Annual precipitation",
       unit = "per SD increase"),
  list(f = "Surv(start,stop,event) ~ hurs",
       t = "hurs",
       covariate = "Relative humidity",
       unit = "per 1 pct increase"),
  list(f = "Surv(start,stop,event) ~ factor(URBN_TYPE)",
       t = c("factor(URBN_TYPE)2", "factor(URBN_TYPE)3"),
       covariate = c("Intermediate urbanisation", "Rural areas"),
       unit = c("vs urban", "vs urban"))
)

tab <- bind_rows(lapply(spec, function(s) {
  res <- fit_univar(s$f, s$t)
  data.frame(covariate = s$covariate, unit = s$unit,
             HR = res$HR, lower95 = res$lower95,
             upper95 = res$upper95, p = res$p)
}))

## ---- format ------------------------------------------------------------------

dec <- function(x, digits) sub(".", "·", formatC(x, format = "f", digits = digits), fixed = TRUE)

fmt_p <- function(p) {
  ifelse(p < 0.0001, "<0·0001",
         ifelse(p < 0.001, dec(p, 5),
                ifelse(p < 0.01, dec(p, 4), dec(p, 3))))
}

tab_fmt <- tab %>%
  mutate(
    `HR (95% CI)` = paste0(dec(HR, 2), " (", dec(lower95, 2), "–", dec(upper95, 2), ")"),
    `p value`     = fmt_p(p)
  ) %>%
  select(Covariate = covariate, `Unit / comparison` = unit,
         `HR (95% CI)`, `p value`)

print(tab, digits = 4)
cat("\n")
print(tab_fmt, right = FALSE)

## ---- write -------------------------------------------------------------------

write.csv(tab, file.path(REV, "rev_table1_univariable.csv"), row.names = FALSE)

caption <- sprintf(
  paste0("Table 1: Univariable Cox proportional-hazards analyses of factors associated ",
         "with first human WNV invasion across European NUTS3 regions, %d–%d"),
  yr_min, yr_max %% 100)

md <- c(
  caption, "",
  paste0("| ", paste(names(tab_fmt), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(tab_fmt)), collapse = "|"), "|"),
  apply(tab_fmt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")),
  "",
  sprintf(paste0("Each row is a separate model with a single covariate. Cohort: %s region-years, ",
                 "%s regions in %s countries, %s first invasions. Regions are followed from 2006 ",
                 "and the models are fitted from %d, because distance to the nearest invaded region ",
                 "is undefined in the first year of follow-up. Efron handling of ties."),
          format(n_ry, big.mark = " "), format(n_reg, big.mark = " "),
          n_cnt, format(n_events, big.mark = " "), yr_min)
)

writeLines(md, file.path(REV, "rev_table1_univariable.md"))

message("\nWritten:")
message("  ", file.path(REV, "rev_table1_univariable.csv"))
message("  ", file.path(REV, "rev_table1_univariable.md"))

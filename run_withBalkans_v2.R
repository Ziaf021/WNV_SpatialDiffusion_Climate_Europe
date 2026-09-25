###############################################################################
# run_withBalkans_v2.R
#
# Same pipeline as run_withBalkans.R, refitted on the repaired covariate file
# produced by 10_climate_nuts_crosswalk.R, which restores the 44 regions that
# the NUTS release mismatch had removed from the cohort.
#
# The run happens inside run_withBalkans_v2/, with its own data/ and output/,
# so nothing produced by the earlier run is touched and the two can be compared
# directly.
#
# Usage, from the project folder:
#   Rscript run_withBalkans_v2.R
#
# To point it at a different covariate file:
#   WNV_CLIMATE_FILE=some_other.csv Rscript run_withBalkans_v2.R
###############################################################################

ROOT <- normalizePath(getwd())
CODE <- file.path(ROOT, "code")
# Configurable so a run at a different decay length lands in its own directory
# and cannot overwrite the primary results:
#   WNV_LAMBDA=75 WNV_RUN_DIR=run_lambda75 bash run_revision_v2.sh
RUN  <- file.path(ROOT, Sys.getenv("WNV_RUN_DIR", "run_withBalkans_v2"))

variant <- "full_follow_up_climatefix"

## ---- 1. Preconditions -------------------------------------------------------

need <- c("tidyverse", "sf", "eurostat", "units", "survival", "broom",
          "ggplot2", "ggthemes", "dplyr")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

prep_src  <- file.path(CODE, "00_data_preparation_withBalkans.R")
script_02 <- file.path(CODE, "02_survival_analyses.R")
script_03 <- file.path(CODE, "03_interaction_models.R")
case_csv  <- file.path(ROOT, "WNV_NUTS3cases_2006to2025_withBalkans.csv")

default_clim <- file.path(ROOT, "run_withBalkans", "output", "revision",
                          "NUTS3ObsClimate_1991to2024_NUTS2021.csv")
clim_csv <- Sys.getenv("WNV_CLIMATE_FILE", default_clim)
if (!file.exists(clim_csv) && file.exists(file.path(ROOT, clim_csv)))
  clim_csv <- file.path(ROOT, clim_csv)

for (f in c(prep_src, script_02, script_03, case_csv, clim_csv)) {
  if (!file.exists(f)) {
    stop("Not found: ", f,
         if (identical(f, clim_csv))
           "\nRun 10_climate_nuts_crosswalk.R first." else "", call. = FALSE)
  }
}
message("covariate file: ", clim_csv)

## ---- 2. Run directory -------------------------------------------------------

dir.create(file.path(RUN, "data"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RUN, "output"), recursive = TRUE, showWarnings = FALSE)

link_or_copy <- function(from, to) {
  if (file.exists(to) || !is.na(Sys.readlink(to)) && nzchar(Sys.readlink(to)))
    file.remove(to)
  ok <- suppressWarnings(file.symlink(from, to))
  if (!isTRUE(ok)) file.copy(from, to, overwrite = TRUE)
  invisible(TRUE)
}
link_or_copy(case_csv, file.path(RUN, "data", basename(case_csv)))
link_or_copy(clim_csv, file.path(RUN, "data", basename(clim_csv)))

# The variant and the covariate filename are set by rewriting three lines of the
# preparation script into the run directory, so the script in the codes folder
# stays exactly as it is.
prep_txt <- readLines(prep_src, warn = FALSE)

set_line <- function(txt, pattern, replacement, what) {
  i <- grep(pattern, txt)
  if (length(i) != 1L) stop("Could not find the ", what, " line in ", prep_src, call. = FALSE)
  txt[i] <- replacement
  txt
}

prep_txt <- set_line(prep_txt, "^sources_only <- ", "sources_only <- FALSE", "sources_only")
prep_txt <- set_line(prep_txt, "^late_entry_year\\s*<- ", "late_entry_year      <- 2006", "late_entry_year")
prep_txt <- set_line(
  prep_txt, "^climate_file <- ",
  paste0('climate_file <- file.path(data_dir, "', basename(clim_csv), '")'),
  "climate_file")

prep_run <- file.path(RUN, "00_data_preparation_run.R")
writeLines(prep_txt, prep_run)

setwd(RUN)
log_path <- file.path("output", paste0("run_log_", variant, ".txt"))
con <- file(log_path, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")

message("run_withBalkans_v2.R started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("variant: ", variant)
set.seed(1)

## ---- 3. Pipeline ------------------------------------------------------------

t0 <- Sys.time()
message("\n", strrep("-", 70), "\nRUN  : data preparation")
source(prep_run, echo = FALSE, local = FALSE)

message("\n", strrep("-", 70), "\nRUN  : 02_survival_analyses.R")
source(script_02, echo = FALSE, local = FALSE)

message("\n", strrep("-", 70), "\nRUN  : 03_interaction_models.R")

if (!exists("theme_classic2")) {
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    theme_classic2 <- ggpubr::theme_classic2
  } else {
    theme_classic2 <- function(base_size = 12, ...) ggplot2::theme_classic(base_size = base_size, ...)
    message("NOTE: ggpubr not installed; using theme_classic() for Figure 4A.")
  }
}

tryCatch(
  source(script_03, echo = FALSE, local = FALSE),
  error = function(e) {
    message("NOTE: 03_interaction_models.R stopped at: ", conditionMessage(e))
    message("      Model objects fitted before that point are kept and summarised below.")
  }
)

if (!exists("cox_interaction")) {
  stop("cox_interaction was not created; 03 failed before the interaction model.", call. = FALSE)
}

message("\nPipeline finished in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")

## ---- 4. Comparison table ----------------------------------------------------

library(dplyr)

get_row <- function(fit, term) {
  s <- summary(fit)
  if (!term %in% rownames(s$conf.int)) return(c(NA, NA, NA, NA))
  c(s$conf.int[term, "exp(coef)"], s$conf.int[term, "lower .95"],
    s$conf.int[term, "upper .95"], s$coefficients[term, "Pr(>|z|)"])
}

int_term <- grep("^distance100:tas$|^tas:distance100$",
                 rownames(summary(cox_interaction)$conf.int), value = TRUE)[1]

# Submitted values, then the withBalkans refit before the covariate repair.
published <- list(
  "Distance per 100 km"                = c(0.56, 0.50, 0.62),
  "Invasion pressure per SD"           = c(1.26, 1.17, 1.35),
  "Source-region intensity (log)"      = c(1.38, 1.18, 1.62),
  "Temperature per 1 C"                = c(1.25, 1.19, 1.31),
  "Relative humidity per 1 pct"        = c(0.96, 0.93, 0.99),
  "Distance x temperature interaction" = c(1.07, 1.04, 1.09)
)
withbalkans <- list(
  "Distance per 100 km"                = c(0.547, 0.494, 0.605),
  "Invasion pressure per SD"           = c(1.274, 1.183, 1.372),
  "Source-region intensity (log)"      = c(1.490, 1.295, 1.716),
  "Temperature per 1 C"                = c(1.212, 1.154, 1.272),
  "Relative humidity per 1 pct"        = c(0.960, 0.932, 0.989),
  "Distance x temperature interaction" = c(1.072, 1.046, 1.099)
)

new_vals <- list(
  "Distance per 100 km"                = get_row(cox_source_full, "distance100"),
  "Invasion pressure per SD"           = get_row(cox_source_full, "pressure_z"),
  "Source-region intensity (log)"      = get_row(cox_source_full, "log_cases_nearest"),
  "Temperature per 1 C"                = get_row(cox_source_full, "tas"),
  "Relative humidity per 1 pct"        = get_row(cox_source_full, "hurs"),
  "Distance x temperature interaction" = if (!is.na(int_term)) get_row(cox_interaction, int_term) else rep(NA_real_, 4)
)

comparison <- bind_rows(lapply(names(published), function(k) {
  p <- published[[k]]; w <- withbalkans[[k]]; n <- new_vals[[k]]
  data.frame(term = k,
             submitted_HR = p[1], submitted_lo = p[2], submitted_hi = p[3],
             balkans_HR = w[1], balkans_lo = w[2], balkans_hi = w[3],
             new_HR = round(n[1], 3), new_lo = round(n[2], 3),
             new_hi = round(n[3], 3), new_p = signif(n[4], 3))
}))

conc_new <- unname(summary(cox_source_full)$concordance["C"])
ev <- analysis_cases %>% filter(event == 1L)
d  <- ev$nearest_invaded_km
share <- function(x) round(100 * mean(d <= x, na.rm = TRUE), 1)

headline <- data.frame(
  quantity = c("First invasion events", "Regions followed", "Countries",
               "Harrell's C", "Invasions within 75 km (pct)",
               "Invasions within 150 km (pct)", "Invasions within 300 km (pct)",
               "Invasions beyond 500 km (n)", "Median invasion distance (km)"),
  submitted = c(269, 1011, 18, 0.875, 45, 71, 88, 6, NA),
  withBalkans = c(305, 1049, 21, 0.878, 47.9, 73.8, 89.2, 6, 76.8),
  new = c(nrow(ev), dplyr::n_distinct(analysis_cases$NUTS_ID),
          dplyr::n_distinct(analysis_cases$CNTR_CODE), round(conc_new, 3),
          share(75), share(150), share(300), sum(d > 500, na.rm = TRUE),
          round(median(d, na.rm = TRUE), 1))
)

by_country <- ev %>% count(CNTR_CODE, name = "first_invasions") %>% arrange(desc(first_invasions))
by_year    <- ev %>% count(year, name = "first_invasions") %>% arrange(year)

suffix <- paste0("_", variant, ".csv")
write.csv(comparison, file.path("output", paste0("comparison_submitted_vs_new", suffix)), row.names = FALSE)
write.csv(headline,   file.path("output", paste0("headline_numbers", suffix)),            row.names = FALSE)
write.csv(univar_results, file.path("output", paste0("univariable_models", suffix)),      row.names = FALSE)
write.csv(coef_table, file.path("output", paste0("multivariable_model", suffix)),         row.names = FALSE)
write.csv(by_country, file.path("output", paste0("events_by_country", suffix)),           row.names = FALSE)
write.csv(by_year,    file.path("output", paste0("events_by_year", suffix)),              row.names = FALSE)

message("\n", strrep("=", 70))
message("SUBMITTED / WITH BALKANS / WITH COVARIATE REPAIR")
print(comparison)
message("")
print(headline)
message("")
print(by_country)

writeLines(capture.output(sessionInfo()), file.path("output", paste0("sessionInfo_", variant, ".txt")))
message("\nWrote CSV outputs to ", file.path(RUN, "output"))

sink(type = "message")
sink()
close(con)

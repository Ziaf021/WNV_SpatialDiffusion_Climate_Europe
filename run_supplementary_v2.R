###############################################################################
# run_supplementary_v2.R
#
# Regenerates the supplementary analyses on the REFITTED cohort.
#
# Why this exists. run_revision_v2.sh runs steps 10, 00/02/03, 08 and 09 only. It
# has never run the sensitivity analyses, the proportional-hazards diagnostics or
# the supplementary figures, so every artefact behind the appendix still reports
# the submitted cohort: 1011 regions, 18 countries, 269 events. The main text now
# reports 1084, 21 and 306. This script closes that gap.
#
#   bash run_supplementary_v2.sh
#
# What it runs, in order:
#   04_sensitivity_prespecified_02.R   lambda profile, kernel comparison,
#                                      leave-country-out, frailty
#   05_ph_diagnostics_2.R              Schoenfeld residuals, PH tests
#   06_supplementary_figures.R         supplementary figures
#
# Two things are rewritten into run-directory copies rather than edited in the
# codes folder, so the originals stay as they are:
#   - climate_file, because the repaired covariate file in run_withBalkans_v2/data
#     carries the _NUTS2021 suffix and 04 looks for the unsuffixed name
#   - output_dir, so supplementary outputs land in output/supplementary and cannot
#     overwrite anything the main chain produced
#
# The leave-country-out analysis will now cover 21 countries rather than 18, so
# that figure and its table change shape as well as their numbers.
###############################################################################

suppressPackageStartupMessages(library(dplyr))

ROOT <- normalizePath(getwd())
RUN  <- file.path(ROOT, Sys.getenv("WNV_RUN_DIR", "run_withBalkans_v2"))
CODE <- file.path(ROOT, "code")
OUT  <- file.path(RUN, "output", "supplementary")

if (!dir.exists(RUN)) stop("Not found: ", RUN, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

scripts <- c("04_sensitivity_prespecified_02.R",
             "05_ph_diagnostics_2.R",
             "06_supplementary_figures.R")
missing <- scripts[!file.exists(file.path(CODE, scripts))]
if (length(missing)) stop("Not found: ", paste(missing, collapse = ", "), call. = FALSE)

## ---- Inputs, and a guard on the cohort --------------------------------------

ac_path <- file.path(RUN, "analysis_cases_with_cats.rds")
fit_path <- file.path(RUN, "cox_source_full.rds")
for (f in c(ac_path, fit_path)) {
  if (!file.exists(f)) stop("Not found: ", f, "\nRun run_revision_v2.sh first.",
                            call. = FALSE)
}

ac_check <- readRDS(ac_path)
message("cohort in ", basename(ac_path), ": ",
        nrow(ac_check), " region-years, ", sum(ac_check$event), " events, ",
        dplyr::n_distinct(ac_check$NUTS_ID), " regions, ",
        dplyr::n_distinct(ac_check$CNTR_CODE), " countries")

if (sum(ac_check$event) != as.integer(Sys.getenv("WNV_EXPECTED_EVENTS", "306"))) {
  stop("Expected 306 events on the refitted cohort, found ", sum(ac_check$event),
       ". The run directory is stale; rerun run_revision_v2.sh.", call. = FALSE)
}
if (!"pr_z" %in% names(ac_check)) {
  stop("pr_z is missing from ", basename(ac_path),
       ". The run directory predates the pr_z change; rerun run_revision_v2.sh.",
       call. = FALSE)
}
rm(ac_check)

# The climate file actually present, whatever it is called.
clim <- list.files(file.path(RUN, "data"), pattern = "^NUTS3ObsClimate.*\\.csv$",
                   full.names = FALSE)
if (length(clim) != 1L) {
  stop("Expected exactly one NUTS3ObsClimate file in ", file.path(RUN, "data"),
       ", found ", length(clim), call. = FALSE)
}
message("climate file: ", clim)

## ---- Objects the scripts expect to find in session ---------------------------
## 05 defines build_pressure_cases inside the else branch of a guard that tests
## whether the spatial objects are in session, then calls it at the end of the
## script, outside that branch. With the objects absent the definition is skipped
## and the call fails, so Table S7, the case-weighted kernel diagnostic the
## Methods cite, is silently unreachable when the script runs on its own. Loading
## the objects both fixes the failure and makes S7 run.

.sp_path  <- file.path(RUN, "output", "spatial_objects_withBalkans.rds")
.csv_path <- file.path(RUN, "output", "WNV_SurvivalDataset_withBalkans.csv")
for (f in c(.sp_path, .csv_path)) {
  if (!file.exists(f)) stop("Not found: ", f, call. = FALSE)
}

# 04 has no guard at all: it expects analysis_cases to be in session and fails
# outright without it. 05 guards and would read this file itself; loading it once
# here serves both and guarantees they use the same object.
analysis_cases <- readRDS(ac_path)
message("preloaded: analysis_cases ", nrow(analysis_cases), " region-years, ",
        sum(analysis_cases$event), " events")
stopifnot("pr_z" %in% names(analysis_cases))

.sp <- readRDS(.sp_path)
dist_mat       <- .sp$dist_mat
region_index   <- .sp$region_index
region_lookup  <- .sp$region_lookup
first_outbreak <- .sp$first_outbreak
dat <- utils::read.csv(.csv_path)
message("preloaded: dist_mat ", paste(dim(dist_mat), collapse = " x "),
        ", region_lookup ", nrow(region_lookup),
        ", case series ", nrow(dat), " rows")
stopifnot(all(c("NUTS_ID", "year", "wnf_sum") %in% names(dat)))

## ---- Patch the scripts into the run directory --------------------------------

patch <- function(name) {
  # Read and write as UTF-8 explicitly. 04 contains a backtick-quoted column name
  # holding a \u escape, and R's parser rejects those unless it knows the input
  # is UTF-8.
  txt <- readLines(file.path(CODE, name), warn = FALSE, encoding = "UTF-8")

  i <- grep('^\\s*output_dir\\s*<-', txt)
  if (length(i)) txt[i] <- 'output_dir <- file.path("output", "supplementary")'

  # R's parser rejects a \uXXXX escape inside a backtick-quoted name unless it is
  # certain the input is UTF-8, and it stays unconvinced even with encoding set.
  # 04 has two such names, the lambda and delta-AIC column headers. Every escape
  # in the file is converted to its literal character, which parses cleanly and is
  # equally valid in the ordinary strings that carry the rest of them.
  .esc <- unique(unlist(regmatches(txt, gregexpr("\\\\u[0-9A-Fa-f]{4}", txt))))
  for (.e in .esc) {
    txt <- gsub(.e, eval(parse(text = paste0('"', .e, '"'))), txt, fixed = TRUE)
  }
  if (length(.esc)) message("  unicode escapes converted: ", length(.esc))

  j <- grep('^\\s*climate_file\\s*<-', txt)
  if (length(j)) {
    txt[j] <- sprintf('climate_file <- file.path(data_dir, "%s")', clim)
  }

  out <- file.path(RUN, sub("\\.R$", "_supp_run.R", name))
  con_out <- file(out, open = "wt", encoding = "UTF-8")
  writeLines(txt, con_out)
  close(con_out)
  message("patched -> ", basename(out),
          "  (output_dir lines: ", length(i), ", climate_file lines: ", length(j), ")")
  out
}

.patched <- vapply(scripts, patch, character(1))

## ---- Run ---------------------------------------------------------------------
## The scripts run in the global environment, in sequence, which is how run_all.R
## has always run them. Sandboxing each in its own environment fails: 05 resolves
## analysis_cases from the global environment at its proportional-hazards section
## and cannot see a local copy. Sharing the environment is safe here because the
## only objects the exists() checks pick up, analysis_cases and cox_source_full,
## are read from the same refitted files either way, and the cohort was verified
## against 306 events before any of this started.
##
## source() rather than sys.source(), because only source() takes an encoding.

## Control names are dot-prefixed. The scripts assign freely into this
## environment, and a script that happened to use `con` or `k` would otherwise
## close the log or break the loop.

.old_wd <- setwd(RUN)
on.exit(setwd(.old_wd), add = TRUE)

.log_path <- file.path("output", "supplementary",
                       paste0("run_log_supplementary_",
                              format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
.con <- file(.log_path, open = "wt")
sink(.con, split = TRUE)
sink(.con, type = "message")
on.exit({ sink(type = "message"); sink(); close(.con) }, add = TRUE)

message("run_supplementary_v2.R started ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("run directory: ", RUN)
set.seed(1)

.results <- character(0)

for (.k in seq_along(.patched)) {
  .nm <- names(.patched)[.k]
  message("\n", strrep("-", 70), "\nRUN  : ", .nm)
  .t0 <- Sys.time()
  .ok <- tryCatch({
    source(.patched[[.k]], local = FALSE, encoding = "UTF-8", echo = FALSE)
    TRUE
  }, error = function(e) {
    message("FAILED: ", .nm, "\n  ", conditionMessage(e))
    FALSE
  })
  .mins <- round(as.numeric(difftime(Sys.time(), .t0, units = "mins")), 1)
  message(if (.ok) "done: " else "failed: ", .nm, "  (", .mins, " min)")
  .results <- c(.results, paste0(if (.ok) "OK     " else "FAILED ", .nm,
                                 "  (", .mins, " min)"))
}

message("\n", strrep("=", 70), "\nSUMMARY")
for (.r in .results) message("  ", .r)
message(strrep("=", 70))

## 06 has no output_dir to patch, so its figures land in the run directory.
## Collect them, otherwise they sit beside the pipeline outputs and the next
## reader cannot tell which run produced them.
.stray <- list.files(".", pattern = "^Supplementary_Figure_S.*\\.pdf$")
if (length(.stray)) {
  .dest <- file.path("output", "supplementary", .stray)
  file.rename(.stray, .dest)
  message("collected from the run directory: ", paste(.stray, collapse = ", "))
}

message("\nOutputs: ", file.path(RUN, "output", "supplementary"))
writeLines(capture.output(sessionInfo()),
           file.path("output", "supplementary", "sessionInfo_supplementary.txt"))
message("run_supplementary_v2.R finished ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

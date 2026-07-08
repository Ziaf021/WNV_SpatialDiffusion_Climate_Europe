###############################################################################
# run_all.R
#
# Master reproduction script for:
#   "Spatial diffusion and climatic suitability shape the West Nile Virus
#    invasion risk across Europe"
#   Author: Dr. Zia Farooq, Umea University
#
# Sources the pipeline in strict dependency order and reproduces every
# reported quantity and figure from the raw inputs in one command.
#
# Settled specification (verified via clean-session sanity check):
#   - kernel decay length      lambda_primary = 50  (AIC-minimising; see 04)
#   - invasion pressure         standardised on the analysis sample
#   - log_cases_nearest         built from the full case series (dat), lag 1 yr
#
# Usage:
#   Set the working directory to the repository root, then:
#     source("run_all.R")
#   or from a shell:
#     Rscript run_all.R
#
# Reproducibility:
#   If using renv, uncomment the renv::restore() line below to install the
#   exact package versions recorded in renv.lock before running.
###############################################################################

## ---- Configuration ----------------------------------------------------------

script_dir <- "code"     # directory holding the numbered/figure scripts
output_dir <- "output"   # derived objects, logs, sessionInfo
log_file   <- file.path(output_dir, "run_all_log.txt")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# renv::restore()   # uncomment to pin package versions

# Global seed. No pipeline step is stochastic; this only fixes visual jitter
# in figures so panels reproduce pixel-for-pixel. It does not affect inference.
set.seed(1)


## ---- Ordered pipeline -------------------------------------------------------
# Each script assumes the objects created by those above it are in the session.
# Filenames below match the code/ folder as staged for the GitHub repository.

pipeline <- c(
  # Foundation
  "00_data_preparation.R",       # -> analysis_cases, dist_mat, region_*, W
  # Primary models
  "02_survival_analyses.R",      # -> cox_source_full, km_fit, spline grid
  "03_interaction_models.R",     # -> cox_interaction, hazard-surface grid
  # Robustness and diagnostics
  "04_sensitivity_prespecified.R",
  "05_ph_diagnostics.R",
  "06_supplementary_figures.R",
  # Main-text and supplementary figures (consume saved model objects from 02/03)
  "fig2_barplot_pressure.R",
  "fig3_forest_spline.R",
  "fig4_tempdist_interaction.R",
  "fig_urbanisation_firstinvaded.R",
  "fig_urbanisation_by_phase.R"
)


## ---- Runner -----------------------------------------------------------------
# Sources one script with timing and a clear failure point. Stops the whole
# run at the first error, since every downstream step inherits upstream state.

run_step <- function(file, dir = script_dir) {
  path <- file.path(dir, file)
  if (!file.exists(path)) {
    stop("Script not found: ", path, call. = FALSE)
  }
  message("\n", strrep("-", 70))
  message("RUN  : ", file)
  t0 <- Sys.time()
  tryCatch(
    source(path, echo = FALSE, local = FALSE),
    error = function(e) {
      stop("FAILED in ", file, " :: ", conditionMessage(e), call. = FALSE)
    }
  )
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message("OK   : ", file, "  (", dt, "s)")
  invisible(dt)
}


## ---- Execute ----------------------------------------------------------------

con <- file(log_file, open = "wt")
sink(con, split = TRUE)          # mirror console to log
sink(con, type = "message")

message("run_all.R started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
t_start <- Sys.time()

for (script in pipeline) run_step(script)

t_total <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 2)
message("\n", strrep("=", 70))
message("Pipeline completed in ", t_total, " min: ",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

# Record the exact environment for the deposit.
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))
message("Wrote: ", file.path(output_dir, "sessionInfo.txt"))

sink(type = "message")
sink()
close(con)

###############################################################################
# End of run_all.R
###############################################################################

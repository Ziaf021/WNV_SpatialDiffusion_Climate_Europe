#!/usr/bin/env bash
###############################################################################
# run_revision_v3.sh   (written 16 September 2026)
#
# Full rerun on the ECDC-patched case file and the 2026 week-37 data
# (to 10 September 2026). Everything goes to run_withBalkans_v3/, so the
# v2 results used in the current drafts are left untouched.
#
# Data changes behind this run (see data_patch_2026-09-16/):
#   - 81 region-years added to WNV_NUTS3cases_2006to2025.csv from ECDC tables
#     (Croatia, Hungary Budapest/Pest, France, Sardinia, Turkey, Bulgaria,
#     Czechia and a few single rows), mostly records lost to NUTS recoding
#   - WNV_2026_ECDC.csv replaced by the ECDC week-37 table
#
# Usage, from the project folder:
#   bash run_revision_v3.sh
# If one step fails, fix it and restart from that step:
#   START=5 bash run_revision_v3.sh
###############################################################################

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
CODE="$ROOT/code"

export WNV_RUN_DIR="${WNV_RUN_DIR:-run_withBalkans_v3}"
export WNV_LAMBDA="${WNV_LAMBDA:-50}"
export WNV_EXPECTED_EVENTS="${WNV_EXPECTED_EVENTS:-324}"
START="${START:-1}"
STOP="${STOP:-999}"

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$ROOT/run_revision_v3_${STAMP}.log"

command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found on PATH." >&2; exit 1; }

SRC_CACHE="$ROOT/run_withBalkans_v2/output/revision/boundary_dist_mat.rds"
[ -f "$SRC_CACHE" ] || SRC_CACHE="$ROOT/run_withBalkans/output/revision/boundary_dist_mat.rds"
DST_CACHE="$ROOT/$WNV_RUN_DIR/output/revision/boundary_dist_mat.rds"
if [ -f "$SRC_CACHE" ] && [ ! -f "$DST_CACHE" ]; then
  mkdir -p "$(dirname "$DST_CACHE")"
  cp "$SRC_CACHE" "$DST_CACHE"
fi

step() {
  local n="$1"; shift
  local label="$1"; shift
  if [ "$n" -lt "$START" ] || [ "$n" -gt "$STOP" ]; then
    echo "=== step $n skipped: $label ===" | tee -a "$LOG"
    return 0
  fi
  echo | tee -a "$LOG"
  echo "=== step $n: $label  ($(date +%H:%M:%S)) ===" | tee -a "$LOG"
  "$@" 2>&1 | tee -a "$LOG"
}

echo "Run directory: $WNV_RUN_DIR   expected events: $WNV_EXPECTED_EVENTS" | tee "$LOG"

step 1  "merge Balkan series into the patched case file" python3 "$ROOT/merge_balkan_cases.py"
step 2  "refit (00, 02, 03)"                 Rscript "$ROOT/run_withBalkans_v2.R"
step 3  "revision analyses (08)"             Rscript "$CODE/08_revision_analyses.R"
step 4  "prospective validation (09)"        Rscript "$CODE/09_prospective_validation_v2.R"
step 5  "supplementary analyses (04 05 06)"  Rscript "$ROOT/run_supplementary_v2.R"
step 6  "Table 1 (13)"                       Rscript "$CODE/13_table1_univariable.R"
step 7  "regional composition (12)"          Rscript "$CODE/12_regional_composition.R"
step 8  "direction and latitude (12b)"       Rscript "$CODE/12b_direction_latitude_4periods.R"

n=9
for key in 1a 1af 1b 1 1f 2a 2b 2 3 4a 4b 4 5 s10; do
  step "$n" "figure $key" bash "$ROOT/run_fig.sh" "$key"
  n=$((n+1))
done

if [ "$n" -ge "$START" ] && [ "$n" -le "$STOP" ]; then
  echo "=== step $n: urbanisation figures (run inside $WNV_RUN_DIR) ===" | tee -a "$LOG"
  ( cd "$ROOT/$WNV_RUN_DIR" && \
    Rscript "$CODE/fig_urbanisation_firstinvaded.R" && \
    Rscript "$CODE/fig_urbanisation_by_phase.R" && \
    Rscript "$CODE/urbanisation_phase_enrichment.R" ) 2>&1 | tee -a "$LOG"
fi

echo | tee -a "$LOG"
echo "Finished $(date). Outputs in $ROOT/$WNV_RUN_DIR/output" | tee -a "$LOG"

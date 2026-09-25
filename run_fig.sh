#!/usr/bin/env bash
###############################################################################
# run_fig.sh
#
# Builds ONE figure, so each can be fixed on its own without rebuilding the set.
#
# Usage, from the project folder:
#   bash run_fig.sh 1a      chronology map, panel A alone
#   bash run_fig.sh 1af     chronology map faceted by period, an experiment
#   bash run_fig.sh 1f      faceted maps over the heatmap, one column
#   bash run_fig.sh 1b      first invasions by country and year, panel B alone
#   bash run_fig.sh 1       panels A and B combined on one canvas
#   bash run_fig.sh 2a      distance at first invasion
#   bash run_fig.sh 2b      invasion rate against neighbourhood pressure
#   bash run_fig.sh 3       forest, distance spline, Kaplan-Meier
#   bash run_fig.sh 4a      temperature at first invasion by distance class
#   bash run_fig.sh 4b      predicted distance by temperature hazard surface
#   bash run_fig.sh 4       panels A and B combined on one canvas
#   bash run_fig.sh 5       out-of-sample test
#   bash run_fig.sh s10     Supplementary Figure S10, direction of invasion (rose)
#   bash run_fig.sh s11     Supplementary Figure S11, latitude trend over time
#   bash run_fig.sh list    show what is available
#
# Output goes to <run dir>/output/revision/figures as PDF, PNG and SVG, with the
# caption as a text file in figures/captions.
#
# The NUTS3 geometry is cached after the first run. To force a refetch:
#   WNV_FIGDATA_REFRESH=1 bash run_fig.sh 1a
###############################################################################

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

CODE="$ROOT/code"

declare -a KEYS=(1a 1af 1b 1 1f 2a 2b 2 3 4a 4b 4 5 s10 s11)
script_for() {
  case "$1" in
    1a) echo "fig01a_chronology_map.R" ;;
    1af) echo "fig01a_facet_periods.R" ;;
    1f) echo "fig01f_combined.R" ;;
    1)  echo "fig01_combined.R" ;;
    1b) echo "fig01b_country_year.R" ;;
    2a) echo "fig02a_distance_distribution.R" ;;
    2b) echo "fig02b_pressure_rate.R" ;;
    2)  echo "fig02_combined.R" ;;
    3)  echo "fig03_forest_spline_km.R" ;;
    4a) echo "fig04a_temperature_by_distance.R" ;;
    4b) echo "fig04b_hazard_surface.R" ;;
    4)  echo "fig04_combined.R" ;;
    5)  echo "fig05_out_of_sample.R" ;;
    s10) echo "figS10_directionality_rose.R" ;;
    s11) echo "figS11_latitude_trend.R" ;;
    *)  echo "" ;;
  esac
}

if [ $# -lt 1 ] || [ "$1" = "list" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  echo "Available figures:"
  for k in "${KEYS[@]}"; do
    s="$(script_for "$k")"
    if [ -f "$CODE/$s" ]; then mark="ready"; else mark="not written yet"; fi
    printf "  %-4s %-36s %s\n" "$k" "$s" "$mark"
  done
  echo
  echo "Usage: bash run_fig.sh <figure>"
  exit 0
fi

KEY="$1"
SCRIPT="$(script_for "$KEY")"

if [ -z "$SCRIPT" ]; then
  echo "Unknown figure '$KEY'. Run 'bash run_fig.sh list' to see the options." >&2
  exit 1
fi

if [ ! -f "$CODE/$SCRIPT" ]; then
  echo "Not written yet: $CODE/$SCRIPT" >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found on PATH." >&2
  exit 1
fi

export WNV_RUN_DIR="${WNV_RUN_DIR:-run_withBalkans_v2}"

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$ROOT/run_fig_${KEY}_${STAMP}.log"

{
  echo "Figure        : $KEY"
  echo "Script        : $SCRIPT"
  echo "Run directory : $WNV_RUN_DIR"
  echo "Started       : $(date)"
  echo
} | tee "$LOG"

Rscript "$CODE/$SCRIPT" 2>&1 | tee -a "$LOG"

{
  echo
  echo "Finished      : $(date)"
  echo "Figures       : $ROOT/$WNV_RUN_DIR/output/revision/figures"
} | tee -a "$LOG"

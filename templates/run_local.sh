#!/usr/bin/env bash
#
# TEAPOT launcher for one machine.
#   cp this, a templates/<case>.yml and its .csv into a directory of your
#   own, fill in the placeholders, then: ./run_local.sh
#
# Use a directory of your own, not the pipeline checkout.
set -euo pipefail

# --- edit these ---------------------------------------------------------------
TEAPOT_DIR="<path/to/teapot>"          # the pipeline checkout
PARAMS=my_run.yml                      # your copy of a templates/*.yml
# ------------------------------------------------------------------------------

[ -d "$TEAPOT_DIR" ] || { echo "TEAPOT_DIR does not exist: $TEAPOT_DIR" >&2; exit 1; }
[ -f "$PARAMS" ]     || { echo "params file not found: $PARAMS" >&2; exit 1; }
grep -q '<' "$PARAMS" && { echo "$PARAMS still has <...> placeholders to fill in" >&2; exit 1; }
command -v nextflow > /dev/null || { echo "nextflow is not on PATH" >&2; exit 1; }

command -v apptainer > /dev/null || {
    echo "apptainer not found: every tool runs in an apptainer container." >&2
    exit 1
}

mkdir -p results/pipeline_info

nextflow -log "$PWD/results/pipeline_info/nextflow.log" \
    run "$TEAPOT_DIR/main.nf" \
    -profile apptainer \
    -params-file "$PARAMS" \
    -resume

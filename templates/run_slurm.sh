#!/usr/bin/env bash
#SBATCH -J teapot
#SBATCH -n 2
#SBATCH --mem 8G
#SBATCH -t 2-00:00:00
#SBATCH -o head.%j.out
#SBATCH -e head.%j.err
#
# TEAPOT launcher for a SLURM cluster.
#   cp this, a templates/<case>.yml and its .csv into a directory of your
#   own, fill in the placeholders, then: sbatch run_slurm.sh
#
# Use a directory of your own, not the pipeline checkout.
set -euo pipefail

# --- edit these ---------------------------------------------------------------
TEAPOT_DIR=<path/to/teapot>          # the pipeline checkout
PARAMS=my_run.yml                    # your copy of a templates/*.yml
SITE_CONFIG="$TEAPOT_DIR/conf/slurm.config"   # resources; copy and adjust per cluster

# Your cluster's module names. Delete if java and nextflow are already on PATH.
module load Java/21.0.11-bdist Nextflow/25.10.6-eb
# ------------------------------------------------------------------------------

[ -d "$TEAPOT_DIR" ]  || { echo "TEAPOT_DIR does not exist: $TEAPOT_DIR" >&2; exit 1; }
[ -f "$PARAMS" ]      || { echo "params file not found: $PARAMS" >&2; exit 1; }
grep -q '<' "$PARAMS" && { echo "$PARAMS still has <...> placeholders to fill in" >&2; exit 1; }

# Keep the head JVM small; the work happens in the child jobs.
export NXF_OPTS='-Xms512m -Xmx4g'

mkdir -p results/pipeline_info

nextflow -log "$PWD/results/pipeline_info/nextflow.log" \
    run "$TEAPOT_DIR/main.nf" \
    -profile hpc \
    -params-file "$PARAMS" \
    -c "$SITE_CONFIG" \
    -resume

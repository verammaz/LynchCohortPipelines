#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if command exits with non-zero status
set -e


script="$LYNCH/star_alignment.sh"
chmod +x $script

echo "script: $script"
echo "project: $project"

module purge
module load star

bsub -J "star" \
    -P ${project} \
    -n ${cores} \
    -M 2000 \
    -R span[hosts=1] \
    -R "rusage[mem=4000]" \
    -W 10:00 \
    -q premium \
    -oo "${LOG_DIR}/star.out" \
    -eo "${LOG_DIR}/star.err" \
    ${script} -r r1, r2


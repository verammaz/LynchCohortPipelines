#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -s <samplesheet.csv> -p <patient_id>

Required Arguments:
    -s <samplesheet.csv>            Configured input file
    -p <patient_id>                 Patient identifier

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --step                          Step to start from (0=alignment, 1=markdup, 2=indelrealign, 3=baserecal)
}

EOF
    exit 1
}


PATIENT=
SAMPLESHEET=
VERBOSE=0
STEP=0

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -s) SAMPLESHEET="$2"; shift ;;
        --step) STEP="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done


script="$LYNCH/star_alignment.sh"
chmod +x $script

echo "script: $script"
echo "project: $project"

module purge
module load star

bsub -J "star" \
    -P ${project} \
    -n ${cores} \
    -M 32000 \
    -R span[hosts=1] \
    -R "rusage[mem=4000]" \
    -W 40:00 \
    -q premium \
    -oo "${LOG_DIR}/star.out" \
    -eo "${LOG_DIR}/star.err" \
    ${script} 


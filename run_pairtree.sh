#!/bin/bash
## This script runs entire pairtree pipeline (pairtree input prep, clustvars, pairtree, plottree, summposterior)

##################################################
#specify path to auxilary scripts
scripts='/hpc/users/mazeev01/matt_lynch/scripts'
##################################################

# Function to print progress with timestamp
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] --patient <patient_id> --data <data_dir>

Required Arguments:
    --patient <patient_id>           Patient id
    --data <data_dir>                Path to data directory, with VCF/ subdirectory

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --samples <s1,s2,...>           List of sample ids to use (comma-separated, no space)
}

EOF
    exit 1
}


VERBOSE=0
PATIENT=
DATA=
SAMPLES=

while [[ "$#" -gt 0 ]]; do
    case "$1" in 
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -patient) PATIENT="$2"; shift ;;
        -data) DATA="$2"; shift ;;
        -samples) SAMPLES="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

module load python

print_progress "Prepping pairtree input."

python $scripts/pairtree_prep.py -patient_id $PATIENT -vcf_dir $DATA/VCF -samples $SAMPLES

chmod +x $scripts/pairtree.sh

$scripts/pairtree.sh $PATIENT $DATA
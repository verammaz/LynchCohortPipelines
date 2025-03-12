#!/bin/bash

source ./config.sh

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -p <patient_id>

Required Arguments:
    -p <patient_id>                         Patient id
Options:
    -h                                      Display this message
    -v                                      Enable verbode mode
    -vcf_dir </path/to/VCF/Patient>         Directory with all sample VCF files for patient
    -patient_sex_file <patient_sex.txt>     File with patient sex info
    --single_sample_trees                   Create single sample trees
    --total_alt                             Flag that VCF file is in total:alt format
}

EOF
    exit 1
}


VERBOSE=0
PATIENT=
SINGLE_TREES=0
VCF_DIR=
SEX_FILE=
TOTAL_ALT=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in 
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -vcf_dir) VCF_DIR="$2"; shift ;;
        -patient_sex_file) SEX_FILE="$2"; shift ;;
        --single_sample_trees) SINGLE_TREES=1 ;;
        --total_alt) TOTAL_ALT=1 ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

module load python

print_progress "Prepping pairtree input."

if [[ -z "$VCF_DIR" ]]; then
    VCF_DIR="${HOME_DIR}/VCF/${PATIENT}"
fi

python "$LYNCH/pairtree_prep.py" -patient_id $PATIENT \
                                 -vcf_dir $VCF_DIR \
                                 -single_sample_trees $SINGLE_TREES \
                                 -total_alt $TOTAL_ALT \
                                 ${SEX_FILE:+-patient_sex_info_file $SEX_FILE}

script="$LYNCH/pairtree_run.sh"
chmod +x $script

job_name="pairtree_${PATIENT}"

print_progress "Submitting pairtree job(s)..."

bsub -J ${job_name} \
     -P ${project} \
     -q premium \
     -n 8 \
     -W 08:00 \
     -M 16000 \
     -R "rusage[mem=4000]" \
     -oo "${LOG_DIR}/${job_name}.out" \
     -eo "${LOG_DIR}/${job_name}.err" \
     ${script} ${PAIRTREE} ${PATIENT} "${HOME_DIR}/Pairtree/${PATIENT}"


if [[ $SINGLE_TREES -eq 1 ]]; then
    for file in "${VCF_DIR}"/*.vcf; do
        sample=$(basename "${file}" .vcf)

        job_name="pairtree_${sample}"

        bsub -J ${job_name} \
            -P ${project} \
            -q premium \
            -n 8 \
            -W 08:00 \
            -M 16000 \
            -R "rusage[mem=4000]" \
            -oo "${LOG_DIR}/${job_name}.out" \
            -eo "${LOG_DIR}/${job_name}.err" \
            ${script} ${PAIRTREE} ${sample} "${HOME_DIR}/Pairtree/${sample}"

    done
fi
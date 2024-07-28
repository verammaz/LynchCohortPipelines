#!/bin/bash

# This is a wrapper script to submit jobs to execute the nextflow 
# nf-core/sarek pipeline

source ./config.sh

# Exit immediately if command exits with non-zero status
set -e


# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -p <patient_id> -s <samplesheet.csv>

Required Arguments:
    -p <patient_id>                 Patient id
    -s <samplesheet.csv>            File with all sample raw data configured

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --variant_call_step             Step of variant call (0 = preprocessing, 1 = variant calling)
}

EOF
    exit 1
}



PATIENT=
SAMPLESHEET=
VERBOSE=0
STEP=1

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -s) SAMPLESHEET="$2"; shift ;;
        --variant_call_step) STEP="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

script="$LYNCH/run_nextflow.sh"
chmod +x $script

RAW_DIR=="${HOME_DIR}/Raw/${PATIENT}"


while IFS=$',' read -r patient sample fastq1 fastq2 status bam bai; do
    
    if [[ $patient == $PATIENT ]] && [ $sample != "Normal" ]; then

        job_name="variant_call_$sample"

        bsub -J ${job_name} \
             -P ${project} \
             -n ${cores} \
             -M 32000 \
             -R span[hosts=1] \
             -R "rusage[mem=4000]" \
             -W 30:00 \
             -q premium \
             -oo "${LOG_DIR}/${job_name}.out" \
             -eo "${LOG_DIR}/${job_name}.err"r \
             ${script} --patient ${PATIENT} --samplesheet ${SAMPLESHEET} --sample ${sample} --step ${STEP}
    
    fi

done < "$SAMPLESHEET"
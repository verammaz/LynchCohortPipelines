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
    --step                          Step of variant call (0 = preprocessing, 1 = variant calling)
    --ref                           Reference genome (default=hg19)

}

EOF
    exit 1
}



PATIENT=
SAMPLESHEET=
VERBOSE=0
STEP=1
REF="hg19"

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -s) SAMPLESHEET="$2"; shift ;;
        --step) STEP="$2"; shift ;;
        --ref) REF="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

GENOME=

if [[ "$REF" == 'hg19' ]]; then
    GENOME=$GENOME_hg19
fi

if [[ "$REF" == 'hg38' ]]; then
    GENOME=$GENOME_hg38
fi

if [[ -z "$GENOME" ]]; then
    echo "Error: GENOME is not set. Check that REF is either 'hg19' or 'hg38', and corresponding variables (GENOME_hg19 / GENOME_hg38) are defined in your config."
    exit 1
fi

script1="$LYNCH/variant_call.sh"
chmod +x $script1

RAW_DIR=="${HOME_DIR}/Raw/${PATIENT}"

job_names=()
sample_names=()

{
read  # Skip the header line
while IFS= read -r line || [[ -n "$line" ]]; do
    
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Split the line into fields using awk to handle potential edge cases
    IFS=',' read -r patient sample fastq1 fastq2 status bam bai ref <<< "$(awk -F',' '{print $1,$2,$3,$4,$5,$6,$7}' OFS=',' <<< "$line")"
    
    if [[ $patient == $PATIENT ]] && [ $sample != "Normal" ]; then

        sample_names+=($sample)
        job_name="variant_call_$sample"
        job_names+=("$job_name")

        bsub -J ${job_name} \
             -P ${project} \
             -n ${cores} \
             -M 32000 \
             -R span[hosts=1] \
             -R "rusage[mem=4000]" \
             -R select[osmajor!=CENT7] \
             -W 30:00 \
             -q premium \
             -oo "${LOG_DIR}/${job_name}.out" \
             -eo "${LOG_DIR}/${job_name}.err" \
             ${script1} --patient ${PATIENT} --samplesheet ${SAMPLESHEET} --sample ${sample} --step ${STEP} --genome ${GENOME}
    
    fi

done 
} < "$SAMPLESHEET"


# Create the wait condition string
wait_condition=""
for job in "${job_names[@]}"; do
    wait_condition+="done(${job}) && "
done
wait_condition=${wait_condition%" && "}  # Remove the trailing ' && '

script2="$LYNCH/variant_call_cleanup.sh"
chmod +x $script2

sample_names_str=$(IFS=','; echo "${sample_names[*]}")

bsub -w "$wait_condition" \
     -q long \
     -n 1 \
     -P ${project} \
     -W 00:20 \
     ${script2} ${PATIENT} ${SAMPLESHEET} ${sample_names_str}


#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -s <samplesheet.csv>

Required Arguments:
    -s <samplesheet.csv>            Configured input file
    -p <patient_id>                 Patient identifier

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --data_dir                      Directory for patient raw data; Sample/ and Normal/ subdirectories with .bam files
}

EOF
    exit 1
}


PATIENT=
REFERENCE=
SAMPLESHEET=
VERBOSE=0

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -s) SAMPLESHEET="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

RAW_DIR="$HOME_DIR/Raw/$PATIENT"

script="$LYNCH/fastq_to_bam.sh"
chmod +x $script

# Create a temporary file to store the modified sample sheet
TEMP_SAMPLESHEET=$(mktemp)

echo "patient,sample,fastq_1,fastq_2,status,bam,bai" >> "$TEMP_SAMPLESHEET"


module purge
module load singularity/3.6.4

while IFS=$',' read -r patient sample fastq1 fastq2 status; do
    
    if [[ $patient == $PATIENT ]] || [[ -z $PATIENT ]]; then

        job_name="fastq_to_bam_$sample"
        create_directory "$RAW_DIR/$sample"
        output_prefix="$RAW_DIR/$sample/$sample"

        bsub -J ${job_name} \
                -P ${project} \
                -n ${cores} \
                -M 32000 \
                -R span[hosts=1] \
                -R "rusage[mem=4000]" \
                -W 40:00 \
                -q premium \
                -oo "${LOG_DIR}/${job_name}.out" \
                -eo "${LOG_DIR}/${job_name}.err"r \
                singularity exec ${CONTAINER_FASTQ2BAM} ${script} -r ${fastq1},${fastq2} -o ${output_prefix} --patient ${PATIENT} --threads ${cores} --post_process


        # Overwrite bam and bai with new paths
        bam="${output_prefix}.bam"
        bai="${output_prefix}.bai"

        # Write the updated line to the temporary file
        echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai" >> "$TEMP_SAMPLESHEET"
    
    fi

done < "$SAMPLESHEET"

# Replace the original sample sheet with the modified one
mv "$TEMP_SAMPLESHEET" "$SAMPLESHEET"

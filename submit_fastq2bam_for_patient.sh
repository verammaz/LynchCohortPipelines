#!/bin/bash

####################### Change the following ! #########################
#specify project
project='acc_FLAI'
#specify cores
cores=24
#specify path to main fastq_to_bam script
script='/hpc/users/mazeev01/matt_lynch/scripts/fastq_to_bam.sh'
########################################################################

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -p <patient_id> -s <samplesheet.csv>

Required Arguments:
    -p <patient_id>                 Patient id
    -s <samplesheet.csv>            Configured input file

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --data_dir                      Directory with patient raw data; Sample/ and Normal/ subdirectories with .bam files
}

EOF
    exit 1
}


PATIENT=
RAW_DIR=
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
        --data_dir) RAW_DIR="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

if [ -z "$RAW_DIR" ]; then
   RAW_DIR=$(dirname "$(realpath "$SAMPLESHEET")") # assume samplesheet in Raw/Patient directory
fi

chmod +x $script

# Create a temporary file to store the modified sample sheet
TEMP_SAMPLESHEET=$(mktemp)

echo "patient,sample,fastq1,fastq2,bam,bai,status" >> "$TEMP_SAMPLESHEET"


while IFS=$',' read -r patient sample fastq1 fastq2 status; do
    
    if [[ $patient == $PATIENT ]]; then

        job_name="fastq_to_bam_$sample"
        mkdir -p $RAW_DIR/$sample
        output_prefix=$RAW_DIR/$sample/$sample

        bsub -J ${job_name} \
                -P ${project} \
                -n ${cores} \
                -M 32000 \
                -R span[hosts=1] \
                -R "rusage[mem=4000]" \
                -W 12:00 \
                -q premium \
                -oo ${output_prefix}.out \
                -eo ${output_prefix}.err \
                ${script} -v -r ${fastq1},${fastq2} -o ${output_prefix} -p ${PATIENT} --threads ${cores} --post_process --step 3

        # Overwrite bam and bai with new paths
        bam="${output_prefix}.bam"
        bai="${output_prefix}.bai"

        # Write the updated line to the temporary file
        echo "$patient,$sample,$fastq1,$fastq2,$bam,$bai,$status" >> "$TEMP_SAMPLESHEET"
    
    fi

done < "$SAMPLESHEET"

# Replace the original sample sheet with the modified one
mv "$TEMP_SAMPLESHEET" "$SAMPLESHEET"

#!/bin/bash

#############################################
#specify path to auxilary scripts
scripts='/hpc/users/mazeev01/matt_lynch/scripts'
############################################

# Function to print progress with timestamp
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] --patient <patient_id> --data <data_info_file> --reference <reference.fasta> --samples <s1,s2,...>

Required Arguments:
    --patient <patient_id>           Patient id
    --data <data_info_file>          File with patient raw data information
    --reference <reference.fasta>    Refence genome 
    --samples <s1,s2,...>            List of sample ids

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --vcf_overwrite                 Overwrite counts in sample vcf file with counts in sample bam file
    --report                        Report vcf / bam count discrepancies
    --keep_intermediate             Keep intermediate file          
}

EOF
    exit 1
}


# Variables
PATIENT=
DATA_FILE=
REF_FASTA=
SAMPLES=
VERBOSE=0
VCF_OVERWRITE=0
REPORT=0
KEEP_INTERMEDIATE=0

# Parse command-line argument
while [[ "$#" -gt 0 ]]; do
    case "$1" in 
        -h) usage ;;
        -v) VERBOSE=1 ;;
        --keep_intermediate) KEEP_INTERMEDIATE=1 ;;
        --report) REPORT=1 ;;
        --vcf_overwrite) VCF_OVERWRITE=1 ;;
        --patient) PATIENT="$2"; shift ;;
        --data) DATA_FILE="$2"; shift ;;
        --reference) REF_FASTA="$2"; shift ;;
        --samples) SAMPLES="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

# Check that mandatory arguments are provided
if [ -z "$PATIENT" ] || [ -z "$DATA_FILE" ] || [ -z "$REF_FASTA" ] || [ -z "$SAMPLES" ]; then
    echo "Error: Not all required arguments provided."
    usage
fi

# Split SAMPLES into array
IFS=',' read -r -a SAMPLE_ARRAY <<< "$SAMPLES"


# Check data file
if [ ! -f "$DATA_FILE" ]; then
    echo "Error: Data info file not found."
    exit 1
fi



#####################################################################################
# Run the pipeline

# Load required module
print_progress "Loading required modules..."
echo ""
module --force purge
module load python
module load bam-readcount >/dev/null 2>&1

print_progress "Runing preprocessing python script." 
echo ""

#python ${scripts}/union_variants_pre.py -patient_id ${PATIENT} -file_info ${DATA_FILE} -samples ${SAMPLE_ARRAY}


print_progress "Prepping for bam-readcounts"
echo ""

# Extract BAM files, BAM counts files, and regions files for the specified patient ID
while IFS=$'\t' read -r patient directory sample fastq vcf bam bamcounts variants regions; do
    if [[ $patient == "$PATIENT" ]]; then  # Check if current line matches the patient ID

        if [[ ",${SAMPLES}," == *",${sample},"* ]]; then
        
            sample_dir=$directory/$sample
            patient_dir=$directory

            if [ $VERBOSE -eq 1 ]; then
                echo "--------------------------"
                echo "Sample: $sample"
                echo "BAM file: $bam"
                echo "BAM counts file: $bamcounts"
                echo "Regions file: $regions"
                echo "Patient directory: $patient_dir"
                echo "Sample directory: $sample_dir"
                echo "--------------------------"
            fi

            # Run the bam-readcounts job in the background
            print_progress "Running bam-readcounts for ${sample}"

            bam-readcount -w1 -i -f $REF_FASTA $sample_dir/$bam -l $patient_dir/$regions > $sample_dir/$bamcounts &
            echo ""
    
        fi
    fi
done < "$DATA_FILE"


# Wait for all jobs to complete
print_progress "Waiting for all bam-readcounts jobs to complete..."
wait

# Continue 
print_progress "All bam-readcounts jobs are complete. Running post-processing python script..."

python ${scripts}/union_variants_post.py -patient_id ${PATIENT} -file_info ${DATA_FILE} -samples ${SAMPLE_ARRAY}

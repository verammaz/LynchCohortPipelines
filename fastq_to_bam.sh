#!/bin/bash

##############################################
#specify file path for sites_of_variation.vcf 
SITES_OF_VARIATION="sites_of_variation.vcf" 
##############################################

# Function to print progress with timestamp
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

# Load required modules
module purge
print_progress "Loading required modules..."
module load samtools >/dev/null 2>&1
module load bwa >/dev/null 2>&1
module load gatk >/dev/null 2>&1
module load java/1.8.0_66 >/dev/null 2>&1
module load picard/2.2.4 >/dev/null 2>&1

# Exit immediately if a command exits with a non-zero status.
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [options] -f <reference.fa> -r <read1.fastq,read2.fastq> -o <output_prefix>

Required Arguments:
  -f <reference.fa>                Reference genome file in FASTA format.
  -r <read1.fastq,read2.fastq>     Comma-separated list of two FASTQ files (no space).
  -o <output_prefix>               Prefix for output files.

Options:
  -h                               Display this message.
  -v                               Enable verbose mode.
  --run_index                      Run bwa_index step (required for bwa_mem alignment, once per reference genome)
  --keep_intermediate              Keep intermediate files generated during the pipeline.
  --view_final                     Index refernece and final BAM file for IGV viewing.
  --post_process                   Carry out post processing of BAM file after alignment and sorting (mark duplicates, base quality recallibration)


EOF
    exit 1
}

# Variables
READS=()
REFERENCE=
OUTPUT_PREFIX=
VIEW=0
VERBOSE=0
INDEX=0
POST_PROCESS=0
KEEP_INTERMEDIATE=0

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        --keep_intermediate) KEEP_INTERMEDIATE=1 ;;
        --view_final) VIEW=1 ;;
        --post_process) POST_PROCESS=1 ;;
        --run_index) INDEX=1 ;;
        -f) REFERENCE="$2"; shift ;;
        -r) READS="$2"; shift ;;
        -o) OUTPUT_PREFIX="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

# Check that mandatory arguments are provided
if [ -z "$REFERENCE" ] || [ -z "$READS" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Error: Not all required arguemnts provided."
    usage
fi

# Split the READS argument into an array
IFS=',' read -r -a READ_ARRAY <<< "$READS"

# Ensure exactly two read files are provided
if [ ${#READ_ARRAY[@]} -ne 2 ]; then
    echo "Error: Exactly two read files must be specified, separated by a comma."
    usage
fi

READS_1=${READ_ARRAY[0]}
READS_2=${READ_ARRAY[1]}


# If verbose mode is enabled, print the parameters
if [ $VERBOSE -eq 1 ]; then

    # Convert 1/0 to true/false for printing
    INDEX_STR=$( [ $INDEX -eq 1 ] && echo "true" || echo "false" )
    KEEP_INTERMEDIATE_STR=$( [ $KEEP_INTERMEDIATE -eq 1 ] && echo "true" || echo "false" )
    POST_PROCESS_STR=$( [ $POST_PROCESS -eq 1 ] && echo "true" || echo "false" )
    VIEW_STR=$( [ $VIEW -eq 1 ] && echo "true" || echo "false" )

    echo "---------------------------------------------"
    echo "Reference: $REFERENCE"
    echo "Read1: $READS_1"
    echo "Read2: $READS_2"
    echo "Output Prefix: $OUTPUT_PREFIX"
    echo "Run Indexing: $INDEX_STR"
    echo "Post Processing: $POST_PROCESS_STR"
    echo "Keep Intermediate Files: $KEEP_INTERMEDIATE_STR"
    echo "View Final File: $VIEW_STR"
    echo "---------------------------------------------"
fi


# Run fastq-->bam pipeline

echo "---------------------------------------------"
print_progress "Starting the pipeline..."
echo "---------------------------------------------"

set -o pipefail 

if [ $INDEX -eq 1 ]; then
    print_progress "Indexing the reference..."
    bwa index -a bwtsw $REFERENCE
fi

# Step 1: BWA MEM and Samtools sort
print_progress "Aligning and sorting..."
bwa mem -t 8 $REFERENCE $READS_1 $READS_2 | samtools sort -l 1 -@8 -o "${OUTPUT_PREFIX}_sorted.bam" -
wait  

if [ $POST_PROCESS -eq 1 ]; then

    # Step 2: Picard MarkDuplicates
    print_progress "Marking duplictaes (Picard MarkDuplicates)..."
    java -jar $PICARD MarkDuplicates \
        I="${OUTPUT_PREFIX}_sorted.bam" \
        O="${OUTPUT_PREFIX}_marked_dup.bam" \
        M="${OUTPUT_PREFIX}_marked_dup_metrics.txt"
    wait  

    if [ -f "$SITES_OF_VARIATION" ]; then

        # Step 3: GATK BaseRecalibrator
        print_progress "Recalibrating bases (GATK BaseRecalibrator)..."
        gatk BaseRecalibrator \
            -I "${OUTPUT_PREFIX}_marked_dup.bam" \
            -R $REFERENCE \
            --known-sites $SITES_OF_VARIATION \
            -O "${OUTPUT_PREFIX}_recal_data.table"
        wait  

        # Step 4: GATK ApplyBQSR
        print_progress "Applying base recalibaration (GATK ApplyBQSR)..."
        gatk ApplyBQSR \
            -R $REFERENCE \
            -I "${OUTPUT_PREFIX}_marked_dup.bam" \
            --bqsr-recal-file "${OUTPUT_PREFIX}_recal_data.table" \
            -O "${OUTPUT_PREFIX}_final.bam"
    
    else
        echo "Sites of variation file $SITES_OF_VARIATION not found."
        echo "Cannot run GATK BaseRecalibrator without sites_of_variation.vcf file."
        echo "Terminating post processing after Picard MarkDuplicates..."
        mv "${OUTPUT_PREFIX}_marked_dup.bam" "${OUTPUT_PREFIX}_final.bam"
    fi
else
    mv "${OUTPUT_PREFIX}_sorted.bam" "${OUTPUT_PREFIX}_final.bam"
fi

print_progress "Pipeline completed successfully. Output is in ${OUTPUT_PREFIX}_final.bam."


# Remove intermediate files if necessary
if [ $KEEP_INTERMEDIATE -eq 0 ] && [ $POST_PROCESS -eq 1 ]; then
    print_progress "Removing intermediate files..."
    rm -f "${OUTPUT_PREFIX}_sorted.bam", "${OUTPUT_PREFIX}_marked_dup.bam", "${OUTPUT_PREFIX}_recal_table.table", "${OUTPUT_PREFIX}_marked_dup_metrics.txt"
fi

# View the final BAM file if the -v flag is set
if [ $VIEW -eq 1 ]; then
    print_progress "Producing .fai files required for IGV viewing..."
    samtools index "${OUTPUT_PREFIX}_final.bam" 
    if [-f "$REFERENCE.fai" ]; then
        samtools faidx $REFERENCE
    fi
fi

exit 1

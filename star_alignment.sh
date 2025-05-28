#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if a command exits with a non-zero status.
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [options] -f <reference.fa> -r <read1.fastq,read2.fastq> -o <output_prefix>

Required Arguments:
  -r <read1.fastq,read2.fastq>     Comma-separated list of two FASTQ files.
  -o <output_prefix>               Prefix for output files.

Options:
  -h                               Display this message.
  -v                               Enable verbose mode.
  --patient                        Patient id.
  --index_ref                      Run bwa_index step.
  --threads                        Number of threads to use.
  --sample                         Sample id


EOF
    exit 1
}

# Variables
PATIENT=
READS=()
OUTPUT_PREFIX=
VERBOSE=0
INDEX=0
THREADS=8
SAMPLE=

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        --index_ref) INDEX=1 ;;
        --threads) THREADS="$2"; shift ;;
        --patient) PATIENT="$2"; shift ;;
        --sample) SAMPLE="$2"; shift ;;
        -r) READS="$2"; shift ;;
        -o) OUTPUT_PREFIX="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done


# Check reference genome file
if [ ! -f "${REF_FASTA}" ]; then
    echo "Reference file ${REF_FASTA} not found!"
    exit 1
fi


GENOME_DIR='/sc/arion/scratch/mazeev01/STAR'
create_directory "$GENOME_DIR"

# If verbose mode is enabled, print the parameters
if [ $VERBOSE -eq 1 ]; then

    # Convert 1/0 to true/false for printing
    INDEX_STR=$( [ $INDEX -eq 1 ] && echo "true" || echo "false" )
    
    echo "---------------------------------------------"
    echo "Reference: $REF_FASTA"
    echo "Read1: $READS_1"
    echo "Read2: $READS_2"
    echo "Output Prefix: $OUTPUT_PREFIX"
    echo "Run Indexing: $INDEX_STR"
    echo "Genome Directory: $GENOME_DIR"
    echo "---------------------------------------------"
fi


#############################################################################################
#Run STAR pipeline


echo "---------------------------------------------"
print_progress "Starting the pipeline..."
echo "---------------------------------------------"


if [ $INDEX -eq 1 ]; then
# As indexing reference takes time and is only required once per reference, double check that indexing step is required:)
    if [ $INDEX -eq 1 ]; then
        print_progress "Indexing the reference..."
        STAR --runMode genomeGenerate \
             --runThreadN ${THREADS} \
             --genomeDir ${GENOME_DIR} \
             --genomeFastaFiles ${REF_FASTA} \
             --sjdbGTFfile ${GTF_FILE} \
             --sjdbOverhang 100  # Adjust depending on read length - 1
    echo "Genome index generation complete."
    else
        echo "STAR genome index already exists. Skipping generation!"
    fi
fi

exit 0
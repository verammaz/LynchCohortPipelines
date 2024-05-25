#!/bin/bash

# Load required modules
module load samtools
module load minimap2

# Exit immediately if a command exits with a non-zero status.
set -e

# Function to print usage
usage() {
    printf "Usage: $0 [-v] -f <reference.fa> -r <read1.fastq> -r <read2.fastq> -o <output_prefix>\n
    \t [-v] flag for viewing alignment immediately \n
    \t -f reference genome file \n 
    \t -r read fastq files can be zipped or unzipped (two files required) \n
    \t -o output file prefix\n"
    exit 1
}

# Variables
READS=()
REFERENCE=
OUTPUT_PREFIX=
VIEW="false"
BAM=

# Parse command line arguments
while getopts 'f:r:o:vh' option; do
    case "$option" in
        f)
            REFERENCE="$OPTARG"
            ;;
        r)
            READS+=("$OPTARG")
            ;;
        o)
            OUTPUT_PREFIX="$OPTARG"
            ;;
        v)
            VIEW=true
            ;;
        h | *)
            usage
            ;;
    esac
done

# Ensure all required arguments are provided
if [[ -z $REFERENCE || ${#READS[@]} -ne 2 || -z $OUTPUT_PREFIX ]]; then
    echo "Not all required arguments provided"
    usage
fi


# Check if output prefix includes .bam
if [[ "$OUTPUT_PREFIX" == *.bam ]]; then
   BAM=$OUTPUT_PREFIX
else
    BAM="${OUTPUT_PREFIX}.bam"
fi

# Function to print the current progress
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

# Function to unzip fastq file
handle_unzipping() {
    local read_file=$1
    local unzipped_file=$2

    if [[ "$read_file" == *.gz ]]; then
            print_progress "Unzipping $read_file..."
            gunzip -c "$read_file" > "$unzipped_file"
    else
        unzipped_file=$read_file
        print_progress "$unzipped_file already exists. Skipping unzipping."
    fi
}


# Check if the final BAM file already exists
if [ -f "$FINAL_BAM" ]; then
    print_progress "Final BAM file $BAM already exists."

# Run fastq-->bam pipline
else
    READS_1=${READS[0]}
    READS_2=${READS[1]}

    print_progress "Starting the pipeline: Aligning reads, fixing mates, sorting, marking duplicates, and creating final BAM..."

    set -o pipefail 
    
    minimap2 -t 8 -a -x sr $REFERENCE $READS_1 $READS_2 | \
    samtools fixmate -u -m - - | \
    samtools sort -u -@8 -T /tmp/example_prefix - -o $BAM | \
    samtools markdup -@8 - -O bam,level=1 $BAM

    print_progress "Pipeline completed successfully. Output is in $BAM."
fi

# View the final BAM file if the -v flag is set
if [ "$VIEW" = true ]; then
    print_progress "Viewing the final BAM file..."
    samtools view -h $BAM | less -S
fi

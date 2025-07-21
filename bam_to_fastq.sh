#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if a command exits with a non-zero status.
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [options] -b <input.bam> -r1 <read1.fastq> -r2 <read2.fastq> 

Required Arguments:                
  -r1 <read1.fastq>                 First end of paired FASTQ output.
  -r2 <read2.fastq>                 Second end of paired FASTQ output.
  -b <input.bam>                    Input bam file from which (paired-end) reads will be extracted.



EOF
    exit 1
}

BAM=
READS1=
READS2=

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -r1) READS1="$2"; shift ;;
        r2) READS2="$2"; shift ;; 
        -b) BAM="$2"; shift ;;
        *) echo "Error: Unknown argument/option: $1" ; usage ;;
    esac
    shift
done

# Check that mandatory arguments are provided
if [ -z "$READS1" ] || [ -z "$READS2" ] || [ -z "$BAM" ] ; then
    echo "Error: Not all required arguments provided."
    usage
fi

module load samtools v1.17
module load java v1.8
module load picard v2.2.4

java -jar $PICARD_JAR SamToFastq \
     I=$BAM \
     FASTQ=$READ1 \ 
     SECOND_END_FASTQ=$READ2 




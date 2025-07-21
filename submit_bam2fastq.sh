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
        -r2) READS2="$2"; shift ;; 
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


job_name="bam2fastq"
script="$LYNCH/bam_to_fastq.sh"

echo "Script: $script"
echo "BAM: $BAM"
echo "Reads: $READS1, $READS2"

chmod +x $script

module purge
module load singularity/3.6.4

bsub -J "$job_name" \
    -P "$project" \
    -n 8 \
    -R "span[hosts=1]" \
    -R "rusage[mem=40000]" \
    -W 6:00 \
    -q premium \
    -oo "${LOG_DIR}/${job_name}.out" \
    -eo "${LOG_DIR}/${job_name}.err" \
    singularity exec "${CONTAINER_FASTQ2BAM}" ${script} -b ${BAM} -r1 ${READS1} -r2 ${READS2}
#!/bin/bash

REFERENCE='/sc/arion/projects/FLAI/vera/matt_lynch/References/ref_b37/human_g1k_v37_decoy.fasta'
scripts='/hpc/users/mazeev01/matt_lynch/scripts'

# Function to print progress with timestamp
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

# Load required module
module load python

id=$1
file_info=$2

print_progress "Runing preprocessing python script." 
echo ""

python ${scripts}/union_variants_pre.py -patient_id ${id} -file_info ${file_info}


print_progress "Prepping for bam-readcounts"
echo ""

# Extract BAM files, BAM counts files, and regions files for the specified patient ID
job_ids=()  # Array to hold job IDs

while IFS=$'\t' read -r Patient Directory Sample FASTQ VCF BAM BAM_COUNTS VARIANTS REGIONS; do
    if [[ $Patient == "$id" ]]; then  # Check if current line matches the patient ID
        echo "--------------------------"
        echo "Sample: $Sample"
        echo "BAM file: $BAM"
        echo "BAM counts file: $BAM_COUNTS"
        echo "Regions file: $REGIONS"

        DIR_SAMPLE=$Directory/$Sample
        DIR_PATIENT=$Directory

        echo "Patient directory: $DIR_PATIENT"
        echo "Sample directory: $DIR_SAMPLE"
        echo ""

        # Run the bam-readcounts job in the background
        print_progress "Starting bam-readcounts for ${Sample}"

        $scripts/run_bamcounts.sh $DIR_SAMPLE/$BAM $DIR_PATIENT/$REGIONS $DIR_SAMPLE/$BAM_COUNTS $REFERENCE &
    fi
done < "$file_info"


# Wait for all jobs to complete
print_progress "Waiting for all bam-readcounts jobs to complete..."

wait

print_progress "All bam-readcounts jobs are complete. Running post-processing python script..."

python ${scripts}/union_variants_post.py -patient_id ${id} -file_info ${file_info}

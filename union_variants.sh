#!/bin/bash

source ./config.sh

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] --patient <patient_id> --samplesheet <samplesheet.csv> --reference <reference.fasta>

Required Arguments:
    --patient <patient_id>           Patient id
    --samples <s1,s2,...>            List of sample ids

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --use_vcf                       Use VCF counts for called sample variants
    --report                        Report vcf / bam count discrepancies
    --keep_intermediate             Keep intermediate files
    --filter_variants               Apply additional filters to variants (other than 'PASS')
    --strelka_mutect_snv_intersect  Only consider snv variants at intersection of strelka mutect callers
    --single_output_file            Write single output vcf file
    --zero_coverage_ok              Include variants even if some sample data has zero coverage at that position
    --reference <reference.fasta>   Refence genome 
}

EOF
    exit 1
}


# Variables
PATIENT=
SAMPLES=
VERBOSE=0
VCF=0
REPORT=0
KEEP_INTERMEDIATE=0
FILTER_VARIANTS=0
SNV_INTERSECT=0
INDEL_INTERSECT=0
SINGLE_FILE=0
ZERO_COVERAGE=0

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in 
        -h) usage ;;
        -v) VERBOSE=1 ;;
        --keep_intermediate) KEEP_INTERMEDIATE=1 ;;
        --report) REPORT=1 ;;
        --vcf_overwrite) VCF=1 ;;
        --filter_variants) FILTER_VARIANTS=1 ;;
        --strelka_mutect_snv_intersect) SNV_INTERSECT=1 ;;
        --strelka_mutect_indel_instersect) INDEL_INTERSECT=1 ;;
        --single_output_file) SINGLE_FILE=1 ;;
        --zero_coverage_ok) ZERO_COVERAGE=1 ;;
        --patient) PATIENT="$2"; shift ;;
        --samples) SAMPLES="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

# Check that mandatory arguments are provided
if [ -z "$PATIENT" ] || [ -z "$SAMPLES" ]; then
    echo "Error: Not all required arguments provided."
    usage
fi

RAW_DIR="$HOME_DIR/Raw/$PATIENT"


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

python ${scripts}/union_variants_pre.py -patient_id ${PATIENT} -data_dir ${RAW_DIR} -samples ${SAMPLES} \
                                        --additional_filter ${FILTER_VARIANTS} \
                                        --strelka_mutect_snv_intersect ${SNV_INTERSECT} \
                                        --strelka_mutect_indel_intersect ${INDEL_INTERSECT}

regions = "${RAW_DIR}/${PATIENT}_regions.txt"

# Split the SAMPLES argument into an array
IFS=',' read -r -a SAMPLES_ARRAY <<< "$SAMPLES"

for sample_id in $SAMPLE_ARRAY; do

    sample_samplesheet="${RAW_DIR}/samplesheet_${sample_id}.csv"

    while IFS=$',' read -r patient sample fastq1 fastq2 bam bai status vcf; do
   
        if [[ "${sample}" != "Normal" ]]; then

            bamcounts="${RAW_DIR}/${sample}/${sample}_bamcounts.txt"

            if [ $VERBOSE -eq 1 ]; then
                echo "--------------------------"
                echo "Sample: $sample"
                echo "BAM file: $bam"
                echo "BAM counts file: $bamcounts"
                echo "Regions file: $regions"
                echo "--------------------------"
            fi

            # Run the bam-readcounts job in the background
            print_progress "Running bam-readcounts for ${sample}"
            echo ""
            bam-readcount -w1 -i -f $REF_FASTA $bam -l $regions > $bamcounts -b 15 &
    
        fi
    done < "$sample_samplesheet"

done


# Wait for all jobs to complete
print_progress "Waiting for all bam-readcounts jobs to complete..."
echo ""
wait

# Continue 
print_progress "All bam-readcounts jobs are complete. Running post-processing python script..."
echo ""

python ${scripts}/union_variants_post.py -patient_id ${PATIENT} -samples ${SAMPLES} -data_dir {$RAW_DIR} --zero_coverage_ok ${ZERO_COVERAGE} --single_combined_vcf ${SINGLE_FILE}

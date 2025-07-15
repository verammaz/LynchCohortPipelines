#!/bin/bash

source ./config.sh

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -p <patient_id> -s <samplesheet.csv>

Required Arguments:
    -p <patient_id>                     Patient id
    -s <samplesheet.csv>                Sample info sheet

Options:
    -h                                  Display this message
    -v                                  Enable verbode mode
    --keep_intermediate                 Keep intermediate files
    --filter_variants                   Apply additional filters to variants (other than 'PASS')
    --strelka_mutect_snv_intersect      Only consider snv variants at intersection of strelka mutect callers
    --strelka_mutect_indel_intersect    Only consider snv variants at intersection of strelka mutect callers
    --single_output_file                Write single output vcf file
    --zero_coverage_ok                  Include variants even if some sample data has zero coverage at that position
    --step                              Step to start (0=pre, 1=bam-readcount, 2=post)
}

EOF
    exit 1
}


# Variables
PATIENT=
SAMPLES=
VERBOSE=0
VCF=0
KEEP_INTERMEDIATE=0
FILTER_VARIANTS=0
SNV_INTERSECT=0
INDEL_INTERSECT=0
SINGLE_FILE=0
ZERO_COVERAGE=0
STEP=0

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in 
        -v) VERBOSE=1 ;;
        --keep_intermediate) KEEP_INTERMEDIATE=1 ;;
        --report) REPORT=1 ;;
        --vcf_overwrite) VCF=1 ;;
        --filter_variants) FILTER_VARIANTS=1 ;;
        --strelka_mutect_snv_intersect) SNV_INTERSECT=1 ;;
        --strelka_mutect_indel_instersect) INDEL_INTERSECT=1 ;;
        --single_output_file) SINGLE_FILE=1 ;;
        --zero_coverage_ok) ZERO_COVERAGE=1 ;;
        --step) STEP="$2"; shift ;;
        -p) PATIENT="$2"; shift ;;
        -s) SAMPLESHEET="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

# Check that mandatory arguments are provided
if [ -z "$PATIENT" ] || [ -z "$SAMPLESHEET" ]; then
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


if [[ "$STEP" -eq 0  ]]; then
    print_progress "Running preprocessing python script." 
    echo ""

    python ${LYNCH}/union_variants_pre.py -patient_id ${PATIENT} -samplesheet ${SAMPLESHEET} -data_dir ${RAW_DIR} \
                                            --additional_filter ${FILTER_VARIANTS} \
                                            --strelka_mutect_snv_intersect ${SNV_INTERSECT} \
                                            --strelka_mutect_indel_intersect ${INDEL_INTERSECT}
    STEP=1
fi


SAMPLE_ARRAY=()

if [[ "$STEP" -ge 1 ]]; then
    regions="${RAW_DIR}/${PATIENT}_regions.txt"


    {
    read  # Skip the header line
    while IFS= read -r line || [[ -n "$line" ]]; do
        
        # Skip empty lines
        [[ -z "$line" ]] && continue

    # Split the line into fields using awk to handle potential edge cases
    IFS=',' read -r patient sample fastq1 fastq2 status bam bai ref <<< "$(awk -F',' '{print $1,$2,$3,$4,$5,$6,$7}' OFS=',' <<< "$line")"

        if [[ $patient == $PATIENT ]] && [[ "${sample}" != "Normal" ]]; then
            
            SAMPLE_ARRAY+=($sample)

            bamcounts="${RAW_DIR}/${sample}/${sample}_bamcounts.txt"

            if [ $VERBOSE -eq 1 ]; then
                echo "--------------------------"
                echo "Sample: $sample"
                echo "BAM file: $bam"
                echo "BAM counts file: $bamcounts"
                echo "Regions file: $regions"
                echo "--------------------------"
            fi

            if [[ "$STEP" -eq 1 ]]; then

                # Run the bam-readcounts job in the background
                print_progress "Running bam-readcounts for ${sample}"
                echo ""
                bam-readcount -w1 -i -f $REF_FASTA $bam -l $regions -b 15 > $bamcounts &
            fi

        fi

    done
    } < "$SAMPLESHEET"



    # Wait for all jobs to complete
    print_progress "Waiting for all bam-readcounts jobs to complete..."
    echo ""
    wait

    STEP=2
fi
    # Continue 
    print_progress "All bam-readcounts jobs are complete. Running post-processing python script..."
    echo ""


if [[ "$STEP" -eq 2 ]]; then
    python ${LYNCH}/union_variants_post.py -patient_id ${PATIENT} -samples ${SAMPLE_ARRAY[@]} -data_dir ${RAW_DIR} \
                                           --zero_coverage_ok ${ZERO_COVERAGE} \
                                           --single_combined_vcf ${SINGLE_FILE}
fi
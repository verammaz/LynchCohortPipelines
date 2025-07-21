#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -s <samplesheet.csv> -p <patient_id>

Required Arguments:
    -s <samplesheet.csv>            Configured input file
    -p <patient_id>                 Patient identifier

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --step                          Step to start from (0=alignment, 1=markdup, 2=indelrealign, 3=baserecal)
    --ref                           Reference genome (default=hg19)

EOF
    exit 1
}


PATIENT=
SAMPLESHEET=
VERBOSE=0
STEP=0
REF="hg19"

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -s) SAMPLESHEET="$2"; shift ;;
        --step) STEP="$2"; shift ;;
        --ref) REF="$2"; shift ;; 
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done


REF_FASTA=

if [[ "$REF" == 'hg19' ]]; then
    REF_FASTA=$REF_FASTA_hg19
fi

if [[ "$REF" == 'hg38' ]]; then
    REF_FASTA=$REF_FASTA_hg38
fi

if [[ -z "$REF_FASTA" ]]; then
    echo "Error: REF_FASTA is not set. Check that REF is either 'hg19' or 'hg38', and corresponding variables (REF_FASTA_hg19 / REF_FASTA_hg38) are defined in your config."
    exit 1
fi

POST_PROCESS_FLAG=""
if [[ "$REF" == "hg19" ]]; then
    POST_PROCESS_FLAG="--post_process"
fi


RAW_DIR="$HOME_DIR/Raw/$PATIENT"

script="$LYNCH/fastq_to_bam.sh"
chmod +x $script

# Create a temporary file to store the modified sample sheet
TEMP_SAMPLESHEET=$(mktemp)

echo "patient,sample,fastq_1,fastq_2,status,bam,bai,ref" >> "$TEMP_SAMPLESHEET"


module purge
module load singularity/3.6.4

{
read  # Skip the header line
while IFS= read -r line || [[ -n "$line" ]]; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Split the line into fields using awk to handle potential edge cases
    IFS=',' read -r patient sample fastq1 fastq2 status <<< "$(awk -F',' '{print $1,$2,$3,$4,$5}' OFS=',' <<< "$line")"
    
    echo $sample
    if [[ $patient == $PATIENT ]] || [[ -z $PATIENT ]]; then

        job_name="fastq_to_bam_$patient_$sample"
        create_directory "$RAW_DIR/$sample"
        output_prefix="$RAW_DIR/$sample/$sample"

        bsub -J "${job_name}" \
             -P "${project}" \
             -n "${cores}" \
             -M 32000 \
             -R span[hosts=1] \
             -R "rusage[mem=4000]" \
             -W 40:00 \
             -q premium \
             -oo "${LOG_DIR}/${job_name}.out" \
             -eo "${LOG_DIR}/${job_name}.err" \
             singularity exec "${CONTAINER_FASTQ2BAM}" "${script}" \
                    -r "${fastq1},${fastq2}" \
                    -o "${output_prefix}" \
                    --patient "${PATIENT}" \
                    --sample "${sample}" \
                    --step "${STEP}" \
                    --threads "${cores}" \
                    $POST_PROCESS_FLAG \
                    --index_ref \
                    -f "${REF_FASTA}"

        # Overwrite bam and bai with new paths
        bam="${output_prefix}.bam"
        bai="${output_prefix}.bai"

        # Write the updated line to the temporary file
        echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,$REF" >> "$TEMP_SAMPLESHEET"
    
    fi

done
} < "$SAMPLESHEET"

# Replace the original sample sheet with the modified one
mv "$TEMP_SAMPLESHEET" "$SAMPLESHEET"

#!/bin/bash

# This script assumes that you have processed bam files, and runs
# the nextflow/nf-core/sarek pipeline from the variant calling step.

######################### Change the following! ########################################
# nextflow executable
NEXTFLOW='/hpc/users/mazeev01/nextflow'
# pon file for mutect variant filtering
PON='/sc/arion/projects/FLAI/vera/References/Mutect2-exome-panel.vcf.gz'
# work directory
WORK_DIR='/sc/arion/scratch/mazeev01/variant_call'
# output directory
OUT_DIR='/sc/arion/work/mazeev01/variant_call'
######################################################################################

# Function to print progress with timestamp
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] --patient <patient_id> --samplesheet <samplesheet.csv>

Required Arguments:
    --patient <patient_id>           Patient id
    --samplesheet <samplesheet.csv>  File with raw data configured

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --out_dir                       Output directory for nextflow vcfs
    --raw_data_dir                  Path to patient raw data; Sample/ and Normal/ subdirectories with bam files; vcf files will be put in respective sample subdirectories
}

EOF
    exit 1
}



PATIENT=
SAMPLESHEET=
DATA_DIR=
VERBOSE=0

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        --patient) PATIENT="$2"; shift ;;
        --samplesheet) SAMPLESHEET="$2"; shift ;;
        --out_dir) OUT_DIR="$2"; shift ;;
        --raw_data_dir) DATA_DIR="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

if [ -z "$DATA_DIR" ]; then
    DATA_DIR=$(dirname "$(realpath "$SAMPLESHEET")")
fi


################################################################################################
# Nextflow variant calling
################################################################################################

print_progress "Launching nextflow to perform variant calling..."

module add java/11.0.2
module add singularity/3.2.1
module add proxies

$NEXTFLOW run nf-core/sarek \
-profile singularity \
-workDir $WORK_DIR \
--wes \
--input $SAMPLESHEET \
--outdir $OUT_DIR  \
--step variant_calling \
--tools strelka,mutect2 \
--pon $PON \
--only_paired_variant_calling \
--genome GATK.GRCh37 \
--max_cpus 24

# Move VCF outputs from nextflow to Raw/Patient directory
print_progress "Variant calilng complete. Getting VCF files."

strelka_dir="${OUT_DIR}/variant_calling/strelka"
mutect_dir="${OUT_DIR}/variant_calling/mutect2"

# Create a temporary file to store the modified sample sheet
TEMP_SAMPLESHEET=$(mktemp)

echo "patient,sample,fastq1,fastq2,bam,bai,status,vcf" >> "$TEMP_SAMPLESHEET"

while IFS=$',' read -r patient sample fastq1 fastq2 bam bai status; do
    if [[ $patient == $PATIENT ]]; then

        if [[ $sample != "Normal" ]] && [[ $status == "1" ]]; then #IMPORTANT: assumes non-tumor sample is named 'Normal'

            strelka_files="${strelka_dir}/${sample}_vs_Normal/*.vcf.gz"
            mutect_files="${mutext_dir}/${sample}_vs_Normal/*.vcf.gz"
            
            mv $strelka_files $DATA_DIR/*$sample
            mv $mutect_files $DATA_DIR/*$sample

            vcf_files="${strelka_files}|${mutect_files}"

            echo "$patient,$sample,$fastq1,$fastq2,$bam,$bai,$status,$vcf_files" >> "$TEMP_SAMPLESHEET"    

        else
            echo "$patient,$sample,$fastq1,$fastq2,$bam,$bai,$status,na" >> "$TEMP_SAMPLESHEET"    
        fi
    fi
done < $SAMPLESHEET


# Replace the original sample sheet with the modified one
mv "$TEMP_SAMPLESHEET" "$SAMPLESHEET"
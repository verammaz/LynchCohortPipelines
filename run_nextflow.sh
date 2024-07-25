#!/bin/bash

######################### Specify some paths ########################################
# nextflow executable
NEXTFLOW='/hpc/users/mazeev01/nextflow'
# pon file for mutect variant filtering
PON='/sc/arion/projects/FLAI/vera/References/Mutect2-exome-panel.vcf.gz'
# output directory
OUT_DIR='/sc/arion/work/mazeev01/variant_call'
######################################################################################

# Function to print progress with timestamp
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

set -e

PATIENT=
SAMPLESHEET=
SAMPLE=
STEP=

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --patient) PATIENT="$2"; shift ;;
        --samplesheet) SAMPLESHEET="$2"; shift ;;
        --sample) SAMPLE="$2"; shift ;; 
        --step) STEP="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; exit ;;
    esac
    shift
done

DATA_DIR=$(dirname "$(realpath "$SAMPLESHEET")")

SAMPLE_SAMPLESHEET="${DATA_DIR}/samplesheet_${SAMPLE}.csv"

> "$SAMPLE_SAMPLESHEET"

if [[ $STEP -eq 0 ]]; then
    echo "patient,sample,fastq_1,fastq_2,status" >> "$SAMPLE_SAMPLESHEET"
else 
    echo "patient,sample,fastq_1,fastq_2,status,bam,bai" >> "$SAMPLE_SAMPLESHEET"
fi

while IFS=$',' read -r patient sample fastq1 fastq2 bam bai status; do
    
    if [[[ $patient == $PATIENT ]] && [[ $sample == $SAMPLE ]]] ||  [[[ $patient == $PATIENT ]] && [[ "$sample" == "Normal" ]]]; then

        if [[ $STEP -eq 0 ]]; then
             echo "$patient,$sample,$fastq1,$fastq2,$status" >> "$SAMPLE_SAMPLESHEET"
        else
            echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai" >> "$SAMPLE_SAMPLESHEET"
        fi
    fi

done < "$SAMPLESHEET"

OUT_DIR="$OUT_DIR/$SAMPLE"
mkdir -p $OUT_DIR

print_progress "Launching nextflow to perform variant calling..."

module load java/11.0.2
module load singularity/3.2.1
module load proxies

step="variant_calling"

if [[ $STEP -eq 0 ]]; then
    step="mapping"
fi

$NEXTFLOW run nf-core/sarek \
            -profile singularity \
            --wes \
            --input $SAMPLE_SAMPLESHEET \
            --outdir $OUT_DIR  \
            --step $step \
            --skip_tools baserecalibrator_report,markduplicates_report \
            --save_mapped \
            --save_output_as_bam \
            --tools strelka,mutect2 \
            --pon $PON \
            --only_paired_variant_calling \
            --genome GATK.GRCh37 \
            --max_cpus 48 \

print_progress "Variant calilng complete."

# Create a temporary file to store the modified sample sheet
TEMP_SAMPLESHEET=$(mktemp)

echo "patient,sample,fastq1,fastq2,status,bam,bai,vcf" >> "$TEMP_SAMPLESHEET"

strelka_dir="${OUT_DIR}/variant_calling/strelka"
mutect_dir="${OUT_DIR}/variant_calling/mutect2"

while IFS=$',' read -r patient sample fastq1 fastq2 status bam bai; do

    if [[ $STEP -eq 0 ]] || [[ -z "$bam" ]] || [[ -z "$bai" ]]; then
       
        print_progress "Getting BAM files for ${sample}."

        # Move BAM outputs from nextflow to Raw/Patient directory (downstream analysis scripts assume specific directory structure)

        recal_dir="${OUT_DIR}/preprocessing/recalibrated/${sample}"
        bam_file="${recal_dir}/${sample}.recal.bam"
        bai_file="${recal_dir}/${sample}.recal.bam.bai"

        mv "$bam_file" "${DATA_DIR}/${sample}"
        mv "$bai_file" "${DATA_DIR}/${sample}"

        bam="${DATA_DIR}/${sample}/$(basename "$bam_file")"
        bai="${DATA_DIR}/${sample}/$(basename "$bai_file")"
    fi

    
    if [[ $sample != "Normal" ]] && [[ $status == "1" ]]; then #IMPORTANT: assumes non-tumor sample is named 'Normal'

        print_progress "Getting VCF files for ${sample}."

        strelka_indel="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_indels.vcf.gz"
        strelka_snv="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_snvs.vcf.gz"
        mutect="${mutect_dir}/${sample}_vs_Normal/${sample}_vs_Normal.mutect2.filtered.vcf.gz"

        strelka_indel_tbi="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_indels.vcf.gz.tbi"
        strelka_snv_tbi="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_snvs.vcf.gz.tbi"
        mutect_tbi="${mutect_dir}/${sample}_vs_Normal/${sample}_vs_Normal.mutect2.filtered.vcf.gz.tbi"
        
        mv "$strelka_indel" "${DATA_DIR}/${sample}"
        mv "$strelka_snv" "${DATA_DIR}/${sample}"
        mv "$mutect" "${DATA_DIR}/${sample}"

        mv "$strelka_indel_tbi" "${DATA_DIR}/${sample}"
        mv "$strelka_snv_tbi" "${DATA_DIR}/${sample}"
        mv "$mutect_tbi" "${DATA_DIR}/${sample}"

        vcf_files="$strelka_indel|$strelka_snv|$mutect"

        echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,$vcf_files" >> "$TEMP_SAMPLESHEET"    

    else
        echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,na" >> "$TEMP_SAMPLESHEET"    
    fi

done < $SAMPLE_SAMPLESHEET


# Replace the original sample sheet with the modified one
mv "$TEMP_SAMPLESHEET" "$SAMPLE_SAMPLESHEET"
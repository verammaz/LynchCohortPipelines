#!/bin/bash

# Get access to global variables
source ./config.sh

set -e

PATIENT=
SAMPLESHEET=
SAMPLE=
STEP=
GENOME="GATK.GRCh37"

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --patient) PATIENT="$2"; shift ;;
        --samplesheet) SAMPLESHEET="$2"; shift ;;
        --sample) SAMPLE="$2"; shift ;; 
        --step) STEP="$2"; shift ;;
        --genome) GENOME="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; exit ;;
    esac
    shift
done


RAW_DIR="$HOME_DIR/Raw/$PATIENT"

SAMPLE_SAMPLESHEET="${RAW_DIR}/samplesheet_${SAMPLE}.csv"

if [[ -f $SAMPLE_SAMPLESHEET ]]; then
    > "$SAMPLE_SAMPLESHEET"
fi

if [[ $STEP -eq 0 ]]; then
    echo "patient,sample,fastq_1,fastq_2,status,ref" >> "$SAMPLE_SAMPLESHEET"
else 
    echo "patient,sample,fastq_1,fastq_2,status,bam,bai,ref" >> "$SAMPLE_SAMPLESHEET"
fi


{
read  # Skip the header line
while IFS= read -r line || [[ -n "$line" ]]; do
    
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Split the line into fields using awk to handle potential edge cases
    IFS=',' read -r patient sample fastq1 fastq2 status bam bai ref <<< "$(awk -F',' '{print $1,$2,$3,$4,$5,$6,$7}' OFS=',' <<< "$line")"

    if [[ "$patient" == "$PATIENT" && ("$sample" == "$SAMPLE" || "$sample" == "Normal") ]]; then

        if [[ $STEP -eq 0 ]]; then
            echo "$patient,$sample,$fastq1,$fastq2,$status,$ref" >> "$SAMPLE_SAMPLESHEET"
        else
            echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,$ref" >> "$SAMPLE_SAMPLESHEET"
        fi
    fi

done 
} < "$SAMPLESHEET"


OUT_DIR="$NEXTFLOW_OUT/$SAMPLE"
mkdir -p $OUT_DIR

print_progress "Launching nextflow to perform variant calling..."

cd $NEXTFLOW_WORK

module load java/11.0.2
module load singularity-ce/4.1.1
module load proxies

step="variant_calling"

if [[ $STEP -eq 0 ]]; then
    print_progress "Mapping step currently not supported, sorry. Exiting..."
    exit 0
    step="mapping" #TODO: mapping step requires lane column in samplesheet
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
            --genome $GENOME \
            --max_cpus 48 \

print_progress "Variant calling complete."

# Create a temporary file to store the modified sample sheet
TEMP_SAMPLESHEET=$(mktemp)

echo "patient,sample,fastq_1,fastq_2,status,bam,bai,vcf,ref" >> "$TEMP_SAMPLESHEET"

strelka_dir="${OUT_DIR}/variant_calling/strelka"
mutect_dir="${OUT_DIR}/variant_calling/mutect2"

{
read  # Skip the header line
while IFS= read -r line || [[ -n "$line" ]]; do
    
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Split the line into fields using awk to handle potential edge cases
    IFS=',' read -r patient sample fastq1 fastq2 status bam bai ref<<< "$(awk -F',' '{print $1,$2,$3,$4,$5,$6,$7}' OFS=',' <<< "$line")"

    if [[ $STEP -eq 0 ]] || [[ -z "$bam" ]] || [[ -z "$bai" ]]; then
       
        print_progress "Getting BAM files for ${sample}."

        # Move BAM outputs from nextflow to Raw/Patient directory (downstream analysis scripts assume specific directory structure)

        recal_dir="${OUT_DIR}/preprocessing/recalibrated/${sample}"
        bam_file="${recal_dir}/${sample}.recal.bam"
        bai_file="${recal_dir}/${sample}.recal.bam.bai"

        mv "$bam_file" "${RAW_DIR}/${sample}"
        mv "$bai_file" "${RAW_DIR}/${sample}"

        bam="${RAW_DIR}/${sample}/$(basename "$bam_file")"
        bai="${RAW_DIR}/${sample}/$(basename "$bai_file")"
    fi

    
    if [[ $sample != "Normal" ]] && [[ $status == "1" ]]; then #IMPORTANT: assumes non-tumor sample is named 'Normal'

        print_progress "Getting VCF files for ${sample}."

        strelka_indel="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_indels.vcf.gz"
        strelka_snv="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_snvs.vcf.gz"
        mutect="${mutect_dir}/${sample}_vs_Normal/${sample}_vs_Normal.mutect2.filtered.vcf.gz"

        strelka_indel_tbi="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_indels.vcf.gz.tbi"
        strelka_snv_tbi="${strelka_dir}/${sample}_vs_Normal/${sample}_vs_Normal.strelka.somatic_snvs.vcf.gz.tbi"
        mutect_tbi="${mutect_dir}/${sample}_vs_Normal/${sample}_vs_Normal.mutect2.filtered.vcf.gz.tbi"
        
        mv "$strelka_indel" "${RAW_DIR}/${sample}"
        mv "$strelka_snv" "${RAW_DIR}/${sample}"
        mv "$mutect" "${RAW_DIR}/${sample}"

        mv "$strelka_indel_tbi" "${RAW_DIR}/${sample}"
        mv "$strelka_snv_tbi" "${RAW_DIR}/${sample}"
        mv "$mutect_tbi" "${RAW_DIR}/${sample}"

        stelka_indel_file="${RAW_DIR}/${sample}/$(basename "$strelka_indel")"
        strelka_snv_file="${RAW_DIR}/${sample}/$(basename "$strelka_snv")"
        mutect_file="${RAW_DIR}/${sample}/$(basename "$mutect")"

        vcf_files="$strelka_indel_file|$strelka_snv_file|$mutect_file"

        echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,$vcf_files,$ref" >> "$TEMP_SAMPLESHEET"    

    elif [[ $sample == "Normal" ]]; then
        echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,na,$ref" >> "$TEMP_SAMPLESHEET"    
    fi

done
} < $SAMPLE_SAMPLESHEET


# Replace the original sample sheet with the modified one
mv "$TEMP_SAMPLESHEET" "$SAMPLE_SAMPLESHEET"

# Remove nextflow output dir
#rm -r $OUT_DIR

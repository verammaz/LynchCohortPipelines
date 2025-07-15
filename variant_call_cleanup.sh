#!/bin/bash

source ./config.sh

patient=$1
samplesheet=$2
samples=$3

IFS=',' read -r -a sample_array <<< "$samples"

RAW_DIR="${HOME_DIR}/Raw/${patient}"

> "$samplesheet"

echo "patient,sample,fastq_1,fastq_2,status,bam,bai,vcf,ref" >> "$samplesheet"

wrote_normal=0

for s in "${sample_array[@]}"; do

    sample_samplesheet="${RAW_DIR}/samplesheet_${s}.csv"

    {
    read
    while IFS= read -r line || [[ -n "$line" ]]; do
    
        [[ -z "$line" ]] && continue

        IFS=',' read -r patient sample fastq1 fastq2 status bam bai vcf ref <<< "$(awk -F',' '{print $1,$2,$3,$4,$5,$6,$7,$8}' OFS=',' <<< "$line")"

        if [[ "$sample" == 'Normal' ]] && [[ $wrote_normal -eq 0 ]]; then
            echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,na,$ref" >> "$samplesheet"
            wrote_normal=1
        
        elif [[ "$sample" == "$s" ]]; then
            echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,$vcf,$ref" >> "$samplesheet"
        fi

    done
    } < "$sample_samplesheet"

    rm $sample_samplesheet
done


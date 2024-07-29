#!/bin/bash

source ./config.sh

patient=$1
samplesheet=$2
samples=$3

RAW_DIR="$HOME_DIR/Raw/$PATIENT"

> "$samplesheet"

echo "patient,sample,fastq_1,fastq_2,status,bam,bai,vcf" >> "$samplesheet"

wrote_normal=0

for s in "${samples[@]}"; do

    sample_samplesheet="${RAW_DIR}/samplesheet_${s}.csv"

    {
    read
    while IFS= read -r line || [[ -n "$line" ]]; do
    
        [[ -z "$line" ]] && continue

        IFS=',' read -r patient sample fastq1 fastq2 status bam bai vcf <<< "$(awk -F',' '{print $1,$2,$3,$4,$5,$6,$7,$8}' OFS=',' <<< "$line")"

        if [[ "$sample" == 'Normal' ]] && [[ $wrote_normal -eq 0 ]]; then
            echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,na" >> "$samplesheet"
            wrote_normal=1
        
        else if [[ "$sample" == "$s" ]]; then
            echo "$patient,$sample,$fastq1,$fastq2,$status,$bam,$bai,$vcf" >> "$samplesheet"
        fi

    done 
    } < "$sample_samplesheet"


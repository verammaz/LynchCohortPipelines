#!/bin/bash

bam_file=$1
regions_file=$2
readcounts_file=$3
ref_file=$4

echo "------------------------"
echo "Running bam-readcount"
echo "BAM File: ${bam_file}"

bam-readcount -w1 -i -f ${ref_file} ${bam_file} -l ${regions_file} > ${readcounts_file}

echo "Done."
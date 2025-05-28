#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if a command exits with a non-zero status.
set -e


#############################################################################################
#Run STAR pipeline


echo "---------------------------------------------"
print_progress "Starting the pipeline..."
echo "---------------------------------------------"

# Define required STAR index files
REQUIRED_FILES=(
    "chrLength.txt"
    "chrName.txt"
    "chrNameLength.txt"
    "chrStart.txt"
    "exonGeTrInfo.tab"
    "exonInfo.tab"
    "geneInfo.tab"
    genomeParameters.txt
    "sjdbInfo.txt"
    "sjdbList.fromGTF.out.tab"
    "sjdbList.out.tab"
    "transcriptInfo.tab"
    "Genome"
    "SA"
    "SAindex"
)

# Function to check if all required STAR files exist
check_star_index_complete() {
    local missing=0
    for file in "${REQUIRED_FILES[@]}"; do
        if [ ! -f "${STAR_DIR}/${file}" ]; then
            echo "[STAR index check] Missing file: ${file}"
            missing=1
        fi
    done
    return $missing
}

# Decide whether to run indexing
RUN_INDEX=0
if ! check_star_index_complete; then
    echo "Indexing required: STAR index is incomplete."
    RUN_INDEX=1
fi

# Actually run STAR indexing if flagged
if [ $RUN_INDEX -eq 1 ]; then
    create_directory "$STAR_DIR"
    
    print_progress "Indexing the reference..."

    STAR --runMode genomeGenerate \
         --runThreadN 8 \
         --genomeDir "$STAR_DIR" \
         --genomeFastaFiles "$REF_FASTA" \
         --sjdbGTFfile "$GTF_FILE" \
         --sjdbOverhang 100 \
         --outFileNamePrefix "${STAR_DIR}/star"
fi

# Run STAR alignment

READ1="$1"
READ2="$2"
OUTPUT="$3"

print_progress "Running STAR alignment..."

STAR --runThreadN 8 \
     --genomeDir "$STAR_DIR" \
     --readFilesIn "$READ1" "$READ2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${OUTPUT}_" \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts

print_progress "Pipeline completed successfully. Output is in  Output is in ${OUTPUT}_*.bam and logs."

exit 0
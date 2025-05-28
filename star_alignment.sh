#!/bin/bash

# Get access to global variables
source ./config.sh

# Exit immediately if a command exits with a non-zero status.
#set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [options] -f <reference.fa> -r <read1.fastq,read2.fastq> -o <output_prefix>

Required Arguments:
  -r <read1.fastq,read2.fastq>     Comma-separated list of two FASTQ files.
  -o <output_prefix>               Prefix for output files.

Options:
  -h                               Display this message.
  -v                               Enable verbose mode.
  --threads                        Number of threads to use.

EOF
    exit 1
}

# Variables
READS=()
OUTPUT_PREFIX=
THREADS=
VERBOSE=


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        --threads) THREADS="$2"; shift ;;
        -r) READS="$2"; shift ;;
        -o) OUTPUT_PREFIX="$2"; shift ;;
        *) echo "Error: Unknown argument/option: $1" ; usage ;;
    esac
    shift
done

# Split the READS argument into an array
IFS=',' read -r -a READ_ARRAY <<< "$READS"

# Ensure exactly two read files are provided
if [ ${#READ_ARRAY[@]} -ne 2 ]; then
    echo "Error: Exactly two read files must be specified, separated by a comma."
    usage
fi

READS_1=${READ_ARRAY[0]}
READS_2=${READ_ARRAY[1]}

# Check fastq read files
#if [ ! -f "$READS_1" ] || [ ! -f "$READS_2" ]; then
    #echo "Error: One or both fastq read files do not exist or are not readable."
    #exit 1
#fi

# Check reference genome file
if [ ! -f "${REF_FASTA}" ]; then
    echo "Reference file ${REF_FASTA} not found!"
    exit 1
fi

# Check gtf annotation file
if [ ! -f "${GTF_FILE}" ]; then
    echo "GTF annotation file ${GTF_FILE} not found!"
    exit 1
fi

# Check STAR directory name not empty
if [ -z "$STAR_DIR" ]; then
    echo "STAR_DIR is not set. Please define the STAR index output directory in config.sh file."
    exit 1
fi

check_gtf_fasta_compatibility() {
    local gtf_chr=$(zgrep -v '^#' "$GTF_FILE" 2>/dev/null | awk '$3 == "exon" {print $1}' | sort | uniq | head -n 1)
    local fasta_chr=$(zgrep '^>' "$REF_FASTA" 2>/dev/null | head -n 1 | sed 's/^>//' | awk '{print $1}')

    echo "[compatibility check] First GTF chromosome: $gtf_chr"
    echo "[compatibility check] First FASTA chromosome: $fasta_chr"

    # Normalize chrM/MT
    gtf_chr="${gtf_chr/chrM/chrMT}"
    fasta_chr="${fasta_chr/MT/chrMT}"

    if [[ "$gtf_chr" == "$fasta_chr" ]]; then
        echo "GTF and FASTA chromosome naming appears compatible."
        return 0
    elif [[ "$gtf_chr" =~ ^chr && ! "$fasta_chr" =~ ^chr ]]; then
        echo "GTF uses 'chr' prefix but FASTA does not. Consider stripping 'chr' from GTF."
        return 1
    elif [[ ! "$gtf_chr" =~ ^chr && "$fasta_chr" =~ ^chr ]]; then
        echo "FASTA uses 'chr' prefix but GTF does not. Consider adding 'chr' to GTF."
        return 1
    else
        echo "Chromosome format mismatch between GTF and FASTA. Manual inspection recommended."
        return 1
    fi
}

check_gtf_fasta_compatibility
if [ $? -ne 0 ]; then
    echo "ERROR: GTF and FASTA chromosome formats do not match. Exiting."
    exit 1
fi




# If verbose mode is enabled, print the parameters
if [ $VERBOSE -eq 1 ]; then
    
    echo "---------------------------------------------"
    echo "Reference: $REF_FASTA"
    echo "Read1: $READS_1"
    echo "Read2: $READS_2"
    echo "Output Prefix: $OUTPUT_PREFIX"
    echo "STAR Genome Directory: $STAR_DIR"
    echo "---------------------------------------------"
fi


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
check_star_index_complete
    if [ $? -ne 0 ]; then
        echo "STAR genome index is incomplete. Running STAR indexing."
        RUN_INDEX=1
    else
        echo "STAR index exists and is complete."
    fi


# Actually run STAR indexing if flagged
if [ $RUN_INDEX -eq 1 ]; then
    create_directory "$STAR_DIR"
    
    print_progress "Indexing the reference..."

    STAR --runMode genomeGenerate \
         --runThreadN "$THREADS" \
         --genomeDir "$STAR_DIR" \
         --genomeFastaFiles "$REF_FASTA" \
         --sjdbGTFfile "$GTF_FILE" \
         --sjdbOverhang 100 \
         --outFileNamePrefix "${STAR_DIR}/star"
fi

# Run STAR alignment

exit 0
print_progress "Running STAR alignment..."

STAR --runThreadN "$THREADS" \
     --genomeDir "$STAR_DIR" \
     --readFilesIn "$READ1" "$READ2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${OUTPUT_PREFIX}_" \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts

print_progress "Pipeline completed successfully. Output is in ${OUTPUT+PREFIX}."

exit 0
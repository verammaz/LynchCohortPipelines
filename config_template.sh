# config.sh

###################################### General ##################################
# specify project
project=''

# specify cores
cores=24

# specify path to the fastq_to_bam singularity image file
CONTAINER_FASTQ2BAM=''

# specify path to OpiType config.ini file
OPTITYE_CONFIG=''

# specify path to conda pairtree environment
PAIRTREE_ENV='' 

#################################################################################


####################################### Reference hg19 ###############################
# reference genome
GENOME_hg19="GATK.GRCh37"

# specify file path for reference fasta file
REF_FASTA_hg19="human_g1k_v37_decoy.fasta"

# specify file path for sites_of_variation.vcf 
SITES_OF_VARIATION_hg19="dbsnp_138.b37.vcf"

# specify file paths for reference files for indel realignment
INDELS_1_hg19="1000G_phase1.indels.b37.vcf"
INDELS_2_hg19="Mills_and_1000G_gold_standard.indels.b37.vcf"

# specify path to exome target intervals file (this one not required!)
EXOME_INTERVALS_hg19="Broad.human.exome.b37.interval_list"

# pon file for mutect variant filtering
PON_hg19="Mutect2-exome-panel.vcf.gz"

# specify file path for COMPATIBLE gtf annotations (required for STAR indexing)
GTF_hg19=""

##################################################################################

####################################### Reference hg38 ###############################
# reference genome
GENOME_hg38="GATK.GRCh38"

# specify file path for reference fasta file
REF_FASTA_hg38=""

# specify file path for sites_of_variation.vcf 
SITES_OF_VARIATION_hg38=""

# specify file paths for reference files for indel realignment
INDELS_1_hg38=""
INDELS_2_hg38=""

# specify path to exome target intervals file (this one not required!)
EXOME_INTERVALS_hg38=""

# pon file for mutect variant filtering
PON_hg38=""

# specify file path for COMPATIBLE gtf annotations (required for STAR indexing)
GTF_hg38=""

##################################################################################

############################### Scripts / Executables ############################
# nextflow executable
NEXTFLOW=''

# NeoPipe directory
NEOPIPE=''

# CFIT directory
CFIT=''

# Pairtree directory
PAIRTREE=''

#specify path to LynchCohortPipelines folder
LYNCH=''

#################################################################################

######### Directories ( make sure these have enough disk space ! ) ###############

# home directory for data and output files
HOME_DIR=''

# nextflow output directory
NEXTFLOW_OUT=''

# nextflow work directory
NEXTFLOW_WORK=''

# specify temporary directory (this is for fastq->bam intermediate files during sorting)
TEMP_DIR=""

# specify directory for log files 
LOG_DIR=''

# specify directory for STAR indexing files
STAR_DIR=""

##################################################################################


#################################### Functions ###################################

# Function to print progress with timestamp
print_progress() {
    echo "[`date +%Y-%m-%dT%H:%M:%S`] $1"
}

# Function to create directory and necessary parent directories
function create_directory() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Directory created: $dir"
    else
        echo "Directory already exists: $dir"
    fi
}


#####################################################################################

# Ensure all global directories exist 
#create_directory "$NEXTFLOW_OUT"
#create_directory "$NEXTFLOW_WORK"
#create_directory "$TEMP_DIR"
#create_directory "$LOG_DIR"
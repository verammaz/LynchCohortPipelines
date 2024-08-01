# config.sh

###################################### General ##################################
# specify project
project='acc_FLAI'

# specify cores
cores=24

#################################################################################


####################################### Reference ###############################
##### Ensure that reference files are compatible! 
#################################################################################

# specify reference genome (GATK.GRCh38, GATK.GRCh37, Ensembl.GRCh37, NCBI.GRCh38)
GENOME="GATK.GRCh37"

# specify file path for reference fasta file
REF_FASTA="/path/to/human_g1k_v37_decoy.fasta"

# specify file path for sites_of_variation.vcf 
SITES_OF_VARIATION="/path/to/dbsnp_138.b37.vcf"

# specify file paths for reference files for indel realignment
INDELS_1="/path/to/1000G_phase1.indels.b37.vcf"
INDELS_2="/path/to/Mills_and_1000G_gold_standard.indels.b37.vcf"

# specify path to exome target intervals file (this one not required!)
EXOME_INTERVALS="/path/to/Broad.human.exome.b37.interval_list"

# pon file for mutect variant filtering
PON='/sc/arion/projects/FLAI/vera/References/Mutect2-exome-panel.vcf.gz'

##################################################################################


############################### Scripts / Executables ############################
# nextflow executable
NEXTFLOW='/hpc/users/mazeev01/nextflow'

# NeoPipe directory
NEOPIPE='/hpc/users/mazeev01/NeoPipe'

# CFIT directory
CFIT='/hpc/users/mazeev01/CFIT'

# Pairtree directory
PAIRTREE='/hpc/users/mazeev01/pairtree'

#specify path to LynchCohortPipelines folder
LYNCH='/hpc/users/mazeev01/LynchCohortPipelines'

# specify path to the fastq_to_bam singularity image file
CONTAINER_FASTQ2BAM='~/fastq2bam_0.3.sif'

# specify path to snpeff
SNPEFF='~/snpEff.v4.3t/'

# specify path to netmhc
NETMHC='~/netMHCpan-4.1/'

# specify path to OpiType config.ini file
OPTITYE_CONFIG='/sc/arion/projects/MSIH-seq/data/MattBrown/Whole_Exome_Sequencing/MattBrownWESAug2023/optitypeResults/config.ini'

#################################################################################


######### Directories ( make sure these have enough disk space ! ) ###############

# home directory for data and output files
HOME_DIR='sc/arion/projects/FLAI/vera/matt_lynch'

# nextflow output directory
NEXTFLOW_OUT='/sc/arion/scratch/mazeev01/nextflow_out'

# nextflow work directory
NEXTFLOW_WORK='/sc/arion/scratch/mazeev01'

# specify temporary directory (this is for fastq->bam intermediate files during sorting)
TEMP_DIR="/sc/arion/scratch/mazeev01/temp"

# specify directory for log files 
LOG_DIR='/sc/arion/work/mazeev01/logs'

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
    fi
}


#####################################################################################

# Ensure all global directories exist 
create_directory "$NEXTFLOW_OUT"
create_directory "$NEXTFLOW_WORK"
create_directory "$TEMP_DIR"
create_directory "$LOG_DIR"
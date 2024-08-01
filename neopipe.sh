#!/bin/bash

source ./config.sh

export PATH=$PATH:$SNPEFF
export PATH=$PATH:$NETMHC
export R_HOME=/opt/hpc/packages/minerva-centos7/modulefiles/R/4.0.2
module load R/4.0.2
module load python/3.7.3


# Exit immediately if command exits with non-zero status
set -e

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -p <patient_id> -sample_info <sample_info.txt>

Required Arguments:
    -p <patient_id>                         Patient id
Options:
    -h                                      Display this message
    -v                                      Enable verbode mode
    -ns <i-j>                               Range of neoantigen peptide lengths
    -kd <i>                                 Kd threshold for neoantigen computation.
    -genome <genome>                        Reference genome
    -hla <hla.txt>                          File with patient hla calls
    --single_sample_trees       

EOF
    exit 1
}

VERBOSE=0
PATIENT=
NS=9
HLA=$HOME_DIR/HLA/HLA_calls.txt
REF=$GENOME
KD=50000
SAMPLE_INFO=
SINGLE_TREES=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in 
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -ns) NS="$2"; shift ;;
        -genome) REF="$2"; shift ;;
        -hla) HLA="$2"; shift ;;
        -kd) KD="$2"; shift ;;
        -sample_info) SAMPLE_INFO="$2"; shift ;;
        --single_sample_trees) SINGLE_TREES=1 ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done



### Run snpeff ####################

genome=
if [[ "$REF" == *"37"* ]]; then
    genome="GRCh37.75"
elif [[ "$REF" == *"38"* ]]; then
    genome="GRCh38.86"
else 
    echo "Error: Unkown genome $GENOME. Available options are grch37 and grch38."
    exit 1
fi


python $NEOPIPE/bin/run_SnpEff_per_patient.py -patient $PATIENT \
                                               -genome $genome \
                                               -dir $HOME_DIR \
                                               -ns $NS -force

### Run netMHC ####################


if [[ ! -f $HLA ]]; then
    echo "Error: Cannont find file with HLA calls."
    exit 1
fi

python $NEOPIPE/bin/compute_neoantigens_per_patient.py -patient $PATIENT \
                                                       -hla $HLA \
                                                       -dir $HOME_DIR \
                                                       -kd_thr 50000 \
                                                       -ns $NS -netmhc pan4.1


### Prep CFIT config and mapping files

python $LYNCH/prep_cfit_configs.py -patient_id $PATIENT \
                                   -sample_info $SAMPLE_INFO \
                                   -hdir $HOME_DIR \
                                   -single_sample_trees $SINGLE_TREES \
                                   -hla $HLA
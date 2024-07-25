#!/bin/bash

###########################################################
#specify project
project='acc_FLAI'
#specify cores
cores=24
#specify path to auxilary script
script='/hpc/users/mazeev01/matt_lynch/run_nextflow.sh'
############################################################


# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] --patient <patient_id> --samplesheet <samplesheet.csv>

Required Arguments:
    --patient <patient_id>           Patient id
    --samplesheet <samplesheet.csv>  File with all sample raw data configured

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --variant_call_step             Step of variant call (0 = preprocessing, 1 = variant calling)
}

EOF
    exit 1
}



PATIENT=
SAMPLESHEET=
VERBOSE=0
STEP=1

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        --patient) PATIENT="$2"; shift ;;
        --samplesheet) SAMPLESHEET="$2"; shift ;;
        --variant_call_step) STEP="$2"; shift ;;
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

chmod +x $script

while IFS=$',' read -r patient sample fastq1 fastq2 status bam bai; do
    
    if [[ $patient == $PATIENT ]] && [ $sample != "Normal" ]; then

        job_name="variant_call_$sample"

        bsub -J ${job_name} \
             -P ${project} \
             -n ${cores} \
             -M 32000 \
             -R span[hosts=1] \
             -R "rusage[mem=4000]" \
             -W 30:00 \
             -q premium \
             -oo ${sample}_variantcall.out \
             -eo ${sample}_variantcall.err \
             ${script} --patient ${PATIENT} --samplesheet ${SAMPLESHEET} --sample ${sample} --step ${STEP}
    
    fi

done < "$SAMPLESHEET"
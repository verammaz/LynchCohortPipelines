#!/bin/bash

source ./config.sh

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS] -s <samplesheet.csv> -p <patient_id>

Required Arguments:
    -s <samplesheet.csv>            Configured input file
    -p <patient_id>                 Patient identifier

Options:
    -h                              Display this message
    -v                              Enable verbode mode
    --normal_only                   Run hla typing only on the noraml fastq
}

EOF
    exit 1
}


PATIENT=
SAMPLESHEET=
VERBOSE=0
NORMAL_ONLY=0

# parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h) usage ;;
        -v) VERBOSE=1 ;;
        -p) PATIENT="$2"; shift ;;
        -s) SAMPLESHEET="$2"; shift ;;
        --normal_only) NORMAL_ONLY=1 ;; 
        *) echo "Error: Unkown argument/option: $1" ; usage ;;
    esac
    shift
done

script=$LYNCH/opitype_run.sh

chmod +x $script

{
read  # Skip the header line
while IFS= read -r line || [[ -n "$line" ]]; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Split the line into fields using awk to handle potential edge cases
    IFS=',' read -r patient sample fastq1 fastq2 status <<< "$(awk -F',' '{print $1,$2,$3,$4,$5}' OFS=',' <<< "$line")"
    
    echo $sample
    if [[ $patient == $PATIENT ]]; then

        if [[ "$sample" == 'Normal' ]] || [[ $NORMAL_ONLY -eq 0 ]]; then
            
            job_name="optitype_${patient}_${sample}"

            bsub -J ${job_name} \
                    -P ${project} \
                    -n ${cores} \
                    -M 32000 \
                    -R span[hosts=1] \
                    -R "rusage[mem=4000]" \
                    -W 40:00 \
                    -q premium \
                    -oo "${LOG_DIR}/${job_name}.out" \
                    -eo "${LOG_DIR}/${job_name}.err" \
                    ${script} -r1 ${fastq1} -r2 ${fastq2} $patient "${patient}_${sample}"
       
       fi
    fi

done
} < "$SAMPLESHEET"
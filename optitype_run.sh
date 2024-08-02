#!/bin/bash

source ./config.sh

module load optitype

r1=$1
r2=$2 
patient=$3
prefix=$4

### Run optitype

OptiTypePipeline.py -i $r1 $2 --dna -v --prefix $prefix -o "${HOME_DIR}/Raw/$patient" -c $OPTITYE_CONFIG


### Cleanup and move output to $HOME_DIR/HLA/HLA_calls.txt file

if [[ $prefix == "${patient}_Normal" ]]; then 
    
    optitype_file="${HOME_DIR}/Raw/$patient/${prefix}_result.tsv"

    if [[ ! -d "${HOME_DIR}/HLA" ]]; then
        mkdir "${HOME_DIR}/HLA"
    fi

    hla_calls_file="${HOME_DIR}/HLA/HLA_calls.txt"

    alleles=$(awk 'NR==2 {print $2","$3","$4","$5","$6","$7}' "$optitype_file")

    new_line="$patient\t$alleles"

    if [ ! -e "$hla_calls_file" ]; then
    touch "$hla_calls_file"
    fi

    if [ -s "$hla_calls_file" ]; then
    if [ "$(tail -c 1 "$hla_calls_file")" != "" ]; then
        echo "" >> "$hla_calls_file"
    fi
    fi

    echo -e "$new_line" >> "$hla_calls_file"
fi
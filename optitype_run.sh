#!/bin/bash

source ./config.sh

module load optitype

r1=$1
r2=$2 
prefix=$3
out=$4

### Run optitype

OptiTypePipeline.py -i $r1 $2 --dna -v --prefix $prefix -o $out -c $OPTITYE_CONFIG


### Cleanup and move output to $HOME_DIR/HLA/HLA_calls.txt file

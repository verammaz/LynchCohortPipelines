#!/bin/bash
#BSUB -J mySTARjob                   # Job name
#BSUB -P acc_FLAI                    # allocation account
#BSUB -q premium                        # queue
#BSUB -n 8                                    # number of compute cores
#BSUB -W 12:00                            # walltime in HH:MM
#BSUB -R rusage[mem=4000]      # 32 GB of memory (4 GB per core)
#BSUB -R span[hosts=1]            # all cores from the same node
#BSUB -o %J.stdout                     # output log (%J : JobID)
#BSUB -eo %J.stderr                   # error log
#BSUB -L /bin/bash                      # Initialize the execution environment

 

source ./config.sh

module load star                           # load star module
GENOME_DIR='/sc/arion/scratch/mazeev01/STAR'
create_directory "$GENOME_DIR"

STAR --runMode genomeGenerate \
             --runThreadN ${THREADS} \
             --genomeDir ${GENOME_DIR} \
             --genomeFastaFiles ${REF_FASTA} \
             --sjdbGTFfile ${GTF_FILE} \
             --sjdbOverhang 100  # Adjust depending on read length - 1
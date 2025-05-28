#!/bin/bash
#BSUB -J mySTARjob                   # Job name
#BSUB -P acc_FLAI                    # allocation account
#BSUB -q premium                        # queue
#BSUB -n 8                                    # number of compute cores
#BSUB -W 12:00                            # walltime in HH:MM
#BSUB -R rusage[mem=40000]      # 32 GB of memory (4 GB per core)
#BSUB -R span[hosts=1]            # all cores from the same node
#BSUB -o %J.stdout                     # output log (%J : JobID)
#BSUB -eo %J.stderr                   # error log
#BSUB -L /bin/bash                      # Initialize the execution environment

 

source ./config.sh

module load star                           # load star module
GENOME_DIR='/sc/arion/scratch/mazeev01/STAR'

READ1=
READ2=

STAR --runThreadN 8 \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn "$READ1" "$READ2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "/sc/arion/scratch/mazeev01/star_out_" \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts
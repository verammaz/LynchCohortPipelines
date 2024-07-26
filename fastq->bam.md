# FASTQ to BAM
> This pipeline calls a sequence mapper / aligner and then carries out optional post processing.

## Pipeline Overview:

0. [BWA Index](https://bio-bwa.sourceforge.net/bwa.shtml)
1. [BWA MEM](https://bio-bwa.sourceforge.net/bwa.shtml) and conversion to BAM format with [samtools](https://www.htslib.org/)
2. [Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)

Optional post-processing:

3. [GATK IndelRealigner](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md)
    - Note: vcf files with known indel sites required for this step.
    - *Note*: for WES data, can speed up this step by using a exome_targets.interval_list file.
4. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - *Note:* vcf file with sites of variation required for this step.
    - *Note*: for WES data, can speed up this step by using a exome_targets.interval_list file.


## Running on Minerva:

### Singularity Container 

Modules:
- bwa v0.7.15
- samtools v1.17
- java v1.8
- picard v2.2.4
- gatk v3.8

```bash
module load singularity/3.6.4
module load proxies 
cd ~/
singularity pull --arch amd64 library://verammaz/bioinformatics/fastq2bam:0.3
```
Set the path to this image file in the `CONTAINER_FASTQ2BAM` variable in the [config](config.sh) file.


### Single sample pipeline

#### Usage

All NGS aligners need the reference sequences to be indexed. On the very first use of the pipeline with a reference genome or if you don't have `genome.fasta.{amb, ann, btw, fai, pac, sa}` files, run

```bash
module singularity/3.6.4
cd LynchCohortPipelines
source ./config.sh
singularity exec $CONTAINER_FASTQ2BAM fastq_to_bam.sh -r <fastq1,fastq2> -o <output_prefix> --index_ref
```
On subsequent execution using the same reference, run without the `--index_ref` option. 

> *Important*: Keep reference fasta and reference index files in the same directory, and make sure the file prefix names are consistent!

> Need at least the following set in the [config](config.sh): `REF_FASTA`, `SITES_OF_VARIATION` (optional), `INDELS_{1,2}` (optional), `EXOME_INTERVALS` (optional) 

#### Required arguments:
```
-r          Two read fastq files, separated by a comma (no space)
-o          Output files prefix
```

> *Note:* current version assumes all reads in fastq files come from same sample and are run on same Illumina sequencing lane, and uses this to set the @RG tag in the BAM files. Read about the @RG tag [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups). 

#### Optional arguments:

| Parameter                 | Description   |	
| :----------------------------------------: | :------: |
| `-v` | Enable verbose mode. |
| `-h` | Display usage message. |
| `--patient` |  Patient id for sample data, required for downstream compatability during nextflow variant calling. Default none.
| `--index_ref` | Flag for running the indexing step for reference.
| `--keep_intermediate` | Flag to keep intermediate files produced during pipeline execution.
| `--post_process` | Flag to carry out post processing of initial alignment bam file (indel local realignment, base recalibration).
| `--threads` | Number of threads to use. Default 8. |
| `--step` | Starting step for pipeline execution. Default 0. Options are 0=alignment, 1=markdup, 2=indelrealign, 3=baserecal. |

#### Submitting a Job (bsub)

Submit to LSF job scheduler with the following header:

```bash
#!/bin/bash --login

#BSUB -J fastq2bam
#BSUB -P acc_ProjectName
#BSUB -q queue_name 
#BSUB -n 8
#BSUB -M 32000 
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=4000]
#BSUB -W 30:00
#BSUB -o fastq2bam_%J.out
#BSUB -e fastq2bam_%J.err
```
*Note:* Can reduce computational time by speeding up alignment. To do this, increase number of cores and call the script with `--threads n`, where `n` matches the number of cores requested for the job. All cores must be on the same compute node. 

### Patient Bulk Pipeline

There is a script to automatically submit fastq->bam jobs, executing the `fastq_to_bam.sh` script, in bulk *per patient.* This is useful when a single patient has more than one sample and you don't want to retype the bsub options and command. The sample processing will be run in parallel. 

#### Usage 

Before running, make sure you have all reference index files (i.e. `genome.fasta.{amb,ann,bwt,fai,pac,sai}`). Keep reference fasta and reference index files in the same directory, and make sure the file prefix names are consistent! If you don't have reference index files, run the following:

```bash
module load bwa/0.7.15
cd LynchCohortPipelines
source ./config.sh
bwa index -a bwtsw $REF_FASTA
```

Now, you can run:

```bash
submit_fastq2bam_for_patient.sh -p <patient_id>  -s <samplesheet.csv>
```

Make sure HOME_DIR is set in the [config](config.sh) file. The script will create `$HOME_DIR/Raw/Patient` directory if it doesn't exist, and place all output .bam and .bai files in 
Sample/ and Normal/ subdirectories.

#### Required arguments:
```
-p         Patient identifier.
-s         CSV file with raw input data files configuration.
```
#### Input .csv file
`samplesheet.csv` needs to have the columns patient, sample, fastq_1, fastq_2, status (0=Normal, 1=Tumor). This file will change to include bam and bai columns before the fastq->bam jobs per sample are submitted. Note, that downstream analysis scripts assume normal sample is named 'Normal', so ensure this is the case for your data. Example:

```csv
patient,sample,fastq_1,fastq_2,status
Patient1,Normal,full/path/to/Normal_R1_001.fastq.gz,full/path/to/Normal_R2_001.fastq.gz,0
Patient1,S1,full/path/to/S1_R1_001.fastq.gz,full/path/to/S1_R2_001.fastq.gz,1
```
> *Important*: Current version assumes all samples are run on single lane!



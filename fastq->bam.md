# fastq_to_bam.sh
This script runs a FASTQ --> BAM pipeline by calling a sequence mapper / aligner and then carrying out optional post processing.

## Pipeline Overview:

0. [BWA Index](https://bio-bwa.sourceforge.net/bwa.shtml)
1. [BWA MEM](https://bio-bwa.sourceforge.net/bwa.shtml) and conversion to BAM format with [samtools](https://www.htslib.org/)
2. [Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)

Optional post-processing:

3. [GATK IndelRealigner](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md)
    - Note: vcf files with known indel sites required for this step. Specify path to these files in top of script in the `INDELS_{1,2}` variables.
    - *Note*: for WES data, can speed up this step by using a exome_targets.interval_list file (can specify path to file in top of script in `EXOME_INTERVALS` variable).
4. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - *Note:* vcf file with sites of variation required for this step. Specify path to this file in top of script in the `SITES_OF_VARIATION` variable.
    - *Note*: for WES data, can speed up this step by using a exome_targets.interval_list file.


## Running on Minerva:

First, specify file paths in the top of the script! Also, specify a temporary directory with sufficient storage.

All NGS aligners need the reference sequences to be indexed. On the very first use of the pipeline with a reference genome or if you don't have ` genome.fasta.{amb, ann, btw, fai, pac, sa}` files, run

```bash
path/to/file/fastq_to_bam.sh -r <fastq1,fastq2> -o <output_prefix> --index_ref
```
On  subsequent execution using the same reference, run without the `--index_ref` option. *IMPORTANT:* Keep reference fasta and reference index files in the same directory, and make sure the file prefix names are consistent!

### Usage 

Required modules:
- bwa
- samtools
- java v1.8
- picard v2.2.4
- gatk v3.6

*Note:* required modules are already available on Minerva. This script loads them in, so no need to load any software modules before running:)

Required arguments:
```
-r          Two read fastq files, separated by a comma (no space)
-o          Output files prefix
```

*Note:* current version assumes all reads in fastq files come from same sample and are run on same Illumina sequencing lane, and uses this to set the @RG tag in the BAM files. Read about the @RG tag [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups). 

Optional arguments:

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

### Submitting a Job (bsub)

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

### Submitting jobs in bulk per patient

There is a script to automatically submit fastq->bam jobs, executing the `fastq_to_bam.sh` script discussed earlier, in bulk *per patient.* This is useful when a single patient has more than one sample and you don't want to retype the bsub options and command. The sample processing will be run in parallel. 

First, change the project name, and script path in the very top of [submit_fastq2bam_for_patient.sh](submit_fastq2bam_for_patient.sh). To run, use the following command:

```bash
submit_fast2bam_for_patient.sh -p <patient_id>  -s <samplesheet.csv>
```

#### Usage 

Required arguments:
```
-p         Patient identifier.
-s         CSV file with raw input data files configuration.
```

`samplesheet.csv` needs to have the columns patient, sample, fastq1, fastq2, status (0=Normal, 1=Tumor). This file will change to include bam and bai columns before the fastq->bam jobs per sample are submitted. Note, that downstream analysis scripts assume normal sample is named 'Normal', so ensure this is the case for your data. Example:

```csv
patient,sample,fastq1,fastq2,status
Patient1,Normal,full/path/to/Normal_R1_001.fastq.gz,full/path/to/Normal_R2_001.fastq.gz,0
Patient1,S1,full/path/to/S1_R1_001.fastq.gz,full/path/to/S1_R2_001.fastq.gz,1
```

Optional arguments:

| Parameter                 | Description   |	
| :----------------------------------------: | :------: |
| `-v` | Enable verbose mode. |
| `-h` | Display usage message. |
| `--data_dir` |  Directory with Sample/ and Normal/ subdirectories that will have the .bam and .bai output files. By default, the script will take the parent directory of the samplesheet.csv file. 



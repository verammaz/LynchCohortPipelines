# fastq_to_bam.sh
This script runs a FASTQ --> BAM pipeline by calling a sequence mapper / aligner and then carrying out optional post processing.

## Pipeline Overview:

[Burrows-Wheeler Alignment tool](https://bio-bwa.sourceforge.net/bwa.shtml)

0. Indexing the reference (bwa_index)
1. Alignment (bwa_mem) and sorting (samtools sort)

Optional post-processing:

2. [Picard MarkDuplictaes](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)
3. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - *Note:* vcf file with sites of variation required for this step. Specify path to this file in top of script in the 'SITES_OF_VARIATION' variable.
    - Pre step: [GATK AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard)
4. [Picard ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR)

## Running on Minerva (HPC cluster):
*Note:* required modules (bwa, samtools, picard, java, gatk) are already available on Minerva. This script loads them in, so no need to load any software modules before running:)

All NGS aligners need the reference sequences to be indexed. On the very first use of the pipeline with a reference genome, run

```bash
path/to/file/fastq_to_bam.sh -f <reference.fa> -r <fastq1,fastq2> -o <output_prefix> --run_index
```
On  subsequent execution using the same reference, run without the `--run_index` option. Keep reference and index files in the same directory.

### Usage

Required arguments:
```
-f          Reference genome fasta file
-r          Two read fastq files, separated by a comma (no space)
-o          Output files prefix
```

Optional arguments:

| Parameter                 | Description   |	
| :----------------------------------------: | :------: |
| `--view_final` |  Flag to create .fai files for reference genome and final aligned bam file, which is required for viewing results in IGV.
| `--run_index` | Flag for running the indexing step.
| `--keep_intermediate` | Flag to keep intermediate files produced during pipeline execution.
| `--post_process` | Flag to carry out post processing of aligned and sorted bam file (mark duplicates, base recalibration).

### Submitting a job (bsub)

Submit to LSF job scheduler with the following header:

```bash
#!/bin/bash --login

#BSUB -J fastq_to_bam
#BSUB -P acc_Project Name
#BSUB -q queue_name
#BSUB -n 8
#BSUB -M 163840 
#BSUB -R "rusage[mem=24576]"
#BSUB -W 10:00
#BSUB -oo fastq_to_bam.out
#BSUB -eo fastq_to_bam.err
```




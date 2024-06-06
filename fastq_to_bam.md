# fastq_to_bam.sh
This script runs a FASTQ --> BAM pipeline by calling a sequence mapper / aligner and then carrying out optional post processing.

## Pipeline Overview:

[Burrows-Wheeler Alignment tool](https://bio-bwa.sourceforge.net/bwa.shtml)

0. Indexing the reference (bwa_index)
1. Alignment (bwa_mem) and conversion to BAM format (samtools)

Optional post-processing:

2. [GATK MarkDuplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/360037224932-MarkDuplicatesSpark)
3. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - *Note:* vcf file with sites of variation required for this step. Specify path to this file in top of script in the 'SITES_OF_VARIATION' variable.
4. [GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR)

## Running on Minerva (HPC cluster):
*Note:* required modules (bwa, samtools, gatk) are already available on Minerva. This script loads them in, so no need to load any software modules before running:)

All NGS aligners need the reference sequences to be indexed. On the very first use of the pipeline with a reference genome, run

```bash
path/to/file/fastq_to_bam.sh -f <reference.fa> -r <fastq1,fastq2> -o <output_prefix> --index_ref
```
On  subsequent execution using the same reference, run without the `--index_ref` option. Keep reference and index files in the same directory.

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
| `-v` | Enable verbose mode. |
| `-h` | Display usage message. |
| `--igv` |  Flag to create .bai files for BAM files, which is required for viewing alignments in IGV.
| `--index_ref` | Flag for running the indexing step for reference.
| `--keep_intermediate` | Flag to keep intermediate files produced during pipeline execution.
| `--post_process` | Flag to carry out post processing of initial alignment bam file (mark duplicates, base recalibration).

### Submitting a Job (bsub)

Submit to LSF job scheduler with the following header:

```bash
#!/bin/bash --login

#BSUB -J fastq_to_bam
#BSUB -P acc_ProjectName
#BSUB -q express (or premium) 
#BSUB -n 48
#BSUB -M 32000 
#BSUB -R span[hosts=1]
#BSUB -W 03:00
#BSUB -oo fastq_to_bam.out
#BSUB -eo fastq_to_bam.err
```




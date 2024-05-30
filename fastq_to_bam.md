# fastq_to_bam.sh
This script runs a FASTQ --> BAM pipeline by calling a sequence mapper and then carrying out alignment and optional post processing.

## Pipeline Overview:

[Burrows-Wheeler Alignment tool](https://bio-bwa.sourceforge.net/bwa.shtml)

0. Indexing the reference (bwa_index)
1. Alignment (bwa_mem) and sorting (samtools sort)

Optional post-processing:

2. [Picard MarkDuplictaes](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)
3. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - *Note:* vcf file with sites of variation required for this step. Specify path to this file in top of script in the 'SITES_OF_VARIATION' variable.
4. [GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR)

## Running on Minerva (HPC cluster):
*Note:* required modules (bwa, samtools, picard, java, gatk) are already available on Minerva 

All NGS aligners need the reference sequences to be indexed. On the very first use of the pipeline with a reference genome, run

```bash
path/to/file/fastq_to_bam.sh -f <reference.fa> -r <fastq1,fastq2> -o <output_prefix> --index
```
On  subsequent execution using the same reference, run without the `--index` option. Keep reference and index files in the same directory.

### Usage

Required arguments:
```
-f          Reference genome fasta file
-r          Two read fastq files, separated by a comma (no space)
-o          Output files prefix
```

Optional arguments:
| Parameter                 | Description   |	
| :------------------------ | :-------------|
| --view_final |  Flag to create .fai files for reference genome and final aligned bam file, which is required for viewing results in IGV.
| --run_index | Flag for running the indexing step.
--keep_intermediate | Flag to keep intermediate files produced during pipeline execution.
--post_process | Flag to carry out post processing of aligned and sorted bam file (mark duplicates, base recalibration).




# fastq_to_bam.sh
This script runs a FASTQ --> BAM pipeline by calling a sequence mapper / aligner and then carrying out optional post processing.

## Pipeline Overview:

[Burrows-Wheeler Alignment tool](https://bio-bwa.sourceforge.net/bwa.shtml)

0. Indexing the reference (bwa_index)
1. Alignment (bwa_mem) and conversion to BAM format (samtools)

Optional post-processing:

2. [Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)
3. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - *Note:* vcf file with sites of variation required for this step. Specify path to this file in top of script in the 'SITES_OF_VARIATION' variable.
4. [GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR)

## Reference File(s):

Recommended directory structure for reference files:
```bash
Reference/
    ├── genome.fasta
    ├── genome.fasta.{amb, ann, btw, fai, pac, sa}  #Output of bwa-index 
    ├── sites_of_variation.vcf
```


## Running on Minerva (HPC cluster):

All NGS aligners need the reference sequences to be indexed. On the very first use of the pipeline with a reference genome, run

```bash
path/to/file/fastq_to_bam.sh -f <reference.fa> -r <fastq1,fastq2> -o <output_prefix> --index_ref
```
On  subsequent execution using the same reference, run without the `--index_ref` option. Keep reference and index files in the same directory.

### Usage 

Required modules:
- [bwa](https://sourceforge.net/projects/bio-bwa/files/)
- [samtools](https://github.com/samtools/samtools)
- java (version 1.8)
- [picard](https://broadinstitute.github.io/picard/)
- [gatk](https://gatk.broadinstitute.org/hc/en-us)

*Note:* required modules are already available on Minerva. This script loads them in, so no need to load any software modules before running:)

Required arguments:
```
-f          Reference genome fasta file
-r          Two read fastq files, separated by a comma (no space)
-o          Output files prefix
```

*Note:* current version assumes all reads in fastq files come from same sample and are run on same Illumina sequencing lane, and uses this to set the @RG tag in the BAM files. Read about the @RG tag [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups). 

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

#BSUB -J fastq2bam
#BSUB -P acc_ProjectName
#BSUB -q queue_name 
#BSUB -n 8
#BSUB -M 32000 
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=4000]
#BSUB -W 07:00
#BSUB -o fastq2bam_%J.out
#BSUB -e fastq2bam_%J.err
```
*Note:* Can reduce computational time (alignment is main bottleneck) by increasing number of cores. To take full advantage of cores, change the number of threads specified in lines 192-193 in the [fastq_to_bam script](https://github.com/verammaz/bioinformatics/blob/main/fastq_to_bam.sh): 
```bash
bwa mem -M -t 8 $REF_FASTA $READS_1 $READS_2 \
        -R "@RG\tID:${id}\tSM:${sample}\tPL:ILLUMINA" $(get_verbosity_flag bwa) | samtools sort -@8 - -o $RAW_BAM
```
to match the number of cores requested with `#BSUB -n`. All cores must be on same host node. 




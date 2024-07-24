# fastq_to_bam.sh
This script runs a FASTQ --> BAM pipeline by calling a sequence mapper / aligner and then carrying out optional post processing.

## Pipeline Overview:

[Burrows-Wheeler Alignment tool](https://bio-bwa.sourceforge.net/bwa.shtml)

0. Indexing the reference (bwa_index)
1. Alignment (bwa_mem) and conversion to BAM format (samtools)
2. [Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard)

Optional post-processing:

3. [GATK IndelRealigner](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md)
    - *Note*: for WES data, can speed up this step by using a exome_targets.interval_list file
4. [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - *Note:* vcf file with sites of variation required for this step. Specify path to this file in top of script in the 'SITES_OF_VARIATION' variable.
    - *Note*: for WES data, can speed up this step by using a exome_targets.interval_list file

## Reference File(s):

Recommended directory structure for reference files:
```bash
Reference/
    ├── genome.fasta
    ├── genome.fasta.{amb, ann, btw, fai, pac, sa}  #Output of bwa-index 
    ├── sites_of_variation.vcf
    ├── known_indels.vcf
    ├── exome_targets.interval_list
```


## Running on Minerva (HPC cluster):

First, specify file paths in the top of the script! 
All NGS aligners need the reference sequences to be indexed. On the very first use of the pipeline with a reference genome, run

```bash
path/to/file/fastq_to_bam.sh -r <fastq1,fastq2> -o <output_prefix> --index_ref
```
On  subsequent execution using the same reference, run without the `--index_ref` option. Keep reference and index files in the same directory.

### Usage 

Required modules:
- [bwa](https://sourceforge.net/projects/bio-bwa/files/)
- [samtools](https://github.com/samtools/samtools)
- java (version 1.8)
- [picard](https://broadinstitute.github.io/picard/)
- [gatk v3.6](https://gatk.broadinstitute.org/hc/en-us)

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
| `--patient` |  Patient id for sample data, required for downstream compatability during nextflow varianrt calling. Defualt none.
| `--index_ref` | Flag for running the indexing step for reference.
| `--keep_intermediate` | Flag to keep intermediate files produced during pipeline execution.
| `--post_process` | Flag to carry out post processing of initial alignment bam file (mark duplicates, indel local realignment, base recalibration).
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
#BSUB -W 20:00
#BSUB -o fastq2bam_%J.out
#BSUB -e fastq2bam_%J.err
```
*Note:* Can reduce computational time (alignment is main bottleneck) by increasing number of cores and calling the script with `--threads n`, where `n` matches the number of cores requested for the job. All cores must be on the same compute node. 




# union_variants.sh
This script takes a union of variants reported in VCF file(s),across multiple samples, and looks them up in each sample's BAM file. The result is a single list of variants, with total count and alt count data for each sample.

## Pipeline Overview:

0. Pre-Process: Extract all variants from all sample vcf files and construct a `regions.txt` file with lines "chrom start end".
1. Bam-readcount: Run the bam-readcount software on each sample bam file, passing in `regions.txt` file as input.
2. Post-Process: Extract bam-readcount outputs and write final VCF file(s).


## Running of Minerva:

Required modules (python, [bam-readcount](https://github.com/genome/bam-readcount)) are available on Minerva and are automatically loaded in by the script.

### Usage

### Submitting a Job (bsub)
Run the following in the terminal:

```bash
cd LynchCohortPipeline
source ./config
patient_id= #specify patient id
samplesheet= #specify path to samplesheet
job_name="union_variants_${patient_id}"
bsub -J ${job_name} \
     -P ${project} \
     -q ${queue} \
     -W 20:00 \
     -n 8 \
     -oo ${LOG_DIR}/{$job_name}.out \
     -eo ${LOG_DIR}/{$job_name}.err \
     ./union_variants.sh -p ${patient_id} -s ${samplesheet} --filter_variants -v
```

#### Required arguments:
```
-p         Patient identifier.
-s         CSV file with raw input data files configuration.
```

#### Optional arguments:

| Parameter                 | Description   |	
| :----------------------------------------: | :------: |
| `-v` | Enable verbose mode. |
| `-h` | Display usage message. |
| `--use_vcf` |  Flag to use VCF counts for called sample variants. 
| `--keep_intermediate` | Flag to keep intermediate files produced during pipeline execution.
| `--filter_variants ` | Flag to add filter to variants (in addition to 'PASS').
|  `--strelka_mutect_snv_intersect` | Only consider snv variants at intersection of strelka mutect callers.
| `--strelka_mutect_indel_intersect` | Only consider snv variants at intersection of strelka mutect callers.
|  `--single_output_file` | Write single output vcf file.
| `--step` | Starting step for pipeline execution. Default 0. Options are 0=preprocess, 1=bam-readcount, 2=postprocess). |

>Intersect options not yet implemented.

#### Input .csv file
`samplesheet.csv` needs to have the columns patient, sample, fastq_1, fastq_2, status (0=Normal, 1=Tumor), bam, bai, vcf. The vcf column should have file paths separated by a '|'. Example:

```csv
patient,sample,fastq_1,fastq_2,status,bam,bai,vcf
Patient1,Normal,full/path/to/Normal_R1_001.fastq.gz,full/path/to/Normal_R2_001.fastq.gz,0,full/path/to/Normal.bam,full/path/to/Normal.bai,na
Patient1,S1,full/path/to/S1_R1_001.fastq.gz,full/path/to/S1_R2_001.fastq.gz,1,full/path/to/S1.bam,full/path/to/S1.bai,full/path/to/S1_1.vcf|full/path/to/S1_2.vcf|full/path/to/S1_3.vcf
```



# Config File

Specify paths and global variables in the [config](config.sh) file!

Run the following:
```bash
cd LynchCohortPipelines
git update-index --assume-unchanged config.sh
```

# Home Directory Structure

    .
    ├── ...
    ├── Raw/                                  # Directory with patient raw data subdirectories
    │   ├── Patient/                          # Patient directory with sample subdirectories
    │   │    ├── samplesheet.csv              # File with information about the samples. Required.   
    │   │    ├── Sample1/                     # Sample directory
    │   │    │    ├── <Sample1>.bam
    │   │    │    ├── <Sample1>.bai
    │   │    │    ├── <Sample1_vs_Tumor>.strelka.somatic_indels.vcf.gz{.tbi}
    │   │    │    ├── <Sample1_vs_Tumor>.strelka.somatic_snvs.vcf.gz{.tbi}
    │   │    │    ├── <Sample1_vs_Tumor>.mutect2.filtered.vcf.gz{.tbi}
    │   │    │    └── ...
    │   │    ├── Normal/                      # Normal directory
    │   │    │    ├── <Normal>.bam
    │   │    │    ├── <Normal>.bai
    │   │    │    └── ...
    │   │    └── ...
    │   └── ...                 
    ├── VCF/                                   # Directory for processed VCF files
    │   ├── Patient/                           # Patient directory with VCF files for each sample
    │   │    ├── <Sample1>.vcf                 # Sample VCF files have same set of variants and in same order
    │   │    ├── <Sample2>.vcf          
    │   │    ├── <Sample1>_ann.vcf             # Annotated vcf files (snpeff, varcode)
    │   │    └── ...
    │   └── ... 
    ├── PairTrees/                             # Directory for pairtree outputs                 
    │   ├── Patient/                          
    │   │    ├── <Patient>.ssm   
    │   │    ├── <Patient>.params.json
    │   │    ├── <Patient>.plottree
    │   │    ├── <Patient>.npz  
    │   │    └── ... 
    │   └── ... 
    ├── Neoantigens/                           # Directory for NoePipe neoantigen calling outputs
    │   ├── neoantigens_<Patient>.txt
    │   ├── neoantigens_other_<Patient>.txt 
    │   └── ... 
    ├── Plots/                                 # Directory for visuals
    │   ├── TreePlots/                         # Directory for patient tree plots
    │   ├── VariantPlot/                       # Directory for patient variant plots
    │
    ├── cfit_config.json                       # Config file for CFIT
    ├── cfit_mapping.json                      # Mapping file for CFIT
    └── ... 

  
# Pipelines in this repository

Note that you should run the pipelines *per patient*.

### 1: Convert FASTQ (paired-end reads) to BAM file
[fastq --> bam readme](fastq->bam.md)

### 2: Variant calling (via nextflow nf-core/sarek)
[variant calling readme](variant_call.md)

### 3: Get union of variants across samples
[union variants readme](union_variants.md)

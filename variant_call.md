# nf-core/sarek 

Visit the official site [here](https://nf-co.re/sarek/3.4.2/) for full details about the nextflow nf-core/sarek pipeline

## Getting nextflow

1. Navigate to the directory where you want your nextflow executable to be located. Specify the path to this executable in the `NEXTFLOW` variable in your [config](config.sh) file.
2. Run the following:
```bash
wget -qO- https://get.nextflow.io | bash
```
3. Check the version:
```bash
/path/to/nextflow -v
```
4. Edit the file ~/.nextflow/assets/nf-core/sarek/nextflow.config by adding the following chunk:
```bash
executor {
    name   = 'local'
    cpus   = 48
    memory = '64GB'
}
```

# Variant Calling Wrapper Script

The script `variant_call.sh` submits nextflow jobs, executing the nf-core/sarek pipeline for multiple samples, per patient, in parallel. It reads input file information from `samplesheet.csv`. It then creates a sample specific `samplesheet_<Sample>.csv` and passes this as input to nextflow, executes nextflow, either from the mapping step or from the variant calling step, and finally moves VCF (and BAM) outputs to `$HOME_DIR/Raw/Patient/Sample`. File information in `samplesheet.csv` will be updated appropriately.

## Usage

#### Arguments:
```
-p         Patient identifier. Required.
-s         CSV file with raw input data files configuration. Required.
--step     Step to start nextflow nf-core/sarek from. 
           Options are 0=mapping, 1=variant_calling (default)
```


#### Input .csv file
`samplesheet.csv` needs to have the columns patient, sample, fastq_1, fastq_2, status (0=Normal, 1=Tumor), bam, bai. Example:

```csv
patient,sample,fastq_1,fastq_2,status,bam,bai
Patient1,Normal,full/path/to/Normal_R1_001.fastq.gz,full/path/to/Normal_R2_001.fastq.gz,0,full/path/to/Normal.bam,full/path/to/Normal.bai
Patient1,S1,full/path/to/S1_R1_001.fastq.gz,full/path/to/S1_R2_001.fastq.gz,1,full/path/to/S1.bam,full/path/to/S1.bai,
```

> Note that you don't need the bam and bai columns if running from the mapping step (step=0).
# nf-core/sarek 

Visit the official site [here](https://nf-co.re/sarek/3.4.2/) for full details about the nextflow nf-core/sarek pipeline

## Usage on Minerva

1. Navigate to the directory where you want your nextflow executable to be located.
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

The script `variant_call.sh` submits nextflow jobs, executing the nf-core/sarek pipeline for multiple samples, per patient, in parallel. It reads input file information from `samplesheet.csv`, passes this as input to nextflow, executes nextflow, either from the mapping step or from the variant calling step, and finally moves VCF (and BAM) outputs to analysis directory.

## Usage

First, change the file project, cores, and script path in the top of [variant_call.sh](variant_call.sh).

Next, change the paths in the top of [run_nextflow.sh](run_nextflow.sh).
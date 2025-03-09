# Pairtree method 

Developed by the Morris Lab.

For more details and full information about pairtree, visit:
 - [github repo](https://github.com/morrislab/pairtree)
 - [paper describing the algorithms](https://aacrjournals.org/bloodcancerdiscov/article/3/3/208/694689/Reconstructing-Complex-Cancer-Evolutionary)
 - [STAR protocol paper](https://pubmed.ncbi.nlm.nih.gov/36129821/)

### Pairtree pipeline overview
1. Cluster variants into subclones suitable for building clone trees.
2. Build clone trees.


## Installation on Minerva

```bash
git clone https://github.com/jwintersinger/pairtree
cd pairtree/lib
git clone https://github.com/ethanumn/projectppm
cd projectppm
bash make.sh
```

### To install the dependencies:
```bash
cd /path/to/pairtree
module load anaconda3
conda create -n pairtree --file requirements.txt --yes
```


## Running pairtree

There is a wrapper script [pairtree.sh](pairtree.sh) that prepares the input files for pairtree and runs the entire pairtree pipeline.

First, set the `PAIRTREE_ENV` variable in your config.sh.

> Note: you can get the full path my running `conda env list` and copying the path corresponding to 'pairtree'.


### Usage

- Input: Single sample VCF files
- Output: Combined sample tree (and optionally, single sample trees)

Make sure your analysis data home folder has a VCF/Patient subfolder with processed VCF files for all samples of the patient. Run the pairtree wrapper script with the following command:
```bash
cd LynchCohortPipelines
./pairtree.sh -p <patient_id> [-patient_sex_file <sex.txt> --single_sample_trees --total_alt]
```

#### Options:
`-patient_sex_file`: Tab-deliminated file with 'patient' and 'sex' columns (multiple patients can be in same file). If this file is not provided, sex will be inferred based on whether there are any Y chromosome variants. \
`--single_sample_trees`: Flag to produce single sample trees.
`--total_alt`: Flag that VCF file(s) in total:alt format. Default assumes total:ref format.

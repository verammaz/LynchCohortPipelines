# pairtree method 

Developed by the Morris Lab.

For more details and full information about pairtree, visit:
 - [github repo](https://github.com/morrislab/pairtree)
 - [paper describing the algorithms](https://aacrjournals.org/bloodcancerdiscov/article/3/3/208/694689/Reconstructing-Complex-Cancer-Evolutionary)
 - [STAR protocol paper](https://pubmed.ncbi.nlm.nih.gov/36129821/)

### pairtree pipeline overview
1. clustervars: cluster variants into subclones suitable for building clone trees.
2. pairtree: build clone trees.


## Installation on Minerva

1. Clone the repository:
```bash
cd ~/
git clone https://github.com/jwintersinger/pairtree  
```

2. Download and build the C code:
```bash
cd pairtree/lib
git clone https://github.com/ethanumn/projectppm
cd projectppm
bash make.sh
```

## Running pairtree

There is a wrapper script [pairtree.sh](pairtree.sh) that prepares the input files for pairtree and runs the entire pairtree pipeline.

- Input: Single sample VCF files
- Output: Combined sample tree (and optionally, single sample trees)

### Usage

Make sure your analysis data home folder has a VCF/Patient subfolder with processed VCF files for all samples of the patient. Run pairtree wrapper script with the following command:
```bash
cd LynchCohortPipelines
./pairtree -p <patient_id> [-patient_sex_file <sex.txt> --single_sample_trees]
```

#### Options:
`-patient_sex_file` Tab-deliminated file with 'patient' and 'sex' columns (multiple patients can be in same file). If this file is not provided, sex will be inferred based on whether there are any Y chromosome variants.
`--single_sample_trees` Flag to produce single sample trees.
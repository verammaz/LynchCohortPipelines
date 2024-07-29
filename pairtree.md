# pairtree method 

Developed by the Morris Lab.

For more details and full information about pairtree, visit:
 - [github repo](https://github.com/morrislab/pairtree)
 - [paper describing the algorithms](https://aacrjournals.org/bloodcancerdiscov/article/3/3/208/694689/Reconstructing-Complex-Cancer-Evolutionary)
 - [STAR protocol paper](https://pubmed.ncbi.nlm.nih.gov/36129821/)


## Installation on Minerva

1. Install requirements:
```bash
cd ~/
git clone https://github.com/jwintersinger/pairtree   
conda create -n pairtree --file requirements.txt --yes
conda activate pairtree
```

2. Download and build the C code:
```bash
cd ~/
git clone https://github.com/jwintersinger/pairtree
cd pairtree/lib
git clone https://github.com/ethanumn/projectppm
cd projectppm
bash make.sh
```

## Running pairtree

### Pipeline
1. clustervars: cluster variants into subclones suitable for building clone trees.
2. pairtree: build clone trees.

There is a wrapper script that prepares the input files for pairtree and runs the entire pairtree pipeline.

- Input: Singe sample VCF files
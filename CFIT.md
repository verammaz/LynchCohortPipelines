# CFIT

### Installation

> Note that this requires an invitation from the git repo owner Marta Luksza, PhD. 
> Git will ask for your username and access token.

After cloning the git repository, run:
```bash
module load python/3.7.3
cd CFIT
pip install . --user
```
To update, run:
```bash
module load python/3.7.3
module load git
cd CFIT
git pull
pip install . --user
```

### Usage
See git repo [here](https://github.com/LukszaLab/CFIT/).

There is a script [import_to_cfit.py](import_to_cfit.py) that loads pairtree tree outputs into CFIT and creates various tree plots. 


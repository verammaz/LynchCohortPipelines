
## netMHCpan-4.1 Installation

Follow the instructions to get the download here: https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/.


1. Once netMHCpan-4.1 is downloaded locally, you can `scp` the .tar folder to your Minerva home directory and untar it. This will produce a directory 'netMHCpan-4.1'.

2. Now, run: 
```bash
cd /path/to/netMHCpan-4.1/
wget https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz
tar -xvf data.tar.gz
rm data.tar.gz
```
This will produce a directory 'data'.  It is necessary  for the NetMHCpan 4.1 software to operate.

3. In the 'netMHCpan-4.1' directory edit the script 'netMHCpan':
   
    a. At the top of the file  locate the part labelled  "GENERAL SETTINGS:
        CUSTOMIZE TO YOUR SITE"  and set  the 'NMHOME' variable  to the full
	    path to the 'netMHCpan-4.1' directory on your system;

    b. Set TMPDIR to the full path to the temporary directory of your choice. It must
        be user-writable. You may for example set it to $NMHOME/tmp (and create
        the tmp folder in the netMHCpan-4.1 directory).

Add the path to this netmhc folder in your config.sh file in the variable `NETMHC`.

## snpEff v4.3t Installation

Run the following:

```bash
wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
unzip snpEff_v4_3t_core.zip
rm snpEff_v4_3t_core.zip
```

Add the path to this snpeff folder in your config.sh file in the variable `SNPEFF`.

# NeoPipe

### Installation

> Note that this requires an invitation from the git repo owner Marta Luksza, PhD. 
> Git will ask for your username and access token.

After cloning the git repository, run:
```bash
module load python/3.7.3
cd NeoPipe
pip install . --user
```
To update, run:
```bash
module load python/3.7.3
cd NeoPipe
git pull
pip install . --user
```


There is a wrapper script that calls various NeoPipe functionalities.
1. Annotate VCF files with snpEff and extract wildtype and mutant peptides.
2. Compute neoantigens using netMHC.
3. Prepare config files required for CFIT.

You can run it with the command    
```bash
cd LynchCohortPipelines
./neopipe.sh [OPTIONS] -p <patient_id> -sample_info <sample_info.txt>
```

#### Required arguments:
```
-p                  Patient identifier.
-sample_info        Tab-separated file with columns 'Sample' 'Patient' 'TimePoint' 'Tissue'
```

Below are the `sample_info.txt` file column meanings (do not include row for 'Normal'):
- Sample: sample identifier (should be consitent with the identifier in `samplesheet.csv`).
- Patient: patient identifier (should be consitent with the identifier in `samplesheet.csv`).
- TimePoint: one of Adenoma, Precancerous, Tumor, Normal, Normal(hysterectomy). If a patient has multiple samples with the same TimePoint label, add a number (e.g. Adenoma1, Adenoma2).
- Tissue: sample tissue origin.


#### Optional arguments:

| Parameter                 | Description   |	
| :----------------------------------------: | :------: |
| `-h` | Display usage message. |
| `-ns`| Neonatigen peptide lengths to consider. Default is 9. Can also be specified as a range (e.g. 8-14).
| `-kd` | Kd (in nM) threshold to use for neoantigen computations. Defualt is 50000 nM. |
| `-hla` | Path to patient HLA calls file. By default, assumes this file is $HOME_DIR/HLA/HLA_calls.txt.
| `--single_sample_trees` | Flag to add single sample tree config infos to CFIT config files. |





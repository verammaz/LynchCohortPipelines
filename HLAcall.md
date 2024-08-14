There is a script that submits an hla typing job for a patient by calling [optitype](https://github.com/FRED-2/OptiType).

> Optitype is installed and available to use on Minerva. The script loads in the required module(s).

# Usage
First, make sure to specify the path to the optitype config file in your `config.sh` in the OPTITYE_CONFIG variable. Example optitype config file for use on Minerva:

```bash
[mapping]

# Absolute path to RazerS3 binary, and number of threads to use for mapping

razers3=/hpc/packages/minerva-common/razers/3.5.4/razers3-3.5.4-Linux-x86_64_sse4/bin/razers3
threads=4

[ilp]

# A Pyomo-supported ILP solver. The solver must be globally accessible in the
# environment OptiType is run, so make sure to include it in PATH.
# Note: this is NOT a path to the solver binary, but a keyword argument for
# Pyomo. Examples: glpk, cplex, cbc.

solver=glpk
threads=1

[behavior]

# tempdir=/path/to/tempdir  # we may enable this setting later. Not used now.

# Delete intermediate bam files produced by RazerS3 after OptiType finished
# loading them. If you plan to re-analyze your samples with different settings
# disabling this option can be a time-saver, as you'll be able to pass the bam
# files to OptiType directly as input and spare the expensive read mapping
# step.

deletebam=true

# In paired-end mode one might want to use reads with just one mapped end (e.g.,
# the other end falls outside the reference region). This setting allows the
# user to keep them with an optionally reduced weight. A value of 0 means they
# are discarded for typing, 0.2 means single reads are "worth" 20% of paired
# reads, and a value of 1 means they are treated as valuable as properly mapped
# read pairs. Note: unpaired reads will be reported on the result coverage plots
# for completeness, regardless of this setting.

unpaired_weight=0

# We call a read pair discordant if its two ends best-map to two disjoint sets
# of alleles. Such reads can be either omitted or either of their ends treated
# as unpaired hits. Note: discordant read pairs are reported on the coverage
# plots as unpaired reads, regardless of this setting.

use_discordant=false
```
## Command: 
```bash
cd LynchCohortPipelines
./submit_hlatyping_for_patient.sh [--normal_only] -p <patient_id> -s <samplesheet.csv>
```

## Arguments:
```
-p         Patient identifier. Required.
-s         CSV file with raw input data files configuration. Required.
--normal_only    (Optional) Run hla typing on patient normal sample only.
```

If a file `${HOME_DIR}/HLA/HLA_calls.txt` doesn't exist, it will be created. The file is in tab-separated format with columns patient_id and HLA_alleles. If the file exists, a new line will be added for the patient (hla types taken from normal sample).
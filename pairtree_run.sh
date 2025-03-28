#!/bin/bash
source ./config.sh

module purge  

module load anaconda3

# manually source Conda (instead of running 'conda init')
export CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Debug: Check if Conda is set up
conda info --envs

# Activate the environment
conda activate $PAIRTREE_ENV
echo $PAIRTREE_ENV

# Debug: Verify activation
conda info --envs
conda list


########################################################################################################################
# input:
#   patient id 
#   data_dir
#
# output:
#    $dir/${patient}.ssm
#    $dir/${patient}.params.json
#    $dir/${patient}_results.npz
#    $dir/${patient}_plottree.html
#    $dir/${patient}_summposterior.html
########################################################################################################################

# pairtree scripts
scripts_path=$1

# patient id
ID=$2

# folder path for outputs 
OUT_DIR=$3

#######################################################################
# --------------------- Constants: file names ----------------------- #
#######################################################################

FN_SSM=$ID".ssm"

FN_PARAMS=$ID".params.json"

FN_RESULTS_NPZ=$ID"_results.npz"

FN_HTML_PLOTTREE=$ID"_plottree.html"

FN_HTML_SUMMPOSTERIOR=$ID"_summposterior.html"

#######################################################################
# --------------------- Start of pipeline  -------------------------- #
#######################################################################


# ------------- Print out welcome message :) -------------- #
echo ""
echo "------------------------------------------------------"
echo "Running pairtree pipeline script for $ID"
echo "------------------------------------------------------"


# -------------------- Create output directory ---------------#
if [ ! -d "$OUT_DIR" ]; then
    # Directory does not exist, so create it
    mkdir -p "$OUT_DIR"
    echo && echo "Creating output directory"
else
    echo && echo "Output directory already exists"
fi


# -------------------- Run clustervars ---------------------- #
echo && echo "Running bin/clustervars"
python $scripts_path/bin/clustervars $OUT_DIR/$FN_SSM $OUT_DIR/$FN_PARAMS $OUT_DIR/$FN_PARAMS
                         
# -------------------- Run pairtree -------------------- #
echo && echo "Running bin/pairtree"
python $scripts_path/bin/pairtree --params $OUT_DIR/$FN_PARAMS $OUT_DIR/$FN_SSM $OUT_DIR/$FN_RESULTS_NPZ --seed 5555
                         
# -------------------- Run plottree -------------------- #
echo && echo "Running bin/plottree"
python $scripts_path/bin/plottree --runid $ID $OUT_DIR/$FN_SSM $OUT_DIR/$FN_PARAMS $OUT_DIR/$FN_RESULTS_NPZ $OUT_DIR/$FN_HTML_PLOTTREE

# -------------------- Run summposterior ---------------- #
echo && echo "Running bin/summposterior"
python $scripts_path/bin/summposterior --runid $ID $OUT_DIR/$FN_SSM $OUT_DIR/$FN_PARAMS $OUT_DIR/$FN_RESULTS_NPZ $OUT_DIR/$FN_HTML_SUMMPOSTERIOR

# ----------------------- Completion message ------------------------ #
echo && echo "Completed."

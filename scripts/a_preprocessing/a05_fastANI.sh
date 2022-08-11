#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/analysis/01_fastANI"

fastANI --ql ${PROJECT_DIR}/data/reference_data/genome_fna.list \
--rl ${PROJECT_DIR}/data/reference_data/genome_fna.list \
--threads 30 --matrix -o ${ANALYSIS_DIR}/ANI_results

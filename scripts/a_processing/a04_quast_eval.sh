#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate omics_py37

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/busco_eval"


sampleId=$1
file_fna="$PROJECT_DIR/data/prokka_annotation/$sampleId/$sampleId.fna"
file_gff="$PROJECT_DIR/data/prokka_annotation/$sampleId/$sampleId.gff"

if [ -f ${file_fna} ] || [ -f ${file_gff} ]
then
    ls "${file_fna} ${file_gff}"
    error_exit $?
fi


process_start "Quast on sample $sampleId"

quast --output-dir data/quast/${sampleId} -t 8 --silent \
--labels ${sampleId} --features ${file_gff} ${file_fna}

error_exit $?

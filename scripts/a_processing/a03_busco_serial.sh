#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/busco_eval"

busco_lineage="data/busco_downloads/lineages/enterobacterales_odb10"

## process file one by one
for sampleId in `cat data/reference_data/data/reference_data/assembly_ids.txt`
do
    file_faa="$PROJECT_DIR/data/prokka_annotation/$sampleId/$sampleId.faa"

    if [ ! -f ${file_faa} ]
    then
        ls "${file_faa}"
        error_exit $?
    fi

    process_start "BUSCO on file $sampleId"
    
    busco -i ${file_faa} -o ${sampleId} -l ${busco_lineage} -m proteins \
    --offline --download_path $PROJECT_DIR/data/busco_downloads \
    --cpu 8 --tar --out_path ${ANALYSIS_DIR} -f --quiet
    
    error_exit $?

done


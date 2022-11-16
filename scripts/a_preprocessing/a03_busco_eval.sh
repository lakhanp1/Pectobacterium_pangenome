#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
BUSCOP_DIR="$PROJECT_DIR/data/busco.prot"
BUSCOG_DIR="$PROJECT_DIR/data/busco.geno"

busco_lineage="data/busco_downloads/lineages/enterobacterales_odb10"

sampleId=$1
file_faa="$PROJECT_DIR/data/prokka_annotation/$sampleId/$sampleId.faa"
file_fna="$PROJECT_DIR/data/prokka_annotation/$sampleId/$sampleId.fna"

if [ ! -f ${file_faa} ]
then
    ls "${file_faa}"
    error_exit $?
fi

## BUSCO protein
file_buscop_log=${BUSCOP_DIR}/${sampleId}/logs/busco.log

if [ ! -f ${file_buscop_log} ] || ! grep -q 'BUSCO analysis done' ${file_busco_log}
then
    conda activate pantools_master
    process_start "BUSCO protein on file $sampleId"

    busco -i ${file_faa} -o ${sampleId} -l ${busco_lineage} -m proteins \
    --offline --download_path $PROJECT_DIR/data/busco_downloads \
    --cpu 8 --tar --out_path ${BUSCOP_DIR} -f --quiet

    error_exit $?
else
    echo "BUSCO results exists for $sampleId"
fi


## BUSCO genome
file_buscog_log=${BUSCOG_DIR}/${sampleId}/logs/busco.log

if [ ! -f ${file_buscog_log} ] || ! grep -q 'BUSCO analysis done' ${file_busco_log}
then
    conda activate pantools
    process_start "BUSCO genome on file $sampleId"

    busco -i ${file_fna} -o ${sampleId} -l ${busco_lineage} -m genome \
    --offline --download_path $PROJECT_DIR/data/busco_downloads \
    --cpu 8 --tar --out_path ${BUSCOG_DIR} -f --quiet

    error_exit $?
else
    echo "BUSCO results exists for $sampleId"
fi


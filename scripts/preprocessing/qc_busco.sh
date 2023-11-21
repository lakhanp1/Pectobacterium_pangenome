#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh

## Setup
PROJECT_DIR="."
BUSCOG_DIR=$(yq '.data.dir , .data.busco_genome.dir' project_config.yaml | paste -d "/" - -)
BUSCOG_DIR="$PROJECT_DIR/$BUSCOG_DIR"

BUSCOP_DIR=$(yq '.data.dir , .data.busco_protein.dir' project_config.yaml | paste -d "/" - -)
BUSCOP_DIR="$PROJECT_DIR/$BUSCOP_DIR"

PROKKA_DIR=$(yq '.data.dir , .data.prokka.dir' project_config.yaml | paste -d "/" - -)
PROKKA_DIR="$PROJECT_DIR/$ANALYSIS_DIR"

busco_lineage=$(yq '.data.busco_db.lineage' project_config.yaml)

sampleId=$1
file_faa="$PROKKA_DIR/$sampleId/${sampleId}.faa"
file_fna="$PROKKA_DIR/$sampleId/${sampleId}.fna"

if [ ! -f ${file_faa} ]; then
    ls "${file_faa}"
    error_exit $?
fi

## BUSCO protein
file_buscop_log=${BUSCOP_DIR}/${sampleId}/logs/busco.log

if [ ! -f ${file_buscop_log} ] || ! grep -q 'BUSCO analysis done' ${file_busco_log}; then
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

if [ ! -f ${file_buscog_log} ] || ! grep -q 'BUSCO analysis done' ${file_busco_log}; then
    conda activate pantools
    process_start "BUSCO genome on file $sampleId"

    busco -i ${file_fna} -o ${sampleId} -l ${busco_lineage} -m genome \
        --offline --download_path $PROJECT_DIR/data/busco_downloads \
        --cpu 8 --tar --out_path ${BUSCOG_DIR} -f --quiet

    error_exit $?
else
    echo "BUSCO results exists for $sampleId"
fi

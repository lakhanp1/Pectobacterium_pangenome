#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh

## Setup
PROJECT_DIR="."

QUAST_DIR=$(yq '.data.dir , .data.quast.dir' project_config.yaml | paste -d "/" - -)
QUAST_DIR="$PROJECT_DIR/$QUAST_DIR"

PROKKA_DIR=$(yq '.data.dir , .data.prokka.dir' project_config.yaml | paste -d "/" - -)
PROKKA_DIR="$PROJECT_DIR/$QUAST_DIR"

sampleId=$1
file_fna="$PROKKA_DIR/$sampleId/${sampleId}.fna"
file_gff="$PROKKA_DIR/$sampleId/${sampleId}.gff"

if [ -f ${file_fna} ] || [ -f ${file_gff} ]; then
    ls "${file_fna}" "${file_gff}"
    error_exit $?
fi

file_report=${QUAST_DIR}/${sampleId}/report.txt

if [ ! -f ${file_report} ]; then
    conda activate omics_py37
    process_start "Quast on sample $sampleId"

    quast --output-dir ${QUAST_DIR}/${sampleId} -t 8 --silent \
        --labels ${sampleId} --features ${file_gff} ${file_fna}

    error_exit $?
else
    echo "Quast results exists for $sampleId"
fi

#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
# Task: Build functional annotation using InterProScan

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
INTERPRO_DIR=$(yq '.data.dir , .data.interproscan.dir' project_config.yaml | paste -d "/" - -)
INTERPRO_DIR="$PROJECT_DIR/$INTERPRO_DIR"

sampleId=$1
file_faa="$PROJECT_DIR/data/prokka_annotation/$sampleId/$sampleId.faa"

if [ ! -f ${file_faa} ]; then
    ls "${file_faa}"
    error_exit $?
fi

## replace * in the FASTA file
sed -i 's/*//g' ${file_faa}

prefix=${sampleId}.interProScan

if [ ! -f $INTERPRO_DIR/${prefix}.gff3 ]; then
    conda activate interproscan
    process_start "InterProScan on file $prefix"

    interproscan.sh -cpu 8 -goterms -f GFF3 -iprlookup \
        --output-file-base $INTERPRO_DIR/${prefix} -i $file_faa

    error_exit $?
else
    echo "${prefix}.gff3 InterProScan results exists"
fi

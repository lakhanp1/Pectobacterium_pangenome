#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
# Task: Build functional annotation using InterProScan

shopt -s expand_aliases
source ~/.bash_aliases

source /mnt/scratch/parde001/tools/miniconda3/etc/profile.d/conda.sh
conda activate interproscan

## Setup
PROJECT_DIR='/mnt/scratch/parde001/projects/03_Pectobacterium'
ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_50g/interProScan"
pan_db="$PROJECT_DIR/analysis/04_pangenome_pecto_50g/pectobacterium.50g.DB"


function run_interproscan(){
    ## replace * in the FASTA file
    sed -i 's/*//g' ${1}
    # sed '/>/{H;$!d} ; x ; s/^/\nSTART-->/ ; s/$/\n<--END/'

    file_outPrefix=`basename ${1} | sed 's/.fasta//'`.interproscan
    interproscan.sh -verbose -cpu 12 -goterms -f GFF3 -iprlookup -dp  --output-file-base $2/${file_outPrefix} -i $1
    error_exit $?
}

export -f run_interproscan

## Run InterProScan on all files using GNU parallel
ls ${pan_db}/proteins/*.fasta | parallel --keep-order --jobs 4 --halt now,fail=1 --results $PROJECT_DIR/logs/interproscan/{/.} --joblog $PROJECT_DIR/logs/interproscan/parallel.log run_interproscan {} ${ANALYSIS_DIR}



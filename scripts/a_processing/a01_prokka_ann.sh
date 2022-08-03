#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
shopt -s extglob
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate prokka

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/prokka_annotation"

file_faa=$1

if [ ! -f ${file_faa} ]
then
    ls "${file_faa}"
    error_exit $?
fi

prefix=`basename -z ${file_faa} | sed -r 's/(.*)\..*$/\1/' | sed -r 's/_genomic//'`
out_dir=${ANALYSIS_DIR}/${prefix}

if [ ! -f ${out_dir}/${prefix}.log ] || ! grep -q 'Annotation finished successfully.' ${out_dir}/${prefix}.log
then
    process_start "prokka $prefix"
    prokka --cpus 8 --outdir ${out_dir} --prefix ${prefix} ${file_faa}
    error_exit $?
else
    echo "prokka results exists"
fi

#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
shopt -s extglob
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh

## Setup
PROJECT_DIR="."
PROKKA_DIR=$(yq '.data.dir , .data.prokka.dir' project_config.yaml | paste -d "/" - -)
PROKKA_DIR="$PROJECT_DIR/$PROKKA_DIR"

file_fna=$1

if [ ! -f ${file_fna} ]; then
    ls "${file_fna}"
    error_exit $?
fi

prefix=$(basename -z ${file_fna} | sed -r 's/(.*)\..*$/\1/' | sed -r 's/_genomic//')
out_dir=${PROKKA_DIR}/${prefix}

if [ ! -f ${out_dir}/${prefix}.log ] || ! grep -q 'Annotation finished successfully.' ${out_dir}/${prefix}.log; then
    conda activate prokka
    process_start "prokka $prefix"
    prokka --cpus 8 --outdir ${out_dir} --prefix ${prefix} ${file_fna}
    error_exit $?
else
    echo "prokka results exists"
fi

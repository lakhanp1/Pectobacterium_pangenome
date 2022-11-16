#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master

## Setup
file_fna=$1
outPrefix=$2

if [ ! -f ${file_fna} ]
then
    ls "${file_fna}"
    error_exit $?
fi

if [ ! -d `dirname ${outPrefix}` ]
then
    ls `dirname ${outPrefix}`
    error_exit $?
fi


fastANI --ql ${file_fna} --rl ${file_fna} --threads 30 --matrix -o ${outPrefix}


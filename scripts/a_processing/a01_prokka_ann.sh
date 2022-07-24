#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
shopt -s extglob
source ~/.bash_aliases

source /mnt/scratch/parde001/tools/miniconda3/etc/profile.d/conda.sh
conda activate prokka

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/prokka_annotation"


function run_prokka(){
    prefix=`basename -z ${1} | sed -r 's/(.*)\..*$/\1/' | sed -r 's/_genomic//'`

    if [ ! -f ${2}/${prefix}/${prefix}.log ] || ! grep -q 'Annotation finished successfully.' ${2}/${prefix}/${prefix}.log
    then
        process_start "prokka $prefix"
        prokka --cpus 8 --outdir ${2}/${prefix} --prefix ${prefix} ${1}
        error_exit $?
    else
        echo "prokka results exists"
    fi
}
 
export -f run_prokka

## Run Prokka on all files using GNU parallel
#ls $PROJECT_DIR/data/genomes/?(ncbi|local)/*.?(fa|fna|fasta) | \
#parallel --gnu --resume --keep-order --jobs 4 --halt now,fail=1 --results $PROJECT_DIR/logs/prokka/{/.} --joblog $PROJECT_DIR/logs/prokka/parallel.log run_prokka {} ${ANALYSIS_DIR}

for fa in $PROJECT_DIR/data/genomes/?(ncbi|local)/*.?(fa|fna|fasta)
do
run_prokka ${fa} ${ANALYSIS_DIR}
done


#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

conda activate pantools

export PANTOOLS="$PANTOOLS_MAIN"

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/analysis/06_pangenome_pecto_v2"
analysis_prefix='pectobacterium.v2'
pan_db="$ANALYSIS_DIR/${analysis_prefix}.DB"

[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR

## Download the raw data
esearch -db assembly -query "pectobacterium [ORGN]" | \
efetch -format docsum -email "lakhansing.pardeshi@wur.nl" | \
xtract -pattern DocumentSummary \
-tab " || " -element FtpPath_RefSeq | \
sed -r 's/((.*\/)([^\/]+))/\1\/\3_genomic.fna.gz/' | \
sed -r 's/^ftp:/rsync:/' > $PROJECT_DIR/data/reference_data/ncbi_genomes_ftp.txt

# data/genomes/ncbi_genomes
for file in `cat $PROJECT_DIR/data/reference_data/ncbi_genomes_ftp.txt`
do
rsync --copy-links --times --verbose ${file} $PROJECT_DIR/data/genomes/ncbi_genomes/
sleep 5
done



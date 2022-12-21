#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

######################################################################
usage='bash b04_subset_optimal_grouping.sh <setId> <genomes>
setId    : string identifier for random set
genomes  : COMMA separated list of genome IDs in the pangenome to be
           used for optimal_grouping
'

if [ $# -ne 2 ] ; then printf "Error: Require 2 args\n${usage}" >&2 ; exit 1 ; fi

setId=$1
genomes=$2

######################################################################
source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master

export PANTOOLS="$PANTOOLS_MASTER"

PANGENOME_NAME='pectobacterium.v2'

## Setup
PROJECT_DIR="/lustre/BIF/nobackup/$USER/projects/03_Pectobacterium"
PANGENOME_DIR="$PROJECT_DIR/data/pangenomes/$PANGENOME_NAME"

## setup for local disk processing on handelsman
PAN_BUILD_DIR="/local_scratch/$USER/03_Pectobacterium"
# PAN_BUILD_DIR="/dev/shm/$USER/03_Pectobacterium"


[ ! -d $PAN_BUILD_DIR ] && mkdir -p $PAN_BUILD_DIR

pan_db="$PAN_BUILD_DIR/${PANGENOME_NAME}.${setId}.DB"
cp -r $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_bsk ${pan_db}

## optimized grouping
process_start optimal_grouping
nice $PANTOOLS optimal_grouping -t 8 --include $genomes ${pan_db} ${pan_db}/busco/enterobacterales_odb10
error_exit $?

cp -r ${pan_db}/optimal_grouping $PROJECT_DIR/analysis/03_pangenome_pecto_v2/subset_optimal_group/optimal_grouping_${setId}
rm -r ${pan_db}



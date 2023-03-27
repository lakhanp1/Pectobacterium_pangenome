#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

scriptName=`basename ${BASH_SOURCE[0]}`

usage="USAGE:
-------------------------------------------------------------
bash ${scriptName} <db_name> <exists> <suffix>
db_name   : STRING name for pangenome database
-------------------------------------------------------------
"

if [ $# -ne 1 ]; then
    printf "Error: Require pangenome name\n${usage}" >&2 ; exit 1
fi

######################################################################

PANGENOME_NAME=$1
# PANGENOME_NAME='pectobacterium.v2'

## Setup
PROJECT_DIR="/lustre/BIF/nobackup/$USER/projects/03_Pectobacterium"
PANGENOME_DIR="$PROJECT_DIR/data/pangenomes/$PANGENOME_NAME"

## setup for local disk processing on handelsman
LOCAL_DIR_PATH="/local_scratch/$USER"
# LOCAL_DIR_PATH="/local/$USER"
# LOCAL_DIR_PATH="/dev/shm/$USER"

# PAN_BUILD_DIR="$LOCAL_DIR_PATH/03_Pectobacterium"
PAN_BUILD_DIR=$PANGENOME_DIR
pan_db="$PAN_BUILD_DIR/${PANGENOME_NAME}.DB${DB_SUFFIX}"

printf "\$PANGENOME_DIR: ${PANGENOME_DIR}
\$PAN_BUILD_DIR: ${PAN_BUILD_DIR}
\$pan_db: ${pan_db}
"

hg_aln_dir="${pan_db}/alignments/msa_per_group/grouping_v1"


#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

conda activate pantools_master

export PANTOOLS="$PANTOOLS_OPT"

PANGENOME_NAME='pectobacterium.v2'

## Setup
PROJECT_DIR="/lustre/BIF/nobackup/$USER/projects/03_Pectobacterium"
PANGENOME_DIR="$PROJECT_DIR/data/pangenomes/$PANGENOME_NAME"

## setup for local disk processing on handelsman
BUILD_DIR="/local/$USER/03_Pectobacterium"

pan_db="$BUILD_DIR/${PANGENOME_NAME}.DB"


[ ! -d $PANGENOME_DIR ] && mkdir $PANGENOME_DIR
[ ! -d $PANGENOME_DIR/backup ] && mkdir $PANGENOME_DIR/backup
[ ! -d $BUILD_DIR ] && mkdir -p $BUILD_DIR



TMPDIR="/local/$USER/tmp"
[ ! -d ${TMPDIR} ] && mkdir -p ${TMPDIR}
[ -d ${TMPDIR}/pantools ] && rm -rd ${TMPDIR}/pantools
[ ! -d ${TMPDIR}/pantools ] && mkdir -p ${TMPDIR}/pantools
[ -d ${TMPDIR}/spark ] && rm -rd ${TMPDIR}/spark
[ ! -d ${TMPDIR}/spark ] && mkdir -p ${TMPDIR}/spark

cache=25000000
PANTOOLS_LOCALIZATION_OUTPUT=${TMPDIR}/pantools
SPARK_LOCAL_DIRS=${TMPDIR}/spark
PANTOOLS_ADDRESS_CACHE_SIZE=${cache}
PANTOOLS_NUM_DB_WRITER_THREADS=10

## Construct pangenome using build_pangenome_parallel
process_start build_pangenome_parallel
$PANTOOLS build_pangenome_parallel --database-path ${pan_db} \
--genomes-file $PANGENOME_DIR/genomes_fa.list -tn 40
error_exit $?

cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_raw





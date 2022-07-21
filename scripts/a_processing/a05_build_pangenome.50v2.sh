#!/bin/bash

shopt -s expand_aliases
source ~/.bash_aliases

#set -e
#set -u
set -o pipefail

conda activate pantools_dev
export PANTOOLS="$PANTOOLS_DEV"

## Setup
PROJECT_DIR='/mnt/scratch/parde001/projects/03_Pectobacterium'
ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_50g"
analysis_prefix='pectobacterium.50gv2'
pan_db="$ANALYSIS_DIR/${analysis_prefix}.DB"


#count=1
#rm $ANALYSIS_DIR/genomes_gff3.list $ANALYSIS_DIR/genomes_fa.list
#touch $ANALYSIS_DIR/genomes_gff3.list $ANALYSIS_DIR/genomes_fa.list

#for i in `cat $ANALYSIS_DIR/genomes.list`
#do
#ls $PROJECT_DIR/analysis/01_prokka_annotation/${i}/*.{fna,gff3}
#printf "$PROJECT_DIR/analysis/01_prokka_annotation/${i}/${i}.fna\n" >> $ANALYSIS_DIR/genomes_fa.list
#printf "$count\t$PROJECT_DIR/analysis/01_prokka_annotation/${i}/${i}.gff3\n" >> $ANALYSIS_DIR/genomes_gff3.list
#((count=count+1))
#done

### Construct pangenome
#process_start build_pangenome
#$PANTOOLS build_pangenome --database-path ${pan_db} --genomes-file $ANALYSIS_DIR/genomes_fa.list
#error_exit $?

## add annotations
process_start add_annotations
$PANTOOLS add_annotations --database-path ${pan_db} -af $ANALYSIS_DIR/genomes_gff3.list -ca
error_exit $?





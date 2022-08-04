#!/bin/bash


shopt -s expand_aliases
source ~/.bash_aliases

source /mnt/scratch/parde001/tools/miniconda3/etc/profile.d/conda.sh
conda deactivate & conda activate pantools

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/analysis/03_pangenome_pecto_6g"
analysis_prefix='pectobacterium.6g'
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
#pantools build_pangenome --database-path ${pan_db} --genomes-file $ANALYSIS_DIR/genomes_fa.list
#error_exit $?

### add annotations
#process_start add_annotations
#pantools add_annotations --database-path ${pan_db} -af $ANALYSIS_DIR/genomes_gff3.list -ca
#error_exit $?

### BUSCO
#process_start busco_protein
#pantools busco_protein -dp ${pan_db} -if enterobacterales_odb10 -tn 20 
#error_exit $?

### optimized grouping
#process_start optimal_grouping
#pantools optimal_grouping -dp ${pan_db} -if ${pan_db}/busco/enterobacterales_odb10 -tn 20
#error_exit $?

Rscript ${pan_db}/optimal_grouping/optimal_grouping.R

## add homology group
process_start add_homology_group
pantools group --database-path ${pan_db} -tn 20
error_exit $?


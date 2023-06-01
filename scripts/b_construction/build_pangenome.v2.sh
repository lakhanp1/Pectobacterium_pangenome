#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

scriptName=`basename $0`

usage="USAGE:
-------------------------------------------------------------
bash ${scriptName} <db_name> <exists> <suffix>
db_name   : STRING name for pangenome database
exists    : INT (0, 1) flag whether the db exists or not
suffix    : STRING suffix for pangenome db_name
-------------------------------------------------------------
"

if [ $# -ne 2 ] && [ $# -ne 3 ]; then
    printf "Error: Require at least 2 args\n${usage}" >&2 ; exit 1
fi

if [ $2 -ne 1 ] && [ $2 -ne 0 ]; then
    printf "Error: exists should be either 0 or 1\n${usage}" >&2
    exit 1
fi

# DB_SUFFIX=""
if [ $# -eq 3 ]; then
    DB_SUFFIX="_"$3
fi

######################################################################

PANGENOME_NAME=$1
# PANGENOME_NAME='pectobacterium.10g'
# PANGENOME_NAME='pectobacterium.v2'

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master

export PANTOOLS="$PANTOOLS_MASTER"

## Setup
PROJECT_DIR="/lustre/BIF/nobackup/$USER/projects/03_Pectobacterium"
PANGENOME_DIR="$PROJECT_DIR/data/pangenomes/$PANGENOME_NAME"

## setup for local disk processing on handelsman
LOCAL_DIR_PATH="/local_scratch/$USER"
# LOCAL_DIR_PATH="/local/$USER"
# LOCAL_DIR_PATH="/dev/shm/$USER"

PAN_BUILD_DIR="$LOCAL_DIR_PATH/03_Pectobacterium"
pan_db="$PAN_BUILD_DIR/${PANGENOME_NAME}.DB${DB_SUFFIX}"

printf "PANGENOME_DIR: ${PANGENOME_DIR}
BUILD_DIR: ${PAN_BUILD_DIR}
pangenome: ${pan_db}
"

if [ $2 -eq 1 ]; then
    printf "Using existing pangenome ${PANGENOME_NAME}\n"
    if [ ! -d ${pan_db} ]; then
        printf "Could not find pangenome db: ${pan_db}\n" >&2
        exit 1
    fi
elif [ $2 -eq 0 ]; then
    printf "Building pangenome ${PANGENOME_NAME}\n"
    [ ! -d $PANGENOME_DIR ] && mkdir $PANGENOME_DIR
    [ ! -d $PANGENOME_DIR/backup ] && mkdir $PANGENOME_DIR/backup
    [ ! -d $PAN_BUILD_DIR ] && mkdir -p $PAN_BUILD_DIR
fi

######################################################################

#export PANTOOLS="$PANTOOLS_OPT"
#TMPDIR="/local/$USER/tmp"
#[ ! -d ${TMPDIR} ] && mkdir -p ${TMPDIR}
#[ -d ${TMPDIR}/pantools ] && rm -rd ${TMPDIR}/pantools
#[ ! -d ${TMPDIR}/pantools ] && mkdir -p ${TMPDIR}/pantools
#[ -d ${TMPDIR}/spark ] && rm -rd ${TMPDIR}/spark
#[ ! -d ${TMPDIR}/spark ] && mkdir -p ${TMPDIR}/spark

#cache=25000000
#PANTOOLS_LOCALIZATION_OUTPUT=${TMPDIR}/pantools
#SPARK_LOCAL_DIRS=${TMPDIR}/spark
#PANTOOLS_ADDRESS_CACHE_SIZE=${cache}
#PANTOOLS_NUM_DB_WRITER_THREADS=10

### Construct pangenome using build_pangenome_parallel
#process_start build_pangenome_parallel
#$PANTOOLS build_pangenome_parallel --database-path ${pan_db} \
#--genomes-file $PANGENOME_DIR/genomes_fa.list -tn 40
#error_exit $?

#cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_raw

### add annotations
#process_start add_annotations
#$PANTOOLS add_annotations --connect ${pan_db} $PANGENOME_DIR/genomes_gff3.list
#error_exit $?

#cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_ann

## add phenotypes
#process_start add_phenotypes
#printf "y\n" | $PANTOOLS remove_phenotypes ${pan_db}
#$PANTOOLS remove_phenotypes ${pan_db}
#$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/genomes_metadata.csv
#error_exit $?

#cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_meta

### add_functions
#process_start add_InterProScan_annotations
#$PANTOOLS add_functions ${pan_db} $PANGENOME_DIR/functional_annotations.txt
#error_exit $?

#cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_fn

### BUSCO
#process_start busco_protein
#$PANTOOLS busco_protein -t 30 --busco10 enterobacterales_odb10 ${pan_db}
#error_exit $?

#cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_bsk

## BUSCO multiqc
#process_start multiqc_busco
#multiqc --filename multiqc_busco -outdir analysis/01_multiqc --ignore "run_enterobacterales_odb10" --module busco -d -dd 1 ${pan_db}/busco/enterobacterales_odb10/protein/results/
#error_exit $?

## optimized grouping
#process_start optimal_grouping
#$PANTOOLS optimal_grouping -t 30 ${pan_db} ${pan_db}/busco/enterobacterales_odb10
#$PANTOOLS optimal_grouping -dp ${pan_db} -if ${pan_db}/busco/enterobacterales_odb10 -tn 25
#error_exit $?

#cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_grp

#cp -r $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_grp ${pan_db}

#Rscript ${pan_db}/optimal_grouping/optimal_grouping.R
#$PANTOOLS grouping_overview ${pan_db}

#$PANTOOLS change_grouping -v 4 ${pan_db}
#$PANTOOLS grouping_overview ${pan_db}

#cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_grp

##################
### extract the metrics from pangenome
#process_start extract_pangenome_metrics
#$PANTOOLS metrics ${pan_db}
#error_exit $?

######################################################################

### Core unique thresholds
#process_start core_unique_thresholds
#$PANTOOLS core_unique_thresholds ${pan_db}
#error_exit $?
#Rscript ${pan_db}/core_unique_thresholds/core_unique_thresholds.R

### gene classification: core and unique
#process_start gene_classification_core_unique
#$PANTOOLS gene_classification ${pan_db} 
#error_exit $?

### Gene distance tree
#Rscript ${pan_db}/gene_classification/gene_distance_tree.R
#mv ${pan_db}/gene_classification ${pan_db}/gene_classification.100.0


### gene classification: soft core and cloud
#process_start gene_classification_softcore_cloud
#$PANTOOLS gene_classification --core-threshold 95 --unique-threshold 5 ${pan_db} 
#error_exit $?

### Gene distance tree
#Rscript ${pan_db}/gene_classification/gene_distance_tree.R
#mv ${pan_db}/gene_classification ${pan_db}/gene_classification.95.5

#######################################################################

### K-mer classification: soft core and cloud
#process_start kmer_classification_core_unique
#$PANTOOLS k_mer_classification ${pan_db}
#error_exit $?

#Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
#mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.100.0


### K-mer classification: soft core and cloud
#process_start kmer_classification_softcore_cloud
#$PANTOOLS k_mer_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
#error_exit $?

#Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
#mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.95.5

#######################################################################

### Pangenome structure: genes
#process_start pangenome_structure_gene
#$PANTOOLS pangenome_structure -t 20 ${pan_db}
#error_exit $?

#Rscript ${pan_db}/pangenome_size/gene/pangenome_growth.R
#Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_or_average.R
#Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_and_average.R
#Rscript ${pan_db}/pangenome_size/gene/heaps_law.R
#Rscript ${pan_db}/pangenome_size/gene/core_access_unique.R

#running
### Pangenome structure: kmer
#process_start pangenome_structure_kmer
#$PANTOOLS pangenome_structure -t 20 --kmer ${pan_db}
#error_exit $?

#Rscript ${pan_db}/pangenome_size/kmer/pangenome_growth.R

######################################################################

### functional_classification: core and unique
#process_start functional_classification_core_unique
#$PANTOOLS functional_classification ${pan_db}
#error_exit $?

#mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.100.0


### functional_classification: softcore and cloud
#process_start functional_classification_softcore_cloud
#$PANTOOLS functional_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
#error_exit $?

#mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.95.5


### function_overview
#process_start function_overview
#$PANTOOLS function_overview ${pan_db} 
#error_exit $?

#Rscript ${pan_db}/cog_per_class.R


### MSA for homology groups
#process_start "msa for homology groups"
#$PANTOOLS msa -t 12 --method per-group --mode nucleotide ${pan_db}
#error_exit $?

#######################################################################

### SNP tree using core gene SNPs
#process_start core_phylogeny
#cp -r ${pan_db}/gene_classification.100.0 ${pan_db}/gene_classification
#$PANTOOLS core_phylogeny -t 20  --clustering-mode ML ${pan_db}
#error_exit $?

#rm ${pan_db}//core_snp_tree/informative.fasta.*
#iqtree -nt 20 -s ${pan_db}/core_snp_tree/informative.fasta -redo -bb 1000

######################################################################

## add the updated phenotypes for association analysis
process_start add_phenotypes
$PANTOOLS remove_phenotype ${pan_db}
$PANTOOLS add_phenotypes ${pan_db} analysis/04_pangenome_pecto_v2/phylogeny/clade_compare_phenotypes.csv
error_exit $?


[ -d ${pan_db}/gene_classification ] && rm -r ${pan_db}/gene_classification
mkdir ${pan_db}/gene_classification.pheno

## Gene classification for each phenotype
phenotypes=(`awk -F "\t" '{ print $1 }' data/analysis_configs/pheno_association_config.tab`)
for phn in ${phenotypes[@]}
do
    process_start "gene_classification for phenotype $phn"
    pheno_arg=`grep "^${phn}\b" data/analysis_configs/pheno_association_config.tab | cut -f2`
    $PANTOOLS gene_classification ${pheno_arg} ${pan_db}
    error_exit $?

    ## move results to a folder
    pheno_dir=${pan_db}/gene_classification.pheno/${phn}
    [ -d ${pheno_dir} ] && rm -r ${pheno_dir}
    mkdir ${pheno_dir}
    mv ${pan_db}/gene_classification/{phenotype_*,gene_classification_phenotype_overview.txt} ${pheno_dir}/
done

rm -r ${pan_db}/gene_classification

######################################################################

### homology group information 
#process_start "extracting homology group information"
#$PANTOOLS group_info ${pan_db} ${pan_db}/gene_classification.100.0/all_homology_groups.csv
#error_exit $?

#conda activate omics_py37
#while IFS=$'\t' read -r phn groups ; do
#    # printf "%b\n" "column1<${phn}>"
#    # printf "%b\n" "column2<${groups}>"
#    process_start "extracting mRNA sequence for homology group specific to phenotype: ${phn}"
#    printf "Homology groups: ${groups}\n"

#    genome=`grep "${phn}" data/analysis_configs/pheno_association_config.tab | cut -f3 | sed 's/,.*//'`

#    file_info="analysis/04_pangenome_pecto_v2/pheno_association/phenotype_${phn}.seq_info.txt"
#    file_fasta="analysis/04_pangenome_pecto_v2/pheno_association/phenotype_${phn}.specific_mRNA.fasta"

#    printf "homology_group_id, genome, mRNA name, mRNA identifier, node identifier, address, strand\n"  > ${file_info}
#    printf "" > ${file_fasta}

#    ## read groups into array
#    IFS=',' read -ra hg_array <<< "${groups}"
#    for hg in "${hg_array[@]}"
#    do
#        # printf "${hg} "
#        seq_info=`grep "^${genome}" ${pan_db}/alignments/msa_per_group/grouping_v4/${hg}/input/sequences.info`
#        mrna_id=`echo ${seq_info} | cut -d"," -f3`
#        samtools faidx ${pan_db}/alignments/msa_per_group/grouping_v4/${hg}/input/nuc.fasta

#        printf "${hg}, ${seq_info}\n" >> ${file_info}
#        samtools faidx ${pan_db}/alignments/msa_per_group/grouping_v4/${hg}/input/nuc.fasta ${mrna_id} >> ${file_fasta}
#    done
#done < analysis/04_pangenome_pecto_v2/pheno_association/phenotype_specific_groups.txt
#conda activate pantools_master

######################################################################

# Pectobacterium pangenome

## Setup

``` bash
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

export PANTOOLS="$PANTOOLS_3_4"

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
```

## Construct pangenome

``` bash
TMPDIR="$LOCAL_DIR_PATH/tmp"
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
export PANTOOLS="$PANTOOLS_OPT"

process_start build_pangenome_parallel
$PANTOOLS build_pangenome_parallel --database-path ${pan_db} \
--genomes-file $PANGENOME_DIR/genomes_fa.list -tn 40
error_exit $?

cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_raw

export PANTOOLS="$PANTOOLS_MASTER"

## add annotations
process_start add_annotations
$PANTOOLS add_annotations --connect ${pan_db} $PANGENOME_DIR/genomes_gff3.list
error_exit $?

cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_ann

## add phenotypes
process_start add_phenotypes
$PANTOOLS remove_phenotype ${pan_db}
$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/genomes_metadata.csv
error_exit $?

```

## Add InterProScan annotations

``` bash
## add_functions
process_start add_InterProScan_annotations
$PANTOOLS add_functions ${pan_db} $PANGENOME_DIR/functional_annotations.txt
error_exit $?

cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_fn
######################################################################
```

## Grouping

### Grouping at the pangenome level

``` bash
cp -r $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_fn ${pan_db}
## BUSCO
process_start busco_protein
$PANTOOLS busco_protein -t 20 --busco10 enterobacterales_odb10 ${pan_db}
error_exit $?

cp -r ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_bsk

## BUSCO multiqc
process_start multiqc_busco
multiqc --filename multiqc_busco -outdir analysis/01_multiqc --ignore "run_enterobacterales_odb10" --module busco -d -dd 1 ${pan_db}/busco/enterobacterales_odb10/protein/results/
error_exit $?

## optimized grouping
process_start optimal_grouping
$PANTOOLS optimal_grouping -t 30 ${pan_db} ${pan_db}/busco/enterobacterales_odb10
# $PANTOOLS optimal_grouping -dp ${pan_db} -if ${pan_db}/busco/enterobacterales_odb10 -tn 20
error_exit $?

Rscript ${pan_db}/optimal_grouping/optimal_grouping.R
$PANTOOLS grouping_overview ${pan_db}

```

### `optimal_grouping` using subset of genomes or type strains

`optimal_grouping` at the Pangenome level is time consuming. An alternative approch using subset of genomes or only type strains was tried and the results were compared with the pangenome scale `optimal_grouping`.

```bash
## type strains
bash ./scripts/b_construction/b04_subset_optimal_grouping.sh "typeStrain" "374,96,385,386,379,388,116,375,377,347,250,335,338,265,256,32,269,387,266"

## generate random subsets to run `optimal_grouping`
Rscript scripts/c_analysis/c05_subset_optimal_grouping.R

## run `optimal_grouping` on these random subsets
nohup \
# sed -n '1,6!p' analysis/03_pangenome_pecto_v2/subset_optimal_group/subsets_conf.tab | \
cat analysis/03_pangenome_pecto_v2/subset_optimal_group/subsets_conf.tab | \
parallel --colsep '\t' --jobs 10  --keep-order --workdir $PWD --halt soon,fail=1 \
--load 100% --results logs/pantools/sub_opt_group/{1} \
--joblog logs/pantools/sub_opt_group/sub_opt_group.log \
bash ./scripts/b_construction/b04_subset_optimal_grouping.sh {1} {2} \
>>logs/pantools/sub_opt_group/nohup.out 2>&1 &

## summarize the results
Rscript scripts/c_analysis/c05_subset_grouping_summary.R

```

### fix a grouping setting

```bash
$PANTOOLS change_grouping -v 4 ${pan_db}
$PANTOOLS grouping_overview ${pan_db}
######################################################################
```

## Pangenome exploration

``` bash
## extract the metrics from pangenome
process_start extract_pangenome_metrics
$PANTOOLS metrics ${pan_db}
error_exit $?
```

### gene classification

#### core-unique variation w.r.t. cutoffs

```bash
## Core unique thresholds
process_start core_unique_thresholds
$PANTOOLS core_unique_thresholds ${pan_db}
error_exit $?
Rscript ${pan_db}/core_unique_thresholds/core_unique_thresholds.R
```

#### core and unique genes

```bash
## gene classification: core and unique
process_start gene_classification_core_unique
$PANTOOLS gene_classification ${pan_db} 
error_exit $?

## Gene distance tree
Rscript ${pan_db}/gene_classification/gene_distance_tree.R

mv ${pan_db}/gene_classification ${pan_db}/gene_classification.100.0
```

#### soft core (95%) and cloud (5%) genes

```bash
## gene classification: soft core and cloud
process_start gene_classification_softcore_cloud
$PANTOOLS gene_classification --core-threshold 95 --unique-threshold 5 ${pan_db} 
error_exit $?

## Gene distance tree
Rscript ${pan_db}/gene_classification/gene_distance_tree.R
mv ${pan_db}/gene_classification ${pan_db}/gene_classification.95.5
```

### Extract information for all homology groups

```bash
## homology group information 
process_start "extracting homology group information"
$PANTOOLS group_info ${pan_db} ${pan_db}/gene_classification.100.0/all_homology_groups.csv
error_exit $?
```

### k-mer classification

#### core and unique kmers

```bash
## K-mer classification: soft core and cloud
process_start kmer_classification_core_unique
$PANTOOLS k_mer_classification ${pan_db}
error_exit $?

Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.100.0
```

#### soft core (95%) and cloud (5%) kmers

```bash
## K-mer classification: soft core and cloud
process_start kmer_classification_softcore_cloud
$PANTOOLS k_mer_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
error_exit $?

Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.95.5
```

### Pangenome structure

#### Use all genomes to determine pangenome structure

```bash
## Pangenome structure: genes
process_start pangenome_structure_gene
$PANTOOLS pangenome_structure -t 20 ${pan_db}
error_exit $?

Rscript ${pan_db}/pangenome_size/gene/pangenome_growth.R
Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_or_average.R
Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_and_average.R
Rscript ${pan_db}/pangenome_size/gene/heaps_law.R

mv ${pan_db}/pangenome_size/gene ${pan_db}/pangenome_size/gene.pangenome

## Pangenome structure: kmer
process_start pangenome_structure_kmer
$PANTOOLS pangenome_structure -t 20 --kmer ${pan_db}
error_exit $?

Rscript ${pan_db}/pangenome_size/kmer/pangenome_growth.R

```

#### Use subset of genomes to determine pangenome structure

```bash
## Pangenome structure for subset of genomes
genomeSets=(`awk -F "\t" '{ print $1 }' $PANGENOME_DIR/genome_sets.tab`)

for gs in ${genomeSets[@]}
do
    process_start "gene_classification for genome subset $gs"
    str_arg=`grep "^${phn}\b" $PANGENOME_DIR/genome_sets.tab | cut -f2`
    $PANTOOLS pangenome_structure -t 20 ${str_arg} ${pan_db}
    error_exit $?

    ## move results to a folder
    gs_dir=${pan_db}/pangenome_size/gene.${gs}
    [ -d ${gs_dir} ] && rm -r ${gs_dir}
    mv ${pan_db}/pangenome_size/gene ${gs_dir}
done
```

### Functional classification

#### core and unique functional classification

```bash
## functional_classification: core and unique
process_start functional_classification_core_unique
$PANTOOLS functional_classification ${pan_db}
error_exit $?

mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.100.0
```

#### softcore (95%) and cloud (5%) functional classification

```bash
## functional_classification: softcore and cloud
process_start functional_classification_softcore_cloud
$PANTOOLS functional_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
error_exit $?

mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.95.5
```

#### function overview
```bash
## function_overview
process_start function_overview
$PANTOOLS function_overview ${pan_db} 
error_exit $?

Rscript ${pan_db}/cog_per_class.R
```

### MSA for all homology groups

```bash
## MSA for homology groups
process_start "msa for homology groups"
$PANTOOLS msa -t 12 --method per-group --mode nucleotide ${pan_db}
error_exit $?
######################################################################
```

## Phylogeny analysis

``` bash

## SNP tree using core gene SNPs
process_start core_phylogeny
cp -r ${pan_db}/gene_classification.100.0 ${pan_db}/gene_classification
$PANTOOLS core_phylogeny -t 20  --clustering-mode ML ${pan_db}
error_exit $?

rm ${pan_db}//core_snp_tree/informative.fasta.*
iqtree -nt 20 -s ${pan_db}/core_snp_tree/informative.fasta -redo -bb 1000
######################################################################
```

### Phenotype association

```bash
## add the updated phenotypes for association analysis
process_start add_phenotypes
$PANTOOLS remove_phenotype ${pan_db}
$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/genomes_metadata.csv
$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/analysis_configs/clade_association_phenotypes.csv
error_exit $?


[ -d ${pan_db}/gene_classification ] && rm -r ${pan_db}/gene_classification
mkdir ${pan_db}/gene_classification.pheno

## Gene classification for each phenotype
phenotypes=(`awk -F "\t" '{ if (NR!=1) {print $1} }' $PANGENOME_DIR/analysis_configs/pheno_association_config.tab`)
for phn in ${phenotypes[@]}
do
    process_start "gene_classification for phenotype $phn"
    pheno_arg=`grep "^${phn}\b" $PANGENOME_DIR/analysis_configs/pheno_association_config.tab | cut -f2`
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
```

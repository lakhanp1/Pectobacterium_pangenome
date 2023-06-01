# Pectobacterium pangenome

## Setup

``` bash
#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

# source scripts/utils/setup_analysis.sh 'pectobacterium.v2'
source scripts/utils/setup_analysis.sh $@

if [ -z ${pan_db+x} ];
then
    echo "\$pan_db is unset"
    error_exit 1
fi

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master
export PANTOOLS="$PANTOOLS_4_1"
#export PANTOOLS="$PANTOOLS_DEV"
######################################################################
```

## Construct pangenome

``` bash
TMPDIR="$LOCAL_DIR_PATH/tmp"
[ ! -d ${TMPDIR} ] && mkdir -p ${TMPDIR}
[ -d ${TMPDIR}/pantools ] && rm -rd ${TMPDIR}/pantools
[ ! -d ${TMPDIR}/pantools ] && mkdir -p ${TMPDIR}/pantools

process_start build_pangenome_parallel
$PANTOOLS build_pangenome --threads 40 --scratch-directory ${TMPDIR}/pantools \
--num-db-writer-threads 10 --cache-size 25000000 \
${pan_db} $PANGENOME_DIR/genomes_fa.list
error_exit $?

cp -rp ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.raw


## add annotations
export PANTOOLS="$PANTOOLS_DEV"
process_start add_annotations
$PANTOOLS add_annotations --connect ${pan_db} $PANGENOME_DIR/genomes_gff3.list
error_exit $?
export PANTOOLS="$PANTOOLS_4_1"

## add phenotypes
process_start add_phenotypes
$PANTOOLS remove_phenotype ${pan_db}
$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/genomes_metadata.csv
error_exit $?

```

## Add InterProScan and COG annotations

``` bash
## InterProScan
process_start add_InterProScan_annotations
$PANTOOLS add_functions ${pan_db} $PANGENOME_DIR/functional_annotations.txt
error_exit $?

## COG
process_start add_COG_annotations 
$PANTOOLS add_functions ${pan_db} $PANGENOME_DIR/deggnog_annotations.txt
error_exit $?

cp -rp ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.fn
######################################################################
```

## Grouping

``` bash
cp -rp $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.backup_fn ${pan_db}
## BUSCO
process_start busco_protein
$PANTOOLS busco_protein -t 20 --busco10 enterobacterales_odb10 ${pan_db}
error_exit $?
```

### `optimal_grouping` using subset of genomes or type strains

`optimal_grouping` at the Pangenome level is time consuming. An alternative approch using subset of genomes or only type strains was tried and the results were compared with the pangenome scale `optimal_grouping`.

``` bash
## generate random subsets to run `optimal_grouping`
Rscript scripts/c_analysis/grouping_subsets_build.R

```

This script will generate a two column TSV file where first column is a set identifier and second column include comma separated genome identifers.

```
typeStrain	374,96,385,386,379,388,116,375,377,347,250,335,338,265,256,32,269,387,266
rand_001	451,375,343,346,220,403,251,228,431,248,116,440,383,322,306,351,430,15,14,314
rand_002	451,375,374,224,32,360,252,226,442,248,116,266,386,385,306,387,205,15,265,92
rand_003	451,375,374,334,279,391,103,445,442,248,116,266,386,27,305,351,393,15,265,21
rand_004	451,375,374,88,96,193,251,235,431,248,116,440,386,124,335,387,294,15,423,123
.
.
.
```

An example to run optimal grouping on subset of genomes in the pangenome.

```bash
## type strains
bash ./scripts/b_construction/grouping_subsets_process.sh "typeStrain" "374,96,385,386,379,388,116,375,377,347,250,335,338,265,256,32,269,387,266"
```

Use GNU parallel to run optimal grouping on all random sets.

```bash
## run `optimal_grouping` on these random subsets
nohup \
# sed -n '1,6!p' analysis/03_pangenome_pecto_v2/subset_optimal_group/subsets_conf.tab | \
cat analysis/03_pangenome_pecto_v2/subset_optimal_group/subsets_conf.tab | \
parallel --colsep '\t' --jobs 10  --keep-order --workdir $PWD --halt soon,fail=1 \
--load 100% --results logs/pantools/sub_opt_group/{1} \
--joblog logs/pantools/sub_opt_group/sub_opt_group.log \
bash ./scripts/b_construction/grouping_subsets_process.sh {1} {2} \
>>logs/pantools/sub_opt_group/nohup.out 2>&1 &
```

Summarize the results by plotting in R

```bash
## summarize the results
Rscript scripts/c_analysis/grouping_subsets_analyze.R
```

### Grouping at the pangenome level

```bash
## grouping with relaxation 4 setting
process_start group_v4
nice $PANTOOLS group -t 30 --relaxation 4 ${pan_db}
error_exit $?

cp -rp ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.grp

## optimized grouping
process_start optimal_grouping
nice $PANTOOLS optimal_grouping -t 30 ${pan_db} ${pan_db}/busco/enterobacterales_odb10
error_exit $?


Rscript ${pan_db}/optimal_grouping/optimal_grouping.R
$PANTOOLS grouping_overview ${pan_db}

## fix a grouping setting
$PANTOOLS change_grouping -v 4 ${pan_db}
$PANTOOLS grouping_overview ${pan_db}

cp -rp ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.opt_grp
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

``` bash
## Core unique thresholds
process_start core_unique_thresholds
$PANTOOLS core_unique_thresholds ${pan_db}
error_exit $?
Rscript ${pan_db}/core_unique_thresholds/core_unique_thresholds.R
```

#### core and unique genes

``` bash
## gene classification: core and unique
process_start gene_classification_core_unique
$PANTOOLS gene_classification ${pan_db} 
error_exit $?

## Gene distance tree
Rscript ${pan_db}/gene_classification/gene_distance_tree.R

mv ${pan_db}/gene_classification ${pan_db}/gene_classification.100.0
```

#### soft core (95%) and cloud (5%) genes

``` bash
## gene classification: soft core and cloud
process_start gene_classification_softcore_cloud
$PANTOOLS gene_classification --core-threshold 95 --unique-threshold 5 ${pan_db} 
error_exit $?

## Gene distance tree
Rscript ${pan_db}/gene_classification/gene_distance_tree.R
mv ${pan_db}/gene_classification ${pan_db}/gene_classification.95.5
```

### Extract information for all homology groups

``` bash
## homology group information
sed -i.bak -r -n '/^#/! p' ${pan_db}/gene_classification.100.0/all_homology_groups.csv
process_start "extracting homology group information"
$PANTOOLS group_info -H ${pan_db}/gene_classification.100.0/all_homology_groups.csv ${pan_db}
error_exit $?
```

### k-mer classification

#### core and unique kmers

``` bash
## K-mer classification: soft core and cloud
process_start kmer_classification_core_unique
$PANTOOLS k_mer_classification ${pan_db}
error_exit $?

Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.100.0
```

#### soft core (95%) and cloud (5%) kmers

``` bash
## K-mer classification: soft core and cloud
process_start kmer_classification_softcore_cloud
$PANTOOLS k_mer_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
error_exit $?

Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.95.5
```

### Pangenome structure

#### Use all genomes to determine pangenome structure

``` bash
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
######################################################################
```

### Functional classification

#### core and unique functional classification

``` bash
## functional_classification: core and unique
process_start functional_classification_core_unique
$PANTOOLS functional_classification ${pan_db}
error_exit $?

mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.100.0
```

#### softcore (95%) and cloud (5%) functional classification

``` bash
## functional_classification: softcore and cloud
process_start functional_classification_softcore_cloud
$PANTOOLS functional_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
error_exit $?

mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.95.5
```

#### function overview

``` bash
## function_overview
process_start function_overview
$PANTOOLS function_overview ${pan_db} 
error_exit $?

Rscript ${pan_db}/cog_per_class.R
```

``` bash
## GO enrichment for core, accessory and unique homology groups
for grp in core accessory unique
do
    process_start GO_enrichment:$grp
    $PANTOOLS go_enrichment -H ${pan_db}/gene_classification.100.0/${grp}_groups.csv  ${pan_db}
    error_exit $?
    mv ${pan_db}/function/go_enrichment ${pan_db}/function/go_enrichment.100.0.${grp} 
done

## GO enrichment for soft-core, accessory and cloud homology groups
for grp in core accessory unique
do
    process_start GO_enrichment:$grp
    $PANTOOLS go_enrichment -H ${pan_db}/gene_classification.95.5/${grp}_groups.csv  ${pan_db}
    error_exit $?
    mv ${pan_db}/function/go_enrichment ${pan_db}/function/go_enrichment.95.5.${grp} 
done

#######################################################################
```

### MSA for all homology groups

``` bash
## MSA for homology groups
process_start "msa for homology groups"
$PANTOOLS msa -t 12 --method per-group --mode nucleotide ${pan_db}
error_exit $?

cp -rp ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.msa
```

## Phylogeny analysis

``` bash

## SNP tree using core gene SNPs
process_start core_phylogeny
cp -rp ${pan_db}/gene_classification.100.0 ${pan_db}/gene_classification
$PANTOOLS core_phylogeny -t 20  --clustering-mode ML ${pan_db}
error_exit $?
rm -r ${pan_db}/gene_classification

## running_job: leunissen
mkdir ${pan_db}/core_snp_tree/ML_tree 
nohup nice iqtree -T 40 -s data/pangenomes/pectobacterium.v2/backup/pectobacterium.v2.DB.msa/core_snp_tree/informative.fasta -redo -B 1000 \
--prefix data/pangenomes/pectobacterium.v2/backup/pectobacterium.v2.DB.msa/core_snp_tree/ML_tree/informative.fasta \
>> logs/v2_pecto/core_tree.log 2>&1 &

nohup nice iqtree -T 40 -s informative.fasta -B 1000 --prefix ML_tree/informative.fasta \
>> core_tree.log 2>&1 &

## run IQ-tree with a specific model
nohup iqtree -T 30 -s ${pan_db}/core_snp_tree/informative.fasta -redo -B 1000 \
-m GTR+F+ASC --prefix ${pan_db}/core_snp_tree/informative.fasta.GTR_F_ASC \
> logs/v2_pecto/iqtree_GTR_F_ASC.log 2>&1 &

######################################################################
```

## Extract pangenome data from Neo4j database

``` bash
python3 scripts/b_construction/neo4j_extract_go.py
```

# Downstream analysis of pangenome data

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
-------------------------------------------------------------
"

if [ $# -ne 1 ]; then
    printf "Error: Require pangenome name\n${usage}" >&2 ; exit 1
fi

######################################################################
source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh

PANGENOME_NAME=$1
# PANGENOME_NAME='pectobacterium.ts'
# PANGENOME_NAME='pectobacterium.v2'

conda activate pantools_master
export PANTOOLS="$PANTOOLS_4_1"

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

printf "PANGENOME_DIR: ${PANGENOME_DIR}
BUILD_DIR: ${PAN_BUILD_DIR}
pangenome: ${pan_db}
"

hg_aln_dir="${pan_db}/alignments/msa_per_group/grouping_v1"

######################################################################
```

## Analysis

### Generate various trees with default layout

```bash
## UPGMA tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./analysis/04_pangenome_pecto_v2/phylogeny/ani/ANI_UPGMA.pangenome.newick \
--name "ani_upgma" --config project_config.yaml

## Neighbor joining tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./analysis/04_pangenome_pecto_v2/phylogeny/ani/ANI_NJ.pangenome.newick \
--name "ani_nj" --config project_config.yaml

## gene distance NJ tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/gene_classification.100.0/gene_distance.tree \
--name "gene_nj" --config project_config.yaml

## kmer distance NJ tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/kmer_classification.100.0/genome_kmer_distance.tree \
--name "kmer_nj" --config project_config.yaml

## core-SNP ML tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/core_snp_tree/informative.fasta.treefile \
--name "core_snp_ml" --config project_config.yaml

## process homology groups
Rscript scripts/c_analysis/homology_group_process.R 

```

### PCR probe BLAST

#### combine multiple genome FASTA into a single file

```bash
## PCR assay blast
ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_v2/insilico_assay"
[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR

## combine genomes and prepare blastn index
python3 scripts/a_preprocessing/a10_merge_genome_faa.py --geno $PANGENOME_DIR/genomes_fa.list \
--out $PANGENOME_DIR/blastdb/genomes_combined.fa

## index the sequences
makeblastdb -dbtype nucl -in $PANGENOME_DIR/blastdb/genomes_combined.fa -parse_seqids \
-title $PANGENOME_NAME -logfile $PANGENOME_DIR/blastdb/makeblastdb.log

```

#### run BLASTN

**NOTE:** 
`blastn` was not showing hits of some primer sequences when searched against all the genomes. Hence, database sequences are split into batches of 100 genomes and `blastdb_aliastool` is run to preprocess the sequence ids. These preprocessed sequence ids are provided as input to `-seqidlist ` argument while running `blastn`. Somehow, when the search space is limited using `-seqidlist` option, missing hits are shown. This is run for each set of sequence ids.

```bash

for chrset in $PANGENOME_DIR/blastdb/chrset_*.acc
do
    ##
    blastdb_aliastool -seqid_file_in ${chrset}

    setName=`basename ${chrset} | sed 's/.acc//'`

    ## run blastn in blastn-short mode
    blastn -db $PANGENOME_DIR/blastdb/genomes_combined.fa \
    -query data/other/pecto_assay/pectobacterium_assay_probes.fasta \
    -task blastn-short -num_threads 12 -max_target_seqs 5000 -seqidlist ${chrset}.bsl \
    -outfmt "6 qseqid qstart qend qlen sseqid sstart send sstrand slen pident length mismatch qcovs gapopen evalue bitscore" \
    -out $ANALYSIS_DIR/${blastnOut}

done

```

#### parse & summarize blastn result table

```bash
Rscript scripts/c_analysis/c02_blastn_parse.R
######################################################################
```

### GO enrichment for homology groups of interest

```bash
## go_enrichment for assay_FN enriched homology groups
process_start "go_enrichment for assay_FN enriched homology groups"
$PANTOOLS go_enrichment -H analysis/04_pangenome_pecto_v2/pheno_association/specific_hgs.assay_FN.txt \
--include=429,439,369,149,29,97,155,366,373,178,181,159,345,371,180,316,414,360,166,243,147,152,173,170,416,433,157,417,191,390,136,419,142,410,146,317,145,194,426,240,340,367,357,364,359,372,358,342,370,196,24,353,52,134,179,187,188,195,192,401,402,413,154,214,153,144,165,176,140,168,156,13,148,163,164,162,418,172,297,302,63,190,415,169,171,167,174,189,193,411,397,398,405,409,412,403,408,399,404,407,175,400,406,158,161,138,60,337,242,368,74,427,308,438,299,391,182,185,236,177,42,43,263,307,379,356,380,141,341,64,352,111,115,114,108,109,99,137 ${pan_db} 
error_exit $?
######################################################################
```


### Use subset of genomes to determine pangenome structure

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
######################################################################
```

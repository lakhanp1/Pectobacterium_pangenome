# Downstream analysis of pangenome data

## setup

```bash
#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master

# PANGENOME_NAME=$1
# PANGENOME_NAME='pectobacterium.10g'
PANGENOME_NAME='pectobacterium.v2'

## Setup
PROJECT_DIR="/lustre/BIF/nobackup/$USER/projects/03_Pectobacterium"
PANGENOME_DIR="$PROJECT_DIR/data/pangenomes/$PANGENOME_NAME"
######################################################################
```

## Phylogeny

### Generate various trees with default layout

```bash
## UPGMA tree
./scripts/c_analysis/c02_plot_evolutionary_trees.R \
--tree ./analysis/04_pangenome_pecto_v2/phylogeny/pangenome_ANI_UPGMA.newick \
--name "ANI UPGMA" --out ANI_UPGMA_rooted --config project_config.yaml

## Neighbor joining tree
./scripts/c_analysis/c02_plot_evolutionary_trees.R \
--tree ./analysis/04_pangenome_pecto_v2/phylogeny/pangenome_ANI_NJ.newick \
--name "Neighbor Joining" --out ANI_NJ_rooted --config project_config.yaml

## gene distance NJ tree
./scripts/c_analysis/c02_plot_evolutionary_trees.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/gene_classification.100.0/gene_distance.tree \
--name "gene distance NJ" --out gene_NJ_rooted --config project_config.yaml

## kmer distance NJ tree
./scripts/c_analysis/c02_plot_evolutionary_trees.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/kmer_classification.100.0/genome_kmer_distance.tree \
--name "kmer distance NJ" --out kmer_NJ_rooted --config project_config.yaml


```

## PCR probe BLAST

#### combine multiple genome FASTA into a single file

```bash
ANALYSIS_DIR="$PROJECT_DIR/analysis/03_pangenome_pecto_v2/insilico_assay"
[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR

## combine genomes and prepare blastn index
python3 scripts/a_preprocessing/a10_merge_genome_faa.py --geno $PANGENOME_DIR/genomes_fa.list \
--out $PANGENOME_DIR/blastdb/genomes_combined.fa

## index the sequences
makeblastdb -dbtype nucl -in $PANGENOME_DIR/blastdb/genomes_combined.fa \
-title $PANGENOME_NAME -logfile $PANGENOME_DIR/blastdb/makeblastdb.log

```

#### run BLASTN

```bash
## run blastn in blastn-short mode
blastn -db $PANGENOME_DIR/blastdb/genomes_combined.fa \
-query data/other/pecto_assay/pectobacterium_assay_probes.fasta \
-task blastn-short -max_target_seqs 2000 -num_threads 12 \
-outfmt "6 qseqid qstart qend qlen sseqid sstart send sstrand slen pident length mismatch qcovs gapopen evalue bitscore" \
-out $ANALYSIS_DIR/probe_blastn.out

```

#### parse & summarize blastn result table

```bash
Rscript scripts/c_analysis/c02_blastn_parse.R

######################################################################
```

### UML 

![Class Diagram](http://www.plantuml.com/plantuml/proxy?src=https://raw.githubusercontent.com/Zingam/Markdown-Document-UML-Use-Test/master/UML/Instance.puml)


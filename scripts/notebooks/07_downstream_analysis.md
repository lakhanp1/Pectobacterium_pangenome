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
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./analysis/04_pangenome_pecto_v2/phylogeny/ANI_UPGMA.pangenome.newick \
--name "ani_upgma" --config project_config.yaml

## Neighbor joining tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./analysis/04_pangenome_pecto_v2/phylogeny/ANI_NJ.pangenome.newick \
--name "ani_nj" --config project_config.yaml

## gene distance NJ tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/gene_classification.100.0/gene_distance.tree \
--name "gene_nj" --config project_config.yaml

## kmer distance NJ tree
./scripts/c_analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/kmer_classification.100.0/genome_kmer_distance.tree \
--name "kmer_nj" --config project_config.yaml


```

## PCR probe BLAST

#### combine multiple genome FASTA into a single file

```bash
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

`blastn` was not showing hits of some primer sequences when searched against all the genomes.
Hence, database sequences are split into batches of 100 genomes and `blastdb_aliastool` is 
run to preprocess the sequence ids. These preprocessed sequence ids are provided
as input to `-seqidlist ` argument while running `blastn`. Somehow, when the search
space is limited using `-seqidlist` option, missing hits are shown. This is run for
each set of sequence ids.

```bash

for chrset in $PANGENOME_DIR/blastdb/chrset_*.acc
do
    ##
    blastdb_aliastool -seqid_file_in ${chrset}

    setName=`basename ${chrset} | sed 's/.acc//'`
    blastnOut="probe_blastn.${setName}.out"

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

### UML 

![Class Diagram](http://www.plantuml.com/plantuml/proxy?src=https://raw.githubusercontent.com/Zingam/Markdown-Document-UML-Use-Test/master/UML/Instance.puml)


# Downstream analysis of pangenome data

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

######################################################################
```

## Use subset of genomes to determine pangenome structure

``` bash
## Pangenome structure for species with more than 20 genomes
for sp in `awk -F "\t" '{ if (NR != 1 && $2 >= 20) {print $1} }' $PANGENOME_DIR/analysis_configs/species_wise_genomes.tab`
do
    genomes=`grep "^${sp}\b" $PANGENOME_DIR/analysis_configs/species_wise_genomes.tab | cut -f3`
    
    process_start "pangenome_structure for genome of species ${sp}: $genomes"
    $PANTOOLS pangenome_structure -t 20 --include ${genomes} ${pan_db}
    error_exit $?

    Rscript ${pan_db}/pangenome_size/gene/pangenome_growth.R
    Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_or_average.R
    Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_and_average.R
    # Rscript ${pan_db}/pangenome_size/gene/heaps_law.R

    ## move results to a folder
    psDir=${pan_db}/pangenome_size/gene.${sp}
    [ -d ${psDir} ] && rm -r ${psDir}
    mv ${pan_db}/pangenome_size/gene ${psDir}
done

## update the path in heaps_law.R script
cd ${pan_db}/pangenome_size
for i in gene.*
do
    sed -i.bak "s/pangenome_size\/gene/pangenome_size\/gene.$i/g" $i/heaps_law.R
done
cd $PROJECT_DIR

```

## Visualize pangenome metadata

```bash
quarto render scripts/analysis/pangenome_metadata.viz.qmd
```

## Phylogeny analysis 

``` bash
## UPGMA tree
./scripts/analysis/phylo.process_newicks.R \
--tree ./analysis/pangenome_v2/phylogeny/ani/ANI_UPGMA.pangenome.newick \
--name "ani_upgma" --config project_config.yaml --save_leaf_order

## Neighbor joining tree
./scripts/analysis/phylo.process_newicks.R \
--tree ./analysis/pangenome_v2/phylogeny/ani/ANI_NJ.pangenome.newick \
--name "ani_nj" --config project_config.yaml --save_leaf_order

## gene distance NJ tree
./scripts/analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/gene_classification.100.0/gene_distance.tree \
--name "gene_nj" --config project_config.yaml 

## kmer distance UPGMA tree
./scripts/analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/kmer_classification.100.0/kmer.upgma.newick \
--name "kmer_upgma" --config project_config.yaml --save_leaf_order

## kmer distance NJ tree
./scripts/analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/kmer_classification.100.0/kmer.nj.newick \
--name "kmer_nj" --config project_config.yaml --save_leaf_order

## core-SNP ML tree
./scripts/analysis/phylo.process_newicks.R \
--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/core_snp_tree/informative.fasta.treefile \
--name "core_snp_ml" --config project_config.yaml --save_leaf_order --bootstrap

## process homology groups
Rscript scripts/analysis/homology_group_process.R 

```

Visualize core-SNP phylogenetic tree with metadata.

```bash
quarto render scripts/analysis/phylo_viz_tree.qmd
```

## Homology groups analysis

### GO enrichment of homology groups categories

This script will use the pangenome `org.db` object (built [here](pangenome_construction.md#build-pangenome-org-db-object))
and perform GO enrichment using `BioConductor::topGO` package.

```bash
Rscript scripts/analysis/HG_GO_enrichment.R
```

### Summarize homology groups

```bash
quarto render scripts/analysis/HG_summary.qmd 
```

### Visualize tandemly mapped homology groups sets on the pangenome

```bash
quarto render scripts/analysis/HG_tandem_match_viz.qmd --execute-dir ./
```

### Visualize homology groups across pagenome for GO terms of interest

to-do: convert this to script with arguments

```bash
quarto render scripts/analysis/homology_groups_go_viz.qmd --execute-dir ./
```

#### viral DNA integration related biological processes

```r
analysisName <- "virulence_GO"

goIds <- c(
  "GO:0046718", "GO:0044826", "GO:0019069", "GO:0044659", "GO:0075713",
  "GO:0044660", "GO:0019076", "GO:0019058"
)
```

#### secretion systems in *Pectobacterium* spp.

```r
analysisName <- "secretion_sys"

goIds <- c(
  "GO:0030253", "GO:0015628", "GO:0030254", "GO:0030255", "GO:0046819",
  "GO:0033103", "GO:0044315"
)
```

#### Homology groups for `lipopolysaccharide biosynthetic process`

```r
analysisName <- "lps_syn"

goIds <- c("GO:0009103", "GO:0009244")
```


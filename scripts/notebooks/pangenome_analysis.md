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

## Analysis

### Generate various trees with default layout

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

### GO enrichment for homology groups of interest

``` bash
## go_enrichment for assay_FN enriched homology groups
process_start "go_enrichment for assay_FN enriched homology groups"
$PANTOOLS go_enrichment -H analysis/04_pangenome_pecto_v2/pheno_association/specific_hgs.assay_FN.txt \
--include=429,439,369,149,29,97,155,366,373,178,181,159,345,371,180,316,414,360,166,243,147,152,173,170,416,433,157,417,191,390,136,419,142,410,146,317,145,194,426,240,340,367,357,364,359,372,358,342,370,196,24,353,52,134,179,187,188,195,192,401,402,413,154,214,153,144,165,176,140,168,156,13,148,163,164,162,418,172,297,302,63,190,415,169,171,167,174,189,193,411,397,398,405,409,412,403,408,399,404,407,175,400,406,158,161,138,60,337,242,368,74,427,308,438,299,391,182,185,236,177,42,43,263,307,379,356,380,141,341,64,352,111,115,114,108,109,99,137 ${pan_db} 
error_exit $?
######################################################################
```

### Use subset of genomes to determine pangenome structure

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
######################################################################
```

## Visualize homology groups of interest

### viral DNA integration related biological processes

```r
analysisName <- "virulence_GO"

goIds <- c(
  "GO:0046718", "GO:0044826", "GO:0019069", "GO:0044659", "GO:0075713",
  "GO:0044660", "GO:0019076", "GO:0019058"
)
```

### secretion systems in *Pectobacterium* spp.

```r
analysisName <- "secretion_sys"

goIds <- c(
  "GO:0030253", "GO:0015628", "GO:0030254", "GO:0030255", "GO:0046819",
  "GO:0033103", "GO:0044315"
)
```


### Prophage region JGAR-100719-1_contig_3_consensus:442,933-461,352

```r
analysisName <- "JGAR-100719-1_contig_3_consensus_442933_461352"

hgs <- c(
  "22426826", "22426802", "22426819", "22427604", "22427607", "22427609", "22427610",
  "22427612", "22427614", "22427616", "22427618", "22427622", "22427623", "22427625",
  "22427626", "22427627", "22427629", "22427634", "22427630", "22427633", "22427636",
  "22427640"
)
```


### Homology groups for `lipopolysaccharide biosynthetic process`

```r
analysisName <- "lps_syn"

goIds <- c("GO:0009103", "GO:0009244")
```
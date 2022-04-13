# Build Pectobacterium pangenome

## Setup

``` bash
shopt -s expand_aliases
source ~/.bash_aliases

set -e
set -u
set -o pipefail

conda activate pantools

## Setup
PROJECT_DIR='/mnt/scratch/parde001/projects/03_Pectobacterium'
ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_50g"
analysis_prefix='pectobacterium.50g'
pan_db="$ANALYSIS_DIR/${analysis_prefix}.DB"


count=1
rm $ANALYSIS_DIR/genomes_gff3.list $ANALYSIS_DIR/genomes_fa.list
touch $ANALYSIS_DIR/genomes_gff3.list $ANALYSIS_DIR/genomes_fa.list

for i in `cat $ANALYSIS_DIR/genomes.list`
do
ls $PROJECT_DIR/analysis/01_prokka_annotation/${i}/*.{fna,gff3}
printf "$PROJECT_DIR/analysis/01_prokka_annotation/${i}/${i}.fna\n" >> $ANALYSIS_DIR/genomes_fa.list
printf "$count\t$PROJECT_DIR/analysis/01_prokka_annotation/${i}/${i}.gff3\n" >> $ANALYSIS_DIR/genomes_gff3.list
((count=count+1))
done

```

## Construct pangenome

``` bash
## Construct pangenome
process_start build_pangenome
pantools build_pangenome --database-path ${pan_db} --genomes-file $ANALYSIS_DIR/genomes_fa.list
error_exit $?

## add annotations
process_start add_annotations
pantools add_annotations --database-path ${pan_db} -af $ANALYSIS_DIR/genomes_gff3.list -ca
error_exit $?

## add phenotypes
process_start add_phenotypes
pantools add_phenotypes --database-path ${pan_db} -ph $ANALYSIS_DIR/phenotypes.txt
error_exit $?

```

## Build InterProScan functional annotations

Once the protein sequences are generated in the `${pan_db}/proteins/` directory, run [InterProScan annotation pipeline](01_interproscan.md) independently. Once InterProScan annotation is generated (as GFF files), add these annotations to the pangenome using `add_function` command.

## Add InterProScan annotations

``` bash
## add_functions
process_start add_InterProScan_annotations
pantools add_functions  -tn 20 -dp ${pan_db}  -if $ANALYSIS_DIR/functional_annotations.txt
error_exit $?
```

## Grouping

``` bash
## BUSCO
process_start busco_protein
pantools busco_protein -dp ${pan_db} -if enterobacterales_odb10 -tn 20 
error_exit $?

## BUSCO multiqc
process_start multiqc_busco
multiqc --filename multiqc_busco -outdir analysis/01_multiqc --ignore "run_enterobacterales_odb10" --module busco -d -dd 1 ${pan_db}/busco/enterobacterales_odb10/protein/results/
error_exit $?

## optimized grouping
process_start optimal_grouping
pantools optimal_grouping -dp ${pan_db} -if ${pan_db}/busco/enterobacterales_odb10 -tn 20
error_exit $?

Rscript ${pan_db}/optimal_grouping/optimal_grouping.R

pantools grouping_overview -dp ${pan_db}

pantools change_grouping -dp ${pan_db} -ref 4

pantools grouping_overview -dp ${pan_db} --mode fast

```

## Pangenome exploration

``` bash
## extract the metrics from pangenome
process_start extract_pangenome_metrics
pantools metrics -dp ${pan_db}
error_exit $?

## Core unique thresholds
process_start core_unique_thresholds
pantools core_unique_thresholds -dp ${pan_db} -tn 20
error_exit $?
Rscript ${pan_db}/core_unique_thresholds/core_unique_thresholds.R

## Gene classification
process_start gene_classification
pantools gene_classification -dp ${pan_db} --core-threshold 95
error_exit $?

## Gene classification for each phenotype
phenotypes=("species" "virulance")
for phn in ${phenotypes[@]}
do
    process_start "gene_classification for phenotype $phn"
    [ -d ${pan_db}/gene_classification/${phn} ] && rm -r ${pan_db}/gene_classification/${phn}
    mkdir ${pan_db}/gene_classification/${phn}
    pantools gene_classification -dp ${pan_db} --core-threshold 95 --phenotype ${phn}
    mv ${pan_db}/gene_classification/{phenotype_*,gene_classification_phenotype_overview.txt} ${pan_db}/gene_classification/${phn}/
    error_exit $?
done


## Pangenome structure
process_start pangenome_size_genes
pantools pangenome_size_genes -dp ${pan_db} -tn 20
error_exit $?

Rscript ${pan_db}/pangenome_size/gene/pangenome_growth.R
Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_or_average.R
Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_and_average.R
Rscript ${pan_db}/pangenome_size/gene/heaps_law.R

# ## Pangenome size k-mers
# process_start pangenome_size_kmer
# pantools pangenome_size_kmer -dp ${pan_db} -tn 20
# error_exit $?

# Rscript ${pan_db}/pangenome_size/kmer/core_accessory_size.R
# Rscript ${pan_db}/pangenome_size/kmer/core_accessory_unique_size.R 

# ## K-mer classification
# pantools kmer_classification -dp ${pan_db} -tn 20
# error_exit $?

# Rscript ${pan_db}/kmer_classification/kmer_distance_tree.R

## functional_classification
process_start functional_classification
pantools functional_classification -dp ${pan_db} --core-threshold 95 --phenotype species
error_exit $?

## function_overview
process_start function_overview
pantools function_overview -dp ${pan_db} 
error_exit $?


## group_info: P_odoriferum specific: FAILED
process_start group_info
pantools group_info -dp ${pan_db}  -hm ${ANALYSIS_DIR}/P_odoriferum.specific_group.csv
error_exit $?

## GO enrichment of P_odoriferum specific: FAILED
process_start go_enrichment
pantools go_enrichment -dp ${pan_db}  -hm ${ANALYSIS_DIR}/P_odoriferum.specific_group.csv
error_exit $?


```

## Phylogeny analysis

``` bash
## Gene distance tree
Rscript ${pan_db}/gene_classification/gene_distance_tree.R 

## SNP tree using core gene SNPs
process_start core_snp_tree
pantools core_snp_tree -dp ${pan_db} --mode ML -tn 20 --phenotype virulance
# pantools core_snp_tree -dp ${pan_db} --mode ML -tn 20 -hm $ANALYSIS_DIR/sco_groups.csv --phenotype species
error_exit $?


## IQtree variable
process_start iqtree_variable
rm ${pan_db}/core_snp_tree/variable.fasta.*
iqtree -nt 20 -s ${pan_db}/core_snp_tree/variable.fasta -redo -bb 1000
error_exit $?

## IQtree informative
process_start iqtree_informative
rm ${pan_db}/core_snp_tree/informative.fasta.*
iqtree -nt 20 -s ${pan_db}/core_snp_tree/informative.fasta -redo -bb 1000
error_exit $?

## ANI tree
process_start ANI_tree
pantools ani -dp ${pan_db} --phenotype species --mode fastani -tn 20 
error_exit $?

## generate tree from ANI distances
Rscript ${pan_db}/ANI/fastANI/ANI_tree.R

## Rename trees
pantools rename_phylogeny -dp ${pan_db} --phenotype strain -if ${pan_db}/gene_classification/gene_distance.tree
pantools rename_phylogeny -dp ${pan_db} --phenotype strain -if ${pan_db}/core_snp_tree/variable.fasta.treefile
pantools rename_phylogeny -dp ${pan_db} --phenotype strain -if ${pan_db}/core_snp_tree/informative.fasta.treefile
pantools rename_phylogeny -dp ${pan_db} --phenotype strain -if ${pan_db}/ANI/fastANI/ANI.newick


## Create iTOL templates
# pantools create_tree_template -dp ${pan_db}
pantools create_tree_template -dp ${pan_db} -ph strain

```

## Explore specific homology groups

``` bash
## go_enrichment for virulent enriched homology groups
process_start "go_enrichment for virulent enriched homology groups"
pantools go_enrichment -dp ${pan_db} -hm ${ANALYSIS_DIR}/virulent_enriched_hm.txt
error_exit $?
```

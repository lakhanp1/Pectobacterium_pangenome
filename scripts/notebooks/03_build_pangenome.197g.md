# Build Pectobacterium pangenome

## Setup
``` bash
shopt -s expand_aliases
source ~/.bash_aliases

source /mnt/scratch/${USER}/tools/miniconda3/etc/profile.d/conda.sh
conda deactivate & conda activate pantools

## Setup
PROJECT_DIR=/mnt/scratch/${USER}/projects/03_Pectobacterium/
ANALYSIS_DIR="$PROJECT_DIR/analysis/02_pangenome_pecto_197g"
analysis_prefix='pectobacterium.197g'
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

### add phenotypes
#process_start add_phenotypes
#pantools add_phenotypes --database-path ${pan_db} -ph $ANALYSIS_DIR/phenotypes.txt
#error_exit $?

```

## Grouping
``` bash
## BUSCO
process_start busco_protein
pantools busco_protein -dp ${pan_db} -if enterobacterales_odb10 -tn 20 
error_exit $?


## optimized grouping
process_start optimal_grouping
pantools optimal_grouping -dp ${pan_db} -if ${pan_db}/busco/enterobacterales_odb10 -tn 20
error_exit $?

Rscript ${pan_db}/optimal_grouping/optimal_grouping.R

pantools grouping_overview -dp ${pan_db} --mode fast

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
pantools gene_classification -dp ${pan_db} --core-threshold 95 --phenotype species
error_exit $?


## Pangenome structure
process_start pangenome_size_genes
pantools pangenome_size_genes -dp ${pan_db} -tn 20
error_exit $?

# Rscript ${pan_db}/pangenome_size/gene/core_accessory_unique_size.R 
# Rscript ${pan_db}/pangenome_size/gene/core_accessory_size.R 
# Rscript ${pan_db}/pangenome_size/gene/gains_losses.R
Rscript ${pan_db}/pangenome_size/gene/pangenome_growth.R
Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_or_average.R
Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_and_average.R
Rscript ${pan_db}/pangenome_size/gene/heaps_law.R

## Pangenome size k-mers
process_start pangenome_size_kmer
pantools pangenome_size_kmer -dp ${pan_db} -tn 20
error_exit $?

Rscript ${pan_db}/pangenome_size/kmer/core_accessory_size.R
Rscript ${pan_db}/pangenome_size/kmer/core_accessory_unique_size.R 

## K-mer classification
pantools kmer_classification -dp ${pan_db} -tn 20
error_exit $?

Rscript ${pan_db}/kmer_classification/kmer_distance_tree.R

## 
```

``` bash
## SNP tree using core gene SNPs
process_start core_snp_tree
pantools core_snp_tree -dp ${pan_db} --mode ML -tn 20 -hm ${pan_db}/gene_classification/single_copy_orthologs.csv --phenotype species
error_exit $?

## IQtree variable
process_start iqtree_variable
rm /mnt/scratch/${USER}/projects/03_Pectobacterium/analysis/04_pangenome3/pectobacterium.3.DB/alignments/grouping_v4/core_snp_tree/variable.fasta.*
iqtree -nt 20 -s /mnt/scratch/${USER}/projects/03_Pectobacterium/analysis/04_pangenome3/pectobacterium.3.DB/alignments/grouping_v4/core_snp_tree/variable.fasta -redo -bb 1000
error_exit $?

## IQtree informative
process_start iqtree_informative
rm /mnt/scratch/${USER}/projects/03_Pectobacterium//analysis/04_pangenome3/pectobacterium.3.DB/alignments/grouping_v4/core_snp_tree/informative.fasta.*
iqtree -nt 20 -s /mnt/scratch/${USER}/projects/03_Pectobacterium//analysis/04_pangenome3/pectobacterium.3.DB/alignments/grouping_v4/core_snp_tree/informative.fasta -redo -bb 1000
error_exit $?

## BUSCO multiqc
process_start multiqc_busco
multiqc --filename multiqc_busco -outdir analysis/01_multiqc --ignore "run_enterobacterales_odb10" --module busco -d -dd 1 ${pan_db}/busco/enterobacterales_odb10/protein/results/
error_exit $?

Rscript ${pan_db}/gene_classification/gene_distance_tree.R 

pantools rename_phylogeny -dp ${pan_db} --phenotype strain_name -if ${pan_db}/alignments/grouping_v1/core_snp_tree/informative.fasta.treefile

cat /mnt/scratch/${USER}/projects/01_pantools_explore/data/pecto_dickeya_DB/alignments/grouping_v1/core_snp_tree/informative.fasta.treefile_RENAMED

```
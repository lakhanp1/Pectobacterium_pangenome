# Build Pectobacterium pangenome

## Setup
``` bash
shopt -s expand_aliases
source ~/.bash_aliases

source /mnt/scratch/${USER}/tools/miniconda3/etc/profile.d/conda.sh
conda deactivate & conda activate pantools

## Setup
PROJECT_DIR=/mnt/scratch/${USER}/projects/03_Pectobacterium/
ANALYSIS_DIR="$PROJECT_DIR/analysis/03_pangenome_pecto_6g"
analysis_prefix='pectobacterium.6g'
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

## add homology group
process_start add_homology_group
pantools group --database-path ${pan_db} -tn 20
error_exit $?

## add phenotypes
#pantools add_phenotypes --database-path ${pan_db} -ph $ANALYSIS_DIR/phenotypes.txt
#error_exit $?

```


## Pangenome exploration
``` bash
## extract the metrics from pangenome
pantools metrics -dp ${pan_db}
error_exit $?

## Gene classification
pantools gene_classification -dp ${pan_db} --core-threshold 95
error_exit $?

## Core unique thresholds
pantools core_unique_thresholds -dp ${pan_db}
error_exit $?
Rscript ${pan_db}/core_unique_thresholds/core_unique_thresholds.R

## Grouping overview
pantools grouping_overview -dp ${pan_db} 
error_exit $?

## Pangenome structure
pantools pangenome_size_genes -dp ${pan_db} -tn 20
error_exit $?

Rscript ${pan_db}/pangenome_size/gene/core_accessory_unique_size.R 
Rscript ${pan_db}/pangenome_size/gene/core_accessory_size.R 
Rscript ${pan_db}/pangenome_size/gene/gains_losses.R
Rscript ${pan_db}/pangenome_size/gene/heaps_law.R

## Pangenome size k-mers
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

## add phenotypes to the pangenome
pantools add_phenotype -dp ${pan_db} -ph ${pan_db}/phenotypes.txt

pantools gene_classification -dp ${pan_db} -ph species

pantools core_snp_tree -dp ${pan_db} --mode ML -tn 4 -hm sco_groups.txt

iqtree -nt 4 -s ${pan_db}/alignments/grouping_v1/core_snp_tree/informative.fasta -redo -bb 1000

Rscript ${pan_db}/gene_classification/gene_distance_tree.R 

pantools rename_phylogeny -dp ${pan_db} --phenotype strain_name -if ${pan_db}/alignments/grouping_v1/core_snp_tree/informative.fasta.treefile

cat /mnt/scratch/${USER}/projects/01_pantools_explore/data/pecto_dickeya_DB/alignments/grouping_v1/core_snp_tree/informative.fasta.treefile_RENAMED

```
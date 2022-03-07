# Annotate all the genomes using prokka

## Setup

``` bash
#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source /mnt/scratch/parde001/tools/miniconda3/etc/profile.d/conda.sh
conda deactivate & conda activate prokka

## Setup
PROJECT_DIR='/mnt/scratch/parde001/projects/03_Pectobacterium/'
ANALYSIS_DIR="$PROJECT_DIR/analysis/01_prokka_annotation"

```

## Annotate each genome using prokka

``` bash
for fa in /mnt/scratch/parde001/projects/03_Pectobacterium/data/fasta/* /mnt/scratch/parde001/projects/03_Pectobacterium/data/new_genomes/*.fa
do
prefix=`basename ${fa} | sed -r 's/\.(fna|fa)$//'`
echo "Annotating $prefix genome"
prokka --cpus 8 --outdir ${ANALYSIS_DIR}/${prefix} --prefix ${prefix} ${fa}
error_exit $?

# echo "prokka --cpus 8 --outdir ${ANALYSIS_DIR}/${prefix} --prefix ${prefix} ${fa}
# error_exit \$?
# "

done

```

## Remove the FASTA sequence at the end of prokka gff files

``` bash
for i in `cat analysis/03_pangenome/genomes.list`
do
sed -n '1,/##FASTA/ {/##FASTA/!p}' analysis/01_prokka_annotation/${i}/${i}.gff > analysis/01_prokka_annotation/${i}/${i}.gff3
done

```

--------------

## Alternatively, generate commands

``` bash
for fa in /mnt/scratch/parde001/projects/03_Pectobacterium/data/fasta/* /mnt/scratch/parde001/projects/03_Pectobacterium/data/new_genomes/*.fa
do
prefix=`basename ${fa} | sed -r 's/\.(fna|fa)$//'`
#echo "Annotating $prefix genome"
echo "prokka --cpus 8 --outdir ${ANALYSIS_DIR}/${prefix} --prefix ${prefix} ${fa}
error_exit \$?
"
done > scripts/prokka_commands.txt

```

## Split the commands into 5 batches and submit the jobs

``` bash
split --additional-suffix ".sh" -d -l 132 scripts/prokka_commands.txt scripts/prokka_fa_set

nohup bash ./scripts/01_prokka_batch.sh ./scripts/prokka_fa_set00.sh >logs/prokka_set00.log 2>&1 &
nohup bash ./scripts/01_prokka_batch.sh ./scripts/prokka_fa_set01.sh >logs/prokka_set01.log 2>&1 &
nohup bash ./scripts/01_prokka_batch.sh ./scripts/prokka_fa_set02.sh >logs/prokka_set02.log 2>&1 &
nohup bash ./scripts/01_prokka_batch.sh ./scripts/prokka_fa_set03.sh >logs/prokka_set03.log 2>&1 &
nohup bash ./scripts/01_prokka_batch.sh ./scripts/prokka_fa_set04.sh >logs/prokka_set04.log 2>&1 &

```

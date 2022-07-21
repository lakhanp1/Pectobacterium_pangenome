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
conda activate prokka

## Setup
PROJECT_DIR='/mnt/scratch/parde001/projects/03_Pectobacterium'
ANALYSIS_DIR="$PROJECT_DIR/data/prokka_annotation"

```

## Annotate each genome using prokka

``` bash

function run_prokka(){
    prefix=`basename -z ${1} | sed -r 's/(.*)\..*$/\1/' | sed -r 's/_genomic//'`

    if [ ! -f ${2}/${prefix}/${prefix}.log ] || ! grep -q 'Annotation finished successfully.' ${2}/${prefix}/${prefix}.log
    then
        echo "Annotating $prefix genome"
        prokka --cpus 8 --outdir ${2}/${prefix} --prefix ${prefix} ${1}
        error_exit $?
    else
        echo "prokka results exists"
    fi
}

export -f run_prokka

# run_prokka $PROJECT_DIR/data/genomes/ncbi/GCF_000011605.1_ASM1160v1_genomic.fna ${ANALYSIS_DIR}

## serial run
for fa in $PROJECT_DIR/data/genomes/?(ncbi|local)/*.?(fa|fna|fasta)
do
run_prokka ${fa} ${ANALYSIS_DIR}
done

# ## Run Prokka on all files using GNU parallel
# ls $PROJECT_DIR/data/genomes/?(ncbi|local)/*.?(fa|fna|fasta) | \
# parallel --gnu --keep-order --jobs 4 --halt now,fail=1 --results $PROJECT_DIR/logs/prokka/{/.} --joblog $PROJECT_DIR/logs/prokka/parallel.log run_prokka {} ${ANALYSIS_DIR}



```

## Remove the FASTA sequence at the end of prokka gff files

``` bash
for i in `cat analysis/03_pangenome/genomes.list`
do
sed -n '1,/##FASTA/ {/##FASTA/!p}' analysis/01_prokka_annotation/${i}/${i}.gff > analysis/01_prokka_annotation/${i}/${i}.gff3
done

```

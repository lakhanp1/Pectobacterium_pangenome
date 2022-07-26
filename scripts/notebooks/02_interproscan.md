# Run InterProScan on the protein sequences

## Setup

``` bash
#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
# Task: Build functional annotation using InterProScan

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate interproscan

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/interproscan"

```

## Generate functional annotation

``` bash

file_faa=$1

## replace * in the FASTA file
sed -i 's/*//g' ${file_faa}
# sed '/>/{H;$!d} ; x ; s/^/\nSTART-->/ ; s/$/\n<--END/'
prefix=`basename -z ${file_faa} | sed -r 's/(.*)\..*$/\1/'`.interProScan

if [ ! -f $ANALYSIS_DIR/${prefix}.gff3 ]
then
    process_start "InterProScan on file $prefix"
    interproscan.sh -cpu 8 -goterms -f GFF3 -iprlookup  --output-file-base $ANALYSIS_DIR/${prefix} -i $file_faa
    error_exit $?
else
    echo "${prefix}.gff3 InterProScan results exists"
fi


```

## execute the script using GNU `parallel`

``` bash
## GNU parallel on single server
nohup \
cat scripts/sub_commands/interproscan_batch.00 | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{/.} --joblog logs/interproscan/parallel_batch00.log \
$PWD/scripts/a_processing/a02_interproscan.sh {} \
>logs/interproscan/nohup_batch00.out 2>&1 &


## GNU parallel: submit jobs to another server
head -n 10 scripts/sub_commands/interproscan_batch.00 | \
env_parallel --jobs 4 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{/.} --joblog logs/interproscan/parallel.log \
--sshlogin 4/waterman.bioinformatics.nl --cleanup \
--env error_exit --env process_start --env TOOLS_PATH --env LUSTRE_HOME \
$PWD/scripts/a_processing/a02_interproscan.sh {}

```

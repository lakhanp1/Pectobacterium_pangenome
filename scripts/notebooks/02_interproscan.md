# Run InterProScan on the protein sequences

## Setup

``` bash
#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
# Task: Build functional annotation using InterProScan

shopt -s expand_aliases
source ~/.bash_aliases

source /mnt/scratch/parde001/tools/miniconda3/etc/profile.d/conda.sh
conda activate interproscan

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/interproscan"

```

## Generate functional annotation

``` bash

function run_interproscan(){
    ## replace * in the FASTA file
    sed -i 's/*//g' ${1}
    # sed '/>/{H;$!d} ; x ; s/^/\nSTART-->/ ; s/$/\n<--END/'
    prefix=`basename -z ${1} | sed -r 's/(.*)\..*$/\1/'`.interProScan

    if [ ! -f $2/${prefix}.gff3 ]
    then
        process_start "InterProScan on file $prefix"
        interproscan.sh -verbose -cpu 12 -goterms -f GFF3 -iprlookup  --output-file-base $2/${prefix} -i $1
        error_exit $?
    else
        echo "${prefix}.gff3 InterProScan results exists"
    fi
    
}

export -f run_interproscan

export -f run_interproscan

# ## GNU parallel on single server
# cat scripts/sub_commands/interproscan_batch.00 | parallel --keep-order --jobs 4 --halt now,fail=1 --results $PROJECT_DIR/logs/interproscan/{/.} --joblog $PROJECT_DIR/logs/interproscan/parallel.log run_interproscan {} ${ANALYSIS_DIR}

## GNU parallel on multiple servers
cat scripts/sub_commands/interproscan_batch.00 | parallel --keep-order --sshloginfile parallel_job_nodes.txt --halt now,fail=1 --dry-run --results $PROJECT_DIR/logs/interproscan/{/.} --joblog $PROJECT_DIR/logs/interproscan/parallel_across.log run_interproscan {} ${ANALYSIS_DIR}

# cat scripts/sub_commands/interproscan_batch.00 | parallel --sshloginfile parallel_job_nodes.txt echo "File {}: Running on \`hostname\`"




```

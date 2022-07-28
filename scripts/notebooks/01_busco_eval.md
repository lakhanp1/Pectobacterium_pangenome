# BUSCO completeness evaluation

## Setup

``` bash
#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
#

shopt -s expand_aliases
source ~/.bash_aliases

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/busco_eval"
busco_lineage='enterobacterales_odb10'

```

## Run BUSCO in parallel for all protein FAA files

``` bash

file_faa=${1}
prefix=`basename -z ${file_faa} | sed -r 's/(.*)\..*$/\1/'`

process_start "BUSCO on file $prefix"
busco -i ${file_faa} -o ${prefix} -l ${busco_lineage} -m proteins \
--offline --download_path $PROJECT_DIR/data/busco_downloads \
--cpu 8 --tar --out_path ${ANALYSIS_DIR} -f --quiet
error_exit $?

```

## execute the script using GNU `parallel`

``` bash
## !GNU parallel on single server
nohup \
cat scripts/sub_commands/interproscan_files.list | \
parallel --jobs 4 --workdir $PWD --halt now,fail=1 --keep-order \
--results logs/busco/{/.} --joblog logs/busco/parallel_batch00.log \
$PWD/scripts/a_processing/a03_busco_eval.sh {} \
>logs/busco/nohup_batch00.out 2>&1 &

```

This is failing with `parallel`, most likely because of the number of open files exceeding the `ulimit`

## Alternatively run BUSCO in serial mode

``` bash

for file_faa in `cat scripts/sub_commands/interproscan_files.list`
do
    prefix=`basename -z ${file_faa} | sed -r 's/(.*)\..*$/\1/'`

    process_start "BUSCO on file $prefix"
    busco -i ${file_faa} -o ${prefix} -l ${busco_lineage} -m proteins \
    --offline --download_path $PROJECT_DIR/data/busco_downloads \
    --cpu 8 --tar --out_path ${ANALYSIS_DIR} -f --quiet
    error_exit $?

done

```

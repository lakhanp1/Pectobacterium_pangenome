# Preprocessing

## Data download from NCBI

``` bash
## Working on script
```

## Prokka annotation for bacterial genomes

#### Processing using GNU `parallel`

``` bash
ls $PROJECT_DIR/data/genomes/?(ncbi|local)/*.?(fa|fna|fasta) | \
parallel --gnu --keep-order --jobs 4 --halt now,fail=1 \
--results logs/prokka/{/.} --joblog logs/prokka/parallel.log \
run_prokka {} ${ANALYSIS_DIR}

```

Processing in parallel is failing without any specific error. Hence the data will be processed in serial mode. Further debugging is required to understand the parallel processing failure.

#### Serial run

``` bash
bash scripts/a_processing/a01_prokka_serial.sh > logs/prokka/prokka.log 2>&1 
```

## BUSCO assembly evaluation

#### Execute the script using GNU `parallel`

``` bash
## !GNU parallel on single server
nohup \
cat scripts/sub_commands/interproscan_files.list | \
parallel --jobs 4 --workdir $PWD --halt now,fail=1 --keep-order \
--results logs/busco/{/.} --joblog logs/busco/parallel_batch00.log \
$PWD/scripts/a_processing/a03_busco_eval.sh {} \
>logs/busco/nohup_batch00.out 2>&1 &

```

This is failing with `parallel`, most likely because of the number of open files exceeding the `ulimit`. Need to debug further.

#### Alternatively run BUSCO in serial mode

``` bash
bash scripts/a_processing/a03_busco_eval.sh > logs/busco/busco.log 2>&1
```

## QUAST assembly evaluation

``` bash

```

## InterProScan functional annotation

#### Execute the script using GNU `parallel`

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

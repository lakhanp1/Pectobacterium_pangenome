# Preprocessing

## Data download from NCBI using EUtils

``` bash
## get all FPT paths for NCBI assemblies
esearch -db assembly -query "pectobacterium [ORGN]" | \
efetch -format docsum -email "lakhansing.pardeshi@wur.nl" | \
xtract -pattern DocumentSummary \
-tab "," -element AssemblyAccession AssemblyName FtpPath_RefSeq FtpPath_GenBank | \
awk '
    BEGIN { FS=","};
    {
    if($4=="") file=$3; else file=$3;
    gsub(/ +/, "_", $0);
    printf "%s/%s_%s_genomic.fna.gz\n", file, $1, $2 > "data/reference_data/ncbi_assembly_ftp.txt";
    printf "%s_%s\n", $1, $2 > "data/reference_data/ncbi_assembly_ids.txt";
    }
'

cat data/reference_data/ncbi_assembly_ftp.txt | \
xargs -I {} wget --timestamping {} -P data/genomes/ncbi/


```

## Prokka annotation for bacterial genomes

#### Processing using GNU `parallel`

``` bash
nohup \
ls data/genomes/?(ncbi|local)/*.?(fa|fna|fasta) | \
parallel --gnu --keep-order --jobs 6 --halt now,fail=1 \
--results logs/prokka/{/.} --joblog logs/prokka/parallel.log \
scripts/a_processing/a01_prokka_ann.sh {} \
>logs/prokka/nohup_batch00.out 2>&1 &

```

Processing in parallel is failing without any specific error. Hence the
data will be processed in serial mode. Further debugging is required to
understand the parallel processing failure.

#### Serial run

``` bash
nohup bash scripts/a_processing/a01_prokka_serial.sh >> logs/prokka/prokka.log 2>&1 &
```

## BUSCO assembly evaluation

#### Execute the script using GNU `parallel`

``` bash
## !GNU parallel on single server
nohup \
cat data/reference_data/sample_subsets/genebank_ass.txt | \
parallel --jobs 4 --workdir $PWD --halt now,fail=1 --keep-order \
--results logs/busco/{/.} --joblog logs/busco/parallel_batch00.log \
$PWD/scripts/a_processing/a03_busco_eval.sh {} \
>>logs/busco/nohup_batch00.out 2>&1 &

```

This is failing with `parallel`, most likely because of the number of open files
exceeding the `ulimit`. Need to debug further.

#### Alternatively run BUSCO in serial mode

``` bash
nohup bash scripts/a_processing/a03_busco_serial.sh >> logs/busco/busco.log 2>&1 &
```

## QUAST assembly evaluation

``` bash
nohup \
cat data/reference_data/assembly_ids.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/quast/{} --resume \
--joblog logs/quast/parallel.log \
./scripts/a_processing/a04_quast_eval.sh {} \
>>logs/quast/nohup.out 2>&1 &
```

## InterProScan functional annotation

#### Execute the script using GNU `parallel`

``` bash
## GNU parallel on single server
nohup \
cat data/reference_data/sample_subsets/genebank_ass.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{} \
--joblog logs/interproscan/parallel_batch05.log \
./scripts/a_processing/a02_interproscan.sh {} \
>>logs/interproscan/nohup_batch05.out 2>&1 &


## GNU parallel: submit jobs to another server
head -n 10 data/reference_data/sample_subsets/genebank_ass.txt | \
env_parallel --jobs 4 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{/.} --joblog logs/interproscan/parallel.log \
--sshlogin 4/waterman.bioinformatics.nl --cleanup \
--env error_exit --env process_start --env TOOLS_PATH --env LUSTRE_HOME \
./scripts/a_processing/a02_interproscan.sh {}

```

# Preprocessing

## Data download from NCBI using EUtils

``` bash

esearch -db assembly -query "pectobacterium [ORGN]" | \
efetch -format docsum -email "lakhansing.pardeshi@wur.nl" > \
data/reference_data/efetch_docsum.xml

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

```

## Check if there are any new assemblies which recently became available

``` bash
## check and extract new assembly IDs
for sampleId in `cat data/reference_data/assembly_ids.txt`
do
    if ! find data/genomes/ -name ${sampleId}* | grep -q "."
    then
        echo $sampleId
    fi
done > data/reference_data/temp_assembly_ids.txt


## download and unzip the FASTA files
for sampleId in `cat data/reference_data/temp_assembly_ids.txt`
do
    grep $sampleId data/reference_data/ncbi_assembly_ftp.txt
done | \
xargs -I {} wget --timestamping {} -P data/genomes/ncbi/

gzip -d data/genomes/ncbi/*.gz

```

## Prokka annotation for bacterial genomes

``` bash
nohup \
ls data/genomes/?(ncbi|local)/*.?(fa|fna|fasta) | \
parallel --gnu --keep-order --jobs 1 --halt now,fail=1 \
--results logs/prokka/{/.} --joblog logs/prokka/parallel.log \
scripts/a_preprocessing/a01_prokka_ann.sh {} \
>logs/prokka/nohup_batch00.out 2>&1 &

```

Processing in parallel is failing without any specific error. Hence the
data will be processed in serial mode. Further debugging is required to
understand the parallel processing failure. For now, running GNU `parallel`
with `--jobs 1` setting in serial mode.

### Remove the FASTA sequence at the end of prokka gff files

``` bash
for i in `cat cat data/reference_data/assembly_ids.txt`
do
sed -n '1,/##FASTA/ {/##FASTA/!p}' data/prokka_annotation/${i}/${i}.gff > \
data/prokka_annotation/${i}/${i}.gff3
done

```

## BUSCO assembly evaluation

``` bash
## !GNU parallel on single server
nohup \
cat data/reference_data/temp_assembly_ids.txt | \
parallel --jobs 1 --workdir $PWD --halt now,fail=1 --keep-order \
--results logs/busco/{} --joblog logs/busco/parallel.log \
$PWD/scripts/a_preprocessing/a03_busco_eval.sh {} \
>>logs/busco/nohup.out 2>&1 &

```

This is failing with `parallel`, most likely because of the number of open files
exceeding the `ulimit`. Need to debug further. For now, running GNU `parallel`
with `--jobs 1` setting is serial mode.

## QUAST assembly evaluation

``` bash
nohup \
cat data/reference_data/temp_assembly_ids.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/quast/{} \
--joblog logs/quast/parallel.log \
./scripts/a_preprocessing/a04_quast_eval.sh {} \
>>logs/quast/nohup.out 2>&1 &
```

## InterProScan functional annotation

``` bash
## GNU parallel on single server
nohup \
cat data/reference_data/temp_assembly_ids.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{} \
--joblog logs/interproscan/parallel_batch05.log \
./scripts/a_preprocessing/a02_interproscan.sh {} \
>>logs/interproscan/nohup_batch05.out 2>&1 &


## GNU parallel: submit jobs to another server
head -n 10 data/reference_data/sample_subsets/genebank_ass.txt | \
env_parallel --jobs 4 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{} \
--joblog logs/interproscan/parallel.log \
--sshlogin 4/waterman.bioinformatics.nl --cleanup \
--env error_exit --env process_start --env TOOLS_PATH --env LUSTRE_HOME \
./scripts/a_preprocessing/a02_interproscan.sh {}

```

## MultiQC

``` bash
## QUAST MultiQC
nohup \
multiqc -f --filename quast_multiqc --interactive --title "QUAST report" \
--outdir analysis/01_multiqc/ --module quast data/quast/ \
>>nohup.out 2>&1 &

## BUSCO MultiQC
nohup \
multiqc -f --filename busco_multiqc --interactive --title "BUSCO report" \
--outdir analysis/01_multiqc/ --module busco data/busco/ \
>>nohup.out 2>&1 &

```

## fastANI

``` bash
nohup bash scripts/a_preprocessing/a05_fastANI.sh >logs/fastANI.log 2>&1 &
```

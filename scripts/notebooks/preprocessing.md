# Preprocessing

## Data download from NCBI using EUtils

``` bash
## number of records available in NCBI
esearch -db assembly -query "pectobacterium [ORGN]" -email "lakhansing.pardeshi@wur.nl"

## download the Assembly DocumentSummary file
esearch -db assembly -query "pectobacterium [ORGN]" -email "lakhansing.pardeshi@wur.nl" | \
efetch -format docsum > data/reference_data/assembly_docsum.xml

## download associated BioSample DocumentSummary file
esearch -db assembly -query "pectobacterium [ORGN]" -email "lakhansing.pardeshi@wur.nl" | \
elink -target biosample  | \
efetch -format docsum > data/reference_data/biosample_docsum.xml


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

## NCBI prokaryotes ANI file
wget --timestamping https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt -P data/other/


cat data/reference_data/local_assembly_ids.txt data/reference_data/ncbi_assembly_ids.txt \
> data/reference_data/pre_qc_assembly_ids.txt

mkdir data/reference_data/batches
split -a 2 -d -n r/8 data/reference_data/pre_qc_assembly_ids.txt data/reference_data/batches/temp_ids_batch.

```

## Check if there are any new assemblies released on NCBI

``` bash
## check and extract new assembly IDs
for sampleId in `cat data/reference_data/pre_qc_assembly_ids.txt`
do
    if ! find data/genomes/ -name ${sampleId}* | grep -q "."
    then
        echo $sampleId
    fi
done > data/reference_data/temp_assembly_ids.txt


## download and unzip the FASTA files
for sampleId in `cat data/reference_data/temp_assembly_ids.txt`
do
    grep "/$sampleId/" data/reference_data/ncbi_assembly_ftp.txt
done | \
xargs -I {} wget --timestamping {} -P data/genomes/ncbi/

gzip -d data/genomes/ncbi/*.gz

```

## Run standard annotations for all the genomes

### Prokka annotation for bacterial genomes

``` bash
nohup \
ls data/genomes/?(NCBI|WUR|NIVIP|NAK)/*.?(fa|fna|fasta) | \
parallel --gnu --keep-order --jobs 1 --halt now,fail=1 \
--results logs/prokka/{/.} --joblog logs/prokka/parallel.log \
scripts/preprocessing/prokka_annotation.sh {} \
>logs/prokka/nohup_batch00.out 2>&1 &

```

Processing in parallel is failing without any specific error. Hence the
data will be processed in serial mode. Further debugging is required to
understand the parallel processing failure. For now, running GNU `parallel`
with `--jobs 1` setting in serial mode.

#### Remove the FASTA sequence at the end of prokka gff files

``` bash
for i in `cat data/reference_data/pre_qc_assembly_ids.txt`
do
sed -n '1,/##FASTA/ {/##FASTA/!p}' data/prokka_annotation/${i}/${i}.gff > \
data/prokka_annotation/${i}/${i}.gff3
done

```

#### FASTA file indexing

``` bash
printf '' > analysis/QC/assembly_chr_size.txt
printf '' > data/reference_data/pre_qc_genomes_fa.list
conda activate omics_py37

for i in `cat data/reference_data/pre_qc_assembly_ids.txt`
do
    samtools faidx data/prokka_annotation/${i}/${i}.fna
    sort -r -n -k 2,2 data/prokka_annotation/${i}/${i}.fna.fai | \
    awk  -v i=${i} '{print i, "\t", $1, "\t", $2}' >> analysis/QC/assembly_chr_size.txt
    ls data/prokka_annotation/${i}/${i}.fna >> data/reference_data/pre_qc_genomes_fa.list

    faToTwoBit data/prokka_annotation/${i}/${i}.fna data/prokka_annotation/${i}/${i}.2bit

done

```

### BUSCO protein and genome assembly evaluation

``` bash
## !GNU parallel on single server
nohup \
cat data/reference_data/temp_assembly_ids.txt | \
parallel --jobs 1 --workdir $PWD --halt now,fail=1 --keep-order \
--results logs/busco/{} --joblog logs/busco/parallel.log \
$PWD/scripts/preprocessing/qc_busco.sh {} \
>logs/busco/nohup03.out 2>&1 &

```

This is failing with `parallel`, most likely because of the number of open files
exceeding the `ulimit`. Need to debug further. For now, running GNU `parallel`
with `--jobs 1` setting is serial mode.

### QUAST assembly evaluation

``` bash
nohup \
cat data/reference_data/temp_assembly_ids.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/quast/{} \
--joblog logs/quast/parallel.log \
./scripts/preprocessing/qc_quast.sh {} \
>>logs/quast/nohup.out 2>&1 &
```

### InterProScan functional annotation

``` bash
## GNU parallel on single server
nohup \
cat data/reference_data/temp_assembly_ids.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{} \
--joblog logs/interproscan/parallel_batch06.log \
./scripts/preprocessing/interproscan_annotation.sh {} \
>>logs/interproscan/nohup_batch07.out 2>&1 &


## GNU parallel: submit jobs to another server
head -n 10 data/reference_data/sample_subsets/genebank_ass.txt | \
env_parallel --jobs 4 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/interproscan/{} \
--joblog logs/interproscan/parallel.log \
--sshlogin 4/waterman.bioinformatics.nl --cleanup \
--env error_exit --env process_start --env TOOLS_PATH --env LUSTRE_HOME \
./scripts/preprocessing/interproscan_annotation.sh {}

```

### eggNOG annotation

``` bash
conda activate eggnog

## download eggNOG data
download_eggnog_data.py -M -F -H -d 1 --dbname root -f -y
download_eggnog_data.py -D -H -d 2759 --dbname eukaryota -q -y
download_eggnog_data.py -D -H -d 2 --dbname bacteria -q -y
download_eggnog_data.py -D -H -d 10239 --dbname viruses -q -y
download_eggnog_data.py -D -H -d 4751 --dbname fungi -q -y
download_eggnog_data.py -D -H -d 2157 --dbname archea -q -y

mmseqs createindex "$EGGNOG_DATA_DIR"/mmseqs/mmseqs.db /local_scratch/parde001/tmp

## create mmseq database for bacteria, archea and viruses
create_dbs.py -m mmseqs --dbname bact_arch_vir --taxids 2,2157,10239 -y

## annotate for eggNOG clusters
for i in `cat data/reference_data/assembly_ids.txt`
do
    file_aa="data/prokka_annotation/${i}/${i}.faa"
    
    emapper.py --cpu 20 -m mmseqs --dbmem --go_evidence all  \
    -i ${file_aa} -o ${i} --output_dir data/eggnog \
    --scratch_dir /local_scratch/parde001/ --temp_dir /local_scratch/parde001/tmp

    error_exit $?
done

```

### MultiQC

``` bash
## QUAST MultiQC
nohup \
multiqc -f --filename quast_multiqc --interactive --title "QUAST report" \
--outdir analysis/multiqc/ --module quast data/quast/ \
>>nohup.out 2>&1 &

## BUSCO protein MultiQC
nohup \
multiqc -f --filename busco_prot_multiqc --interactive --title "BUSCO report" \
--outdir analysis/multiqc/ --module busco data/busco.prot/ \
>>nohup.out 2>&1 &

## BUSCO genome MultiQC
nohup \
multiqc -f --filename busco_geno_multiqc --interactive --title "BUSCO report" \
--outdir analysis/multiqc/ --module busco data/busco.geno/ \
>>nohup.out 2>&1 &

```

### fastANI

``` bash
process_start "ANI on all genomes"
nohup bash scripts/preprocessing/ANI_processing.sh \
data/reference_data/pre_qc_genomes_fa.list \
analysis/02_fastANI/ANI_results >logs/fastANI.log 2>&1 &
error_exit $?
```

### Defense-finder

```bash
conda activate defensefinder

for i in `cat data/reference_data/assembly_ids.txt`
do
    file_aa="data/prokka_annotation/${i}/${i}.faa"
    out_dir="data/defense_finder/${i}"
    mkdir -p ${out_dir}

    process_start "defense-finder: ${i}"

    defense-finder run --workers 20 -o ${out_dir} \
    --models-dir $TOOLS_PATH/defense-finder-models ${file_aa}

    error_exit $?
done
```

Combine Defense-finder output into single file

```bash
# combine systems
for i in `cat data/reference_data/assembly_ids.txt`
do
    awk -v var="${i}" '!/^sys_id/ {print var "\t" $0}' data/defense_finder/${i}/${i}_defense_finder_systems.tsv
done > data/defense_finder/defense_systems.tsv

# combine genes
for i in `cat data/reference_data/assembly_ids.txt`
do
    awk -v var="${i}" 'NR>1 {print var "\t" $0}' data/defense_finder/${i}/${i}_defense_finder_genes.tsv
done > data/defense_finder/defense_systems_genes.tsv
```

## MMSeq2 database

```
conda activate mmseq2

mmseqs createdb data/pangenomes/pectobacterium.v2/mmseq_panproteome_db
```

## Prepare data for pangenome construction

### Extract metadata from NCBI XML files

```bash
Rscript scripts/preprocessing/raw_prebuild_metadata.R
```

### Correct taxonomy

```bash
Rscript scripts/preprocessing/taxonomy_check_inhouse.R
```

### Prepare final metadata table with updated taxonomy

```bash
Rscript scripts/preprocessing/build_metadata.R
```

### Prepare PanTools input files

```bash
Rscript scripts/preprocessing/prepare_PanTools_input.R
```

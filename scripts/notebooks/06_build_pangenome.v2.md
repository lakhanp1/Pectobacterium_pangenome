# Pectobacterium pangenome

## Setup

``` bash
shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

conda activate pantools

export PANTOOLS="$PANTOOLS_MAIN"

## Setup
PROJECT_DIR='/mnt/scratch/parde001/projects/03_Pectobacterium'
ANALYSIS_DIR="$PROJECT_DIR/analysis/06_pangenome_pecto_v2"
analysis_prefix='pectobacterium.v2'
pan_db="$ANALYSIS_DIR/${analysis_prefix}.DB"

[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR

```

## Download the raw data

```bash
esearch -db assembly -query "pectobacterium [ORGN]" | \
efetch -format docsum -email "lakhansing.pardeshi@wur.nl" | \
xtract -pattern DocumentSummary \
-tab " || " -element FtpPath_RefSeq | \
sed -r 's/((.*\/)([^\/]+))/\1\/\3_genomic.fna.gz/' | \
sed -r 's/^ftp:/rsync:/' > $PROJECT_DIR/data/reference_data/ncbi_genomes_ftp.txt

# data/genomes/ncbi_genomes
for file in `cat $PROJECT_DIR/data/reference_data/ncbi_genomes_ftp.txt`
do
rsync --copy-links --times --verbose ${file} $PROJECT_DIR/data/genomes/ncbi_genomes/
sleep 5
done

wget -i $PROJECT_DIR/data/reference_data/ncbi_genomes_ftp.txt
rsync --copy-links --times --verbose --files-from=$PROJECT_DIR/data/reference_data/ncbi_genomes_ftp.txt $PROJECT_DIR/data/genomes/ncbi_genomes/

rsync --dry-run --copy-links --times --verbose --files-from=temp.txt $PROJECT_DIR/data/genomes/ncbi_genomes/
rsync --dry-run --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/914/135/GCF_023914135.1_ASM2391413v1/GCF_023914135.1_ASM2391413v1_genomic.fna.gz $PROJECT_DIR/data/genomes/ncbi_genomes/

```


```
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

# *In*-silico PCR assay mapping to genomes

## Setup

``` bash
#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

# source scripts/utils/setup_analysis.sh 'pectobacterium.v2'
source scripts/utils/setup_analysis.sh $@

if [ -z ${pan_db+x} ];
then
    echo "\$pan_db is unset"
    error_exit 1
fi

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master
export PANTOOLS="$PANTOOLS_4_1"

######################################################################
```

## LZI and TIR assay genes

### LZI assay

- gene: GOECMFEC_02714
- sample: NAK238
- genome: 401
- location pangenome: 401_10:140994-141296(+)
- location fna: NAK238_contig_10_consensus:140994-141296(+)
- homology group: 22434221
- COG description: 4-amino-4-deoxy-alpha-L-arabinopyranosyl undecaprenyl phosphate biosynthetic process
- This gene seems to be part of prophage. However, this region was not marked as prophage
  by the genomad pipeline.


### TIR assay

- gene: GOECMFEC_02263
- sample: NAK238
- genome: 401
- location pangenome: 401_8:166700-167047(-)
- location fna: NAK238_contig_8_consensus:166700-167047(-)
- homology group: 22434077
- COG description: MTH538 TIR-like domain (DUF1863)
- This is a loss of gene event. No prophage nearby.


## PCR primer sequence BLAST against the genomes

### combine multiple genome FASTA into a single file

``` bash
## PCR assay blast
ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_v2/insilico_assay"
[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR

## combine genomes and prepare blastn index
python3 scripts/a_preprocessing/a10_merge_genome_faa.py --geno $PANGENOME_DIR/genomes_fa.list \
--out $PANGENOME_DIR/blastdb/genomes_combined.fa

## index the sequences
makeblastdb -dbtype nucl -in $PANGENOME_DIR/blastdb/genomes_combined.fa -parse_seqids \
-title $PANGENOME_NAME -logfile $PANGENOME_DIR/blastdb/makeblastdb.log

```

### run BLASTN

**NOTE:**
`blastn` was not showing hits of some primer sequences when searched against all the genomes. Hence, database sequences are split into batches of 100 genomes and `blastdb_aliastool` is run to preprocess the sequence ids. These preprocessed sequence ids are provided as input to `-seqidlist` argument while running `blastn`. Somehow, when the search space is limited using `-seqidlist` option, missing hits are shown. This is run for each set of sequence ids.

``` bash

for chrset in $PANGENOME_DIR/blastdb/chrset_*.acc
do
    ##
    blastdb_aliastool -seqid_file_in ${chrset}

    setName=`basename ${chrset} | sed 's/.acc//'`
    blastnOut="pcr_primer_blastn."${setName}".out"

    ## run blastn in blastn-short mode
    blastn -db $PANGENOME_DIR/blastdb/genomes_combined.fa \
    -query data/other/pecto_assay/pectobacterium_assay_probes.fasta \
    -task blastn-short -num_threads 12 -max_target_seqs 5000 -seqidlist ${chrset}.bsl \
    -outfmt "6 qseqid qstart qend qlen sseqid sstart send sstrand slen pident length mismatch qcovs gapopen evalue bitscore" \
    -out $ANALYSIS_DIR/${blastnOut}

done

cat $ANALYSIS_DIR/pcr_primer_blastn.chrset_*.out > $ANALYSIS_DIR/pcr_primer_blastn.out 
```

## parse & summarize blastn result table

``` bash
quarto render scripts/analysis/insilico_assay.viz.qmd
```

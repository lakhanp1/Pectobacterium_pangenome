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

## PCR assay blast
ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_v2/insilico_assay"
[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR

## combine genomes and prepare blastn index
python3 scripts/preprocessing/merge_genome_faa.py --geno $PANGENOME_DIR/genomes_fa.list \
--out $PANGENOME_DIR/blastdb/genomes_combined.fa

## index the sequences
makeblastdb -dbtype nucl -in $PANGENOME_DIR/blastdb/genomes_combined.fa -parse_seqids \
-title $PANGENOME_NAME -logfile $PANGENOME_DIR/blastdb/makeblastdb.log


## run BLASTN
for chrset in $PANGENOME_DIR/blastdb/chrset_*.acc
do
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


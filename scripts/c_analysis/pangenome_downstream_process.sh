#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

scriptName=`basename $0`

usage="USAGE:
-------------------------------------------------------------
bash ${scriptName} <db_name> <exists> <suffix>
db_name   : STRING name for pangenome database
-------------------------------------------------------------
"

if [ $# -ne 1 ]; then
    printf "Error: Require pangenome name\n${usage}" >&2 ; exit 1
fi

######################################################################

PANGENOME_NAME=$1
# PANGENOME_NAME='pectobacterium.10g'
# PANGENOME_NAME='pectobacterium.v2'

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master

export PANTOOLS="$PANTOOLS_4_1"

## Setup
PROJECT_DIR="/lustre/BIF/nobackup/$USER/projects/03_Pectobacterium"
PANGENOME_DIR="$PROJECT_DIR/data/pangenomes/$PANGENOME_NAME"
pan_db="$PANGENOME_DIR/${PANGENOME_NAME}.DB${DB_SUFFIX}"

printf "PANGENOME_DIR: ${PANGENOME_DIR}
pangenome: ${pan_db}
"

######################################################################
### UPGMA tree
#./scripts/c_analysis/phylo.process_newicks.R \
#--tree ./analysis/04_pangenome_pecto_v2/phylogeny/ani/ANI_UPGMA.pangenome.newick \
#--name "ani_upgma" --config project_config.yaml

### Neighbor joining tree
#./scripts/c_analysis/phylo.process_newicks.R \
#--tree ./analysis/04_pangenome_pecto_v2/phylogeny/ani/ANI_NJ.pangenome.newick \
#--name "ani_nj" --config project_config.yaml

### gene distance NJ tree
#./scripts/c_analysis/phylo.process_newicks.R \
#--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/gene_classification.100.0/gene_distance.tree \
#--name "gene_nj" --config project_config.yaml

### kmer distance NJ tree
#./scripts/c_analysis/phylo.process_newicks.R \
#--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/kmer_classification.100.0/genome_kmer_distance.tree \
#--name "kmer_nj" --config project_config.yaml

### core-SNP ML tree
#./scripts/c_analysis/phylo.process_newicks.R \
#--tree ./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/core_snp_tree/informative.fasta.treefile \
#--name "core_snp_ml" --config project_config.yaml

######################################################################
### Phenotype association
### add the updated phenotypes for association analysis
#process_start add_phenotypes
#$PANTOOLS remove_phenotype ${pan_db}
#$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/genomes_metadata.csv
#$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/analysis_configs/clade_association_phenotypes.csv
#error_exit $?


#[ -d ${pan_db}/gene_classification ] && rm -r ${pan_db}/gene_classification
#mkdir ${pan_db}/gene_classification.pheno

### Gene classification for each phenotype
#phenotypes=(`awk -F "\t" '{ if (NR!=1) {print $1} }' $PANGENOME_DIR/analysis_configs/pheno_association_config.tab`)
#for phn in ${phenotypes[@]}
#do
#    process_start "gene_classification for phenotype $phn"
#    pheno_arg=`grep "^${phn}\b" $PANGENOME_DIR/analysis_configs/pheno_association_config.tab | cut -f2`
#    $PANTOOLS gene_classification ${pheno_arg} ${pan_db}
#    error_exit $?

#    ## move results to a folder
#    pheno_dir=${pan_db}/gene_classification.pheno/${phn}
#    [ -d ${pheno_dir} ] && rm -r ${pheno_dir}
#    mkdir ${pheno_dir}
#    mv ${pan_db}/gene_classification/{phenotype_*,gene_classification_phenotype_overview.txt} ${pheno_dir}/
#done

#rm -r ${pan_db}/gene_classification

######################################################################
### Use subset of genomes to determine pangenome structure
#genomeSets=(`awk -F "\t" '{ print $1 }' $PANGENOME_DIR/genome_sets.tab`)

#for gs in ${genomeSets[@]}
#do
#    process_start "gene_classification for genome subset $gs"
#    str_arg=`grep "^${phn}\b" $PANGENOME_DIR/genome_sets.tab | cut -f2`
#    $PANTOOLS pangenome_structure -t 20 ${str_arg} ${pan_db}
#    error_exit $?

#    ## move results to a folder
#    gs_dir=${pan_db}/pangenome_size/gene.${gs}
#    [ -d ${gs_dir} ] && rm -r ${gs_dir}
#    mv ${pan_db}/pangenome_size/gene ${gs_dir}
#done

######################################################################
### PCR assay blast
#ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_v2/insilico_assay"
#[ ! -d $ANALYSIS_DIR ] && mkdir $ANALYSIS_DIR

### combine genomes and prepare blastn index
#python3 scripts/a_preprocessing/a10_merge_genome_faa.py --geno $PANGENOME_DIR/genomes_fa.list \
#--out $PANGENOME_DIR/blastdb/genomes_combined.fa

### index the sequences
#makeblastdb -dbtype nucl -in $PANGENOME_DIR/blastdb/genomes_combined.fa -parse_seqids \
#-title $PANGENOME_NAME -logfile $PANGENOME_DIR/blastdb/makeblastdb.log

#for chrset in $PANGENOME_DIR/blastdb/chrset_*.acc
#do
#    ##
#    blastdb_aliastool -seqid_file_in ${chrset}

#    setName=`basename ${chrset} | sed 's/.acc//'`

#    ## run blastn in blastn-short mode
#    blastn -db $PANGENOME_DIR/blastdb/genomes_combined.fa \
#    -query data/other/pecto_assay/pectobacterium_assay_probes.fasta \
#    -task blastn-short -num_threads 12 -max_target_seqs 5000 -seqidlist ${chrset}.bsl \
#    -outfmt "6 qseqid qstart qend qlen sseqid sstart send sstrand slen pident length mismatch qcovs gapopen evalue bitscore" \
#    -out $ANALYSIS_DIR/${blastnOut}

#done

#Rscript scripts/c_analysis/c02_blastn_parse.R
#######################################################################

### Extract specific information from pangenome
## gene_classification.pheno/${phn}/phenotype_association.csv file
#conda activate omics_py37
#while IFS=$'\t' read -r phn groups ; do
#    # printf "%b\n" "column1<${phn}>"
#    # printf "%b\n" "column2<${groups}>"
#    process_start "extracting mRNA sequence for homology group specific to phenotype: ${phn}"
#    printf "Homology groups: ${groups}\n"

#    genome=`grep "^${phn}\b" $PANGENOME_DIR/analysis_configs/pheno_association_config.tab | cut -f3 | sed 's/,.*//'`
#    pheno_dir="analysis/04_pangenome_pecto_v2/pheno_association/${phn}"

#    file_info="${pheno_dir}/${phn}.pheno_specific.seq_info.txt"
#    file_info_g1="${pheno_dir}/${phn}.pheno_specific.${genome}g.seq_info.txt"
#    file_fasta="${pheno_dir}/${phn}.pheno_specific.${genome}g.mRNA.fasta"

#    printf "homology_group_id\tGenome\tmRNA_name\tmRNA_identifier\tnode_identifier\tgenome_id\tchr\tstart\tend\tstrand\n"  > ${file_info}
#    printf "homology_group_id\tGenome\tmRNA_name\tmRNA_identifier\tnode_identifier\tgenome_id\tchr\tstart\tend\tstrand\n"  > ${file_info_g1}
#    printf "" > ${file_fasta}

#    ## read groups into array
#    IFS=',' read -ra hg_array <<< "${groups}"
#    for hg in "${hg_array[@]}"
#    do
#        # printf "${hg} "
#        seq_info=`grep "^${genome}" ${pan_db}/alignments/msa_per_group/grouping_v4/${hg}/input/sequences.info`
#        mrna_id=`echo ${seq_info} | cut -d"," -f3`
#        samtools faidx ${pan_db}/alignments/msa_per_group/grouping_v4/${hg}/input/nuc.fasta

#        printf "${hg}, ${seq_info}\n" | sed -r 's/(, | )/\t/g' >> ${file_info_g1}

#        cat ${pan_db}/alignments/msa_per_group/grouping_v4/${hg}/input/sequences.info |
#        sed -n "1,/^#genome/! s/^\(.\)/${hg}, \1/p" |
#        sed -r 's/(, | )/\t/g' >> ${file_info}

#        samtools faidx ${pan_db}/alignments/msa_per_group/grouping_v4/${hg}/input/nuc.fasta ${mrna_id} >> ${file_fasta}
#    done
#done < <(grep '^assay_FN\b' analysis/04_pangenome_pecto_v2/pheno_association/phenotype_specific_groups.txt)

######################################################################

#conda activate omics_py37

### run blastn on phenotype specific sequences against pangenome to verify that
### the homology group based specificity is not because of lack of annotation
#process_start "blastn phenotype specific sequences against pangenome"

#blastn -db $PANGENOME_DIR/blastdb/genomes_combined.fa \
#-query analysis/04_pangenome_pecto_v2/pheno_association/assay_FN/assay_FN.pheno_specific.399g.mRNA.fasta \
#-num_threads 12 -max_target_seqs 5000  \
#-outfmt "6 qseqid qstart qend qlen sseqid sstart send sstrand slen pident length mismatch qcovs gapopen evalue bitscore" \
#-out analysis/04_pangenome_pecto_v2/pheno_association/assay_FN/pheno_specific.399g.pangenome.blastn.out

#error_exit $?
######################################################################



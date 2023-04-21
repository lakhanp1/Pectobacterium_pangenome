#!/usr/bin/env bash

## assoication analysis using pyseer
## genotype data used: PAV, SNP, k-mer data
## association analysis: fixed effect, mixed effect, lineage effect, elastic net

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
conda activate pyseer
######################################################################

## analysis setup
analysis_prefix='pbr_virulence'
analysis_dir="analysis/04_pangenome_pecto_v2/pheno_association/${analysis_prefix}"

## PAV
#genotype='pav'
#genotype_data_arg="--pres <>"

## SNP
#genotype='snp'
#genotype_data_arg="--vcf <>"

## k-mers
genotype='kmer'
unitigs_name='pbr'  # pangenome, pbr
genotype_data_arg="--kmers data/unitigs_${unitigs_name}/unitigs_${unitigs_name}.pyseer.gz"

#####################################################################
## global setup for pyseer
phenotype='virulence'
phn_file="${analysis_dir}/pyseer_pbr.phenotype.tab"
result_dir="${analysis_dir}/${genotype}_association"
[ ! -d $result_dir ] && mkdir ${result_dir}
result_prefix="${genotype}.${phenotype}_asso"

## pyseer: fixed effect model
pyseer ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--distances ${analysis_dir}/phylogeny_dist.tab \
--save-m ${analysis_dir}/phylogeny_dist.MDS \
--output-patterns ${result_dir}/${genotype}_patterns.fixed_eff.txt --cpu 30 \
> ${result_dir}/${result_prefix}.fixed_eff.txt 2>>logs/v2_pecto/pyseer.log

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${result_dir}/${genotype}_patterns.fixed_eff.txt \
> ${result_dir}/threshold.fixed_eff.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${result_dir}/threshold.fixed_eff.txt`

cat <(head -1 ${result_dir}/${result_prefix}.fixed_eff.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${result_dir}/${result_prefix}.fixed_eff.txt) \
> ${result_dir}/${result_prefix}.fixed_eff.significant.txt

#############################

## pyseer: mixed effect model
pyseer --lmm ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--similarity ${analysis_dir}/phylogeny_K.tsv \
--save-lmm ${analysis_dir}/pyseer_pbr.lmm \
--output-patterns ${result_dir}/${genotype}_patterns.lmm.txt --cpu 30 \
> ${result_dir}/${result_prefix}.lmm.txt 2>>logs/v2_pecto/pyseer.log

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${result_dir}/${genotype}_patterns.lmm.txt \
> ${result_dir}/threshold.lmm.txt

cat <(head -1 ${result_dir}/${result_prefix}.lmm.txt) \
<(awk '$4<1.23E-07 {print $0}' ${result_dir}/${result_prefix}.lmm.txt) \
> ${result_dir}/${result_prefix}.lmm.significant.txt

#############################

## pyseer: lineage effect
pyseer --lineage ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--load-m ${analysis_dir}/phylogeny_dist.MDS.pkl \
--lineage-file ${result_dir}/lineage_effects.txt \
--output-patterns ${result_dir}/${genotype}_patterns.lineage_eff.txt --cpu 30 \
> ${result_dir}/${result_prefix}.lineage_eff.txt 2>>logs/v2_pecto/pyseer.log

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${result_dir}/${genotype}_patterns.lineage_eff.txt \
> ${result_dir}/threshold.lineage_eff.txt

cat <(head -1 ${result_dir}/${result_prefix}.lineage_eff.txt) \
<(awk '$4<1.23E-07 {print $0}' ${result_dir}/${result_prefix}.lineage_eff.txt) \
> ${result_dir}/${result_prefix}.lineage_eff.significant.txt

#############################




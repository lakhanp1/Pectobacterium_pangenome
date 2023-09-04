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
genotype='pav'
genotype_data_arg="--pres ${analysis_dir}/pbr.accessory_PAV.tab"

## SNP
genotype='snp'
genotype_data_arg="--vcf <>"

genotype='kmer'
unitigs_name='pbr' #pangenome, pbr
genotype_data_arg="--kmers data/unitigs_${unitigs_name}/unitigs_${unitigs_name}.pyseer.gz"

## global setup
phenotype='virulence'
phn_file="${analysis_dir}/pbr_phenotypes.tab"
result_dir="${analysis_dir}/${phenotype}"
[ ! -d $result_dir ] && mkdir ${result_dir}
out_prefix="${result_dir}/${genotype}.${phenotype}.asso"

#############################

## pyseer: fixed effect model
pyseer ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--distances ${analysis_dir}/phylogeny_dist.tab \
--save-m ${analysis_dir}/phylogeny_dist.MDS \
--output-patterns ${out_prefix}.patterns.fixed_eff.txt --cpu 30 \
> ${out_prefix}.fixed_eff.txt

error_exit $?

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${out_prefix}.patterns.fixed_eff.txt \
> ${out_prefix}.threshold.fixed_eff.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${result_dir}/threshold.fixed_eff.txt`

cat <(head -1 ${out_prefix}.fixed_eff.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${out_prefix}.fixed_eff.txt) \
> ${out_prefix}.fixed_eff.significant.txt

#############################

## similarity
python $TOOLS_PATH/pyseer/scripts/phylogeny_distance.py --lmm \
${analysis_dir}/phylogeny_tree.newick > ${analysis_dir}/phylogeny_K.tsv

## pyseer: mixed effect model
pyseer --lmm ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--similarity ${analysis_dir}/phylogeny_K.tsv \
--save-lmm ${analysis_dir}/pyseer_pbr.lmm \
--output-patterns ${out_prefix}.patterns.lmm.txt --cpu 30 \
> ${out_prefix}.lmm.txt

error_exit $?

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${out_prefix}.patterns.lmm.txt \
> ${out_prefix}.threshold.lmm.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${result_dir}/threshold.lmm.txt`

cat <(head -1 ${out_prefix}.lmm.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${out_prefix}.lmm.txt) \
> ${out_prefix}.lmm.significant.txt

#############################

## pyseer: lineage effect
pyseer --lineage ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--load-m ${analysis_dir}/phylogeny_dist.MDS.pkl \
--lineage-file ${out_prefix}.lineage_effects.txt \
--output-patterns ${out_prefix}.patterns.lineage_eff.txt --cpu 30 \
> ${out_prefix}.lineage_eff.txt

error_exit $?

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${out_prefix}.patterns.lineage_eff.txt \
> ${out_prefix}.threshold.lineage_eff.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${result_dir}/threshold.lineage_eff.txt`

cat <(head -1 ${out_prefix}.lineage_eff.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${out_prefix}.lineage_eff.txt) \
> ${out_prefix}.lineage_eff.significant.txt

#############################



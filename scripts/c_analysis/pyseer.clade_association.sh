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
analysis_prefix='pyseer_clades'
analysis_dir="analysis/04_pangenome_pecto_v2/pheno_association/${analysis_prefix}"

## PAV
genotype='pav'
genotype_data_arg="--pres ${analysis_dir}/clade.accessory_PAV.tab"

#genotype='kmer'
#unitigs_name='pbr' #pangenome, pbr
#genotype_data_arg="--kmers data/unitigs_${unitigs_name}/unitigs_${unitigs_name}.pyseer.gz"

## global setup
phenotype='assay_FN'
phn_file="${analysis_dir}/clade_phenotypes.tab"
result_dir="${analysis_dir}/${phenotype}"
[ ! -d $result_dir ] && mkdir ${result_dir}
out_prefix="${result_dir}/${genotype}.${phenotype}.asso"

#############################

## pyseer: fixed effect model without population structure
pyseer  --no-distances ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--output-patterns ${out_prefix}.patterns.no_dist.txt --cpu 30 \
> ${out_prefix}.no_dist.txt

error_exit $?

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${out_prefix}.patterns.no_dist.txt \
> ${out_prefix}.threshold.no_dist.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${out_prefix}.threshold.no_dist.txt`

cat <(head -1 ${out_prefix}.no_dist.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${out_prefix}.no_dist.txt) \
> ${out_prefix}.no_dist.significant.txt

#############################



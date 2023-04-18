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
conda activate pyseer
######################################################################

## pangenome unitigs
genomes_list="${PANGENOME_DIR}/genomes_fa.list"
unitigs_name='pangenome'

## P. brasiliense unitigs
genomes_list="analysis/04_pangenome_pecto_v2/pheno_association/pbr_virulence/pbr.fa.list"
unitigs_name='pbr'

## generate unitigs
build_dir="$PAN_BUILD_DIR/unitigs_${unitigs_name}"
[ ! -d $build_dir ] && mkdir $build_dir

unitig-caller --call --threads 30 --write-graph --pyseer \
--out ${build_dir}/unitigs_${unitigs_name} --refs ${genomes_list}

error_exit $?

unitigs_dir='data/unitigs'
[ ! -d $unitigs_dir ] && mkdir $unitigs_dir

mv ${build_dir} $unitigs_dir/
gzip ${unitigs_dir}/unitigs_${unitigs_name}/*.pyseer

## store graph
graph_base="${genomes_list%.*}"
mv ${graph_base}.gfa ${unitigs_dir}/unitigs_${unitigs_name}/${unitigs_name}.gfa
mv ${graph_base}.color.bfg ${unitigs_dir}/unitigs_${unitigs_name}/${unitigs_name}.color.bfg
######################################################################


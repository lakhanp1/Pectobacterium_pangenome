# Association analysis using `pyseer`

## Preprocessing

### Setup

```bash
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
```

### Generate unitigs for pangenome

```bash
## pangenome unitigs
genomes_list="${PANGENOME_DIR}/genomes_fa.list"
unitigs_name='pangenome'
```

### generate unitigs for *P. brasiliense* genomes

```bash
## P. brasiliense unitigs
genomes_list="analysis/04_pangenome_pecto_v2/pheno_association/pbr_virulence/pbr.fa.list"
unitigs_name='pbr'
```

Call unitigs

```bash
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
```

## Run `pyseer`

```bash
## analysis setup
analysis_prefix='pbr_virulence'
analysis_dir="analysis/04_pangenome_pecto_v2/pheno_association/${analysis_prefix}"
```

### population structure

```bash
## similarity
python $TOOLS_PATH/pyseer/scripts/phylogeny_distance.py --lmm \
${analysis_dir}/phylogeny_tree.newick > ${analysis_dir}/phylogeny_K.tsv
```

### `pyseer` command setup

#### `pyseer` on homology group PAV

```bash
## PAV
genotype='pav'
genotype_data_arg="--pres ${analysis_dir}/pbr.accessory_PAV.tab"
```

#### `pyseer` on core SNPs

```bash
## SNP
genotype='snp'
genotype_data_arg="--vcf <>"
```

#### `pyseer` on unitigs

```bash
genotype='kmer'
unitigs_name='pbr' #pangenome, pbr
genotype_data_arg="--kmers data/unitigs_${unitigs_name}/unitigs_${unitigs_name}.pyseer.gz"
```

### `pyseer` commands

Setup the paths for analysis files.

```bash
## global setup
phenotype='virulence'
phn_file="${analysis_dir}/pbr_phenotypes.tab"
result_dir="${analysis_dir}/${genotype}_association"
[ ! -d $result_dir ] && mkdir ${result_dir}
result_prefix="${genotype}.${phenotype}_asso"
```

#### Fixed effects model

Running

```bash
## pyseer: fixed effect model
pyseer ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--distances ${analysis_dir}/phylogeny_dist.tab \
--save-m ${analysis_dir}/phylogeny_dist.MDS \
--output-patterns ${result_dir}/${genotype}_patterns.fixed_eff.txt --cpu 30 \
> ${result_dir}/${result_prefix}.fixed_eff.txt 

error_exit $?

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${result_dir}/${genotype}_patterns.fixed_eff.txt \
> ${result_dir}/threshold.fixed_eff.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${result_dir}/threshold.fixed_eff.txt`

cat <(head -1 ${result_dir}/${result_prefix}.fixed_eff.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${result_dir}/${result_prefix}.fixed_eff.txt) \
> ${result_dir}/${result_prefix}.fixed_eff.significant.txt
#############################
```

#### Mixed effect model (FaST-LMM)

Running

```bash
## pyseer: mixed effect model
pyseer --lmm ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--similarity ${analysis_dir}/phylogeny_K.tsv \
--save-lmm ${analysis_dir}/pyseer_pbr.lmm \
--output-patterns ${result_dir}/${genotype}_patterns.lmm.txt --cpu 30 \
> ${result_dir}/${result_prefix}.lmm.txt 

error_exit $?

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${result_dir}/${genotype}_patterns.lmm.txt \
> ${result_dir}/threshold.lmm.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${result_dir}/threshold.lmm.txt`

cat <(head -1 ${result_dir}/${result_prefix}.lmm.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${result_dir}/${result_prefix}.lmm.txt) \
> ${result_dir}/${result_prefix}.lmm.significant.txt
#############################
```

#### Lineage effect: using MDS components as lineage

```bash
## pyseer: lineage effect
pyseer --lineage ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--load-m ${analysis_dir}/phylogeny_dist.MDS.pkl \
--lineage-file ${result_dir}/lineage_effects.txt \
--output-patterns ${result_dir}/${genotype}_patterns.lineage_eff.txt --cpu 30 \
> ${result_dir}/${result_prefix}.lineage_eff.txt 

error_exit $?

## get significance threshold
python $TOOLS_PATH/pyseer/scripts/count_patterns.py ${result_dir}/${genotype}_patterns.lineage_eff.txt \
> ${result_dir}/threshold.lineage_eff.txt

pval_cut=`sed -nr '/^Threshold/ s/Threshold:\s+//p' ${result_dir}/threshold.lineage_eff.txt`

cat <(head -1 ${result_dir}/${result_prefix}.lineage_eff.txt) \
<(awk -v pval_cut="${pval_cut}" '$4<pval_cut {print $0}' ${result_dir}/${result_prefix}.lineage_eff.txt) \
> ${result_dir}/${result_prefix}.lineage_eff.significant.txt
#############################
```

#### Whole genome models (elastic net)

```bash
pyseer --wg enet ${genotype_data_arg} \
--phenotypes ${phn_file} --phenotype-column ${phenotype} \
--load-m ${analysis_dir}/phylogeny_dist.MDS.pkl \
--cpu 30 --n-folds 5 --save-vars ${result_dir}/${result_prefix}.wg_vars \
--save-model ${result_dir}/${result_prefix}.wg_lasso \ 
> ${result_dir}/${result_prefix}.elastic_net.txt 

error_exit $?

```

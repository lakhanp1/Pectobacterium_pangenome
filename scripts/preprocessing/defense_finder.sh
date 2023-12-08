#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate defensefinder

######################################################################

for id in $(cat data/reference_data/assembly_ids.txt); do
    file_aa="data/prokka_annotation/${id}/${id}.faa"
    out_dir="data/defense_finder/${id}"
    mkdir -p ${out_dir}

    process_start "defense-finder: ${id}"

    defense-finder run --workers 20 -o ${out_dir} \
        --models-dir $TOOLS_PATH/defense-finder-models ${file_aa}

    error_exit $?
done

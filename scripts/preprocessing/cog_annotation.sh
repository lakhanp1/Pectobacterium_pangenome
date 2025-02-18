#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate eggnog
######################################################################

for id in $(cat data/reference_data/assembly_ids.txt); do
    file_aa="data/prokka_annotation/${id}/${id}.faa"

    emapper.py --cpu 20 -m mmseqs --dbmem --go_evidence all \
        -i ${file_aa} -o ${id} --output_dir data/prokka_annotation/${id} \
        --scratch_dir /local_scratch/${USER}/ \
        --temp_dir /local_scratch/${USER}/tmp

    error_exit $?
done

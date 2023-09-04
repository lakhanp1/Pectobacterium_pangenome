#!/usr/bin/env bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
# Task: Prophage detection using genomad and assessment using checkv

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh

scriptName=$(basename ${BASH_SOURCE[0]})

usage="USAGE:
-----------------------------------------------------------------------
bash ${scriptName} <sample_id>
sample_id : STRING sample identifier that will be used to get
            genome FASTA file
-----------------------------------------------------------------------
"

if [ $# -ne 1 ]; then
    printf "Error: Wrong args\n${usage}" >&2
    exit 1
fi

## Setup
PROJECT_DIR="$LUSTRE_HOME/projects/03_Pectobacterium"
ANALYSIS_DIR="$PROJECT_DIR/data/interproscan"

sampleId=$1
file_fna="$PROJECT_DIR/data/prokka_annotation/$sampleId/$sampleId.fna"
out_path="data/prophage_genomad/${sampleId}"

[ ! -d "data/prophage_genomad" ] && mkdir "data/prophage_genomad"

if [ ! -f ${file_fna} ]; then
    ls "${file_fna}"
    error_exit $?
fi

## run genomad to annotate prophage and then assessment with checkv
if [ ! -f ${out_path}/${sampleId}_summary/${sampleId}_virus_summary.tsv ]; then
    conda activate genomad
    process_start "geNamad on $sampleId"

    genomad end-to-end --cleanup -t 12 --enable-score-calibration \
        $file_fna $out_path $TOOLS_PATH/genomad_db

    error_exit $?

    n_vir=$(tail +2 $out_path/${sampleId}_summary/${sampleId}_virus_summary.tsv | wc -l)

    if [ ${n_vir} -gt 0 ]; then
        # run checkv
        process_start "checkv on $sampleId"

        checkv end_to_end --remove_tmp -t 12 \
            $out_path/${sampleId}_summary/${sampleId}_virus.fna \
            $out_path/${sampleId}_checkv

        error_exit $?
    fi

else
    echo "genomad output for ${sampleId} already exists"
fi

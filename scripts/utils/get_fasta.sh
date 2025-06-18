#!/usr/bin/env bash

## owner: Lakhansing Pardeshi, lakhansing.pardeshi@wur.nl
## get FASTA sequence from a genome for given locus

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate omics_py37

scriptName=$(basename ${BASH_SOURCE[0]})

usage="USAGE:
-------------------------------------------------------------
${scriptName} --region_file <FILE> 

region_file : FILE   Path to the genomic region file.
-------------------------------------------------------------
"

if [ $# -ne 2 ]; then
    printf "Error: argument missing\n${usage}" >&2
    error_exit 1
fi

## argument parsing
while [[ "$#" -gt 0 ]]; do
    case "$1" in
    --region_file)
        region_file="$2"
        shift # Shift past the argument value
        ;;
    # -o | --out)
    #     out_fasta="$2"
    #     shift # Shift past the argument value
    #     ;;
    -h | --help)
        echo "${usage}"
        exit 0
        ;;
    *)
        echo "Unknown option: $1"
        error_exit 1
        ;;
    esac
    shift # Shift past the current argument (option or value)
done

dir_path=$(dirname "${region_file}")
out_fasta="${region_file%.*}".fasta

true >"${out_fasta}" # Clear the output file

declare -A COL_INDEX

IFS=$'\t' read -r -a HEADER_ARRAY <"${region_file}"

# Map header names to their column index
for i in "${!HEADER_ARRAY[@]}"; do
    COL_INDEX["${HEADER_ARRAY[i]}"]="$i"
done

# echo "--- Header Mapping ---"
# for col_name in "${!COL_INDEX[@]}"; do
#     echo "Column '$col_name' is at index ${COL_INDEX[$col_name]}"
# done
# echo "----------------------"

tail -n +2 "${region_file}" |
    while IFS=$'\t' read -r -a DATA_ARRAY; do
        name="${DATA_ARRAY[${COL_INDEX[name]}]}"
        sampleId="${DATA_ARRAY[${COL_INDEX[sampleId]}]}"
        chr_name="${DATA_ARRAY[${COL_INDEX[chr_name]}]}"
        start="${DATA_ARRAY[${COL_INDEX[start]}]}"
        end="${DATA_ARRAY[${COL_INDEX[end]}]}"
        strand="${DATA_ARRAY[${COL_INDEX[strand]}]}"

        pos="${chr_name}:${start}-${end}"

        rc=""
        if [ ${strand} == "-" ]; then
            rc="--reverse-complement"
        fi

        printf "Extracting sequence for sample ${sampleId} at position ${pos} on strand ${strand}\n"

        seq=">${name}\n"
        seq+=$(samtools faidx ${rc} data/prokka_annotation/${sampleId}/${sampleId}.fna $pos | tail -n +2)

        error_exit $?

        printf "${seq}\n" >>${out_fasta}

    done

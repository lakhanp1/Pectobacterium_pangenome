#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

######################################################################
scriptName=$(basename $0)

usage="bash ${scriptName} <setId> <genomes>
setId       : STRING identifier for random set
genomes     : COMMA separated list of genome IDs in the pangenome to be
              used for optimal_grouping
db_name     : STRING name for pangenome database
"

if [ $# -ne 3 ]; then
    printf "Error: Require 2 args\n${usage}" >&2
    exit 1
fi

setId=$1
genomes=$2

######################################################################
## Setup
source scripts/utils/setup_analysis.sh ${3} .subset_opt_grp

if [ -z ${pan_db+x} ]; then
    echo "\$pan_db is unset"
    error_exit 1
fi

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_master

export PANTOOLS="$PANTOOLS_4_1"

busco_lineage=$(yq '.data.busco_db.lineage' project_config.yaml)

if [ -z ${pan_db+x} ]; then
    echo "\$pan_db is unset"
    error_exit 1
fi

if [ ! -d ${pan_db} ]; then
    cp -r $PANGENOME_DIR/backup/${3}.DB.backup_bsk ${pan_db}
else
    printf "Removing the grouping of the existing pan.db\n" >&2
    nice $PANTOOLS remove_grouping ${pan_db}
fi

## optimized grouping
process_start optimal_grouping
nice $PANTOOLS optimal_grouping -t 8 --include $genomes ${pan_db} ${pan_db}/busco/${busco_lineage}
error_exit $?

cp -r ${pan_db}/optimal_grouping $PROJECT_DIR/analysis/subset_optimal_group/optimal_grouping_${setId}

######################################################################


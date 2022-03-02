# Run InterProScan on the protein sequences

## Setup

``` bash
#!/bin/bash

# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
# Task: Build functional annotation using InterProScan

shopt -s expand_aliases
source ~/.bash_aliases

source /mnt/scratch/parde001/tools/miniconda3/etc/profile.d/conda.sh
conda activate interproscan

## Setup
PROJECT_DIR='/mnt/scratch/parde001/projects/03_Pectobacterium'
ANALYSIS_DIR="$PROJECT_DIR/analysis/04_pangenome_pecto_50g/interProScan"
pan_db="$PROJECT_DIR/analysis/04_pangenome_pecto_50g/pectobacterium.50g.DB"

```

## Generate functional annotation

``` bash

function run_interproscan(){
    ## replace * in the FASTA file
    sed -i 's/*//g' ${1}
    # sed '/>/{H;$!d} ; x ; s/^/\nSTART-->/ ; s/$/\n<--END/'

    file_outPrefix=`basename ${1} | sed 's/.fasta//'`.interproscan
    interproscan.sh -cpu 12 -goterms -f GFF3 -iprlookup -dp  --output-file-base $2/${file_outPrefix} -i $1
    error_exit $?
}

export -f run_interproscan

## Run InterProScan on all files using GNU parallel
ls ${pan_db}/proteins/*.fasta | parallel --keep-order --jobs 4 --halt now,fail=1 --results $PROJECT_DIR/logs/interproscan/{/.} --joblog $PROJECT_DIR/logs/interproscan/parallel.log run_interproscan {} ${ANALYSIS_DIR}

conda deactivate

```

## Explore `sed` to replace `*` in the FASTA records

``` bash
:'
>FEDMPDKE_04655_mRNA
MNKKLPLGLLLVLSFFLTGCDSEPSDADVSKAMQSFTEEMNKGAPSESARVVFNSAKKIACKDKDSDGGY
KCTVEYNAKIPFLGERTSTMELKFFKFEGKWKVSPSK
>FEDMPDKE_04654_mRNA
LFRHRIEAHTHTRGGGIQQVNGFVRQLAAGQVTAGQGDRGTHRIVGDVHAVMFGIAGFQAPEHQTRRVIV
RLVHQHPLEAALQRGIAFEVLLVLGPGGRRDGAQFAAGQRRF
>FEDMPDKE_04652_mRNA
VRELLSQYDFPGDDTPVVRGSALKALEGEAEWEAKIIELAEHLDSYIPEPERAIDKPFLLPIEDVFSISG
RGTVVTGRVERGIVKVGEEVEIVGIKDTAKSTCTGVEMFRKLLDEGRAGENVGVLLRGIKREEIERGQVL
AKPGSIKPHTQFESEVYILSKDEGGRHTPFFKGYRPQFYFRTTDVTGTIELPEGVEMVMPGDNIKMVVTL
IHPIAMDDGLRFAIREGGRTVGAGVVAKVIA
>FEDMPDKE_04651_mRNA
VFLLVGPSGVGKTETALALADVLYGGEKSLITINLSEYQEPHTVSQLKGSPPGYVGYGQGGILTEAVRKR
PYSVVLLDEVEKAHRDVMNLFYQVFDRGFMRDGEGREIDFRNTVILMTKKKTSF



p, Print the pattern space. 
h, hold function saves content from pattern space into hold space
g, get function retrieves content from hold space into pattern space
n, (next) If auto-print is not disabled, print the pattern space, then, regardless, replace the pattern space with the next line of input. If there is no more input then sed exits without processing any more commands. 
x, exchange function swaps the contents of pattern space and hold space
d, Delete the pattern space; immediately start next cycle. 

P, print function prints the first line of pattern space to standard output
H, hold function adds a \n character to hold space, and then appends the content from pattern space
G, get function adds a \n character to pattern space, and then appends the content from hold space
N, next function adds a \n character to pattern space and then appends the next line from the input stream

'

sed '/./{H;$!d} ; x ; s/^/\nSTART-->\n/ ; s/$/\n<--END/' input.txt

sed  'H; {/^>/!G;d}; x; s/^/START-->\n/ ; s/$/\n<--END/' seq.fa
```

# InterProScan installation

## Prepare the conda environment with the required packages

``` bash
mamba create -n interproscan python=3.8

nohup mamba install -y -c conda-forge openjdk=11.0.13 perl=5.32.1 >nohup.out 2>&1 &

mamba activate interproscan
```

## Download the InterProScan and follow the installation instructions

``` bash
python3 /path/to/initial_setup.py
```

## Some of the HMM index files are not generated properly. Hence run the following commands to regerate the HMM index

``` bash
cd  /path/to/InterProScan/folder

bin/hmmer/hmmer3/3.3/hmmpress -f data/gene3d/4.3.0/gene3d_main.hmm
bin/hmmer/hmmer3/3.3/hmmpress -f data/hamap/2020_05/hamap.hmm.lib
bin/hmmer/hmmer3/3.3/hmmpress -f data/panther/15.0/panther.hmm
bin/hmmer/hmmer3/3.3/hmmpress -f data/pfam/34.0/pfam_a.hmm
bin/hmmer/hmmer3/3.1b1/hmmpress -f data/sfld/4/sfld.hmm
bin/hmmer/hmmer3/3.1b1/hmmpress -f data/superfamily/1.75/hmmlib_1.75
bin/hmmer/hmmer3/3.3/hmmpress -f data/tigrfam/15.0/TIGRFAMs_HMM.LIB
```

## Test the InterProScan

``` bash
interproscan.sh -i /path/to/InterProScan/folder/test_all_appl.fasta -f tsv -dp

interproscan.sh -i /path/to/InterProScan/folder/test_all_appl.fasta -f tsv
```

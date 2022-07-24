# Test commands

``` bash

nohup mamba env create -n pantools_test -f conda.yaml &
git checkout master
git pull
git checkout add_gff_parser
git pull
mamba env create -n pantools-dev -f conda.yaml
conda activate pantools-dev
mvn package

java -jar target/pantools-3.3.4.jar


```

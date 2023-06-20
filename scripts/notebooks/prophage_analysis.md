# Prophage analysis in *Pectobacterium* pangenome

## Workflow

```mermaid

```

## Detect prophages in pangenome

### genomad prophage detection

```bash
# run genomad and checkv pipeline
nohup \
cat data/reference_data/assembly_ids.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/genomad/{} \
--joblog logs/genomad/parallel.log \
scripts/a_preprocessing/genomad_prophage_annotation.sh {} \
>>logs/genomad/nohup.out 2>&1 &

# count of prophage detected in each genome
for i in `cat data/reference_data/assembly_ids.txt`
do
  printf "$i\t"
  tail -n +2 data/prophage_genomad/$i/${i}_summary/${i}_virus_summary.tsv | wc -l
done | sort -nr -k2

```

### Post processing of prophage 

Process genomad output and separate prophages into independent FASTA files. Each
prophage is assigned a unique identifier with format `g_{genome}.vir_{#count}`.
Prophages are saved in independent FASTA files for downstream analysis, for e.g.
ANI clustering of prophages.

```bash
# post processing of prophage data
Rscript scripts/a_preprocessing/genomad_checkv_merge.R

# save prophages independently
bash scripts/a_preprocessing/prophage_genomes_extract.sh
```

Run ANI on prophage genomes to identify the representative prophages in the 
pangenome.

```bash
# list prophage FASTA files
ls data/prophage_genomad/phage_seqs/*.fna  > data/prophage_genomad/prophage_fasta.list

# run fastANI
fastANI --ql data/prophage_genomad/prophage_fasta.list \
--rl data/prophage_genomad/prophage_fasta.list \
--threads 30 --matrix --visualize --fragLen 500 \
-o data/prophage_genomad/prophage_ANI
```

## Analysis

### Cluster prophages

Use shared homology groups between two prophages to establish the child-parent
relationship between the prophages in pangenome.

```bash
Rscript scripts/c_analysis/prophage_DAG.R
```

```bash
Rscript scripts/c_analysis/prophage_DAG_viz.R
```

Process prophage ANI results and cluster prophages based on ANI.

```bash
Rscript scripts/c_analysis/prophage_ANI_summary.R
```

Merge the homology groups based prophage DAG with ANI to make final representative
prophages list.

```bash
Rscript scripts/c_analysis/prophages_representatives.R
```

### Visualize

Visualize the representative prophages across pangenome.

```bash
Rscript scripts/c_analysis/prophage_pangenome_viz.R
```

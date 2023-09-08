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

### Post processing of prophages 

Process genomad output and separate prophages into independent FASTA files. Each
prophage is assigned a unique identifier with format `g_{genome}.vir_{#count}`.
Prophages are saved in independent FASTA files for downstream analysis, for e.g.
ANI clustering of prophages.

```bash
# post processing of prophage data
Rscript scripts/preprocessing/genomad_checkv_merge.R

# save prophages independently
bash scripts/preprocessing/prophage_genomes_extract.sh
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

# use mash for distance estimation
mkdir data/prophage_genomad/mash

mash sketch -p 30 -s 15000 -k 17 -o data/prophage_genomad/mash/prophages_ref \
-l data/prophage_genomad/prophage_fasta.list >logs/prophage/mash.log 2>&1

mash info data/prophage_genomad/mash/prophages_ref.msh > \
data/prophage_genomad/mash/prophage_sketch_info.txt

mash dist -p 30 data/prophage_genomad/mash/prophages_ref.msh \
data/prophage_genomad/mash/prophages_ref.msh > \
data/prophage_genomad/mash/prophage_distance.tab

```

## Analysis

### Cluster prophages

Use shared homology groups between two prophages to establish the child-parent
relationship between the prophages in pangenome.

```bash
Rscript scripts/analysis/prophage_DAG.R
Rscript scripts/analysis/prophage_DAG_postprocess.R
```

Visualize the clusters of prophages as directed acyclic graphs.

```r
knitr::knit("scripts/analysis/prophage_DAG_viz.qmd")
```

Process prophage MASH and ANI results and visualize MASH similarity matrix for
the non-redundant prophages.

```bash
Rscript scripts/analysis/prophage_seq_similarity.qmd
```

### Visualize

Summarize the representative prophages in the pangenome.

```bash
Rscript scripts/analysis/prophage_representative_summary.qmd
```

Visualize the representative prophages across pangenome. Use MASH distance to
cluster prophages and ANI distance to cluster the genomes in the pangenome.

```bash
Rscript scripts/analysis/prophage_pangenome_viz.qmd
```

Visualize homology group sets for selected prophages as PAV heatmap across the
pangenome. Replace the code block in `scripts/analysis/homology_group_viz.R` 
script with the following.

```r
analysisName <- "vir_lineage_prophages"

outDir <- file.path(confs$analysis$prophages$dir)
outPrefix <- file.path(outDir, analysisName)

setIds <- c(
  "g_189.vir_2", "g_400.vir_2", "g_399.vir_2", "g_448.vir_2", "g_189.vir_1", 
  "g_3.vir_4", "g_438.vir_3", "g_350.vir_1",
  "g_172.vir_2", "g_393.vir_1"
)
```

Finally, run the script to generate the plot.

```bash
Rscript scripts/analysis/homology_group_viz.R
```

### carotovoricin analysis

#### Incomplete ctv cluster in *P. atrosepticum*

- Species: *P. atrosepticum*
- Genome assembly: GCF_000011605.1
- Pangenome genome: 396
- ECA3695-3742

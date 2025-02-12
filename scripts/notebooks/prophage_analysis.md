# Prophage analysis in *Pectobacterium* pangenome

## Detect prophages in genome 

### genomad prophage detection

```bash
# run genomad and checkv pipeline
nohup \
cat data/reference_data/assembly_ids.txt | \
parallel --jobs 6 --workdir $PWD --halt now,fail=1 \
--keep-order --results logs/genomad/{} \
--joblog logs/genomad/parallel.log \
scripts/preprocessing/genomad_prophage_annotation.sh {} \
>>logs/genomad/nohup.out 2>&1 &

# count of prophage detected in each genome
for i in `cat data/reference_data/assembly_ids.txt`
do
  printf "$i\t"
  tail -n +2 data/prophage_genomad/$i/${i}_summary/${i}_virus_summary.tsv | wc -l
done | sort -nr -k2

```

Check the checkV contamination.tsv file to see if there is any other contamination
than the expected `host,viral`, `viral,host` and `host,viral,host` ones.

```bash
for cn in `find -name '*contamination.tsv'`
do
cat $cn | awk '$9 != "NA" && NR != 1 && $9 != "host,viral" && $9 != "viral,host" && $9 != "host,viral,host"'
done
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

Output:
- yaml: `confs$data$prophages$files$data`
- file: `data/prophage_genomad/prophage_summary.tab`

Extract homology groups for prophages from the pangenome org.db object.

```bash
Rscript scripts/analysis/prophage_homology_groups.R
```

Output:
- yaml: `confs$analysis$prophages$preprocessing$files$raw_prophage_hg`
- file: `analysis/pangenome_v2/prophages/preprocessing/prophage_hgs.tab`

Run ANI and MASH on prophage genomes.

```bash
conda activate pantools_v4_3

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

### QC

Identify the fragmented prophages resulting because of the genome assembly
that can be merged into a complete prophage. Prepare a consolidated list of
prophages which will be used for the subsequent analysis.

```bash
Rscript scripts/analysis/prophage_merge.R
```

Output:
- yaml: `confs$analysis$prophages$preprocessing$files$consolidated`
- file: `analysis/prophages/preprocessing/consolidated_phages.tab`

Summarize the consolidated prophage list to generate various statistics and 
generate summary plots. Prophages that are fragmented or smaller than the length
cutoff are removed and the filtered list of prophages is saved for downstream
analysis.

```bash
Rscript scripts/analysis/prophage_summary.R
```

Output:
- yaml: `confs$analysis$prophages$preprocessing$files$prophage_pool`
- File: `analysis/pangenome_v2/prophages/preprocessing/prophages_pool.tab`

Process prophage MASH and ANI results and visualize MASH similarity matrix for
the filtered prophages.

```bash
quarto render scripts/analysis/prophage_seq_similarity.qmd --execute-dir ./
```

## Cluster prophages

Use shared homology groups between two prophages to establish the child-parent
relationship between the prophages in pangenome. First, all vs all syntenic
Jaccard similarity is calculated between the prophages. This information is later
used to cluster the prophages using two different approaches.

```bash
Rscript scripts/analysis/prophage_pairwise_similarity.R
```

### Approach 1

A heuristic DAG construcion is extracted from the output of the script 
`scripts/analysis/prophage_merge.R`. This output includes a best matching parent
prophage for each prophage. This is further processed to build a DAG and
visualize the clusters.

Post processing of the child-parent relationships to add additional metadata
and save a DAG in tabular format:

```bash
Rscript scripts/analysis/prophage_DAG_postprocess.R
```

Visualize the clusters of prophages as directed acyclic graphs:

```r
knitr::knit("scripts/analysis/prophage_DAG_viz.qmd")
```

### Approach 2

We first tried to use the Affinity propagation (AP) clustering to group the
prophages into clusters and have a set of cluster representatives. However,
we encountered some clusters generated by the AP clustering where the cluster
members did not even share a homology group. Therefore, we decided not to use
AP clustering.

```r
knitr::knit("scripts/analysis/prophage_AP_clustering.qmd")
```

Cluster the all-vs-all syntenic Jaccard distance matrix using "complete" linkage
hierarchical clustering and cut the cluster using carotovoricin cluster as a
reference.

```r
knitr::knit("scripts/analysis/prophage_clustering.qmd")
```

Output:
- yaml: `confs$analysis$prophages$files$clusters`
- file: `analysis/pangenome_v2/prophages/prophages_clusters.tab`

### Summarize prophage clusters

Summarize the representative prophages in the pangenome.

```bash
Rscript scripts/analysis/prophage_cluster_summary.qmd
```

Visualize the representative prophages across pangenome. Use MASH distance to
cluster prophages and ANI distance to cluster the genomes in the pangenome.

:::{.callout-caution}
Need a fix. Do not run!

```bash {eval=FALSE}
quarto render scripts/analysis/prophage_pangenome_viz.qmd --execute-dir ./
```

:::

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

## Individual prophage cluster analysis

Script `scripts/analysis/clustersmap_data_prophages.R` is used to combine the
prophages, homology groups and phylogeny to generate input for `clustermap.js`
visualization tool and generate a report. Below, some configuration is provided
for individual prophage clusters.

### phage_grp_46: highly conserved in all but 5 genomes from n23 clade

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22429500,hg_22429499 \
--inner_region --genomes "g_63,g_190,g_193,g_411,g_415,g_406,g_194,g_263"
```

```r
cluster_title <- "phage_grp_46"
outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$prophages$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- "phage_grp_46"
other_regions <- character(0)

subSample <- TRUE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- TRUE 

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list(
  g_190_reg = list(chr = "NZ_CP059960.1", start = 1997124, end = 2000186, genomeId = "g_190"),
  g_193_reg = list(chr = "NZ_CP059957.1", start = 5005351, end = 5008413, genomeId = "g_193"),
  g_411_reg = list(chr = "NAK682_contig_2_consensus", start = 461100, end = 464162, genomeId = "g_411"),
  g_415_reg = list(chr = "NAK701_contig_5_consensus", start = 415799, end = 418861, genomeId = "g_415"),
  g_63_reg = list(chr = "NZ_CP092039.1", start = 4288654, end = 4291716, genomeId = "g_63"),
  g_406_reg = list(chr = "NAK641_contig_10_consensus", start = 469664, end = 474521, genomeId = "g_406"),
  g_194_reg = list(chr = "NZ_CP059956.1", start = 4163454, end = 4167352, genomeId = "g_194")
  # g_263_reg = list(chr = "NZ_CP047495.1", start = 956237, end = 956907, genomeId = "g_263")
)
```

### phage_grp_45: prophage present in all BL-causing isolates

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22429319,hg_22429327 \
--genomes "g_406,g_194"
```

```r
cluster_title <- "grp_45_mosaicsm"
outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$prophages$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- NA
other_regions <- c("g_302.vir_1", "g_399.vir_1", "g_400.vir_2")

subSample <- TRUE
cutHeight <- 0.5
addFlankingRegions <- TRUE
flankingRegion <- 7000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- FALSE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list(
  g_406_reg = list(chr = "NAK641_contig_10_consensus", start = 671040, end = 674984, genomeId = "g_406"),
  g_194_reg = list(chr = "NZ_CP059956.1", start = 3940785, end = 3958060, genomeId = "g_194")
)
```

### phage_grp_71

```r
cluster_title <- "phage_grp_71"
outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$prophages$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- "phage_grp_71"
other_regions <- c()

subSample <- FALSE
cutHeight <- 0.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "hg_pav" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- FALSE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list()

```

### phage_grp_107

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427004,hg_22426993 \
--inner_region --genomes "g_63,g_194"
```

```r
cluster_title <- "phage_grp_107"
outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz", sep = "")
hg_color_categories <- confs$analysis$prophages$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- "phage_grp_107"
other_regions <- character(0)

subSample <- FALSE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- TRUE 

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list(
  g_63_reg = list(chr = "NZ_CP092039.1", start = 2572138, end = 2573569, genomeId = "g_63"),
  g_194_reg = list(chr = "NZ_CP059956.1", start = 2620331, end = 2621761, genomeId = "g_194")
)

```

### BL-causing clade specific signatures

```r
grpToView <- "Pbr_FN_specific"
subSample <- FALSE
cutHeight <- 0.5
addFlankingRegions <- FALSE
flankingRegion <- 5000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "default"  # completeness, host, hg_pav, cluster_mash, default


# optionally, a custom region list can be provided to generate the plot
grp <- list(
  phage_grp = grpToView,
  members = c("g_302.vir_1", "g_400.vir_2", "g_409.vir_1", "g_302.vir_2")
)
```

### phage_grp_36

```r
grpToView <- "phage_grp_36"
subSample <- FALSE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "host" # completeness, host, hg_pav, cluster_mash
```

### phage_grp_50

```r
grpToView <- "phage_grp_50"
subSample <- FALSE
addFlankingRegions <- TRUE 
flankingRegion <- 10000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "cluster_mash"  # completeness, host, hg_pav, cluster_mash
```

### phage_grp_21

```r
grpToView <- "phage_grp_21"
subSample <- FALSE
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "cluster_mash"  # completeness, host, hg_pav, cluster_mash
```

### T3SS missing from 4 P. brasiliense genomes

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22430875,hg_22430876 --genomes g_187,g_368,g_182,g_185,g_236,g_177
```

```r
grpToView <- "t3ss_regions"
subSample <- FALSE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 6000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "host"  # completeness, host, hg_pav, cluster_mash

# a vector of prophage identifiers that will be included in clustermap plot
appendPhages <- c()

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list(
  g_177_reg  = list(genomeId = "g_177", chr = "NZ_JACGEP010000027.1", start = 56850, end = 59337),
  g_182_reg = list(genomeId = "g_182", chr = "NZ_JACGZZ010000068.1", start = 116773, end = 119260),
  g_185_reg = list(genomeId = "g_185", chr = "NZ_JACGEN010000027.1", start = 117747, end = 120234),
  g_236_reg = list(genomeId = "g_236", chr = "NZ_JACDSF010000026.1", start = 2275, end = 4762),
  g_187_reg = list(genomeId = "g_187", chr = "NZ_CP059963.1", start = 2603432, end = 2605919),
  g_368_reg = list(genomeId = "g_368", chr = "NZ_JUJT01000001.1", start = 566617, end = 569104)
)

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- TRUE
```

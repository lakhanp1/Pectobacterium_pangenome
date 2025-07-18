
# Carotovoricin (CTV) cluster or phage_grp_1 data analysis

## Extract carotovoricin data from prophage clustering output

Prophage cluster `phage_grp_1` belongs to the carotovoricin. Hence, prepare
a table with prophage identifiers, their genomic coordinates, homology group
signatures and other relevant information for all the prophages in the
`phage_grp_1` cluster.

```bash
Rscript scripts/analysis/ctv_prophages_table.R
```

Manually check the output file `analysis/pangenome_v2/carotovoricin/data/ctv_prophages_data.tsv`
and verify carotovoricin presence/absence status in the genomes. The manually
curated table is saved as `analysis/pangenome_v2/carotovoricin/data/ctv_genome_table.tsv`.

## Extract all homology groups representing CTV region across the pangenome

Genes and their respective homology groups flanking CTV:
- tolC_2: hg_22427641
- ybiB: hg_22427603

Extract genomic coordinates for CTV region together with the haplotypes

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427641,hg_22427603 \
--inner_region --haplotypes \
--out analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab
```

Process the output file from the previous step to build CTV homology group
table and then manually annotate the homology groups with broad functionl
categories based on the keywords in GO, COG and PFAM annotations.

CTV present in 428 genomes.

- Intact CTV along with flanking genes *tolC_2* and *ybiB*: 365 genomes
- Present at the end of the contig: 31 genomes
- Fragmented into multiple contigs: 32 genomes

CTV is missing or deleted in the following 26 genomes:

- Complete deletion:

  - *P. brasiliense*: g_177, g_182, g_185, g_236
  - *P. carotovorum*: g_15 (now renamed to *P. araliae*)
  - *P. cacticida*: g_451

- Partial deletion:

  - *P. atrosepticum* (all 16 genomes): g_124, g_200, g_27, g_283, g_311, g_322, g_336, g_348, g_382,
    g_385, g_389, g_392, g_396, g_50, g_51, g_55
  - *P. betavasculorum*: g_383, g_386
  - *P. versatile*: g_313 (deletion of 5 genes including tail sheath, tail tube
    and tape measure genes). Check g_313.vir_2

- Insertion of another prophage:
  - *P. brasiliense*: g_149 CTV cluster has integration of another prophage and   hence it is disrupted. It was detected as g_149.vir_1	and was part of phage_grp_89 with 3 members.

## Clustermap visualization of CTV clusters

### CTV across pangenome

Include *P. atrosepticum* (g_385), *P. betavasculorum* (g_386) and *P. cacticida*
(g_451) genomes that lack CTV cluster.

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427603,hg_22427641 \
--inner_region --genomes "g_385,g_386,g_451"
```

Update the `scripts/analysis/clustersmap_data_prophages.R` script as below to
visualize the CTV cluster in *P. brasiliense*.

```r
cluster_title <- "ctv_typeStrains"
outDir <- paste(confs$analysis$ctv$path, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- NA
other_regions <- c(
  "g_345.vir_1", "g_446.vir_4", "g_66.vir_3", "g_222.vir_2", "g_365.vir_3",
  "g_8.vir_2", "g_38.vir_2", "g_273.vir_2", "g_259.vir_4",
  "g_305.vir_1", "g_378.vir_6", "g_428.vir_1", "g_248.vir_1", "g_449.vir_1",
  "g_54.vir_1", "g_116.vir_3", "g_423.vir_3", "g_375.vir_2", "g_381.vir_2"
)

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
  g_385_reg = list(chr = "NZ_JQHK01000003.1", start = 203963, end = 207120, genomeId = "g_385"),
  g_386_reg = list(chr = "NZ_JQHM01000001.1", start = 553213, end = 555615, genomeId = "g_386"),
  g_451_reg = list(chr = "Contig_2_668.636", start = 191452, end = 191490, genomeId = "g_451")
)

```

### carotovoricin cluster in *P. brasiliense* isolates

```r
cluster_title <- "ctv_pbr"
outDir <- paste(confs$analysis$ctv$path, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- "phage_grp_1"
other_regions <- c()

subSample <- TRUE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- FALSE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list()

regionClusters <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$clusters)
)

# optional filter
regionClusters %<>%
  dplyr::filter(
    SpeciesName == "P. brasiliense",
    !prophage_id %in% c("g_408.vir_3", "g_403.vir_3", "g_399.vir_3")
  )

```

### Visualize CTV deletion

Extract genomic coordinates for the species/genomes in which CTV is deleted.

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427603,hg_22427641 \
--inner_region --genomes "g_149,g_177,g_182,g_185,g_236,g_15,g_313,g_385,g_386"
```

Update the `scripts/analysis/clustersmap_data_prophages.R` script as below to
visualize the CTV clusters.

```r
cluster_title <- "ctv_deletion"
outDir <- paste(confs$analysis$ctv$path, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- NA
other_regions <- c("g_345.vir_1", "g_66.vir_3", "g_8.vir_2", "g_381.vir_2")

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
  g_177_reg = list(chr = "NZ_JACGEP010000002.1", start = 102352, end = 103187, genomeId = "g_177"),
  g_182_reg = list(chr = "NZ_JACGZZ010000050.1", start = 81626, end = 82461, genomeId = "g_182"),
  g_185_reg = list(chr = "NZ_JACGEN010000006.1", start = 81658, end = 82493, genomeId = "g_185"),
  g_236_reg = list(chr = "NZ_JACDSF010000027.1", start = 89639, end = 90475, genomeId = "g_236"),
  g_149_reg = list(chr = "NZ_JACGFD010000001.1", start = 2861612, end = 2913419, genomeId = "g_149"),
  g_385_reg = list(chr = "NZ_JQHK01000003.1", start = 203963, end = 207120, genomeId = "g_385"),
  g_386_reg = list(chr = "NZ_JQHM01000001.1", start = 553213, end = 555615, genomeId = "g_386"),
  g_451_reg = list(chr = "Contig_2_668.636", start = 191452, end = 191490, genomeId = "g_451"),
  g_313_reg = list(chr = "NZ_CP024842.1", start = 2917505, end = 2929851, genomeId = "g_313"),
  g_15_reg = list(chr = "BRCR01000012.1", start = 138656, end = 139632, genomeId = "g_15")
)

```

Prophage clusters found in the ctv-lacking Pbr: phage_grp_30, phage_grp_6, phage_grp_29

### A prophage integrated within the CTV in *P. brasiliense* genome `g_149`

<!--  -->
Genome g_149 CTV cluster has integration of another prophage and hence it is 
disrupted. It was detected as g_149.vir_1	and was part of phage_grp_89 which has
3 members, `r c("g_149.vir_1", "g_116.vir_1", "g_46.vir_2")`

```r
cluster_title <- "prophage_in_ctv"
outDir <- paste(confs$analysis$ctv$path, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$ctv$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- NA
other_regions <- c("g_368.vir_3", "g_149.vir_1", "g_116.vir_1", "g_46.vir_2")

subSample <- FALSE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "default" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- FALSE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list()

```

### CTV in *P. versatile*

Update the `scripts/analysis/clustersmap_data_prophages.R` script as below to
visualize the CTV cluster in *P. versatile*.

```r
cluster_title <- "ctv_pv"
outDir <- paste(confs$analysis$ctv$path, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$ctv$data$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- "phage_grp_1"
other_regions <- c()

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
customRegions <- list()

regionClusters <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$clusters)
)

# optional filter
regionClusters %<>%
  dplyr::filter(
    SpeciesName == "P. versatile"
  )
```

### CTV cluster in *P. versatile* collected from France

Update the `scripts/analysis/clustersmap_data_prophages.R` script as below to
visualize the CTV clusters.

```r
grpToView <- "ctv_pvs_fr"
subSample <- TRUE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

grp <- list(
  phage_grp = grpToView,
  members = dplyr::filter(
    regionClusters,
    SpeciesName == "P. versatile", nFragments == 1, phage_grp == "phage_grp_1",
    geo_loc_country == "France"
  ) %>%
    dplyr::pull(prophage_id)
)
```

## Variable regions in CTV cluster

### Nitric oxide deoxygenase (hg_22427623)

```bash
# all genomes that have CTV
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427625,hg_22427622 \
--haplotypes --max_genes 10 \
--genomes_file analysis/pangenome_v2/carotovoricin/data/genomes_ctv.txt \
--out analysis/pangenome_v2/carotovoricin/ctv_nitric_oxide/hg_regions_all_genomes.tab

# only for intact CTV genomes
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427625,hg_22427622 \
--haplotypes --max_genes 10 \
--genomes_file analysis/pangenome_v2/carotovoricin/data/genomes_ctv_intact.txt \
--out analysis/pangenome_v2/carotovoricin/ctv_nitric_oxide/hg_regions.tab
```

### Tape measure protein

Tape measure genes is located between the homology groups hg_22427622 and hg_22427616.
This region is also represented by multiple combinations of homology groups.

```bash
# all genomes that have CTV
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427622,hg_22427616 \
--haplotypes --inner_region --overlapping --max_genes 10 \
--genomes_file analysis/pangenome_v2/carotovoricin/data/genomes_ctv.txt \
--out analysis/pangenome_v2/carotovoricin/ctv_tape_measure/hg_regions_all_genomes.tab

# only for intact CTV genomes
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427622,hg_22427616 \
--haplotypes --inner_region --overlapping --max_genes 10 \
--genomes_file analysis/pangenome_v2/carotovoricin/data/genomes_ctv_intact.txt \
--out analysis/pangenome_v2/carotovoricin/ctv_tape_measure/hg_regions.tab
```

### Tail fiber region

```bash
# all genomes that have CTV
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--inner_region --haplotypes \
--genomes_file analysis/pangenome_v2/carotovoricin/data/genomes_ctv.txt \
--out analysis/pangenome_v2/carotovoricin/ctv_tfl/hg_regions_all_genomes.tab

# only for intact CTV genomes
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--inner_region --haplotypes \
--genomes_file analysis/pangenome_v2/carotovoricin/data/genomes_ctv_intact.txt \
--out analysis/pangenome_v2/carotovoricin/ctv_tfl/hg_regions.tab
```

### CTV conserved loci

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427640,hg_22427604 \
--genomes_file analysis/pangenome_v2/carotovoricin/data/genomes_ctv_intact.txt \
--out analysis/pangenome_v2/carotovoricin/ctv_conserved/hg_regions.tab
```

### Haplotype analysis and visualization for the tail fiber locus

```bash
quarto render scripts/analysis/ctv_tfl_haplotype_analysis.qmd --execute-dir $PWD
```

Optionally, tandem homology groups can be mapped to the genomes and visualized
on the phylogenetic tree. Need to change the `hgSets` variable in the script
`scripts/analysis/HG_tandem_match_viz.qmd` to include the homology groups
of interest as a `list()` object. 

For example,

```r
hgSets <- list(
    hg_comb_2 = c("hg_22427604", "hg_22426818", "hg_22427603"),
    hg_comb_3 = c("hg_22427604", "hg_22426828", "hg_22431676", "hg_22431685", "hg_22427603"),
    hg_comb_5 = c("hg_22427604", "hg_22426817", "hg_22426799", "hg_22427603")
)
```

```bash
quarto render scripts/analysis/HG_tandem_match_viz.qmd --execute-dir $PWD
```

### Summary of carotovoricin cluster

```bash
quarto render scripts/analysis/ctv_summary.qmd --execute-dir $PWD
```

## Inter-species horizontal gene transfer of CTV

Clustermap visualization for some manually selected candidates. Change the code
blocks in the script `scripts/analysis/clustersmap_data_prophages.R` and run.

```r
grpToView <- "ctv_hgt"
subSample <- FALSE 
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "hg_pav" # completeness, host, hg_pav, cluster_mash


# optionally, a custom region list can be provided to generate the plot
grp <- list(
  phage_grp = grpToView,
  members = c(
    "g_145.vir_1", "g_194.vir_1", "g_429.vir_1", "g_442.vir_1", "g_421.vir_1", 
    "g_150.vir_4", "g_447.vir_1", "g_434.vir_4",
    "g_221.vir_3", "g_53.vir_3", "g_106.vir_2", "g_57.vir_1", "g_125.vir_1"
  )
)
```

Run script `scripts/utils/HG_range_coordinates.R` to extract genomic coordinates
for CTV region in selected genomes.

Complete CTV region along with 3 flanking genes:

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427643,hg_22427599 \
--inner_region --out analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab 

# create GFF3 files for CTV regions of interest
Rscript scripts/utils/HGs_gff_from_regions.R \
--regions analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab \
--genomes g_145,g_194,g_429,g_442,g_421,g_150,g_447,g_434,g_53,g_106,g_57,g_221,g_125

# create CTV regions GFF3 files for one species representative each
Rscript scripts/utils/HGs_gff_from_regions.R \
--regions analysis/pangenome_v2/carotovoricin/ctv_region/hg_regions.tab \
--genomes g_183,g_423,g_277,g_375,g_116,g_449,g_446,g_378,g_296,g_428,g_259,g_66,g_337,g_442,g_305,g_381,g_8,g_273,g_38
```

### DNA sequence distance comparison for conserved and TFL regions

#### Extract DNA sequence for the regions

```bash
bash scripts/utils/get_fasta.sh --region_file analysis/pangenome_v2/carotovoricin/ctv_tfl/hg_regions.tab
bash scripts/utils/get_fasta.sh --region_file analysis/pangenome_v2/carotovoricin/ctv_conserved/hg_regions.tab
bash scripts/utils/get_fasta.sh --region_file analysis/pangenome_v2/carotovoricin/ctv_tape_measure/hg_regions.tab

bash scripts/utils/get_fasta.sh --separate --region_file analysis/pangenome_v2/carotovoricin/ctv_tfl/hg_regions.tab
bash scripts/utils/get_fasta.sh --separate --region_file analysis/pangenome_v2/carotovoricin/ctv_conserved/hg_regions.tab
bash scripts/utils/get_fasta.sh --separate --region_file analysis/pangenome_v2/carotovoricin/ctv_tape_measure/hg_regions.tab
```

#### Calculate the distance using `mash`

```bash
cd analysis/pangenome_v2/carotovoricin/ctv_tfl
mash dist -p 8 -k 12 -s 2000 -i -S 124 hg_regions.fasta hg_regions.fasta > distance_mash.tab

cd analysis/pangenome_v2/carotovoricin/ctv_conserved
mash dist -p 8 -k 12 -s 2000 -i -S 124 hg_regions.fasta hg_regions.fasta > distance_mash.tab

cd analysis/pangenome_v2/carotovoricin/ctv_tape_measure
mash dist -p 8 -k 12 -s 2000 -i -S 124 hg_regions.fasta hg_regions.fasta > distance_mash.tab

grep -e 'g_(53|106|57|221|125).*g_(53|106|57|221|125)' ctv_*/distance_mash.tab
grep -e 'g_(150|447|434).*g_(150|447|434)' ctv_*/distance_mash.tab
grep -e 'g_(53|106|57).*g_(53|106|57)' ctv_*/distance_mash.tab

```

#### Calculate Jaccard distance matrix using `sourmash`

MASH distance is low for haplotype pairs where one is subset of the other, eg,
A-B-C-D and A-B. Therefore, another tool was tried to calculate the Jaccard
index between sequences.

```bash
mamba activate smash

file_fa="ctv_conserved/hg_regions.fasta"
# file_fa="ctv_tfl/hg_regions.fasta"

file_sig="${file_fa%/*}"/sourmash_signatures.zip
file_dist="${file_fa%/*}/"distance_sourmash.csv

# # sketch: scaled=5
# sourmash sketch dna -f -p k=7,k=9,k=11,scaled=5,abund,seed=124 --singleton \
# -o ${file_sig} ${file_fa}

# sketch: num=2000
sourmash sketch dna -f -p k=7,k=9,k=11,k=13,num=2000,abund,seed=124 --singleton \
-o "${file_sig}" "${file_fa}"

sourmash compare --ignore-abundance --distance-matrix -p 12 -k 13 --dna \
--csv "${file_dist}" "${file_sig}"

```

#### Calculate Jaccard distance using `dashing`

```bash
mamba activate smash

cd analysis/pangenome_v2/carotovoricin

# dir_fa="ctv_conserved/region_fasta"
dir_fa="ctv_tfl/region_fasta"

file_dist="${dir_fa%/*}/"distance_dashing.txt

dashing cmp -k 11 --seed 124 --nthreads 12 --full-mash-dist  \
-O "${file_dist}" "${dir_fa}"/*.fasta

zip -m -r "${dir_fa}".zip "${dir_fa}"
```

#### Process the sequence distance comparison results

```bash
quarto render scripts/analysis/ctv_tfl_hgt.qmd --execute-dir $PWD
```

### Tail fiber HGT demonstration

Extract sequences for some example tail fiber loci haplotypes that are present
in multiple *Pectobacterium* species.

```bash
genomeSet="g_279,g_425,g_299,g_377,g_100,g_106,g_331,g_249,g_125,g_221,g_53,g_395,g_108,g_444,g_160,g_352"
setName="genomeset_1"

# tail fiber locus
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--genomes "${genomeSet}" \
--out "analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_tfl_flanking.tab"

Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--genomes "${genomeSet}" \
--inner_region --out "analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_tfl.tab"

# conserved CTV locus
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427640,hg_22427604 \
--genomes "${genomeSet}" \
--out "analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_conserved.tab"

# upstream core
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427642,hg_22427641 \
--genomes "${genomeSet}" \
--out "analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_upstream_core.tab"
```

Extract DNA sequence for these regions.

```bash
for f in analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_*.tab
do
  bash scripts/utils/get_fasta.sh --region_file ${f}
done
```

#### Perform MSA and generate phylogenetic trees

```bash
mamba activate pantools_4_3_3

nohup mafft --reorder --allowshift --unalignlevel 0.2 --leavegappyregion \
--maxiterate 0 --globalpair \
"analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_tfl.fasta" \
> "analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_tfl.msa.fasta" 2>nohup.out &

nohup mafft --reorder --allowshift --unalignlevel 0.2 --leavegappyregion \
--maxiterate 0 --globalpair \
"analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_conserved.fasta" \
> "analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_conserved.msa.fasta" 2>nohup.out &

nohup mafft --reorder --allowshift --unalignlevel 0.2 --leavegappyregion \
--maxiterate 0 --globalpair \
"analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_upstream_core.fasta" \
> "analysis/pangenome_v2/carotovoricin/tfl_hgt/${setName}_upstream_core.msa.fasta" 2>nohup.out &
```

Generate a maximum-likelihood phylogenetic tree for MSAs

```bash
cd analysis/pangenome_v2/carotovoricin/tfl_hgt
mkdir iqtree

nohup nice iqtree2 -T 40 -s "${setName}_tfl.msa.fasta" -B 1000 \
--prefix "iqtree/${setName}_tfl.msa.fasta" >> iqtree_tree.log 2>&1 &

nohup nice iqtree2 -T 40 -s "${setName}_conserved.msa.fasta" -B 1000 \
--prefix "iqtree/${setName}_conserved.msa.fasta" >> iqtree_tree2.log 2>&1 &

nohup nice iqtree2 -T 40 -s "${setName}_upstream_core.msa.fasta" -B 1000 \
--prefix "iqtree/${setName}_upstream_core.msa.fasta" >> iqtree_tree3.log 2>&1 &

```

## TFL region Mauve analysis

Follow the steps below to perform Mauve analysis for genomic regions of interest.

1. Define the genome subset for the analysis.

For example:

```bash
subsetName="tfl_pecto"
genomeSet="g_302,g_364,g_337,g_439,g_403,g_399,g_408,g_175,g_138,g_368,g_345,g_366,g_308,g_191,g_173,g_155,g_166,g_299,g_438,g_263,g_43,g_391"
outDir="analysis/pangenome_v2/carotovoricin/mauve/${subsetName}"
```

2. Extract the genomic coordinates for a region based on the homology group range.

The TFL region is flanked by the homology groups `hg_22427604` and `hg_22427603`.
These homology groups will be used to extract the genomic coordinates for the
region of interest. We include the region encoding homology groups and therefore
do not use the `--inner_region` option.

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--genomes "${genomeSet}" --out "${outDir}/hg_regions.tab"
```

3. Create GFF3 files for these regions.

```bash
Rscript scripts/utils/HGs_gff_from_regions.R --regions "${outDir}/hg_regions.tab"
```

4. Extract DNA sequence for these regions.

```bash
bash scripts/utils/get_fasta.sh --region_file "${outDir}/hg_regions.tab" --separate
```

5. Generate GenBank files from the `fasta` sequences and `gff3` annotations.

Mauve aligns sequences provided in `fasta` and `genbank` format. However, the
sequence annotation is shown only if the alignment input was in `genbank`
format. Therefore, we need to combine the homology group annotation in `gff3`
format with the `fasta` sequences to generate a `genbank` formatted files. To
do this, Emboss tool `seqret` was used.

```bash
conda activate omics_py37

cd "${outDir}"

mkdir gbk_files

for vir in `cut -f 1 hg_regions.tab | tail -n +2`
do
  seqret -sequence region_fasta/${vir}.fasta -feature -fformat1 gff3 \
  -fopenfile1 gff3/${vir}.gff3 -osformat2 genbank -osextension2 gbk \
  -osname2 ${vir} -osdirectory2 gbk_files -auto
done
```

6. Correct GenBank file as per Mauve requirement

```bash
sed -i.bak -r -e 's/(^\s+CDS\s+(complement\()?)[^:]+:([[:digit:]]+\.\.[[:digit:]]+\)?).*/\1\3/' \
  -e 's/\/note="\*([^:]+): /\/\1="/' gbk_files/*.gbk
```

7. Run `progressiveMauve` to align the sequences.

```bash
mauve_out="${subsetName}"

progressiveMauve --output="${mauve_out}.mauve.xmfa" \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
gbk_files/*.gbk > "${mauve_out}".log 2>&1
```

8. Analyze the `xmfa` output files using Mauve GUI.

### Specific Mauve analysis

#### TFL haplotype variation in *P. brasiliense* strains

Extract the region between the two homology groups, `hg_22427604` and `hg_22427603`
for the carotovoricin cluster. Additionally, generate `GFF3` files with the
homology group, COG, PFAM and other metadata information to visualize.

```bash
subsetName="tfl_pbr"
genomeSet="g_302,g_364,g_337,g_439,g_403,g_399,g_408,g_175,g_138,g_368,g_345,g_366,g_308,g_191,g_173,g_155,g_166,g_299,g_438,g_263,g_43,g_391"
```

Specific use cases

```bash
# genomes where Ein gene present in cluster
mauve_out="${subsetName}_ein"

progressiveMauve --output="${mauve_out}.mauve.xmfa" \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
gbk_files/g_345.P_brasiliense.gbk gbk_files/g_366.P_brasiliense.gbk \
gbk_files/g_308.P_brasiliense.gbk gbk_files/g_191.P_brasiliense.gbk \
gbk_files/g_173.P_brasiliense.gbk gbk_files/g_155.P_brasiliense.gbk \
gbk_files/g_166.P_brasiliense.gbk  > "${mauve_out}".log 2>&1

# represetative genomes
mauve_out="${subsetName}_rep"

progressiveMauve --output="${mauve_out}.mauve.xmfa" \
--backbone-output="${mauve_out}.mauve.backbone" \
--output-guide-tree="${mauve_out}.mauve.guide_tree.newick" \
gbk_files/g_302.P_brasiliense.gbk gbk_files/g_337.P_brasiliense.gbk \
gbk_files/g_175.P_brasiliense.gbk gbk_files/g_173.P_brasiliense.gbk \
gbk_files/g_299.P_brasiliense.gbk gbk_files/g_438.P_brasiliense.gbk \
gbk_files/g_391.P_brasiliense.gbk  > "${mauve_out}".log 2>&1

```

#### Variation of TFL haplotype with *ein* across *Pectobacterium* species

```bash
subsetName="tfl_ein_pecto"
genomeSet="g_148,g_110,g_279,g_428,g_53,g_160,g_121,g_218,g_421,g_259,g_206,g_14,g_247,g_344,g_333,g_305,g_351,g_31"
outDir="analysis/pangenome_v2/carotovoricin/mauve/${subsetName}"
```

#### Variation of TFL haplotype (without *ein*) across *Pectobacterium* species

```bash
subsetName="tfl_ein_del_pecto"
genomeSet=""
outDir="analysis/pangenome_v2/carotovoricin/mauve/${subsetName}"
```

#### Mauve analysis for TFL HGT example genomes

```bash
subsetName="tfl_hgt"
genomeSet="g_279,g_425,g_299,g_133,g_377,g_100,g_106,g_331,g_249,g_125,g_221,g_53,g_395,g_108,g_444,g_160,g_352"
outDir="analysis/pangenome_v2/carotovoricin/mauve/${subsetName}"
```

Here, only internal region between the two homology groups is extracted. Therefore,
the `--inner_region` option is used.

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22427604,hg_22427603 \
--inner_region --genomes "${genomeSet}" --out "${outDir}/hg_regions.tab"
```

Rest of the steps are similar to the ones described above.

### Perform MSA using `MAFFT`

```bash
conda activate pantools_v4_3
cd analysis/pangenome_v2/carotovoricin/mauve/tfl_pbr/

mafft --globalpair --quiet --maxiterate 1000 --treeout ctv_pbr.variable_regions.fasta
```

Smash++ pairwise sequence comparison

```bash
conda activate omics_py37

cd analysis/pangenome_v2/carotovoricin/mauve/tfl_pbr/region_fasta

function smashpp_compare
{
  smashpp -r $1 -t $2
  smashpp viz -o $(basename $1 .fasta)$(basename $2 .fasta)".svg" "${1}.${2}.pos"
}

export -f smashpp_compare

smashpp_compare g_345.vir_1.fasta g_345.vir_1.fasta 
smashpp_compare g_345.vir_1.fasta g_366.vir_3.fasta 
smashpp_compare g_345.vir_1.fasta g_302.vir_3.fasta
smashpp_compare g_345.vir_1.fasta g_403.vir_3.fasta
smashpp_compare g_345.vir_1.fasta g_173.vir_2.fasta
smashpp_compare g_345.vir_1.fasta g_138.vir_2.fasta
smashpp_compare g_345.vir_1.fasta g_263.vir_2.fasta
smashpp_compare g_345.vir_1.fasta g_155.vir_1.fasta

```

Left inverted repeat for Ein: CTCCCGCAAACCTCGGTTTTGGGGAC (CTCCCGCAAACCTCGGTTT)  
Left inverted repeat for Ein(rev-com): GTCCCCAAAACCGAGGTTTGCGGGAG (AAACCGAGGTTTGCGGGAG, AAACCGAGGTTTGCG)  

Right inverted repeat for Ein: TTCTCGCAAACCTCGGTTTTGGAGAA  
Right inverted repeat for Ein(rev): AAGAGGTTTTGGCTCCAAACGCTCTT  
Right inverted repeat for Ein(rev-com): TTCTCCAAAACCGAGGTTTGCGAGAA (AAACCGAGGTTTGCGAGAA, AAACCGAGGTTTGCG)  
Right inverted repeat for Ein(comp): AAGAGCGTTTGGAGCCAAAACCTCTT

## CTV gene tree vs species tree comparison

Normalised Robinson-Foulds distance between the gene trees and the species tree.

## O-antigen coevolution

### O-antigen homology groups

#### ABC transporter dependent O-antigen synthesis

- *Wzm* (ABC transporter permease):
- *Wzt* (ABC transporter ATP binding protein):

#### G4C or O-antigen capsule cluster

- *Wza*: hg_22428205
- *Wzb*: hg_22428204
- *Wzc*: hg_22428201
- *WaaL* (ligase):
- *rfbB*: found & conserved, GCF_009873295.1
- *rfbA*: found & conserved
- *rfbC*: found & conserved

- ?: hg_22428182
- *GfcB*: hg_22428181
- *GfcC*: hg_22428167
- *YbjH*: hg_22428166

- *wecA* (hg_22429666)

#### *Wzx/Wzy*-dependent pathway

- *Wzx* (flipase):
- *Wxy* (O-antigen polymerase):
- *Wzz* (chain length regulator):

*Dickeya solani* genome : https://www.ncbi.nlm.nih.gov/nuccore/CP015137.1

GCF_002068115.1: genes in LPS biosynthesis

### O-antigen HGs

hg_22428205,hg_22428166

```bash
Rscript scripts/utils/HG_range_coordinates.R --hgs hg_22428205,hg_22428166 \
--out analysis/pangenome_v2/o_antigen/abc_transporter_hg_regions.tab
```
